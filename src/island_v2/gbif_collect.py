"""Collect succeeded GBIF block downloads into island-by-species occurrences.

A block download is an occurrence dump for a *buffered* regional catchment that
may span up to 125 islands. This module does the two things the v2 design
requires before any flora is trusted:

1. assign every returned occurrence coordinate back to the **original exact
   island polygon** (never the buffered query catchment), dropping records that
   fall in the buffer but on no real island; and
2. retain the observation-process fields (basis of record, year, coordinate
   uncertainty, dataset, establishment means) so the later audit can tell
   specimen-backed island membership from a single iNaturalist point or a
   cultivated record.

The spatial assignment and archive parsing are pure functions so they can be
unit-tested without GBIF; the CLI only adds download + ledger plumbing.
"""

from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

# GBIF SIMPLE_CSV columns we keep. Names are GBIF's; everything else is dropped.
_SOURCE_COLUMNS = {
    "gbifID": "gbif_id",
    "datasetKey": "dataset_key",
    "family": "family",
    "genus": "genus",
    "species": "species",
    "scientificName": "scientific_name",
    "taxonRank": "taxon_rank",
    "taxonKey": "taxon_key",
    "speciesKey": "species_key",
    "decimalLatitude": "decimal_latitude",
    "decimalLongitude": "decimal_longitude",
    "coordinateUncertaintyInMeters": "coordinate_uncertainty_m",
    "year": "year",
    "basisOfRecord": "basis_of_record",
    "occurrenceStatus": "occurrence_status",
    "establishmentMeans": "establishment_means",
}


@app.callback()
def main() -> None:
    """Assign succeeded GBIF block downloads to exact island polygons."""


def _tabular_member(archive: zipfile.ZipFile) -> str:
    names = [n for n in archive.namelist() if n.lower().endswith((".csv", ".txt", ".tsv"))]
    if not names:
        raise RuntimeError("GBIF archive contains no delimited occurrence file")
    # SIMPLE_CSV ships a single big occurrence table; pick the largest candidate.
    return max(names, key=lambda n: archive.getinfo(n).file_size)


def read_block_occurrences(archive_path: Path) -> pd.DataFrame:
    """Read a GBIF SIMPLE_CSV block archive into normalised occurrence rows.

    SIMPLE_CSV is tab-separated. Only the columns in `_SOURCE_COLUMNS` are kept,
    and rows without a finite coordinate are dropped (they cannot be assigned to
    an island and would otherwise silently vanish or misassign).
    """
    with zipfile.ZipFile(archive_path) as archive:
        member = _tabular_member(archive)
        with archive.open(member) as handle:
            text = io.TextIOWrapper(handle, encoding="utf-8", newline="")
            frame = pd.read_csv(
                text,
                sep="\t",
                dtype=str,
                usecols=lambda c: c in _SOURCE_COLUMNS,
                on_bad_lines="skip",
                quoting=3,  # csv.QUOTE_NONE: GBIF SIMPLE_CSV is unquoted tab-delimited
            )
    frame = frame.rename(columns=_SOURCE_COLUMNS)
    for column in _SOURCE_COLUMNS.values():
        if column not in frame.columns:
            frame[column] = pd.NA
    lon = pd.to_numeric(frame["decimal_longitude"], errors="coerce")
    lat = pd.to_numeric(frame["decimal_latitude"], errors="coerce")
    frame["decimal_longitude"] = lon
    frame["decimal_latitude"] = lat
    valid = lon.notna() & lat.notna() & lon.between(-180, 180) & lat.between(-90, 90)
    return frame.loc[valid].reset_index(drop=True)


def assign_occurrences_to_islands(
    occurrences: pd.DataFrame,
    islands: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """Point-in-polygon assign occurrences to the exact island polygons.

    `islands` must carry `island_id` and exact (unbuffered) geometry in EPSG:4326.
    Returns every input occurrence with an `island_id` column that is null when
    the point fell in the query buffer but inside no real island. A point on a
    shared boundary of two islands (which should not happen for disjoint islands)
    keeps only its first match so records are never duplicated.
    """
    if islands.crs is None:
        raise typer.BadParameter("island geometry has no CRS")
    points = gpd.GeoDataFrame(
        occurrences.copy(),
        geometry=gpd.points_from_xy(
            occurrences["decimal_longitude"], occurrences["decimal_latitude"]
        ),
        crs=4326,
    )
    islands_4326 = islands.to_crs(4326)[["island_id", "geometry"]]
    joined = gpd.sjoin(points, islands_4326, how="left", predicate="within")
    joined = joined[~joined.index.duplicated(keep="first")]
    result = occurrences.copy()
    result["island_id"] = joined["island_id"].to_numpy()
    return result


def summarize_observation_effort(assigned: pd.DataFrame) -> pd.DataFrame:
    """Per-island observation-process summary for the later coverage audit."""
    on_island = assigned.loc[assigned["island_id"].notna()].copy()
    if on_island.empty:
        return pd.DataFrame(
            columns=[
                "island_id",
                "n_records",
                "n_species",
                "n_datasets",
                "n_preserved_specimen",
                "year_min",
                "year_max",
            ]
        )
    year = pd.to_numeric(on_island["year"], errors="coerce")
    on_island["_year"] = year
    is_specimen = on_island["basis_of_record"].eq("PRESERVED_SPECIMEN")
    rows = []
    for island_id, group in on_island.groupby("island_id", sort=True):
        species = group["species"].dropna()
        species = species[species.str.strip().astype(bool)]
        years = group["_year"].dropna()
        rows.append(
            {
                "island_id": island_id,
                "n_records": int(len(group)),
                "n_species": int(species.nunique()),
                "n_datasets": int(group["dataset_key"].nunique()),
                "n_preserved_specimen": int(is_specimen.loc[group.index].sum()),
                "year_min": int(years.min()) if not years.empty else pd.NA,
                "year_max": int(years.max()) if not years.empty else pd.NA,
            }
        )
    return pd.DataFrame(rows)


def island_species_table(assigned: pd.DataFrame) -> pd.DataFrame:
    """Deduplicated island x species candidate table with retained provenance."""
    on_island = assigned.loc[assigned["island_id"].notna()].copy()
    named = on_island.loc[on_island["species"].fillna("").str.strip().astype(bool)].copy()
    if named.empty:
        return pd.DataFrame(
            columns=["island_id", "species", "n_records", "basis_of_record_set", "review_status"]
        )
    rows = []
    for (island_id, species), group in named.groupby(["island_id", "species"], sort=True):
        bases = sorted(set(group["basis_of_record"].dropna()))
        rows.append(
            {
                "island_id": island_id,
                "species": species,
                "n_records": int(len(group)),
                "basis_of_record_set": "|".join(bases),
                "review_status": "unresolved_taxonomy",
            }
        )
    return pd.DataFrame(rows)


def island_taxa_table(assigned: pd.DataFrame) -> pd.DataFrame:
    """Unique accepted-species taxa list that feeds `island-v2-traits run`.

    trait_extraction requires exactly `accepted_species, genus, family`, so this
    is the direct hand-off from occurrence collection to trait acquisition. Only
    records assigned to a real island and resolved to species rank are used, and
    genus/family are taken from GBIF's backbone-normalised columns.
    """
    on_island = assigned.loc[assigned["island_id"].notna()].copy()
    named = on_island.loc[on_island["species"].fillna("").str.strip().astype(bool)].copy()
    if named.empty:
        return pd.DataFrame(columns=["accepted_species", "genus", "family", "n_islands", "n_records"])
    rows = []
    for species, group in named.groupby("species", sort=True):
        genus = next((g for g in group["genus"].dropna() if str(g).strip()), "")
        family = next((f for f in group["family"].dropna() if str(f).strip()), "")
        rows.append(
            {
                "accepted_species": species,
                "genus": genus,
                "family": family,
                "n_islands": int(group["island_id"].nunique()),
                "n_records": int(len(group)),
            }
        )
    return pd.DataFrame(rows)


def _succeeded_blocks(campaign: dict[str, Any]) -> list[dict[str, Any]]:
    ready = []
    for entry in campaign.get("ledger", []):
        if str(entry.get("request_status")) == "succeeded" and str(entry.get("download_key", "")).strip():
            ready.append(entry)
    return ready


def _download_archive(
    client: httpx.Client, download_key: str, download_dir: Path, link: str | None
) -> Path:
    destination = download_dir / f"{download_key}.zip"
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    url = link or f"https://api.gbif.org/v1/occurrence/download/request/{download_key}.zip"
    response = client.get(url)
    response.raise_for_status()
    destination.write_bytes(response.content)
    return destination


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True, help="Campaign ledger JSON."),
    block_members_csv: Path = typer.Option(..., exists=True, help="gbif_block_members.csv."),
    islands_gpkg: Path = typer.Option(..., exists=True, help="Exact island polygons GeoPackage."),
    download_dir: Path = typer.Option(..., help="Where block archives are cached."),
    output_dir: Path = typer.Option(..., help="Where island-level outputs are written."),
) -> None:
    """Download succeeded blocks and assign occurrences to exact island polygons."""
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    members = pd.read_csv(block_members_csv, dtype=str).fillna("")
    if "analysis_island_id" not in members.columns:
        members["analysis_island_id"] = members["island_id"]
    islands = gpd.read_file(islands_gpkg, layer="islands")
    islands_by_id = islands.set_index("island_id")

    ready = _succeeded_blocks(campaign)
    download_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    all_assigned: list[pd.DataFrame] = []
    with httpx.Client(timeout=300.0, follow_redirects=True, headers={"User-Agent": "island-floral-v2/0.3"}) as client:
        for entry in ready:
            block_id = str(entry["block_id"])
            key = str(entry["download_key"])
            block_islands = members.loc[members["block_id"].eq(block_id), "analysis_island_id"].unique()
            present = [i for i in block_islands if i in islands_by_id.index]
            if not present:
                continue
            island_subset = islands_by_id.loc[present].reset_index()
            island_subset = gpd.GeoDataFrame(island_subset, geometry="geometry", crs=islands.crs)
            archive = _download_archive(client, key, download_dir, entry.get("last_download_link") or entry.get("download_link"))
            occurrences = read_block_occurrences(archive)
            assigned = assign_occurrences_to_islands(occurrences, island_subset)
            assigned["block_id"] = block_id
            assigned["gbif_download_key"] = key
            all_assigned.append(assigned)

    if not all_assigned:
        typer.echo("No succeeded blocks are ready to collect.")
        (output_dir / "collection_status.json").write_text(
            json.dumps({"n_succeeded_blocks": len(ready), "n_collected_blocks": 0}, indent=2),
            encoding="utf-8",
        )
        return

    assigned = pd.concat(all_assigned, ignore_index=True)
    species_table = island_species_table(assigned)
    effort = summarize_observation_effort(assigned)
    taxa = island_taxa_table(assigned)
    n_unassigned = int(assigned["island_id"].isna().sum())

    species_table.to_csv(output_dir / "island_species_occurrences.csv", index=False)
    effort.to_csv(output_dir / "island_observation_effort.csv", index=False)
    taxa.to_csv(output_dir / "island_taxa.csv", index=False)
    (output_dir / "collection_status.json").write_text(
        json.dumps(
            {
                "n_succeeded_blocks": len(ready),
                "n_collected_blocks": len(all_assigned),
                "n_occurrences_on_island": int(assigned["island_id"].notna().sum()),
                "n_occurrences_in_buffer_only": n_unassigned,
                "n_island_species_pairs": int(len(species_table)),
                "n_islands_with_records": int(effort["island_id"].nunique()) if not effort.empty else 0,
                "n_accepted_species": int(len(taxa)),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    typer.echo(
        f"Collected {len(all_assigned)} blocks: {len(species_table)} island-species pairs, "
        f"{len(taxa)} accepted species across "
        f"{effort['island_id'].nunique() if not effort.empty else 0} islands "
        f"({n_unassigned} occurrences fell in buffer only)."
    )


if __name__ == "__main__":
    app()
