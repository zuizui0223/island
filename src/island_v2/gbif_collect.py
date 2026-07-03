"""Collect succeeded GBIF block downloads into island-by-species occurrences.

A block download is an occurrence dump for a *buffered* regional catchment that
may span up to 125 islands. This module does the two things the v2 design
requires before any flora is trusted:

1. assign every returned occurrence coordinate back to the **original exact
   island polygon** (never the buffered query catchment), dropping records that
   fall in the buffer but on no real island; and
2. retain the observation-process fields needed to audit sampling effort,
   coordinate quality, record basis, and establishment status.

The spatial assignment and archive parsing are pure functions so they can be
unit-tested without GBIF; the CLI only adds download + ledger plumbing.
"""

from __future__ import annotations

import hashlib
import io
import json
import zipfile
from collections.abc import Iterator
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

EFFORT_COLUMNS = [
    "island_id",
    "n_records",
    "n_unique_gbif_ids",
    "n_species",
    "n_datasets",
    "n_preserved_specimen",
    "n_human_observation",
    "n_other_basis_of_record",
    "n_uncertainty_reported",
    "median_coordinate_uncertainty_m",
    "p90_coordinate_uncertainty_m",
    "n_coordinate_uncertainty_gt_1000m",
    "n_establishment_native",
    "n_establishment_introduced",
    "n_establishment_cultivated",
    "year_min",
    "year_max",
]


@app.callback()
def main() -> None:
    """Assign succeeded GBIF block downloads to exact island polygons."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tabular_member(archive: zipfile.ZipFile) -> str:
    names = [name for name in archive.namelist() if name.lower().endswith((".csv", ".txt", ".tsv"))]
    if not names:
        raise RuntimeError("GBIF archive contains no delimited occurrence file")
    # SIMPLE_CSV ships a single big occurrence table; pick the largest candidate.
    return max(names, key=lambda name: archive.getinfo(name).file_size)


def _normalise_occurrence_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Rename GBIF columns and retain only finite WGS84 coordinates."""
    n_source_rows = int(len(frame))
    frame = frame.rename(columns=_SOURCE_COLUMNS)
    for column in _SOURCE_COLUMNS.values():
        if column not in frame.columns:
            frame[column] = pd.NA
    lon = pd.to_numeric(frame["decimal_longitude"], errors="coerce")
    lat = pd.to_numeric(frame["decimal_latitude"], errors="coerce")
    frame["decimal_longitude"] = lon
    frame["decimal_latitude"] = lat
    valid = lon.notna() & lat.notna() & lon.between(-180, 180) & lat.between(-90, 90)
    result = frame.loc[valid].reset_index(drop=True)
    result.attrs["n_source_rows"] = n_source_rows
    return result


def iter_block_occurrences(archive_path: Path, chunksize: int = 250_000) -> Iterator[pd.DataFrame]:
    """Yield normalised GBIF SIMPLE_CSV rows in bounded-memory chunks.

    Every yielded frame has ``attrs['n_source_rows']`` so callers can retain a
    complete audit of source rows versus coordinate-valid rows. Empty frames are
    yielded when a source chunk contains no valid coordinates.
    """
    if chunksize < 1:
        raise typer.BadParameter("chunksize must be at least 1")
    with zipfile.ZipFile(archive_path) as archive:
        member = _tabular_member(archive)
        with archive.open(member) as handle:
            text = io.TextIOWrapper(handle, encoding="utf-8", newline="")
            reader = pd.read_csv(
                text,
                sep="\t",
                dtype=str,
                usecols=lambda column: column in _SOURCE_COLUMNS,
                on_bad_lines="skip",
                quoting=3,  # csv.QUOTE_NONE: GBIF SIMPLE_CSV is unquoted tab-delimited
                chunksize=chunksize,
            )
            for raw_chunk in reader:
                yield _normalise_occurrence_frame(raw_chunk)


def read_block_occurrences(archive_path: Path) -> pd.DataFrame:
    """Read a GBIF SIMPLE_CSV block archive into normalised occurrence rows.

    This convenience function is retained for tests and small pilot archives.
    Production collection uses :func:`iter_block_occurrences` so a single large
    regional archive does not need to fit in memory while it is parsed.
    """
    chunks = [chunk for chunk in iter_block_occurrences(archive_path) if not chunk.empty]
    if chunks:
        return pd.concat(chunks, ignore_index=True)
    return pd.DataFrame(columns=list(_SOURCE_COLUMNS.values()))


def assign_occurrences_to_islands(
    occurrences: pd.DataFrame,
    islands: gpd.GeoDataFrame,
) -> pd.DataFrame:
    """Assign coordinates to original exact island polygons, including boundaries.

    The GBIF request geometry is intentionally buffered, so a returned point is
    accepted only when it intersects an original, unbuffered island polygon.
    ``intersects`` includes a coordinate recorded exactly on the shoreline;
    previous strict ``within`` logic discarded those boundary points. A rare point
    touching two island polygons is assigned deterministically to the first
    lexicographic island ID and is never duplicated.
    """
    if islands.crs is None:
        raise typer.BadParameter("island geometry has no CRS")
    result = occurrences.reset_index(drop=True).copy()
    if result.empty:
        result["island_id"] = pd.Series(dtype="object")
        return result
    points = gpd.GeoDataFrame(
        result.copy(),
        geometry=gpd.points_from_xy(result["decimal_longitude"], result["decimal_latitude"]),
        crs=4326,
    )
    points["_source_row"] = range(len(points))
    islands_4326 = islands.to_crs(4326)[["island_id", "geometry"]]
    joined = gpd.sjoin(points, islands_4326, how="left", predicate="intersects")
    joined["_matched"] = joined["island_id"].notna().astype(int)
    joined = joined.sort_values(
        ["_source_row", "_matched", "island_id"],
        ascending=[True, False, True],
        na_position="last",
    )
    winners = joined.drop_duplicates("_source_row", keep="first").set_index("_source_row")["island_id"]
    result["island_id"] = [winners.get(index, pd.NA) for index in range(len(result))]
    return result


def deduplicate_occurrences(assigned: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Collapse cross-block copies of a GBIF occurrence deterministically.

    Regional request catchments can overlap. The same GBIF occurrence may thus
    occur in more than one completed archive. ``gbifID`` is the only deduplication
    key used; rows with no gbifID are deliberately retained because their identity
    cannot be established safely. When duplicate copies disagree, an exact-island
    assignment is preferred over a buffer-only copy, then lower coordinate
    uncertainty, then block ID for deterministic tie breaking.
    """
    if assigned.empty:
        return assigned.copy(), {
            "n_input_rows": 0,
            "n_unique_gbif_ids": 0,
            "n_duplicate_gbif_ids": 0,
            "n_duplicate_rows_removed": 0,
            "n_cross_block_duplicate_gbif_ids": 0,
            "n_conflicting_island_assignments": 0,
        }

    work = assigned.reset_index(drop=True).copy()
    work["_source_order"] = range(len(work))
    gbif_id = work["gbif_id"].fillna("").astype(str).str.strip()
    work["_has_gbif_id"] = gbif_id.ne("")
    work["_gbif_id_norm"] = gbif_id
    with_id = work.loc[work["_has_gbif_id"]].copy()
    without_id = work.loc[~work["_has_gbif_id"]].copy()

    if with_id.empty:
        result = work.copy()
        result["gbif_duplicate_group_size"] = pd.NA
        result["is_cross_block_duplicate"] = False
        result = result.drop(columns=["_source_order", "_has_gbif_id", "_gbif_id_norm"])
        return result, {
            "n_input_rows": int(len(work)),
            "n_unique_gbif_ids": 0,
            "n_duplicate_gbif_ids": 0,
            "n_duplicate_rows_removed": 0,
            "n_cross_block_duplicate_gbif_ids": 0,
            "n_conflicting_island_assignments": 0,
        }

    group_size = with_id.groupby("_gbif_id_norm")["_gbif_id_norm"].transform("size")
    with_id["gbif_duplicate_group_size"] = group_size.astype(int)
    duplicate_ids = with_id.loc[group_size.gt(1), "_gbif_id_norm"].unique().tolist()
    cross_block_ids = set(
        with_id.groupby("_gbif_id_norm")["block_id"].nunique().loc[lambda values: values.gt(1)].index
    )
    conflicting_islands = 0
    for _, group in with_id.loc[with_id["_gbif_id_norm"].isin(duplicate_ids)].groupby("_gbif_id_norm"):
        if group["island_id"].dropna().nunique() > 1:
            conflicting_islands += 1

    with_id["_has_island"] = with_id["island_id"].notna().astype(int)
    uncertainty = pd.to_numeric(with_id["coordinate_uncertainty_m"], errors="coerce")
    with_id["_uncertainty_sort"] = uncertainty.fillna(float("inf"))
    with_id = with_id.sort_values(
        ["_gbif_id_norm", "_has_island", "_uncertainty_sort", "block_id", "_source_order"],
        ascending=[True, False, True, True, True],
        na_position="last",
    )
    with_id["is_cross_block_duplicate"] = with_id["_gbif_id_norm"].isin(cross_block_ids)
    winners = with_id.drop_duplicates("_gbif_id_norm", keep="first")

    without_id["gbif_duplicate_group_size"] = pd.NA
    without_id["is_cross_block_duplicate"] = False
    result = pd.concat([winners, without_id], ignore_index=True).sort_values("_source_order").reset_index(drop=True)
    result = result.drop(
        columns=[
            "_source_order",
            "_has_gbif_id",
            "_gbif_id_norm",
            "_has_island",
            "_uncertainty_sort",
        ],
        errors="ignore",
    )
    audit = {
        "n_input_rows": int(len(work)),
        "n_unique_gbif_ids": int(with_id["_gbif_id_norm"].nunique()),
        "n_duplicate_gbif_ids": int(len(duplicate_ids)),
        "n_duplicate_rows_removed": int(len(work) - len(result)),
        "n_cross_block_duplicate_gbif_ids": int(len(cross_block_ids)),
        "n_conflicting_island_assignments": int(conflicting_islands),
    }
    return result, audit


def summarize_observation_effort(assigned: pd.DataFrame) -> pd.DataFrame:
    """Summarise island-specific sampling effort and record quality."""
    on_island = assigned.loc[assigned["island_id"].notna()].copy()
    if on_island.empty:
        return pd.DataFrame(columns=EFFORT_COLUMNS)

    on_island["_year"] = pd.to_numeric(on_island["year"], errors="coerce")
    on_island["_uncertainty"] = pd.to_numeric(on_island["coordinate_uncertainty_m"], errors="coerce")
    on_island["_basis"] = on_island["basis_of_record"].fillna("").astype(str).str.strip().str.upper()
    on_island["_establishment"] = on_island["establishment_means"].fillna("").astype(str).str.strip().str.upper()
    rows = []
    for island_id, group in on_island.groupby("island_id", sort=True):
        species = group["species"].dropna().astype(str).str.strip()
        species = species.loc[species.ne("")]
        gbif_ids = group["gbif_id"].dropna().astype(str).str.strip()
        gbif_ids = gbif_ids.loc[gbif_ids.ne("")]
        datasets = group["dataset_key"].dropna().astype(str).str.strip()
        datasets = datasets.loc[datasets.ne("")]
        years = group["_year"].dropna()
        uncertainty = group["_uncertainty"].dropna()
        basis = group["_basis"]
        establishment = group["_establishment"]
        rows.append(
            {
                "island_id": island_id,
                "n_records": int(len(group)),
                "n_unique_gbif_ids": int(gbif_ids.nunique()),
                "n_species": int(species.nunique()),
                "n_datasets": int(datasets.nunique()),
                "n_preserved_specimen": int(basis.eq("PRESERVED_SPECIMEN").sum()),
                "n_human_observation": int(basis.eq("HUMAN_OBSERVATION").sum()),
                "n_other_basis_of_record": int((basis.ne("") & ~basis.isin(["PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"])).sum()),
                "n_uncertainty_reported": int(uncertainty.notna().sum()),
                "median_coordinate_uncertainty_m": float(uncertainty.median()) if not uncertainty.empty else pd.NA,
                "p90_coordinate_uncertainty_m": float(uncertainty.quantile(0.9)) if not uncertainty.empty else pd.NA,
                "n_coordinate_uncertainty_gt_1000m": int(uncertainty.gt(1000).sum()),
                "n_establishment_native": int(establishment.str.contains("NATIVE", regex=False).sum()),
                "n_establishment_introduced": int(establishment.str.contains("INTRODUCED", regex=False).sum()),
                "n_establishment_cultivated": int(establishment.str.contains("CULTIVATED", regex=False).sum()),
                "year_min": int(years.min()) if not years.empty else pd.NA,
                "year_max": int(years.max()) if not years.empty else pd.NA,
            }
        )
    return pd.DataFrame(rows, columns=EFFORT_COLUMNS)


def island_species_table(assigned: pd.DataFrame) -> pd.DataFrame:
    """Deduplicated island x species candidate table with retained provenance."""
    on_island = assigned.loc[assigned["island_id"].notna()].copy()
    named = on_island.loc[on_island["species"].fillna("").astype(str).str.strip().astype(bool)].copy()
    if named.empty:
        return pd.DataFrame(
            columns=["island_id", "species", "n_records", "n_unique_gbif_ids", "basis_of_record_set", "review_status"]
        )
    rows = []
    for (island_id, species), group in named.groupby(["island_id", "species"], sort=True):
        bases = sorted(set(group["basis_of_record"].dropna().astype(str).str.strip()) - {""})
        gbif_ids = group["gbif_id"].dropna().astype(str).str.strip()
        gbif_ids = gbif_ids.loc[gbif_ids.ne("")]
        rows.append(
            {
                "island_id": island_id,
                "species": species,
                "n_records": int(len(group)),
                "n_unique_gbif_ids": int(gbif_ids.nunique()),
                "basis_of_record_set": "|".join(bases),
                "review_status": "unresolved_taxonomy",
            }
        )
    return pd.DataFrame(rows)


def island_taxa_table(assigned: pd.DataFrame) -> pd.DataFrame:
    """Unique accepted-species taxa list that feeds ``island-v2-traits run``."""
    on_island = assigned.loc[assigned["island_id"].notna()].copy()
    named = on_island.loc[on_island["species"].fillna("").astype(str).str.strip().astype(bool)].copy()
    if named.empty:
        return pd.DataFrame(columns=["accepted_species", "genus", "family", "n_islands", "n_records"])
    rows = []
    for species, group in named.groupby("species", sort=True):
        genus = next((value for value in group["genus"].dropna() if str(value).strip()), "")
        family = next((value for value in group["family"].dropna() if str(value).strip()), "")
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
    return [
        entry
        for entry in campaign.get("ledger", [])
        if str(entry.get("request_status")) == "succeeded" and str(entry.get("download_key", "")).strip()
    ]


def _download_archive(
    client: httpx.Client,
    download_key: str,
    download_dir: Path,
    link: str | None,
) -> Path:
    """Download a GBIF archive atomically without loading it all into RAM."""
    destination = download_dir / f"{download_key}.zip"
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    url = link or f"https://api.gbif.org/v1/occurrence/download/request/{download_key}.zip"
    temporary = destination.with_suffix(".zip.part")
    with client.stream("GET", url) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for chunk in response.iter_bytes():
                handle.write(chunk)
    temporary.replace(destination)
    return destination


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True, help="Campaign ledger JSON."),
    block_members_csv: Path = typer.Option(..., exists=True, help="gbif_block_members.csv."),
    islands_gpkg: Path = typer.Option(..., exists=True, help="Exact island polygons GeoPackage."),
    download_dir: Path = typer.Option(..., help="Where block archives are cached."),
    output_dir: Path = typer.Option(..., help="Where island-level outputs are written."),
    chunksize: int = typer.Option(250_000, min=1, help="Rows read per GBIF SIMPLE_CSV chunk."),
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
    manifest_rows: list[dict[str, Any]] = []
    failures: list[str] = []
    with httpx.Client(timeout=300.0, follow_redirects=True, headers={"User-Agent": "island-floral-v2/0.4"}) as client:
        for entry in ready:
            block_id = str(entry["block_id"])
            key = str(entry["download_key"])
            block_islands = members.loc[members["block_id"].eq(block_id), "analysis_island_id"].unique()
            present = [island_id for island_id in block_islands if island_id in islands_by_id.index]
            metrics: dict[str, Any] = {
                "block_id": block_id,
                "download_key": key,
                "doi": str(entry.get("doi", "")),
                "n_block_member_islands": int(len(block_islands)),
                "n_islands_available_for_assignment": int(len(present)),
                "archive_path": "",
                "archive_sha256": "",
                "archive_bytes": 0,
                "n_source_rows": 0,
                "n_valid_coordinate_rows": 0,
                "n_exact_island_rows_before_dedup": 0,
                "n_buffer_only_rows_before_dedup": 0,
                "n_rows_retained_after_global_dedup": 0,
                "collection_status": "pending",
                "error": "",
            }
            if not present:
                metrics["collection_status"] = "skipped_no_matching_islands"
                manifest_rows.append(metrics)
                continue
            try:
                island_subset = islands_by_id.loc[present].reset_index()
                island_subset = gpd.GeoDataFrame(island_subset, geometry="geometry", crs=islands.crs)
                archive = _download_archive(client, key, download_dir, entry.get("last_download_link") or entry.get("download_link"))
                metrics["archive_path"] = str(archive)
                metrics["archive_sha256"] = sha256_file(archive)
                metrics["archive_bytes"] = int(archive.stat().st_size)
                for occurrences in iter_block_occurrences(archive, chunksize=chunksize):
                    metrics["n_source_rows"] += int(occurrences.attrs.get("n_source_rows", len(occurrences)))
                    metrics["n_valid_coordinate_rows"] += int(len(occurrences))
                    if occurrences.empty:
                        continue
                    assigned = assign_occurrences_to_islands(occurrences, island_subset)
                    assigned["block_id"] = block_id
                    assigned["gbif_download_key"] = key
                    metrics["n_exact_island_rows_before_dedup"] += int(assigned["island_id"].notna().sum())
                    metrics["n_buffer_only_rows_before_dedup"] += int(assigned["island_id"].isna().sum())
                    all_assigned.append(assigned)
                metrics["collection_status"] = "collected"
            except Exception as exc:  # retain a block-level audit rather than silently losing it
                metrics["collection_status"] = "collection_failed"
                metrics["error"] = f"{type(exc).__name__}: {exc}"
                failures.append(block_id)
            manifest_rows.append(metrics)

    if all_assigned:
        assigned_before_dedup = pd.concat(all_assigned, ignore_index=True)
        assigned, duplicate_audit = deduplicate_occurrences(assigned_before_dedup)
    else:
        assigned_before_dedup = pd.DataFrame(columns=list(_SOURCE_COLUMNS.values()) + ["island_id", "block_id", "gbif_download_key"])
        assigned = assigned_before_dedup.copy()
        duplicate_audit = {
            "n_input_rows": 0,
            "n_unique_gbif_ids": 0,
            "n_duplicate_gbif_ids": 0,
            "n_duplicate_rows_removed": 0,
            "n_cross_block_duplicate_gbif_ids": 0,
            "n_conflicting_island_assignments": 0,
        }

    if manifest_rows and not assigned.empty:
        retained_by_block = assigned.groupby("block_id").size().to_dict()
        for metrics in manifest_rows:
            metrics["n_rows_retained_after_global_dedup"] = int(retained_by_block.get(metrics["block_id"], 0))

    species_table = island_species_table(assigned)
    effort = summarize_observation_effort(assigned)
    taxa = island_taxa_table(assigned)
    manifest = pd.DataFrame(manifest_rows)
    if manifest.empty:
        manifest = pd.DataFrame(
            columns=[
                "block_id", "download_key", "doi", "n_block_member_islands", "n_islands_available_for_assignment",
                "archive_path", "archive_sha256", "archive_bytes", "n_source_rows", "n_valid_coordinate_rows",
                "n_exact_island_rows_before_dedup", "n_buffer_only_rows_before_dedup",
                "n_rows_retained_after_global_dedup", "collection_status", "error",
            ]
        )

    species_table.to_csv(output_dir / "island_species_occurrences.csv", index=False)
    effort.to_csv(output_dir / "island_observation_effort.csv", index=False)
    taxa.to_csv(output_dir / "island_taxa.csv", index=False)
    manifest.to_csv(output_dir / "collection_manifest.csv", index=False)
    status = {
        "n_succeeded_blocks": len(ready),
        "n_collected_blocks": int(manifest["collection_status"].eq("collected").sum()) if not manifest.empty else 0,
        "n_failed_blocks": int(manifest["collection_status"].eq("collection_failed").sum()) if not manifest.empty else 0,
        "n_occurrences_before_global_dedup": int(len(assigned_before_dedup)),
        "n_occurrences_after_global_dedup": int(len(assigned)),
        "n_occurrences_on_island": int(assigned["island_id"].notna().sum()) if "island_id" in assigned else 0,
        "n_occurrences_in_buffer_only": int(assigned["island_id"].isna().sum()) if "island_id" in assigned else 0,
        "n_island_species_pairs": int(len(species_table)),
        "n_islands_with_records": int(effort["island_id"].nunique()) if not effort.empty else 0,
        "n_accepted_species": int(len(taxa)),
        "duplicate_audit": duplicate_audit,
        "failed_block_ids": failures,
    }
    (output_dir / "collection_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    typer.echo(
        f"Collected {status['n_collected_blocks']} blocks: {len(species_table)} island-species pairs, "
        f"{len(taxa)} accepted species across {status['n_islands_with_records']} islands "
        f"({status['n_occurrences_in_buffer_only']} occurrences fell in buffer only; "
        f"{duplicate_audit['n_duplicate_rows_removed']} duplicate rows removed)."
    )
    if failures:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
