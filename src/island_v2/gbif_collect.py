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
import sys
import time
import zipfile
from collections.abc import Iterator
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer


def _progress(message: str) -> None:
    """Emit a flushed, timestamped progress line to stderr for CI visibility."""
    print(f"[collect {time.strftime('%H:%M:%S')}] {message}", file=sys.stderr, flush=True)

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
    islands_4326 = islands.to_crs(4326)[["island_id", "geometry"]]
    result["island_id"] = pd.NA
    # A point outside the combined bounding box of this block's islands cannot
    # intersect any of them, so it is buffer-only by definition. Skipping the
    # spatial join for those points is exact (not an approximation) and avoids
    # joining the millions of mainland points in a dense northern catchment.
    lon = pd.to_numeric(result["decimal_longitude"], errors="coerce")
    lat = pd.to_numeric(result["decimal_latitude"], errors="coerce")
    minx, miny, maxx, maxy = islands_4326.total_bounds
    in_bounds = lon.between(minx, maxx) & lat.between(miny, maxy)
    candidate_positions = result.index[in_bounds.to_numpy()]
    if len(candidate_positions):
        points = gpd.GeoDataFrame(
            {"_source_row": candidate_positions},
            geometry=gpd.points_from_xy(lon.loc[candidate_positions], lat.loc[candidate_positions]),
            crs=4326,
        )
        joined = gpd.sjoin(points, islands_4326, how="left", predicate="intersects")
        joined["_matched"] = joined["island_id"].notna().astype(int)
        joined = joined.sort_values(
            ["_source_row", "_matched", "island_id"],
            ascending=[True, False, True],
            na_position="last",
        )
        winners = joined.drop_duplicates("_source_row", keep="first").set_index("_source_row")["island_id"]
        matched = winners.dropna()
        if not matched.empty:
            result.loc[matched.index, "island_id"] = matched
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
    duplicate_rows = with_id.loc[with_id["_gbif_id_norm"].isin(duplicate_ids)]
    island_nunique = duplicate_rows.groupby("_gbif_id_norm")["island_id"].nunique()
    conflicting_islands = int((island_nunique > 1).sum())

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

    basis = on_island["basis_of_record"].fillna("").astype(str).str.strip().str.upper()
    establishment = on_island["establishment_means"].fillna("").astype(str).str.strip().str.upper()
    unc = pd.to_numeric(on_island["coordinate_uncertainty_m"], errors="coerce")

    def _nonblank(column: str) -> pd.Series:
        text = on_island[column].fillna("").astype(str).str.strip()
        return text.where(text.ne(""), other=pd.NA)

    work = pd.DataFrame(
        {
            "island_id": on_island["island_id"].to_numpy(),
            "_year": pd.to_numeric(on_island["year"], errors="coerce").to_numpy(),
            "_uncertainty": unc.to_numpy(),
            "_species": _nonblank("species").to_numpy(),
            "_gbif": _nonblank("gbif_id").to_numpy(),
            "_dataset": _nonblank("dataset_key").to_numpy(),
            "_preserved": basis.eq("PRESERVED_SPECIMEN").to_numpy(),
            "_human": basis.eq("HUMAN_OBSERVATION").to_numpy(),
            "_other_basis": (basis.ne("") & ~basis.isin(["PRESERVED_SPECIMEN", "HUMAN_OBSERVATION"])).to_numpy(),
            "_unc_reported": unc.notna().to_numpy(),
            "_unc_gt1000": unc.gt(1000).to_numpy(),
            "_native": establishment.str.contains("NATIVE", regex=False).to_numpy(),
            "_introduced": establishment.str.contains("INTRODUCED", regex=False).to_numpy(),
            "_cultivated": establishment.str.contains("CULTIVATED", regex=False).to_numpy(),
        }
    )
    grouped = work.groupby("island_id", sort=True)
    out = grouped.agg(
        n_records=("_year", "size"),
        n_unique_gbif_ids=("_gbif", "nunique"),
        n_species=("_species", "nunique"),
        n_datasets=("_dataset", "nunique"),
        n_preserved_specimen=("_preserved", "sum"),
        n_human_observation=("_human", "sum"),
        n_other_basis_of_record=("_other_basis", "sum"),
        n_uncertainty_reported=("_unc_reported", "sum"),
        median_coordinate_uncertainty_m=("_uncertainty", "median"),
        p90_coordinate_uncertainty_m=("_uncertainty", lambda values: values.quantile(0.9)),
        n_coordinate_uncertainty_gt_1000m=("_unc_gt1000", "sum"),
        n_establishment_native=("_native", "sum"),
        n_establishment_introduced=("_introduced", "sum"),
        n_establishment_cultivated=("_cultivated", "sum"),
        year_min=("_year", "min"),
        year_max=("_year", "max"),
    ).reset_index()
    int_columns = [
        "n_records", "n_unique_gbif_ids", "n_species", "n_datasets", "n_preserved_specimen",
        "n_human_observation", "n_other_basis_of_record", "n_uncertainty_reported",
        "n_coordinate_uncertainty_gt_1000m", "n_establishment_native",
        "n_establishment_introduced", "n_establishment_cultivated",
    ]
    for column in int_columns:
        out[column] = out[column].astype(int)
    for column in ("year_min", "year_max"):
        out[column] = out[column].map(lambda value: int(value) if pd.notna(value) else pd.NA)
    for column in ("median_coordinate_uncertainty_m", "p90_coordinate_uncertainty_m"):
        out[column] = out[column].map(lambda value: float(value) if pd.notna(value) else pd.NA)
    return out[EFFORT_COLUMNS]


def island_species_table(assigned: pd.DataFrame) -> pd.DataFrame:
    """Deduplicated island x species candidate table with retained provenance."""
    on_island = assigned.loc[assigned["island_id"].notna()]
    species = on_island["species"].fillna("").astype(str).str.strip()
    named = on_island.loc[species.ne("")]
    if named.empty:
        return pd.DataFrame(
            columns=["island_id", "species", "n_records", "n_unique_gbif_ids", "basis_of_record_set", "review_status"]
        )
    gbif = named["gbif_id"].fillna("").astype(str).str.strip()
    work = pd.DataFrame(
        {
            "island_id": named["island_id"].to_numpy(),
            "species": species.loc[named.index].to_numpy(),
            "_gbif": gbif.where(gbif.ne(""), other=pd.NA).to_numpy(),
            "_basis": named["basis_of_record"].fillna("").astype(str).str.strip().to_numpy(),
        }
    )
    out = work.groupby(["island_id", "species"], sort=True).agg(
        n_records=("species", "size"),
        n_unique_gbif_ids=("_gbif", "nunique"),
        basis_of_record_set=("_basis", lambda values: "|".join(sorted(set(values) - {""}))),
    ).reset_index()
    out["n_records"] = out["n_records"].astype(int)
    out["n_unique_gbif_ids"] = out["n_unique_gbif_ids"].astype(int)
    out["review_status"] = "unresolved_taxonomy"
    return out[["island_id", "species", "n_records", "n_unique_gbif_ids", "basis_of_record_set", "review_status"]]


def island_taxa_table(assigned: pd.DataFrame) -> pd.DataFrame:
    """Unique accepted-species taxa list that feeds ``island-v2-traits run``."""
    on_island = assigned.loc[assigned["island_id"].notna()]
    species = on_island["species"].fillna("").astype(str).str.strip()
    named = on_island.loc[species.ne("")]
    if named.empty:
        return pd.DataFrame(columns=["accepted_species", "genus", "family", "n_islands", "n_records"])

    def _first_nonblank(values: pd.Series) -> str:
        return next((value for value in values.dropna() if str(value).strip()), "")

    work = pd.DataFrame(
        {
            "species": species.loc[named.index].to_numpy(),
            "genus": named["genus"].to_numpy(),
            "family": named["family"].to_numpy(),
            "island_id": named["island_id"].to_numpy(),
        }
    )
    out = work.groupby("species", sort=True).agg(
        genus=("genus", _first_nonblank),
        family=("family", _first_nonblank),
        n_islands=("island_id", "nunique"),
        n_records=("species", "size"),
    ).reset_index().rename(columns={"species": "accepted_species"})
    out["n_islands"] = out["n_islands"].astype(int)
    out["n_records"] = out["n_records"].astype(int)
    return out[["accepted_species", "genus", "family", "n_islands", "n_records"]]


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
                block_started = time.monotonic()
                _progress(f"block {block_id}: {len(present)} islands, downloading {key}")
                archive = _download_archive(client, key, download_dir, entry.get("last_download_link") or entry.get("download_link"))
                metrics["archive_path"] = str(archive)
                metrics["archive_sha256"] = sha256_file(archive)
                metrics["archive_bytes"] = int(archive.stat().st_size)
                _progress(f"block {block_id}: archive {metrics['archive_bytes'] / 1e6:.0f} MB; assigning chunks")
                chunk_index = 0
                for occurrences in iter_block_occurrences(archive, chunksize=chunksize):
                    chunk_index += 1
                    metrics["n_source_rows"] += int(occurrences.attrs.get("n_source_rows", len(occurrences)))
                    metrics["n_valid_coordinate_rows"] += int(len(occurrences))
                    if occurrences.empty:
                        continue
                    assigned = assign_occurrences_to_islands(occurrences, island_subset)
                    assigned["block_id"] = block_id
                    assigned["gbif_download_key"] = key
                    on_island_mask = assigned["island_id"].notna()
                    metrics["n_exact_island_rows_before_dedup"] += int(on_island_mask.sum())
                    metrics["n_buffer_only_rows_before_dedup"] += int((~on_island_mask).sum())
                    # Only on-island rows contribute to any island output; buffer-only
                    # rows are counted for provenance but never retained, so a dense
                    # mainland-heavy northern block no longer accumulates millions of
                    # rows in memory.
                    on_island = assigned.loc[on_island_mask]
                    if not on_island.empty:
                        all_assigned.append(on_island)
                    if chunk_index % 10 == 0:
                        _progress(
                            f"block {block_id}: chunk {chunk_index}, "
                            f"source_rows={metrics['n_source_rows']}, "
                            f"on_island={metrics['n_exact_island_rows_before_dedup']}, "
                            f"{time.monotonic() - block_started:.0f}s"
                        )
                metrics["collection_status"] = "collected"
                _progress(
                    f"block {block_id}: assigned {metrics['n_source_rows']} source rows in "
                    f"{chunk_index} chunks, {time.monotonic() - block_started:.0f}s"
                )
            except Exception as exc:  # retain a block-level audit rather than silently losing it
                metrics["collection_status"] = "collection_failed"
                metrics["error"] = f"{type(exc).__name__}: {exc}"
                failures.append(block_id)
            manifest_rows.append(metrics)

    if all_assigned:
        _progress(f"concatenating {sum(len(part) for part in all_assigned)} on-island rows")
        phase_started = time.monotonic()
        assigned_before_dedup = pd.concat(all_assigned, ignore_index=True)
        assigned, duplicate_audit = deduplicate_occurrences(assigned_before_dedup)
        _progress(f"dedup -> {len(assigned)} rows, {time.monotonic() - phase_started:.0f}s")
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

    tables_started = time.monotonic()
    species_table = island_species_table(assigned)
    effort = summarize_observation_effort(assigned)
    taxa = island_taxa_table(assigned)
    _progress(f"built island tables, {time.monotonic() - tables_started:.0f}s")
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
        # Buffer-only rows are dropped before accumulation, so sum their per-block
        # provenance counts rather than counting retained NA-island rows.
        "n_occurrences_in_buffer_only": int(
            sum(int(row.get("n_buffer_only_rows_before_dedup", 0)) for row in manifest_rows)
        ),
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
