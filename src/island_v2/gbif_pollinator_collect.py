"""Collect GBIF pollinator downloads into exact-island occurrence records.

The existing plant collector intentionally writes compact flora products. The
Bombus diagnostic needs a different retained product: one quality-auditable
record table with family, genus, dataset, year, coordinates, record basis, and
exact island assignment. This collector reuses the same exact-polygon and
global-deduplication logic without reinterpreting a missing Bombus record as
absence.

A pollinator campaign should be prepared with the generic ``island-v2-gbif-
blocks prepare`` command for the declared target group (currently Apidae). This
collector operates only on succeeded downloads recorded in that campaign ledger.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import pandas as pd
import typer

from island_v2.gbif_collect import (
    _download_archive,
    _succeeded_blocks,
    assign_occurrences_to_islands,
    deduplicate_occurrences,
    iter_block_occurrences,
    sha256_file,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Collect exact-island pollinator occurrences for Bombus diagnostics."""


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def exact_island_pollinator_records(assigned: pd.DataFrame, expected_family: str) -> pd.DataFrame:
    """Keep exact-island rows and add a transparent expected-family diagnostic.

    Rows are not discarded merely because GBIF's family field is blank or differs
    from the expected campaign family. Keeping them preserves a direct audit of
    taxonomy and download scope; downstream Bombus diagnostics decide which rows
    meet their quality and target-group requirements.
    """
    result = assigned.loc[assigned["island_id"].notna()].copy()
    if result.empty:
        empty = assigned.head(0).copy()
        empty["expected_family"] = pd.Series(dtype="object")
        empty["expected_family_match"] = pd.Series(dtype="bool")
        return empty
    family = _text(result.get("family", pd.Series("", index=result.index))).str.upper()
    expected = str(expected_family).strip().upper()
    if not expected:
        raise typer.BadParameter("expected_family must be non-empty")
    result["expected_family"] = expected
    result["expected_family_match"] = family.eq(expected)
    return result.sort_values(["island_id", "gbif_id", "block_id"], na_position="last").reset_index(drop=True)


def summarize_pollinator_records(records: pd.DataFrame, expected_family: str) -> pd.DataFrame:
    """Summarise retained record coverage without classifying non-detection."""
    if records.empty:
        return pd.DataFrame(
            columns=[
                "island_id",
                "n_exact_island_records",
                "n_expected_family_records",
                "n_bombus_records",
                "n_datasets",
                "year_min",
                "year_max",
            ]
        )
    rows: list[dict[str, Any]] = []
    expected = str(expected_family).strip().upper()
    for island_id, group in records.groupby("island_id", sort=True):
        family = _text(group.get("family", pd.Series("", index=group.index))).str.upper()
        genus = _text(group.get("genus", pd.Series("", index=group.index))).str.upper()
        datasets = _text(group.get("dataset_key", pd.Series("", index=group.index)))
        years = pd.to_numeric(group.get("year", pd.Series(dtype=object)), errors="coerce").dropna()
        rows.append(
            {
                "island_id": island_id,
                "n_exact_island_records": int(len(group)),
                "n_expected_family_records": int(family.eq(expected).sum()),
                "n_bombus_records": int((family.eq(expected) & genus.eq("BOMBUS")).sum()),
                "n_datasets": int(datasets.loc[datasets.ne("")].nunique()),
                "year_min": int(years.min()) if not years.empty else pd.NA,
                "year_max": int(years.max()) if not years.empty else pd.NA,
            }
        )
    return pd.DataFrame(rows)


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True, help="Pollinator campaign ledger JSON."),
    block_members_csv: Path = typer.Option(..., exists=True, help="Pollinator GBIF block membership CSV."),
    islands_gpkg: Path = typer.Option(..., exists=True, help="Frozen exact island polygons GeoPackage."),
    download_dir: Path = typer.Option(..., help="Directory where GBIF archives are cached."),
    output_dir: Path = typer.Option(..., help="Directory for exact-island pollinator records."),
    expected_family: str = typer.Option("Apidae", help="Declared acquisition family for transparent diagnostics."),
    chunksize: int = typer.Option(250_000, min=1, help="Rows parsed per GBIF SIMPLE_CSV chunk."),
) -> None:
    """Download succeeded pollinator blocks and retain exact-island record rows."""
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    members = pd.read_csv(block_members_csv, dtype=str).fillna("")
    if "analysis_island_id" not in members.columns:
        members["analysis_island_id"] = members["island_id"]
    required_members = {"block_id", "analysis_island_id"}
    missing_members = required_members.difference(members.columns)
    if missing_members:
        raise typer.BadParameter(
            f"Pollinator block-members table missing columns: {', '.join(sorted(missing_members))}"
        )
    islands = gpd.read_file(islands_gpkg, layer="islands")
    if "island_id" not in islands.columns:
        raise typer.BadParameter("Exact island GeoPackage must contain island_id")
    islands_by_id = islands.set_index("island_id")
    ready = _succeeded_blocks(campaign)

    download_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_assigned: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, Any]] = []
    failures: list[str] = []

    with httpx.Client(timeout=300.0, follow_redirects=True, headers={"User-Agent": "island-floral-v2/0.6"}) as client:
        for entry in ready:
            block_id = str(entry["block_id"])
            download_key = str(entry["download_key"])
            block_islands = members.loc[members["block_id"].eq(block_id), "analysis_island_id"].unique()
            present = [island_id for island_id in block_islands if island_id in islands_by_id.index]
            metrics: dict[str, Any] = {
                "block_id": block_id,
                "download_key": download_key,
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
                archive = _download_archive(
                    client,
                    download_key,
                    download_dir,
                    entry.get("last_download_link") or entry.get("download_link"),
                )
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
                    assigned["gbif_download_key"] = download_key
                    metrics["n_exact_island_rows_before_dedup"] += int(assigned["island_id"].notna().sum())
                    metrics["n_buffer_only_rows_before_dedup"] += int(assigned["island_id"].isna().sum())
                    all_assigned.append(assigned)
                metrics["collection_status"] = "collected"
            except Exception as exc:
                metrics["collection_status"] = "collection_failed"
                metrics["error"] = f"{type(exc).__name__}: {exc}"
                failures.append(block_id)
            manifest_rows.append(metrics)

    if all_assigned:
        assigned_before_dedup = pd.concat(all_assigned, ignore_index=True)
        assigned, duplicate_audit = deduplicate_occurrences(assigned_before_dedup)
    else:
        assigned_before_dedup = pd.DataFrame()
        assigned = pd.DataFrame(columns=["island_id", "family", "genus"])
        duplicate_audit = {
            "n_input_rows": 0,
            "n_unique_gbif_ids": 0,
            "n_duplicate_gbif_ids": 0,
            "n_duplicate_rows_removed": 0,
            "n_cross_block_duplicate_gbif_ids": 0,
            "n_conflicting_island_assignments": 0,
        }

    records = exact_island_pollinator_records(assigned, expected_family)
    coverage = summarize_pollinator_records(records, expected_family)
    manifest = pd.DataFrame(manifest_rows)
    if not manifest.empty and not assigned.empty and "block_id" in assigned.columns:
        retained_by_block = assigned.loc[assigned["island_id"].notna()].groupby("block_id").size().to_dict()
        for metrics in manifest_rows:
            metrics["n_rows_retained_after_global_dedup"] = int(retained_by_block.get(metrics["block_id"], 0))
        manifest = pd.DataFrame(manifest_rows)

    records.to_csv(output_dir / "island_pollinator_occurrences.csv", index=False)
    coverage.to_csv(output_dir / "island_pollinator_record_coverage.csv", index=False)
    manifest.to_csv(output_dir / "pollinator_collection_manifest.csv", index=False)
    status = {
        "expected_family": expected_family,
        "n_succeeded_blocks": int(len(ready)),
        "n_collected_blocks": int(manifest["collection_status"].eq("collected").sum()) if not manifest.empty else 0,
        "n_failed_blocks": int(manifest["collection_status"].eq("collection_failed").sum()) if not manifest.empty else 0,
        "n_occurrences_before_global_dedup": int(len(assigned_before_dedup)),
        "n_occurrences_after_global_dedup": int(len(assigned)),
        "n_exact_island_pollinator_records": int(len(records)),
        "n_islands_with_exact_pollinator_records": int(records["island_id"].nunique()) if not records.empty else 0,
        "n_expected_family_records": int(records["expected_family_match"].sum()) if not records.empty else 0,
        "n_bombus_records": int(
            (
                _text(records.get("family", pd.Series(dtype=object))).str.upper().eq(str(expected_family).upper())
                & _text(records.get("genus", pd.Series(dtype=object))).str.upper().eq("BOMBUS")
            ).sum()
        ) if not records.empty else 0,
        "duplicate_audit": duplicate_audit,
        "failed_block_ids": failures,
    }
    (output_dir / "pollinator_collection_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Collected {status['n_collected_blocks']} pollinator blocks: "
        f"{status['n_exact_island_pollinator_records']} exact-island records across "
        f"{status['n_islands_with_exact_pollinator_records']} islands."
    )
    if failures:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
