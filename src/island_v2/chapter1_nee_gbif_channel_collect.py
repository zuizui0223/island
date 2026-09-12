"""Collect arbitrary frozen GBIF background-group campaigns for N1 channels.

Unlike the legacy Bombus-specific collector, this module has no family assumption.
It retains quality-auditable exact-island occurrence rows for any predeclared broad
taxon campaign (e.g. Lepidoptera, Aves, Diptera or one of the bee families). Target
functional-channel detections are classified later against the independent GloBI
catalog; this collector never infers pollination function from taxonomy alone.
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


def exact_island_records(assigned: pd.DataFrame) -> pd.DataFrame:
    """Retain only rows assigned to original exact island polygons."""
    if "island_id" not in assigned.columns:
        raise ValueError("assigned occurrence table requires island_id")
    result = assigned.loc[assigned["island_id"].notna()].copy()
    if result.empty:
        return assigned.head(0).copy()
    sort_columns = [column for column in ["island_id", "gbif_id", "block_id"] if column in result.columns]
    return result.sort_values(sort_columns, na_position="last").reset_index(drop=True)


def summarize_background(records: pd.DataFrame, campaign_taxon_name: str) -> pd.DataFrame:
    """Compact audit of raw exact-island background coverage before N1 classification."""
    columns = [
        "island_id",
        "campaign_taxon_name",
        "n_records",
        "n_species",
        "n_datasets",
        "year_min",
        "year_max",
    ]
    if records.empty:
        return pd.DataFrame(columns=columns)
    rows: list[dict[str, Any]] = []
    for island_id, group in records.groupby("island_id", sort=True):
        species = group.get("species", pd.Series("", index=group.index)).fillna("").astype(str).str.strip()
        datasets = group.get("dataset_key", pd.Series("", index=group.index)).fillna("").astype(str).str.strip()
        years = pd.to_numeric(group.get("year", pd.Series(dtype=object)), errors="coerce").dropna()
        rows.append(
            {
                "island_id": str(island_id),
                "campaign_taxon_name": str(campaign_taxon_name),
                "n_records": int(len(group)),
                "n_species": int(species.loc[species.ne("")].nunique()),
                "n_datasets": int(datasets.loc[datasets.ne("")].nunique()),
                "year_min": int(years.min()) if not years.empty else pd.NA,
                "year_max": int(years.max()) if not years.empty else pd.NA,
            }
        )
    return pd.DataFrame(rows, columns=columns)


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True),
    block_members_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    download_dir: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    campaign_taxon_name: str = typer.Option(...),
    chunksize: int = typer.Option(250_000, min=1),
) -> None:
    """Download succeeded broad-group blocks and retain exact-island occurrence rows."""
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    members = pd.read_csv(block_members_csv, dtype=str).fillna("")
    if "analysis_island_id" not in members.columns:
        if "island_id" not in members.columns:
            raise typer.BadParameter("block-members table requires island_id or analysis_island_id")
        members["analysis_island_id"] = members["island_id"]
    if "block_id" not in members.columns:
        raise typer.BadParameter("block-members table requires block_id")

    islands = gpd.read_file(islands_gpkg, layer="islands")
    if "island_id" not in islands.columns:
        raise typer.BadParameter("exact island GeoPackage requires island_id")
    islands_by_id = islands.set_index("island_id")
    ready = _succeeded_blocks(campaign)

    download_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_assigned: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, Any]] = []
    failures: list[str] = []

    with httpx.Client(
        timeout=300.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/nee-channel-collector"},
    ) as client:
        for entry in ready:
            block_id = str(entry["block_id"])
            download_key = str(entry["download_key"])
            block_islands = members.loc[
                members["block_id"].eq(block_id), "analysis_island_id"
            ].unique()
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
                    metrics["n_source_rows"] += int(
                        occurrences.attrs.get("n_source_rows", len(occurrences))
                    )
                    metrics["n_valid_coordinate_rows"] += int(len(occurrences))
                    if occurrences.empty:
                        continue
                    assigned = assign_occurrences_to_islands(occurrences, island_subset)
                    assigned["block_id"] = block_id
                    assigned["gbif_download_key"] = download_key
                    metrics["n_exact_island_rows_before_dedup"] += int(
                        assigned["island_id"].notna().sum()
                    )
                    all_assigned.append(assigned)
                metrics["collection_status"] = "collected"
            except Exception as exc:  # noqa: BLE001
                metrics["collection_status"] = "collection_failed"
                metrics["error"] = f"{type(exc).__name__}: {exc}"
                failures.append(block_id)
            manifest_rows.append(metrics)

    if all_assigned:
        assigned_before_dedup = pd.concat(all_assigned, ignore_index=True)
        assigned, duplicate_audit = deduplicate_occurrences(assigned_before_dedup)
    else:
        assigned_before_dedup = pd.DataFrame()
        assigned = pd.DataFrame(columns=["island_id"])
        duplicate_audit = {
            "n_input_rows": 0,
            "n_unique_gbif_ids": 0,
            "n_duplicate_gbif_ids": 0,
            "n_duplicate_rows_removed": 0,
            "n_cross_block_duplicate_gbif_ids": 0,
            "n_conflicting_island_assignments": 0,
        }

    records = exact_island_records(assigned)
    coverage = summarize_background(records, campaign_taxon_name)
    manifest = pd.DataFrame(manifest_rows)
    if not manifest.empty and not assigned.empty and "block_id" in assigned.columns:
        retained_by_block = (
            assigned.loc[assigned["island_id"].notna()].groupby("block_id").size().to_dict()
        )
        for metrics in manifest_rows:
            metrics["n_rows_retained_after_global_dedup"] = int(
                retained_by_block.get(metrics["block_id"], 0)
            )
        manifest = pd.DataFrame(manifest_rows)

    records.to_csv(
        output_dir / "island_channel_background_occurrences.csv.gz",
        index=False,
        compression="gzip",
    )
    coverage.to_csv(output_dir / "island_channel_background_coverage.csv", index=False)
    manifest.to_csv(output_dir / "channel_collection_manifest.csv", index=False)
    status = {
        "campaign_taxon_name": campaign_taxon_name,
        "n_succeeded_blocks": int(len(ready)),
        "n_collected_blocks": int(manifest["collection_status"].eq("collected").sum())
        if not manifest.empty
        else 0,
        "n_failed_blocks": int(manifest["collection_status"].eq("collection_failed").sum())
        if not manifest.empty
        else 0,
        "n_occurrences_before_global_dedup": int(len(assigned_before_dedup)),
        "n_occurrences_after_global_dedup": int(len(assigned)),
        "n_exact_island_records": int(len(records)),
        "n_islands_with_exact_records": int(records["island_id"].nunique())
        if not records.empty
        else 0,
        "duplicate_audit": duplicate_audit,
        "failed_block_ids": failures,
        "functional_channel_not_inferred_here": True,
    }
    (output_dir / "channel_collection_status.json").write_text(
        json.dumps(status, indent=2), encoding="utf-8"
    )
    typer.echo(json.dumps(status))
    if failures:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
