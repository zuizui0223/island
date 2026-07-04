"""Memory-bounded collection of succeeded GBIF blocks into island flora products.

GBIF requests use buffered catchments to protect against shoreline mismatch.
Buffer-only rows are useful for audit counts but cannot contribute to an island
flora; retaining all of them globally can exhaust memory before outputs exist.
This collector keeps exact-island rows after each join and retains buffer-only
counts in the block audit.
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
    _SOURCE_COLUMNS,
    _download_archive,
    _succeeded_blocks,
    assign_occurrences_to_islands,
    deduplicate_occurrences,
    island_species_table,
    island_taxa_table,
    iter_block_occurrences,
    sha256_file,
    summarize_observation_effort,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Collect succeeded GBIF blocks with exact-island-only in-memory retention."""


def _empty_assigned_frame() -> pd.DataFrame:
    return pd.DataFrame(
        columns=list(_SOURCE_COLUMNS.values()) + ["island_id", "block_id", "gbif_download_key"]
    )


def _empty_manifest() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "block_id",
            "download_key",
            "doi",
            "n_block_member_islands",
            "n_islands_available_for_assignment",
            "archive_path",
            "archive_sha256",
            "archive_bytes",
            "n_source_rows",
            "n_valid_coordinate_rows",
            "n_exact_island_rows_before_dedup",
            "n_buffer_only_rows",
            "n_rows_retained_after_global_dedup",
            "collection_status",
            "error",
        ]
    )


def collect_succeeded_blocks(
    campaign: dict[str, Any],
    members: pd.DataFrame,
    islands: gpd.GeoDataFrame,
    download_dir: Path,
    chunksize: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int], list[str]]:
    """Download, exact-assign, and globally deduplicate succeeded blocks."""
    if "analysis_island_id" not in members.columns:
        members = members.copy()
        members["analysis_island_id"] = members["island_id"]
    required_members = {"block_id", "analysis_island_id"}
    missing_members = required_members.difference(members.columns)
    if missing_members:
        raise typer.BadParameter(
            f"Block-members table missing columns: {', '.join(sorted(missing_members))}"
        )
    if "island_id" not in islands.columns:
        raise typer.BadParameter("Exact island GeoPackage must contain island_id")

    islands_by_id = islands.set_index("island_id")
    ready = _succeeded_blocks(campaign)
    exact_chunks: list[pd.DataFrame] = []
    manifest_rows: list[dict[str, Any]] = []
    failures: list[str] = []

    download_dir.mkdir(parents=True, exist_ok=True)
    with httpx.Client(
        timeout=300.0,
        follow_redirects=True,
        headers={"User-Agent": "island-floral-v2/0.7"},
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
                "n_buffer_only_rows": 0,
                "n_rows_retained_after_global_dedup": 0,
                "collection_status": "pending",
                "error": "",
            }
            if not present:
                metrics["collection_status"] = "skipped_no_matching_islands"
                manifest_rows.append(metrics)
                continue
            try:
                island_subset = gpd.GeoDataFrame(
                    islands_by_id.loc[present].reset_index(), geometry="geometry", crs=islands.crs
                )
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
                    exact = assigned.loc[assigned["island_id"].notna()].copy()
                    metrics["n_exact_island_rows_before_dedup"] += int(len(exact))
                    metrics["n_buffer_only_rows"] += int(len(assigned) - len(exact))
                    if not exact.empty:
                        exact_chunks.append(exact)
                metrics["collection_status"] = "collected"
            except Exception as exc:
                metrics["collection_status"] = "collection_failed"
                metrics["error"] = f"{type(exc).__name__}: {exc}"
                failures.append(block_id)
            manifest_rows.append(metrics)

    exact_before_dedup = pd.concat(exact_chunks, ignore_index=True) if exact_chunks else _empty_assigned_frame()
    assigned, duplicate_audit = deduplicate_occurrences(exact_before_dedup)
    manifest = pd.DataFrame(manifest_rows) if manifest_rows else _empty_manifest()
    if not manifest.empty and not assigned.empty:
        retained_by_block = assigned.groupby("block_id").size().to_dict()
        manifest["n_rows_retained_after_global_dedup"] = (
            manifest["block_id"].map(retained_by_block).fillna(0).astype(int)
        )
    return assigned, manifest, duplicate_audit, failures


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True, help="Campaign ledger JSON."),
    block_members_csv: Path = typer.Option(..., exists=True, help="gbif_block_members.csv."),
    islands_gpkg: Path = typer.Option(..., exists=True, help="Exact island polygons GeoPackage."),
    download_dir: Path = typer.Option(..., help="Where block archives are cached."),
    output_dir: Path = typer.Option(..., help="Where island-level outputs are written."),
    chunksize: int = typer.Option(250_000, min=1, help="Rows read per GBIF SIMPLE_CSV chunk."),
) -> None:
    """Collect succeeded blocks into exact-island flora products and audits."""
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    members = pd.read_csv(block_members_csv, dtype=str).fillna("")
    islands = gpd.read_file(islands_gpkg, layer="islands")
    output_dir.mkdir(parents=True, exist_ok=True)

    assigned, manifest, duplicate_audit, failures = collect_succeeded_blocks(
        campaign, members, islands, download_dir, chunksize
    )
    species_table = island_species_table(assigned)
    effort = summarize_observation_effort(assigned)
    taxa = island_taxa_table(assigned)

    species_table.to_csv(output_dir / "island_species_occurrences.csv", index=False)
    effort.to_csv(output_dir / "island_observation_effort.csv", index=False)
    taxa.to_csv(output_dir / "island_taxa.csv", index=False)
    manifest.to_csv(output_dir / "collection_manifest.csv", index=False)
    status = {
        "collection_mode": "exact_island_rows_only_in_memory",
        "n_succeeded_blocks": int(len(_succeeded_blocks(campaign))),
        "n_collected_blocks": int(manifest["collection_status"].eq("collected").sum()) if not manifest.empty else 0,
        "n_failed_blocks": int(manifest["collection_status"].eq("collection_failed").sum()) if not manifest.empty else 0,
        "n_exact_island_rows_before_global_dedup": int(duplicate_audit["n_input_rows"]),
        "n_occurrences_after_global_dedup": int(len(assigned)),
        "n_occurrences_on_island": int(len(assigned)),
        "n_occurrences_in_buffer_only": int(manifest["n_buffer_only_rows"].sum()) if not manifest.empty else 0,
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
        f"({status['n_occurrences_in_buffer_only']} buffer-only rows discarded before global dedup; "
        f"{duplicate_audit['n_duplicate_rows_removed']} duplicate exact-island rows removed)."
    )
    if failures:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
