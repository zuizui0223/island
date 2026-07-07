"""Resumable, one-block-at-a-time GBIF collection.

The global collection job can be cancelled before an all-block run reaches its
final write step. This wrapper selects a small deterministic batch of succeeded
GBIF blocks, calls the existing exact-island bounded collector for only that
batch, then merges the compact products into committed cumulative outputs.

Raw request geometry remains buffered only for acquisition. Every retained row
was assigned by the underlying collector to an original exact island polygon.
The potentially large island-by-species snapshot is stored as a gzip-compressed
CSV in the repository; pandas reads it directly without altering its schema.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.gbif_collect import _succeeded_blocks
from island_v2.gbif_collect_bounded import collect as collect_bounded

app = typer.Typer(add_completion=False, no_args_is_help=True)

MANIFEST_COLUMNS = [
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
DONE_STATUSES = {"collected", "skipped_no_matching_islands"}
CUMULATIVE_SPECIES_SNAPSHOT = "island_species_occurrences.csv.gz"
LEGACY_SPECIES_SNAPSHOT = "island_species_occurrences.csv"


@app.callback()
def main() -> None:
    """Collect succeeded GBIF blocks incrementally into exact-island products."""


def _read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str).fillna("") if path.exists() else pd.DataFrame()


def _write_csv(table: pd.DataFrame, path: Path) -> None:
    """Atomically write CSV, explicitly compressing gzip snapshot paths."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    compression = "gzip" if path.name.endswith(".gz") else None
    table.to_csv(temporary, index=False, compression=compression)
    temporary.replace(path)


def _read_cumulative_species(output_dir: Path) -> pd.DataFrame:
    """Prefer the compressed snapshot, with a one-time fallback for legacy runs."""
    compressed = output_dir / CUMULATIVE_SPECIES_SNAPSHOT
    if compressed.exists():
        return _read_csv(compressed)
    return _read_csv(output_dir / LEGACY_SPECIES_SNAPSHOT)


def _normalise_manifest(table: pd.DataFrame) -> pd.DataFrame:
    result = table.copy()
    for column in MANIFEST_COLUMNS:
        if column not in result.columns:
            result[column] = ""
    return result[MANIFEST_COLUMNS].fillna("")


def _member_counts(block_members_csv: Path) -> dict[str, int]:
    members = _read_csv(block_members_csv)
    if "block_id" not in members.columns:
        raise typer.BadParameter("block-members table is missing block_id")
    return members.groupby("block_id").size().astype(int).to_dict()


def select_batch(
    campaign: dict[str, Any], previous_manifest: pd.DataFrame, member_counts: dict[str, int], max_blocks: int
) -> list[dict[str, Any]]:
    """Pick uncollected succeeded blocks, starting with smaller island batches."""
    if max_blocks < 1:
        raise typer.BadParameter("max_blocks must be at least 1")
    previous = _normalise_manifest(previous_manifest)
    completed = set(
        previous.loc[previous["collection_status"].isin(DONE_STATUSES), "block_id"].astype(str)
    )
    candidates = [entry for entry in _succeeded_blocks(campaign) if str(entry["block_id"]) not in completed]
    return sorted(
        candidates,
        key=lambda entry: (member_counts.get(str(entry["block_id"]), 10**9), str(entry["block_id"])),
    )[:max_blocks]


def _merge_species(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    columns = ["island_id", "species", "n_records", "n_unique_gbif_ids", "basis_of_record_set", "review_status"]
    table = pd.concat([old, new], ignore_index=True) if not old.empty or not new.empty else pd.DataFrame(columns=columns)
    if table.empty:
        return pd.DataFrame(columns=columns)
    for column in columns:
        if column not in table.columns:
            table[column] = ""
    table["n_records"] = pd.to_numeric(table["n_records"], errors="coerce").fillna(0).astype(int)
    table["n_unique_gbif_ids"] = pd.to_numeric(table["n_unique_gbif_ids"], errors="coerce").fillna(0).astype(int)
    rows: list[dict[str, Any]] = []
    for (island_id, species), group in table.groupby(["island_id", "species"], sort=True):
        bases = {item for value in group["basis_of_record_set"].astype(str) for item in value.split("|") if item}
        rows.append(
            {
                "island_id": island_id,
                "species": species,
                "n_records": int(group["n_records"].sum()),
                "n_unique_gbif_ids": int(group["n_unique_gbif_ids"].sum()),
                "basis_of_record_set": "|".join(sorted(bases)),
                "review_status": "unresolved_taxonomy",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _merge_effort(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    table = pd.concat([old, new], ignore_index=True)
    if table.empty:
        return table
    duplicate_islands = table.loc[table["island_id"].duplicated(), "island_id"].tolist()
    if duplicate_islands:
        raise typer.BadParameter(
            "Incremental block products overlap on island_id; frozen block partition is invalid: "
            f"{duplicate_islands[:5]}"
        )
    return table.sort_values("island_id").reset_index(drop=True)


def _merge_taxa(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    columns = ["accepted_species", "genus", "family", "n_islands", "n_records"]
    table = pd.concat([old, new], ignore_index=True) if not old.empty or not new.empty else pd.DataFrame(columns=columns)
    if table.empty:
        return pd.DataFrame(columns=columns)
    for column in columns:
        if column not in table.columns:
            table[column] = ""
    table["n_islands"] = pd.to_numeric(table["n_islands"], errors="coerce").fillna(0).astype(int)
    table["n_records"] = pd.to_numeric(table["n_records"], errors="coerce").fillna(0).astype(int)
    rows = []
    for species, group in table.groupby("accepted_species", sort=True):
        genus = next((value for value in group["genus"].astype(str) if value.strip()), "")
        family = next((value for value in group["family"].astype(str) if value.strip()), "")
        rows.append(
            {
                "accepted_species": species,
                "genus": genus,
                "family": family,
                "n_islands": int(group["n_islands"].sum()),
                "n_records": int(group["n_records"].sum()),
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _merge_manifest(old: pd.DataFrame, new: pd.DataFrame) -> pd.DataFrame:
    old = _normalise_manifest(old)
    new = _normalise_manifest(new)
    if new.empty:
        return old.sort_values("block_id").reset_index(drop=True)
    retained = old.loc[~old["block_id"].isin(new["block_id"])].copy()
    return pd.concat([retained, new], ignore_index=True).sort_values("block_id").reset_index(drop=True)


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True),
    block_members_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    download_dir: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    max_blocks: int = typer.Option(1, min=1),
    chunksize: int = typer.Option(250_000, min=1),
) -> None:
    """Process a small resumable batch and update cumulative island products."""
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    old_manifest = _read_csv(output_dir / "collection_manifest.csv")
    batch = select_batch(campaign, old_manifest, _member_counts(block_members_csv), max_blocks)
    if not batch:
        typer.echo("No uncollected succeeded GBIF blocks remain.")
        return

    filtered = dict(campaign)
    filtered["ledger"] = batch
    with tempfile.TemporaryDirectory(prefix="island_gbif_batch_") as temp:
        temp_dir = Path(temp)
        filtered_path = temp_dir / "campaign.json"
        filtered_path.write_text(json.dumps(filtered), encoding="utf-8")
        batch_output = temp_dir / "collected"
        batch_exit: int | None = None
        try:
            collect_bounded(
                campaign_json=filtered_path,
                block_members_csv=block_members_csv,
                islands_gpkg=islands_gpkg,
                download_dir=download_dir,
                output_dir=batch_output,
                chunksize=chunksize,
            )
        except typer.Exit as exc:
            batch_exit = int(exc.exit_code or 1)

        new_manifest = _read_csv(batch_output / "collection_manifest.csv")
        merged_manifest = _merge_manifest(old_manifest, new_manifest)
        species = _merge_species(
            _read_cumulative_species(output_dir),
            _read_csv(batch_output / LEGACY_SPECIES_SNAPSHOT),
        )
        effort = _merge_effort(
            _read_csv(output_dir / "island_observation_effort.csv"),
            _read_csv(batch_output / "island_observation_effort.csv"),
        )
        taxa = _merge_taxa(
            _read_csv(output_dir / "island_taxa.csv"),
            _read_csv(batch_output / "island_taxa.csv"),
        )

    _write_csv(merged_manifest, output_dir / "collection_manifest.csv")
    _write_csv(species, output_dir / CUMULATIVE_SPECIES_SNAPSHOT)
    legacy_snapshot = output_dir / LEGACY_SPECIES_SNAPSHOT
    if legacy_snapshot.exists():
        legacy_snapshot.unlink()
    _write_csv(effort, output_dir / "island_observation_effort.csv")
    _write_csv(taxa, output_dir / "island_taxa.csv")
    succeeded = _succeeded_blocks(campaign)
    completed = set(
        merged_manifest.loc[merged_manifest["collection_status"].isin(DONE_STATUSES), "block_id"].astype(str)
    )
    status = {
        "collection_mode": "resumable_exact_island_block_batches",
        "island_species_snapshot": CUMULATIVE_SPECIES_SNAPSHOT,
        "selected_block_ids": [str(entry["block_id"]) for entry in batch],
        "max_blocks_per_run": int(max_blocks),
        "n_succeeded_blocks_in_campaign": int(len(succeeded)),
        "n_collected_succeeded_blocks": int(len(completed.intersection({str(entry["block_id"]) for entry in succeeded}))),
        "n_succeeded_blocks_remaining": int(len({str(entry["block_id"]) for entry in succeeded}.difference(completed))),
        "n_failed_blocks": int(merged_manifest["collection_status"].eq("collection_failed").sum()),
        "n_island_species_pairs": int(len(species)),
        "n_islands_with_records": int(effort["island_id"].nunique()) if not effort.empty else 0,
        "n_accepted_species": int(len(taxa)),
    }
    (output_dir / "collection_status.json").write_text(json.dumps(status, indent=2), encoding="utf-8")
    typer.echo(
        f"Processed {len(batch)} block(s); {status['n_collected_succeeded_blocks']} of "
        f"{status['n_succeeded_blocks_in_campaign']} succeeded blocks are committed."
    )
    if batch_exit:
        raise typer.Exit(code=batch_exit)


if __name__ == "__main__":
    app()
