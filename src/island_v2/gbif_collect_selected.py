"""Collect one explicit already-succeeded GBIF block into cumulative outputs.

This is an operational ordering override only. It filters the existing frozen
campaign to one named succeeded, uncollected block, then delegates all exact
point-in-polygon assignment and cumulative merging to the standard incremental
collector. It never reads traits, pollinator observations, or model outputs.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2 import gbif_collect_incremental as incremental
from island_v2.gbif_collect import _succeeded_blocks

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Collect one named succeeded GBIF block without changing collection semantics."""


def _completed_blocks(manifest: pd.DataFrame) -> set[str]:
    normalised = incremental._normalise_manifest(manifest)
    return set(
        normalised.loc[
            normalised["collection_status"].isin(incremental.DONE_STATUSES), "block_id"
        ].astype(str)
    )


def select_succeeded_block(
    campaign: dict[str, Any], previous_manifest: pd.DataFrame, block_id: str
) -> dict[str, Any]:
    """Return exactly one uncollected succeeded campaign ledger entry."""
    requested = str(block_id or "").strip()
    if not requested:
        raise typer.BadParameter("block_id is required")
    if requested in _completed_blocks(previous_manifest):
        raise typer.BadParameter(f"Requested block_id is already collected: {requested}")
    matches = [entry for entry in _succeeded_blocks(campaign) if str(entry["block_id"]) == requested]
    if not matches:
        raise typer.BadParameter(
            "Requested block_id is not a succeeded uncollected campaign block: " f"{requested}"
        )
    return matches[0]


def _status_after_selected_collection(
    campaign: dict[str, Any], output_dir: Path, block_id: str, selection_rationale: str
) -> None:
    """Restore global campaign counts after the delegated one-block run."""
    status_path = output_dir / "collection_status.json"
    manifest_path = output_dir / "collection_manifest.csv"
    status = json.loads(status_path.read_text(encoding="utf-8")) if status_path.exists() else {}
    manifest = incremental._read_csv(manifest_path)
    succeeded = _succeeded_blocks(campaign)
    succeeded_ids = {str(entry["block_id"]) for entry in succeeded}
    completed = _completed_blocks(manifest)
    status.update(
        {
            "collection_mode": "resumable_exact_island_block_batches",
            "selection_mode": "explicit_succeeded_uncollected_block",
            "requested_block_id": block_id,
            "selection_rationale": selection_rationale,
            "selected_block_ids": [block_id],
            "max_blocks_per_run": 1,
            "n_succeeded_blocks_in_campaign": len(succeeded_ids),
            "n_collected_succeeded_blocks": len(completed.intersection(succeeded_ids)),
            "n_succeeded_blocks_remaining": len(succeeded_ids.difference(completed)),
            "n_failed_blocks": int(manifest["collection_status"].eq("collection_failed").sum())
            if not manifest.empty
            else 0,
        }
    )
    status_path.write_text(json.dumps(status, indent=2), encoding="utf-8")


@app.command("collect")
def collect(
    campaign_json: Path = typer.Option(..., exists=True),
    block_members_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    download_dir: Path = typer.Option(...),
    output_dir: Path = typer.Option(...),
    block_id: str = typer.Option(..., help="Exact already-succeeded, uncollected campaign block ID."),
    selection_rationale: str = typer.Option(
        ..., help="Predeclared geographic/source-region coverage rationale; never trait/outcome based."
    ),
    chunksize: int = typer.Option(250_000, min=1),
) -> None:
    """Collect one explicit succeeded block using the normal exact-island path."""
    rationale = str(selection_rationale or "").strip()
    if not rationale:
        raise typer.BadParameter("selection_rationale is required")
    campaign = json.loads(campaign_json.read_text(encoding="utf-8"))
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = select_succeeded_block(
        campaign, incremental._read_csv(output_dir / "collection_manifest.csv"), block_id
    )
    filtered = dict(campaign)
    filtered["ledger"] = [selected]
    with tempfile.TemporaryDirectory(prefix="island_gbif_selected_") as temporary:
        filtered_path = Path(temporary) / "single_block_campaign.json"
        filtered_path.write_text(json.dumps(filtered), encoding="utf-8")
        incremental.collect(
            campaign_json=filtered_path,
            block_members_csv=block_members_csv,
            islands_gpkg=islands_gpkg,
            download_dir=download_dir,
            output_dir=output_dir,
            max_blocks=1,
            chunksize=chunksize,
        )
    _status_after_selected_collection(campaign, output_dir, str(block_id), rationale)
    typer.echo(f"Collected explicit succeeded block {block_id} through the exact-island pipeline.")


if __name__ == "__main__":
    app()
