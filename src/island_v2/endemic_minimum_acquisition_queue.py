"""Build the smallest ranked acquisition queue needed to reach the island target.

The upstream endemic targeter ranks species by unsupported-island gain,
isolation-edge gain, and GBIF record count. This module adds a stop rule: for
strata already in the 30-49-island targeted-acquisition zone, greedily select
ranked species only until the declared confirmatory island gap is covered.

It never assigns a trait value. A queued species is only a request for new direct
evidence. Below-pilot outcomes are not converted into a 50-island shopping list;
they remain recoverability diagnostics.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def build_minimum_queue(
    decisions: pd.DataFrame,
    unsupported_islands: pd.DataFrame,
    targets: pd.DataFrame,
    status_ledger: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_dec = {"regime", "outcome", "gap_to_50", "decision"}
    missing = required_dec - set(decisions.columns)
    if missing:
        raise typer.BadParameter(f"decision table missing columns: {sorted(missing)}")
    required_islands = {"regime", "outcome", "island_id"}
    missing = required_islands - set(unsupported_islands.columns)
    if missing:
        raise typer.BadParameter(f"unsupported-island table missing columns: {sorted(missing)}")
    required_targets = {"regime", "outcome", "accepted_species", "priority_rank"}
    missing = required_targets - set(targets.columns)
    if missing:
        raise typer.BadParameter(f"target table missing columns: {sorted(missing)}")
    required_status = {"island_id", "accepted_species", "origin_status", "endemic_status"}
    missing = required_status - set(status_ledger.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")

    status = status_ledger.copy()
    status["island_id"] = _text(status["island_id"])
    status["accepted_species"] = _text(status["accepted_species"])
    status = status.loc[
        _text(status["origin_status"]).str.lower().eq("native")
        & _text(status["endemic_status"]).str.lower().eq("endemic")
    ].drop_duplicates(["island_id", "accepted_species"])

    queue_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    targeted = decisions.loc[decisions["decision"].eq("targeted_trait_acquisition")].copy()
    for decision in targeted.itertuples(index=False):
        regime = str(decision.regime)
        outcome = str(decision.outcome)
        gap = int(decision.gap_to_50)
        allowed_islands = set(
            unsupported_islands.loc[
                unsupported_islands["regime"].eq(regime)
                & unsupported_islands["outcome"].eq(outcome),
                "island_id",
            ].astype(str)
        )
        ranked = targets.loc[
            targets["regime"].eq(regime) & targets["outcome"].eq(outcome)
        ].sort_values(["priority_rank", "accepted_species"])
        covered: set[str] = set()
        selected = 0
        for target in ranked.itertuples(index=False):
            species = str(target.accepted_species)
            species_islands = set(
                status.loc[
                    status["accepted_species"].eq(species)
                    & status["island_id"].isin(allowed_islands),
                    "island_id",
                ]
            )
            adds = sorted(species_islands - covered)
            if not adds:
                continue
            selected += 1
            covered.update(adds)
            row = target._asdict()
            row.update(
                {
                    "selection_order": selected,
                    "new_island_ids": "|".join(adds),
                    "new_islands_added": len(adds),
                    "cumulative_new_islands": len(covered),
                    "required_gap": gap,
                }
            )
            queue_rows.append(row)
            if len(covered) >= gap:
                break
        summary_rows.append(
            {
                "regime": regime,
                "outcome": outcome,
                "gap_to_confirmatory_target": gap,
                "n_selected_species": selected,
                "n_new_islands_covered": len(covered),
                "queue_reaches_target": len(covered) >= gap,
                "n_available_unsupported_islands": len(allowed_islands),
                "policy": "stop as soon as ranked direct-evidence targets cover the island gap",
            }
        )

    queue = pd.DataFrame(queue_rows)
    summary = pd.DataFrame(summary_rows)
    return queue, summary


@app.command("run")
def run(
    decisions_csv: Path = typer.Option(..., exists=True),
    unsupported_islands_csv: Path = typer.Option(..., exists=True),
    targets_csv: Path = typer.Option(..., exists=True),
    status_ledger_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    queue, summary = build_minimum_queue(
        pd.read_csv(decisions_csv),
        pd.read_csv(unsupported_islands_csv),
        pd.read_csv(targets_csv),
        pd.read_csv(status_ledger_csv),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    queue.to_csv(output_dir / "minimum_endemic_trait_acquisition_queue.csv", index=False)
    summary.to_csv(output_dir / "minimum_endemic_trait_acquisition_summary.csv", index=False)
    manifest = {
        "contract": "minimum_endemic_acquisition_queue_v1",
        "policy": (
            "Only targeted-acquisition-zone outcomes are queued. Species follow the upstream "
            "island-gain/isolation-edge/recoverability ranking and stop once the island gap is covered."
        ),
        "summary": summary.to_dict("records"),
        "n_queue_rows": int(len(queue)),
    }
    (output_dir / "minimum_endemic_trait_acquisition_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
