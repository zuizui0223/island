"""Diagnose acquisition headroom for below-pilot endemic trait strata.

This module does not infer any trait state and does not treat GBIF records as
literature evidence.  It asks a narrower planning question: among endemic
species on currently unsupported islands, how many islands could in principle
be reached if we only attempted species above progressively looser GBIF-record
thresholds?

The diagnostic is intended for strata classified upstream as
``recoverability_test_before_acquisition``.  It is used to decide whether a
small pilot-to-30 acquisition attempt is plausible before any broader campaign.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

DEFAULT_THRESHOLDS = (200, 100, 50, 20, 10, 5, 1)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def recoverability_headroom(
    decisions: pd.DataFrame,
    unsupported_islands: pd.DataFrame,
    targets: pd.DataFrame,
    status_ledger: pd.DataFrame,
    *,
    thresholds: Iterable[int] = DEFAULT_THRESHOLDS,
    pilot_target: int = 30,
    confirmatory_target: int = 50,
) -> pd.DataFrame:
    required_decisions = {
        "regime",
        "outcome",
        "n_endemic_status_islands",
        "n_direct_supported_islands",
        "decision",
    }
    missing = required_decisions - set(decisions.columns)
    if missing:
        raise typer.BadParameter(f"decision table missing columns: {sorted(missing)}")
    required_unsupported = {"regime", "outcome", "island_id"}
    missing = required_unsupported - set(unsupported_islands.columns)
    if missing:
        raise typer.BadParameter(f"unsupported-island table missing columns: {sorted(missing)}")
    required_targets = {"regime", "outcome", "accepted_species", "n_records"}
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

    target_table = targets.copy()
    target_table["accepted_species"] = _text(target_table["accepted_species"])
    target_table["n_records"] = pd.to_numeric(target_table["n_records"], errors="coerce").fillna(0)
    threshold_values = sorted({int(value) for value in thresholds if int(value) >= 0}, reverse=True)
    if not threshold_values:
        raise typer.BadParameter("at least one non-negative record threshold is required")

    rows: list[dict[str, object]] = []
    strata = decisions.loc[
        decisions["decision"].eq("recoverability_test_before_acquisition")
    ].copy()
    for decision in strata.itertuples(index=False):
        regime = str(decision.regime)
        outcome = str(decision.outcome)
        current = int(decision.n_direct_supported_islands)
        ceiling = int(decision.n_endemic_status_islands)
        allowed_islands = set(
            unsupported_islands.loc[
                unsupported_islands["regime"].eq(regime)
                & unsupported_islands["outcome"].eq(outcome),
                "island_id",
            ].astype(str)
        )
        ranked = target_table.loc[
            target_table["regime"].eq(regime)
            & target_table["outcome"].eq(outcome)
        ].copy()
        for threshold in threshold_values:
            eligible_species = set(
                ranked.loc[ranked["n_records"].ge(threshold), "accepted_species"].astype(str)
            )
            reached = set(
                status.loc[
                    status["accepted_species"].isin(eligible_species)
                    & status["island_id"].isin(allowed_islands),
                    "island_id",
                ]
            )
            max_support = min(ceiling, current + len(reached))
            rows.append(
                {
                    "regime": regime,
                    "outcome": outcome,
                    "current_direct_supported_islands": current,
                    "endemic_status_ceiling_islands": ceiling,
                    "gbif_record_threshold": threshold,
                    "eligible_candidate_species": len(eligible_species),
                    "unsupported_islands_reached": len(reached),
                    "max_support_if_all_candidates_recovered": max_support,
                    "pilot_target": pilot_target,
                    "pilot_shortfall_after_full_recovery": max(0, pilot_target - max_support),
                    "pilot_reachable": max_support >= pilot_target,
                    "confirmatory_target": confirmatory_target,
                    "confirmatory_shortfall_after_full_recovery": max(
                        0, confirmatory_target - max_support
                    ),
                    "confirmatory_reachable": max_support >= confirmatory_target,
                }
            )
    return pd.DataFrame(rows)


def strictest_pilot_threshold(table: pd.DataFrame) -> pd.DataFrame:
    """Return the highest record threshold that can reach the pilot target per stratum."""

    if table.empty:
        return pd.DataFrame(
            columns=[
                "regime",
                "outcome",
                "current_direct_supported_islands",
                "strictest_pilot_reachable_gbif_threshold",
                "eligible_candidate_species_at_threshold",
                "unsupported_islands_reached_at_threshold",
                "max_support_at_threshold",
                "pilot_reachable_under_tested_thresholds",
            ]
        )
    rows: list[dict[str, object]] = []
    for (regime, outcome), group in table.groupby(["regime", "outcome"], sort=True):
        passing = group.loc[group["pilot_reachable"]].sort_values(
            "gbif_record_threshold", ascending=False
        )
        first = passing.iloc[0] if not passing.empty else None
        current = int(group["current_direct_supported_islands"].iloc[0])
        rows.append(
            {
                "regime": regime,
                "outcome": outcome,
                "current_direct_supported_islands": current,
                "strictest_pilot_reachable_gbif_threshold": (
                    int(first["gbif_record_threshold"]) if first is not None else pd.NA
                ),
                "eligible_candidate_species_at_threshold": (
                    int(first["eligible_candidate_species"]) if first is not None else 0
                ),
                "unsupported_islands_reached_at_threshold": (
                    int(first["unsupported_islands_reached"]) if first is not None else 0
                ),
                "max_support_at_threshold": (
                    int(first["max_support_if_all_candidates_recovered"])
                    if first is not None
                    else current
                ),
                "pilot_reachable_under_tested_thresholds": first is not None,
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    decisions_csv: Path = typer.Option(..., exists=True),
    unsupported_islands_csv: Path = typer.Option(..., exists=True),
    targets_csv: Path = typer.Option(..., exists=True),
    status_ledger_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    thresholds: str = typer.Option("200,100,50,20,10,5,1"),
    pilot_target: int = typer.Option(30),
    confirmatory_target: int = typer.Option(50),
) -> None:
    parsed = tuple(int(value.strip()) for value in thresholds.split(",") if value.strip())
    table = recoverability_headroom(
        pd.read_csv(decisions_csv),
        pd.read_csv(unsupported_islands_csv),
        pd.read_csv(targets_csv),
        pd.read_csv(status_ledger_csv),
        thresholds=parsed,
        pilot_target=pilot_target,
        confirmatory_target=confirmatory_target,
    )
    summary = strictest_pilot_threshold(table)
    output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(output_dir / "endemic_trait_recoverability_headroom.csv", index=False)
    summary.to_csv(output_dir / "endemic_trait_recoverability_summary.csv", index=False)
    manifest = {
        "contract": "endemic_trait_recoverability_headroom_v1",
        "interpretation": (
            "GBIF record count is an observation-effort proxy used only to bound acquisition "
            "headroom. It is not evidence that a direct trait source exists."
        ),
        "pilot_target": pilot_target,
        "confirmatory_target": confirmatory_target,
        "thresholds": list(parsed),
        "n_diagnostic_rows": int(len(table)),
        "summary": summary.to_dict("records"),
    }
    (output_dir / "endemic_trait_recoverability_manifest.json").write_text(
        json.dumps(manifest, indent=2, default=str) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2, default=str))


if __name__ == "__main__":
    app()
