"""Decide whether an outcome is model-ready, trait-limited, or status-limited.

The key distinction is between the number of islands that *could* contribute to
an outcome once floristic status is resolved and the number that currently have
at least one direct trait trial. Trait acquisition cannot solve a status ceiling
below the requested analysis target.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)


def classify_support_decision(
    summary: pd.DataFrame,
    *,
    pilot_min_islands: int = 30,
    confirmatory_min_islands: int = 50,
) -> pd.DataFrame:
    required = {
        "regime",
        "stratum",
        "outcome",
        "n_islands_with_stratum",
        "n_direct_eligible_islands",
    }
    missing = required - set(summary.columns)
    if missing:
        raise typer.BadParameter(f"support summary missing columns: {sorted(missing)}")

    work = summary.copy()
    work["n_islands_with_stratum"] = pd.to_numeric(
        work["n_islands_with_stratum"], errors="coerce"
    ).fillna(0).astype(int)
    work["n_direct_eligible_islands"] = pd.to_numeric(
        work["n_direct_eligible_islands"], errors="coerce"
    ).fillna(0).astype(int)

    decisions: list[str] = []
    actions: list[str] = []
    for row in work.itertuples(index=False):
        ceiling = int(row.n_islands_with_stratum)
        direct = int(row.n_direct_eligible_islands)
        if ceiling < pilot_min_islands:
            decision = "status_ceiling_below_pilot"
            action = "improve_floristic_status_or_report_unidentified"
        elif ceiling < confirmatory_min_islands:
            decision = "status_ceiling_below_confirmatory"
            action = "improve_floristic_status_not_trait_acquisition"
        elif direct >= confirmatory_min_islands:
            decision = "confirmatory_count_met"
            action = "model_before_acquiring_more_traits"
        elif direct >= pilot_min_islands:
            decision = "targeted_trait_acquisition_zone"
            action = "acquire_only_species_that_add_supported_islands"
        else:
            decision = "trait_support_below_pilot"
            action = "test_recoverable_trait_head_before_any_large_campaign"
        decisions.append(decision)
        actions.append(action)

    work["support_decision"] = decisions
    work["next_action"] = actions
    work["status_ceiling_gap_to_confirmatory"] = (
        confirmatory_min_islands - work["n_islands_with_stratum"]
    ).clip(lower=0)
    work["direct_gap_to_confirmatory"] = (
        confirmatory_min_islands - work["n_direct_eligible_islands"]
    ).clip(lower=0)
    work["trait_acquisition_allowed"] = (
        work["n_islands_with_stratum"].ge(confirmatory_min_islands)
        & work["n_direct_eligible_islands"].lt(confirmatory_min_islands)
    )
    return work


def threshold_robustness(
    island_support: pd.DataFrame,
    *,
    thresholds: tuple[int, ...] = (1, 3, 5, 10),
    regime_column: str = "analysis_regime",
) -> pd.DataFrame:
    required = {"stratum", "outcome", "island_id", "n_direct_species"}
    missing = required - set(island_support.columns)
    if missing:
        raise typer.BadParameter(f"island support missing columns: {sorted(missing)}")
    work = island_support.copy()
    if regime_column not in work.columns:
        work[regime_column] = "all"
    work[regime_column] = work[regime_column].fillna("unresolved").astype(str)
    work["n_direct_species"] = pd.to_numeric(work["n_direct_species"], errors="coerce").fillna(0)

    rows: list[dict[str, object]] = []
    for (regime, stratum, outcome), group in work.groupby(
        [regime_column, "stratum", "outcome"], sort=True
    ):
        for threshold in thresholds:
            rows.append(
                {
                    "regime": regime,
                    "stratum": stratum,
                    "outcome": outcome,
                    "min_direct_species": int(threshold),
                    "n_islands": int(group["n_direct_species"].ge(threshold).sum()),
                }
            )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    support_summary_csv: Path = typer.Option(..., exists=True),
    island_support_csv: Path | None = typer.Option(None),
    output_dir: Path = typer.Option(...),
    pilot_min_islands: int = typer.Option(30),
    confirmatory_min_islands: int = typer.Option(50),
) -> None:
    summary = pd.read_csv(support_summary_csv)
    decision = classify_support_decision(
        summary,
        pilot_min_islands=pilot_min_islands,
        confirmatory_min_islands=confirmatory_min_islands,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    decision.to_csv(output_dir / "support_decisions.csv", index=False)
    if island_support_csv is not None:
        island_support = pd.read_csv(island_support_csv)
        threshold_robustness(island_support).to_csv(
            output_dir / "direct_species_threshold_robustness.csv", index=False
        )


if __name__ == "__main__":
    app()
