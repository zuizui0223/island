"""Joint source-pool and observation-selection sensitivity for PR138.

This wrapper applies the already-frozen outcome-blind branch-selection sensitivity to
source-adjusted island syndrome scores. Source assignment and source expectation are built
upstream without consulting island trait outcomes; this module never refits or changes those
source pools. Each predeclared source mode is analysed separately so source definitions are
never pooled or selected post hoc.

The result is a robustness layer only. It does not replace the observed-assemblage primary
estimand and does not identify historical ancestry, in-situ evolution, pollinator loss, or
effective pollination service.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_global_branch_selection import run_selection_sensitivity

app = typer.Typer(add_completion=False, no_args_is_help=True)

OUTPUT_NAMES = [
    "joint_selection_weights.csv.gz",
    "joint_selection_model_fits.csv",
    "joint_selection_weight_diagnostics.csv",
    "joint_selection_slopes.csv",
    "joint_selection_within_context.csv",
    "joint_selection_between_context.csv",
    "joint_selection_states.csv",
    "joint_selection_headline_summary.csv",
]


def _validate_source_adjusted_scores(scores: pd.DataFrame) -> pd.DataFrame:
    required = {
        "island_id",
        "stratum",
        "syndrome",
        "syndrome_score",
        "n_species",
        "source_mode",
        "source_expectation_eligible",
    }
    missing = required - set(scores.columns)
    if missing:
        raise typer.BadParameter(
            f"source-adjusted island scores missing columns: {sorted(missing)}"
        )
    work = scores.copy()
    work["island_id"] = work["island_id"].astype(str)
    work["stratum"] = work["stratum"].astype(str)
    work["syndrome"] = work["syndrome"].astype(str)
    work["source_mode"] = work["source_mode"].astype(str)
    eligible = work["source_expectation_eligible"]
    if eligible.dtype != bool:
        eligible = eligible.astype(str).str.lower().isin(["true", "1"])
    work = work.loc[eligible].copy()
    if work.empty:
        raise typer.BadParameter("no source-adjusted scores have eligible source expectations")
    keys = ["source_mode", "island_id", "stratum", "syndrome"]
    if work.duplicated(keys).any():
        examples = work.loc[work.duplicated(keys, keep=False), keys].head(5).to_dict("records")
        raise typer.BadParameter(f"duplicate source-adjusted island scores: {examples}")
    return work


def run_joint_sensitivity(
    source_adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    selection_config: dict[str, Any],
) -> tuple[pd.DataFrame, ...]:
    scores = _validate_source_adjusted_scores(source_adjusted_scores)
    source_modes = sorted(scores["source_mode"].unique().astype(str))
    output_parts: list[list[pd.DataFrame]] = [[] for _ in OUTPUT_NAMES]

    for source_mode in source_modes:
        mode_scores = scores.loc[scores["source_mode"].eq(source_mode)].copy()
        outputs = run_selection_sensitivity(
            mode_scores,
            covariates,
            realm_assignment,
            pattern_config,
            branching_config,
            selection_config,
        )
        for index, frame in enumerate(outputs):
            tagged = frame.copy()
            tagged.insert(0, "source_mode", source_mode)
            output_parts[index].append(tagged)

    combined: list[pd.DataFrame] = []
    for parts in output_parts:
        combined.append(pd.concat(parts, ignore_index=True) if parts else pd.DataFrame())
    return tuple(combined)


@app.command("run")
def run(
    source_adjusted_scores_csv: Annotated[Path, typer.Option(exists=True)],
    covariates_csv: Annotated[Path, typer.Option(exists=True)],
    realm_assignment_csv: Annotated[Path, typer.Option(exists=True)],
    pattern_config_path: Annotated[Path, typer.Option(exists=True)],
    branching_config_path: Annotated[Path, typer.Option(exists=True)],
    selection_config_path: Annotated[Path, typer.Option(exists=True)],
    output_dir: Annotated[Path, typer.Option()],
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(branching_config_path.read_text(encoding="utf-8"))
    selection_config = yaml.safe_load(selection_config_path.read_text(encoding="utf-8"))
    scores = pd.read_csv(source_adjusted_scores_csv)
    outputs = run_joint_sensitivity(
        scores,
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pattern_config,
        branching_config,
        selection_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame, name in zip(outputs, OUTPUT_NAMES, strict=True):
        frame.to_csv(output_dir / name, index=False)

    source_modes = sorted(scores["source_mode"].dropna().astype(str).unique())
    diagnostics = outputs[2]
    manifest = {
        "contract": "chapter1_pr138_joint_source_selection_v1",
        "source_modes": source_modes,
        "source_adjustment_recomputed": False,
        "selection_uses_branch_score_values": False,
        "selection_uses_effect_estimates": False,
        "selection_uses_p_values": False,
        "source_assignment_uses_island_trait_outcomes": False,
        "observed_assemblage_primary_replaced": False,
        "maximum_realized_weight": float(diagnostics["weight_max"].max())
        if not diagnostics.empty
        else None,
        "minimum_effective_sample_fraction": float(
            diagnostics.loc[
                diagnostics["selection_mode"].ne("unweighted"),
                "effective_sample_fraction",
            ].min()
        )
        if not diagnostics.empty
        else None,
        "claim_boundary": (
            "Joint source-pool plus measured observation-selection sensitivity only; "
            "persistence does not identify historical source ancestry or a pollinator mechanism."
        ),
    }
    (output_dir / "joint_source_selection_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
