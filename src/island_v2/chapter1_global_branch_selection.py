"""Outcome-blind observation-selection sensitivity for global plant-side branches.

Each primary branch axis gets its own score-availability model over the complete island
universe. Stabilized inverse-probability weights are then applied to the already-declared
continuous branch responses. Score values, fitted ecological effects and p-values never
enter the selection model.

This is a missing-at-random sensitivity over measured geography, area, climate and broad
context. It is not a correction for an unknown source pool and not evidence about realized
pollination service.
"""

from __future__ import annotations

import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design
from island_v2.chapter1_global_branching import (
    _run_context_layer,
    build_branch_scores,
    classify_branch_states,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _expit(values: np.ndarray) -> np.ndarray:
    out = np.empty_like(values, dtype=float)
    positive = values >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    negative = ~positive
    exp_values = np.exp(values[negative])
    out[negative] = exp_values / (1.0 + exp_values)
    return out


def _context_universe(
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    layer_name: str,
    layer_spec: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    context = str(layer_spec["column"])
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    required = {"island_id", geography, cluster, *baseline}
    missing = required - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    universe = covariates.copy()
    if context not in universe.columns:
        assignment_required = {"island_id", context}
        missing = assignment_required - set(realm_assignment.columns)
        if missing:
            raise typer.BadParameter(
                f"context assignment missing columns: {sorted(missing)}"
            )
        assignment = realm_assignment[["island_id", context]].drop_duplicates("island_id")
        universe = universe.merge(
            assignment,
            on="island_id",
            how="left",
            validate="one_to_one",
        )
    if universe["island_id"].duplicated().any():
        raise typer.BadParameter("selection universe contains duplicate island_id values")
    numeric = [geography, *baseline]
    for column in numeric:
        universe[column] = pd.to_numeric(universe[column], errors="coerce")
    universe[context] = universe[context].fillna("").astype(str)
    universe[cluster] = universe[cluster].fillna("").astype(str)
    declared_contexts = [str(x) for x in layer_spec["contexts"]]
    full_n = int(universe["island_id"].nunique())
    complete = universe.dropna(subset=numeric).copy()
    complete = complete.loc[
        complete[context].isin(declared_contexts) & complete[cluster].ne("")
    ].copy()
    return complete, {
        "context_layer": layer_name,
        "n_full_islands": full_n,
        "n_context_assigned_complete": int(complete["island_id"].nunique()),
        "n_context_unassigned_or_incomplete": full_n
        - int(complete["island_id"].nunique()),
    }


def _selection_design(
    universe: pd.DataFrame,
    pattern_config: dict[str, Any],
    layer_spec: dict[str, Any],
) -> tuple[np.ndarray, list[str]]:
    context = str(layer_spec["column"])
    geography = str(pattern_config["geography_column"])
    continuous = [
        *[str(x) for x in pattern_config["baseline_covariates"]],
        geography,
    ]
    names = ["intercept"]
    columns = [np.ones(len(universe), dtype=float)]
    standardized: dict[str, np.ndarray] = {}
    for predictor in continuous:
        values = pd.to_numeric(universe[predictor], errors="coerce").to_numpy(float)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=0))
        if not math.isfinite(sd) or sd <= 0:
            raise typer.BadParameter(f"constant selection predictor: {predictor}")
        standardized[predictor] = (values - mean) / sd
        names.append(f"z_{predictor}")
        columns.append(standardized[predictor])

    present = set(universe[context].astype(str))
    contexts = [str(x) for x in layer_spec["contexts"] if str(x) in present]
    if not contexts:
        raise typer.BadParameter("selection universe has no declared contexts")
    reference = contexts[0]
    for value in contexts:
        if value == reference:
            continue
        indicator = universe[context].eq(value).to_numpy(float)
        names.append(f"context[{value}]")
        columns.append(indicator)
        names.append(f"z_{geography}:context[{value}]")
        columns.append(standardized[geography] * indicator)
    return np.column_stack(columns), names


def build_selection_weights(
    branch_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    selection_config: dict[str, Any],
    *,
    layer_name: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build bounded, stabilized weights without consulting any branch-score value."""

    layer_spec = branching_config["context_layers"][layer_name]
    context = str(layer_spec["column"])
    cluster = str(pattern_config["cluster_column"])
    primary_set = str(selection_config["primary_axis_set"])
    axes = [str(x) for x in branching_config["axis_sets"][primary_set]["axes"]]
    strata = [str(x) for x in pattern_config["strata"]]
    universe, universe_summary = _context_universe(
        covariates,
        realm_assignment,
        pattern_config,
        layer_name,
        layer_spec,
    )
    design, names = _selection_design(universe, pattern_config, layer_spec)
    minimum_observed = int(selection_config["selection_model"]["minimum_observed"])
    minimum_unobserved = int(selection_config["selection_model"]["minimum_unobserved"])

    weight_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    diagnostic_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for axis in axes:
            observed_ids = set(
                branch_scores.loc[
                    branch_scores["stratum"].eq(stratum)
                    & branch_scores["syndrome"].eq(axis),
                    "island_id",
                ].astype(str)
            )
            outcome = universe["island_id"].astype(str).isin(observed_ids).to_numpy(float)
            n_observed = int(outcome.sum())
            n_unobserved = int(len(outcome) - n_observed)
            fit_base = {
                **universe_summary,
                "axis_set": primary_set,
                "stratum": stratum,
                "axis": axis,
                "n_observed": n_observed,
                "n_unobserved": n_unobserved,
            }
            unweighted = universe.loc[
                outcome.astype(bool), ["island_id", context]
            ].copy()
            unweighted["propensity"] = np.nan
            unweighted["context_layer"] = layer_name
            unweighted["axis_set"] = primary_set
            unweighted["stratum"] = stratum
            unweighted["syndrome"] = axis
            unweighted["selection_mode"] = "unweighted"
            unweighted["analysis_weight"] = 1.0
            weight_parts.append(unweighted)
            diagnostic_rows.append(
                {
                    "context_layer": layer_name,
                    "axis_set": primary_set,
                    "stratum": stratum,
                    "axis": axis,
                    "selection_mode": "unweighted",
                    "n_observed": n_observed,
                    "n_propensity_clipped": 0,
                    "weight_min": 1.0,
                    "weight_median": 1.0,
                    "weight_max": 1.0,
                    "effective_sample_size": float(n_observed),
                    "effective_sample_fraction": 1.0,
                }
            )
            if n_observed < minimum_observed or n_unobserved < minimum_unobserved:
                fit_rows.append({**fit_base, "status": "below_selection_support"})
                continue
            coefficients, fit, _ = _fit_grouped_binomial_design(
                outcome,
                np.ones(len(outcome), dtype=float),
                design,
                names,
                universe[cluster].astype(str).to_numpy(),
            )
            if coefficients.empty or fit.get("status") != "fit":
                fit_rows.append({**fit_base, **fit})
                continue
            beta = coefficients.set_index("predictor").loc[names, "estimate_log_odds"].to_numpy(
                float
            )
            propensity = np.clip(_expit(design @ beta), 1e-8, 1 - 1e-8)
            fit_rows.append(
                {
                    **fit_base,
                    **fit,
                    "propensity_min": float(np.min(propensity)),
                    "propensity_median": float(np.median(propensity)),
                    "propensity_max": float(np.max(propensity)),
                }
            )
            predicted = universe[["island_id", context]].copy()
            predicted["observed"] = outcome.astype(bool)
            predicted["propensity"] = propensity
            prevalence = predicted.groupby(context)["observed"].mean()
            observed = predicted.loc[predicted["observed"]].copy()
            observed["stabilizing_numerator"] = observed[context].map(prevalence)
            for mode, mode_spec in selection_config["weight_modes"].items():
                if mode == "unweighted":
                    continue
                lower = float(mode_spec["minimum_propensity"])
                upper = float(mode_spec["maximum_propensity"])
                raw_propensity = observed["propensity"].to_numpy(float)
                bounded = np.clip(raw_propensity, lower, upper)
                clipped = bounded != raw_propensity
                weights = observed["stabilizing_numerator"].to_numpy(float) / bounded
                if bool(
                    selection_config["stabilization"][
                        "normalize_within_axis_stratum_context"
                    ]
                ):
                    normalizer = (
                        pd.DataFrame({context: observed[context].to_numpy(), "weight": weights})
                        .groupby(context)["weight"]
                        .transform("mean")
                        .to_numpy(float)
                    )
                    weights = weights / normalizer
                weights = np.minimum(weights, float(mode_spec["maximum_weight"]))
                frame = observed[["island_id", context, "propensity"]].copy()
                frame["context_layer"] = layer_name
                frame["axis_set"] = primary_set
                frame["stratum"] = stratum
                frame["syndrome"] = axis
                frame["selection_mode"] = str(mode)
                frame["analysis_weight"] = weights
                weight_parts.append(frame)
                sum_weight = float(np.sum(weights))
                ess = sum_weight**2 / float(np.sum(np.square(weights)))
                diagnostic_rows.append(
                    {
                        "context_layer": layer_name,
                        "axis_set": primary_set,
                        "stratum": stratum,
                        "axis": axis,
                        "selection_mode": str(mode),
                        "n_observed": len(observed),
                        "n_propensity_clipped": int(clipped.sum()),
                        "weight_min": float(np.min(weights)),
                        "weight_median": float(np.median(weights)),
                        "weight_max": float(np.max(weights)),
                        "effective_sample_size": ess,
                        "effective_sample_fraction": ess / len(observed),
                    }
                )
    weights = pd.concat(weight_parts, ignore_index=True) if weight_parts else pd.DataFrame()
    return weights, pd.DataFrame(fit_rows), pd.DataFrame(diagnostic_rows)


def _headline_summary(
    states: pd.DataFrame,
    between: pd.DataFrame,
    selection_config: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    headline_strata = [
        value
        for value in ("all_native", "native_nonendemic")
        if value in set(states["stratum"].astype(str))
    ]
    for mode in states["selection_mode"].drop_duplicates().astype(str):
        mode_states = states.loc[
            states["selection_mode"].eq(mode)
            & states["stratum"].isin(headline_strata)
            & states["support_tier"].eq("confirmatory")
        ]
        mode_between = between.loc[
            between["selection_mode"].eq(mode)
            & between["stratum"].isin(headline_strata)
            & between["support_tier"].eq("confirmatory")
        ]
        for layer, contrasts in selection_config["headline_contrasts"].items():
            for context_a, context_b in contrasts:
                for stratum in headline_strata:
                    subset = mode_states.loc[
                        mode_states["context_layer"].eq(layer)
                        & mode_states["stratum"].eq(stratum)
                        & mode_states["context"].isin([context_a, context_b])
                    ].set_index("context")
                    direct = mode_between.loc[
                        mode_between["context_layer"].eq(layer)
                        & mode_between["stratum"].eq(stratum)
                        & mode_between["context_a"].eq(context_a)
                        & mode_between["context_b"].eq(context_b)
                    ]
                    rows.append(
                        {
                            "selection_mode": mode,
                            "context_layer": layer,
                            "stratum": stratum,
                            "context_a": context_a,
                            "context_b": context_b,
                            "context_a_state": subset.loc[context_a, "plant_side_branch_state"]
                            if context_a in subset.index
                            else "not_testable",
                            "context_b_state": subset.loc[context_b, "plant_side_branch_state"]
                            if context_b in subset.index
                            else "not_testable",
                            "direct_vector_difference_q": float(
                                direct.iloc[0].get("q_between_context_family", np.nan)
                            )
                            if not direct.empty
                            else np.nan,
                            "direct_vector_difference_supported": bool(
                                not direct.empty
                                and direct.iloc[0].get(
                                    "context_vector_difference_supported", False
                                )
                            ),
                        }
                    )
    return pd.DataFrame(rows)


def run_selection_sensitivity(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    selection_config: dict[str, Any],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    branch_scores = build_branch_scores(island_scores, branching_config)
    primary_set = str(selection_config["primary_axis_set"])
    primary_axes = [str(x) for x in branching_config["axis_sets"][primary_set]["axes"]]
    primary_scores = branch_scores.loc[branch_scores["syndrome"].isin(primary_axes)].copy()
    primary_config = deepcopy(branching_config)
    primary_config["axis_sets"] = {primary_set: primary_config["axis_sets"][primary_set]}

    all_weights: list[pd.DataFrame] = []
    all_fits: list[pd.DataFrame] = []
    all_diagnostics: list[pd.DataFrame] = []
    all_slopes: list[pd.DataFrame] = []
    all_within: list[pd.DataFrame] = []
    all_between: list[pd.DataFrame] = []
    all_states: list[pd.DataFrame] = []
    for layer_name, layer_spec in branching_config["context_layers"].items():
        weights, fits, diagnostics = build_selection_weights(
            primary_scores,
            covariates,
            realm_assignment,
            pattern_config,
            branching_config,
            selection_config,
            layer_name=str(layer_name),
        )
        all_weights.append(weights)
        all_fits.append(fits)
        all_diagnostics.append(diagnostics)
        layer_covariates, _ = _context_universe(
            covariates,
            realm_assignment,
            pattern_config,
            str(layer_name),
            layer_spec,
        )
        layer_pattern = deepcopy(pattern_config)
        layer_pattern["context_column"] = str(layer_spec["column"])
        layer_pattern["contexts"] = [str(x) for x in layer_spec["contexts"]]
        layer_config = deepcopy(primary_config)
        layer_config["context_layers"] = {str(layer_name): layer_spec}
        for mode in selection_config["weight_modes"]:
            mode_weights = weights.loc[weights["selection_mode"].eq(str(mode))].copy()
            weighted_scores = primary_scores.merge(
                mode_weights[
                    ["island_id", "stratum", "syndrome", "analysis_weight"]
                ],
                on=["island_id", "stratum", "syndrome"],
                how="inner",
                validate="one_to_one",
            )
            slopes, within, between = _run_context_layer(
                weighted_scores,
                layer_covariates,
                layer_pattern,
                layer_name=str(layer_name),
                alpha=float(branching_config["alpha"]),
                branching_config=layer_config,
                weight_column="analysis_weight",
            )
            states = classify_branch_states(slopes, within, layer_config)
            for frame in (slopes, within, between, states):
                frame.insert(0, "selection_mode", str(mode))
            all_slopes.append(slopes)
            all_within.append(within)
            all_between.append(between)
            all_states.append(states)

    weight_table = pd.concat(all_weights, ignore_index=True)
    fit_table = pd.concat(all_fits, ignore_index=True)
    diagnostic_table = pd.concat(all_diagnostics, ignore_index=True)
    slopes = pd.concat(all_slopes, ignore_index=True)
    within = pd.concat(all_within, ignore_index=True)
    between = pd.concat(all_between, ignore_index=True)
    states = pd.concat(all_states, ignore_index=True)
    summary = _headline_summary(states, between, selection_config)
    return (
        weight_table,
        fit_table,
        diagnostic_table,
        slopes,
        within,
        between,
        states,
        summary,
    )


@app.command("run")
def run(
    island_scores_csv: Annotated[Path, typer.Option(exists=True)],
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
    outputs = run_selection_sensitivity(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pattern_config,
        branching_config,
        selection_config,
    )
    names = [
        "branch_selection_weights.csv.gz",
        "branch_selection_model_fits.csv",
        "branch_selection_weight_diagnostics.csv",
        "branch_selection_slopes.csv",
        "branch_selection_within_context.csv",
        "branch_selection_between_context.csv",
        "branch_selection_states.csv",
        "branch_selection_headline_summary.csv",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame, name in zip(outputs, names, strict=True):
        frame.to_csv(output_dir / name, index=False)
    weights, fits, diagnostics, _, _, _, _, summary = outputs
    manifest = {
        "contract": str(selection_config["contract"]),
        "n_full_islands_regime_selection_universe": int(
            fits.loc[fits["context_layer"].eq("analysis_regime"), "n_full_islands"].max()
        ),
        "selection_uses_branch_score_values": False,
        "selection_uses_effect_estimates": False,
        "selection_uses_p_values": False,
        "weight_modes": [str(x) for x in selection_config["weight_modes"]],
        "maximum_realized_weight": float(diagnostics["weight_max"].max()),
        "minimum_effective_sample_fraction": float(
            diagnostics.loc[
                diagnostics["selection_mode"].ne("unweighted"),
                "effective_sample_fraction",
            ].min()
        ),
        "n_weight_rows": len(weights),
        "n_headline_rows": len(summary),
        "pollination_mechanism_identified": False,
        "source_pool_standardized": False,
    }
    (output_dir / "branch_selection_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
