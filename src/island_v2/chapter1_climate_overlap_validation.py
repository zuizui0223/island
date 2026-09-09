"""Prospective V1 falsification of categorical biogeographic branching.

The baseline Chapter 1 model adjusts for climate main effects. This module asks two
stronger questions without changing that historical result:

1. Does the predeclared distance-by-context response-vector contrast survive
   outcome-blind overlap weighting on island area and frozen climate PCs?
2. Does it survive after distance effects are allowed to vary continuously with the
   frozen climate PCs?

Neither weighting nor model construction uses trait values, fitted ecological effects,
or p-values. Support loss and poor overlap are retained as scientific results. The module
does not identify a pollinator, pollination service, historical loss, or replacement.
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
from island_v2.chapter1_global_branch_selection import _context_universe, _expit
from island_v2.chapter1_global_branching import build_branch_scores
from island_v2.chapter1_pr136_biogeographic_residual import (
    _fit_weighted_clustered_design,
)
from island_v2.chapter1_pr138_syndrome_analysis import (
    _between_contexts,
    _bh,
    _joint_wald,
    _prepare,
    _standardize,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _validation_spec(explanation_config: dict[str, Any]) -> dict[str, Any]:
    try:
        return explanation_config["validations"][
            "V1_H2_environment_vs_biogeography"
        ]["frozen_implementation"]
    except KeyError as exc:
        raise typer.BadParameter("V1 frozen implementation is incomplete") from exc


def _weighted_moments(values: np.ndarray, weights: np.ndarray) -> tuple[float, float]:
    total = float(np.sum(weights))
    if total <= 0:
        return float("nan"), float("nan")
    mean = float(np.sum(weights * values) / total)
    variance = float(np.sum(weights * np.square(values - mean)) / total)
    return mean, max(variance, 0.0)


def _standardized_mean_difference(
    values_a: np.ndarray,
    values_b: np.ndarray,
    weights_a: np.ndarray,
    weights_b: np.ndarray,
) -> float:
    mean_a, variance_a = _weighted_moments(values_a, weights_a)
    mean_b, variance_b = _weighted_moments(values_b, weights_b)
    denominator = math.sqrt(max((variance_a + variance_b) / 2.0, 0.0))
    if not math.isfinite(denominator) or denominator <= 0:
        return float("nan")
    return (mean_b - mean_a) / denominator


def _effective_sample_size(weights: np.ndarray) -> float:
    total = float(np.sum(weights))
    squared = float(np.sum(np.square(weights)))
    return total**2 / squared if squared > 0 else 0.0


def build_context_overlap_weights(
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    explanation_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Fit frozen context-overlap weights on the complete covariate universe."""

    spec = _validation_spec(explanation_config)
    predictors = [
        str(x)
        for x in explanation_config["validations"][
            "V1_H2_environment_vs_biogeography"
        ]["analyses"][0]["balance_predictors"]
    ]
    overlap_spec = spec["overlap_model"]
    minimum_context_islands = int(overlap_spec["minimum_context_islands"])
    clip = float(overlap_spec["propensity_clip"])
    poor_spec = overlap_spec["poor_overlap_if_either"]
    minimum_ess_fraction = float(poor_spec["minimum_effective_sample_fraction"])
    maximum_weighted_smd = float(
        poor_spec["maximum_absolute_weighted_standardized_mean_difference"]
    )
    cluster = str(pattern_config["cluster_column"])

    weight_parts: list[pd.DataFrame] = []
    coefficient_parts: list[pd.DataFrame] = []
    balance_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for layer_name, contrasts in spec["headline_contrasts"].items():
        if layer_name not in branching_config["context_layers"]:
            raise typer.BadParameter(f"unknown V1 context layer: {layer_name}")
        layer_spec = branching_config["context_layers"][layer_name]
        context_column = str(layer_spec["column"])
        universe, universe_summary = _context_universe(
            covariates,
            realm_assignment,
            pattern_config,
            str(layer_name),
            layer_spec,
        )
        missing = set(predictors) - set(universe.columns)
        if missing:
            raise typer.BadParameter(
                f"V1 balance predictors missing: {sorted(missing)}"
            )
        for context_a, context_b in contrasts:
            context_a = str(context_a)
            context_b = str(context_b)
            work = universe.loc[
                universe[context_column].isin([context_a, context_b])
            ].copy()
            counts = work[context_column].value_counts()
            n_a = int(counts.get(context_a, 0))
            n_b = int(counts.get(context_b, 0))
            base = {
                "context_layer": str(layer_name),
                "context_column": context_column,
                "context_a": context_a,
                "context_b": context_b,
                "n_context_a": n_a,
                "n_context_b": n_b,
                **universe_summary,
            }
            if n_a < minimum_context_islands or n_b < minimum_context_islands:
                summary_rows.append(
                    {
                        **base,
                        "status": "below_context_support",
                        "poor_overlap": True,
                    }
                )
                continue

            names = ["intercept"]
            columns = [np.ones(len(work), dtype=float)]
            for predictor in predictors:
                names.append(f"z_{predictor}")
                columns.append(_standardize(work[predictor]))
            design = np.column_stack(columns)
            context_b_indicator = work[context_column].eq(context_b).to_numpy(float)
            coefficients, fit, _ = _fit_grouped_binomial_design(
                context_b_indicator,
                np.ones(len(work), dtype=float),
                design,
                names,
                work[cluster].astype(str).to_numpy(),
            )
            coefficient_frame = coefficients.copy()
            for key, value in base.items():
                coefficient_frame.insert(len(coefficient_frame.columns), key, value)
            coefficient_frame["fit_status"] = str(fit.get("status", ""))
            coefficient_parts.append(coefficient_frame)
            if coefficients.empty or fit.get("status") != "fit":
                summary_rows.append(
                    {
                        **base,
                        **fit,
                        "poor_overlap": True,
                    }
                )
                continue

            beta = (
                coefficients.set_index("predictor")
                .loc[names, "estimate_log_odds"]
                .to_numpy(float)
            )
            propensity = np.clip(_expit(design @ beta), clip, 1.0 - clip)
            raw_weight = np.where(
                context_b_indicator.astype(bool),
                1.0 - propensity,
                propensity,
            )
            analysis_weight = raw_weight.copy()
            if bool(overlap_spec["normalize_mean_weight_within_context"]):
                normalizer = (
                    pd.DataFrame(
                        {
                            "context": work[context_column].to_numpy(str),
                            "weight": analysis_weight,
                        }
                    )
                    .groupby("context")["weight"]
                    .transform("mean")
                    .to_numpy(float)
                )
                analysis_weight = analysis_weight / normalizer

            weights = work[["island_id", context_column]].copy()
            weights["context_layer"] = str(layer_name)
            weights["context_a"] = context_a
            weights["context_b"] = context_b
            weights["propensity_context_b"] = propensity
            weights["raw_overlap_weight"] = raw_weight
            weights["overlap_weight"] = analysis_weight
            weight_parts.append(weights)

            mask_a = work[context_column].eq(context_a).to_numpy()
            mask_b = work[context_column].eq(context_b).to_numpy()
            maximum_absolute_after = 0.0
            for predictor in predictors:
                values = pd.to_numeric(work[predictor], errors="coerce").to_numpy(float)
                before = _standardized_mean_difference(
                    values[mask_a],
                    values[mask_b],
                    np.ones(n_a),
                    np.ones(n_b),
                )
                after = _standardized_mean_difference(
                    values[mask_a],
                    values[mask_b],
                    analysis_weight[mask_a],
                    analysis_weight[mask_b],
                )
                if math.isfinite(after):
                    maximum_absolute_after = max(maximum_absolute_after, abs(after))
                balance_rows.append(
                    {
                        **base,
                        "predictor": predictor,
                        "unweighted_standardized_mean_difference": before,
                        "overlap_weighted_standardized_mean_difference": after,
                    }
                )

            ess_a = _effective_sample_size(analysis_weight[mask_a])
            ess_b = _effective_sample_size(analysis_weight[mask_b])
            ess_fraction_a = ess_a / n_a
            ess_fraction_b = ess_b / n_b
            poor_overlap = bool(
                min(ess_fraction_a, ess_fraction_b) < minimum_ess_fraction
                or maximum_absolute_after > maximum_weighted_smd
            )
            summary_rows.append(
                {
                    **base,
                    **fit,
                    "propensity_min": float(np.min(propensity)),
                    "propensity_median": float(np.median(propensity)),
                    "propensity_max": float(np.max(propensity)),
                    "weight_min": float(np.min(analysis_weight)),
                    "weight_median": float(np.median(analysis_weight)),
                    "weight_max": float(np.max(analysis_weight)),
                    "context_a_effective_sample_size": ess_a,
                    "context_b_effective_sample_size": ess_b,
                    "context_a_effective_sample_fraction": ess_fraction_a,
                    "context_b_effective_sample_fraction": ess_fraction_b,
                    "maximum_absolute_weighted_smd": maximum_absolute_after,
                    "poor_overlap": poor_overlap,
                    "positivity_pass": not poor_overlap,
                }
            )

    weights = pd.concat(weight_parts, ignore_index=True) if weight_parts else pd.DataFrame()
    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True)
        if coefficient_parts
        else pd.DataFrame()
    )
    return weights, coefficients, pd.DataFrame(balance_rows), pd.DataFrame(summary_rows)


def _climate_interaction_between_contexts(
    data: pd.DataFrame,
    *,
    stratum: str,
    context_a: str,
    context_b: str,
    support_tier: str,
    threshold: int,
    pattern_config: dict[str, Any],
    climate_predictors: list[str],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    geography = str(pattern_config["geography_column"])
    context = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    work = data.loc[
        data["stratum"].eq(stratum)
        & data[context].isin([context_a, context_b])
    ].copy()
    counts = work.groupby(["syndrome", context])["island_id"].nunique().unstack(
        fill_value=0
    )
    for value in [context_a, context_b]:
        if value not in counts.columns:
            counts[value] = 0
    syndromes = sorted(
        counts.index[
            (counts[context_a] >= threshold) & (counts[context_b] >= threshold)
        ].astype(str)
    )
    base = {
        "stratum": stratum,
        "support_tier": support_tier,
        "threshold": threshold,
        "context_a": context_a,
        "context_b": context_b,
        "n_retained_syndromes": len(syndromes),
        "retained_syndromes": "|".join(syndromes),
    }
    if len(syndromes) < 2:
        return pd.DataFrame(), {**base, "status": "not_testable"}
    work = work.loc[work["syndrome"].isin(syndromes)].copy()
    context_b_indicator = work[context].eq(context_b).to_numpy(float)
    names: list[str] = []
    columns: list[np.ndarray] = []
    target_names: list[str] = []
    try:
        for syndrome in syndromes:
            mask = work["syndrome"].eq(syndrome).to_numpy()
            indicator = mask.astype(float)
            names.append(f"syndrome[{syndrome}]")
            columns.append(indicator)
            standardized: dict[str, np.ndarray] = {}
            for predictor in baseline:
                z = np.zeros(len(work), dtype=float)
                z[mask] = _standardize(work.loc[mask, predictor])
                standardized[predictor] = z
                names.append(f"syndrome[{syndrome}]:z_{predictor}")
                columns.append(z)
            names.append(f"syndrome[{syndrome}]:context[{context_b}]")
            columns.append(indicator * context_b_indicator)
            z_geography = np.zeros(len(work), dtype=float)
            z_geography[mask] = _standardize(work.loc[mask, geography])
            names.append(f"syndrome[{syndrome}]:z_{geography}")
            columns.append(z_geography)
            for predictor in climate_predictors:
                interaction_name = (
                    f"syndrome[{syndrome}]:z_{geography}:z_{predictor}"
                )
                names.append(interaction_name)
                columns.append(z_geography * standardized[predictor])
            target = (
                f"syndrome[{syndrome}]:z_{geography}:context[{context_b}]"
            )
            names.append(target)
            columns.append(z_geography * context_b_indicator)
            target_names.append(target)
    except ValueError:
        return pd.DataFrame(), {**base, "status": "constant_or_invalid_predictor"}

    coefficients, covariance, fit = _fit_weighted_clustered_design(
        work["syndrome_score"].to_numpy(float),
        np.ones(len(work), dtype=float),
        np.column_stack(columns),
        names,
        work[cluster].to_numpy(str),
    )
    if coefficients.empty:
        return coefficients, {**base, "status": str(fit.get("status", "fit_failed"))}
    indexed = coefficients.set_index("predictor")
    target_indices = [names.index(name) for name in target_names]
    target_vector = indexed.loc[target_names, "estimate"].to_numpy(float)
    target_covariance = covariance[np.ix_(target_indices, target_indices)]
    statistic, degrees_freedom, p_value = _joint_wald(
        target_vector, target_covariance
    )
    target_rows: list[dict[str, Any]] = []
    for syndrome, target_name in zip(syndromes, target_names, strict=True):
        row = indexed.loc[target_name]
        target_rows.append(
            {
                **base,
                "syndrome": syndrome,
                "residual_distance_x_context_estimate": float(row["estimate"]),
                "cluster_robust_se": float(row["cluster_robust_se"]),
                "p_value": float(row["p_value"]),
                "n_context_a": int(counts.loc[syndrome, context_a]),
                "n_context_b": int(counts.loc[syndrome, context_b]),
            }
        )
    return pd.DataFrame(target_rows), {
        **base,
        "status": "fit",
        "n_unique_islands": int(work["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "joint_context_difference_chisq": statistic,
        "joint_context_difference_df": degrees_freedom,
        "p_value": p_value,
    }


def _run_pairwise_tests(
    primary_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    overlap_weights: pd.DataFrame,
    overlap_summary: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    explanation_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    spec = _validation_spec(explanation_config)
    climate_predictors = [
        str(x)
        for x in spec["distance_x_climate_model"][
            "climate_interaction_predictors"
        ]
    ]
    vector_rows: list[dict[str, Any]] = []
    climate_coefficient_parts: list[pd.DataFrame] = []
    strata = [str(x) for x in pattern_config["strata"]]
    support_tiers = {
        str(name): int(value)
        for name, value in pattern_config["support_tiers"].items()
    }
    primary_axis_set = str(spec["primary_axis_set"])

    for layer_name, contrasts in spec["headline_contrasts"].items():
        layer_spec = branching_config["context_layers"][layer_name]
        context_column = str(layer_spec["column"])
        layer_covariates, _ = _context_universe(
            covariates,
            realm_assignment,
            pattern_config,
            str(layer_name),
            layer_spec,
        )
        for context_a, context_b in contrasts:
            context_a = str(context_a)
            context_b = str(context_b)
            pair_pattern = deepcopy(pattern_config)
            pair_pattern["context_column"] = context_column
            pair_pattern["contexts"] = [context_a, context_b]
            pair_scores = primary_scores.copy()
            overlap_available = bool(
                not overlap_weights.empty
                and {
                    "context_layer",
                    "context_a",
                    "context_b",
                    "island_id",
                    "overlap_weight",
                }.issubset(overlap_weights.columns)
            )
            if overlap_available:
                pair_weights = overlap_weights.loc[
                    overlap_weights["context_layer"].eq(str(layer_name))
                    & overlap_weights["context_a"].eq(context_a)
                    & overlap_weights["context_b"].eq(context_b),
                    ["island_id", "overlap_weight"],
                ]
                overlap_available = not pair_weights.empty
            else:
                pair_weights = pd.DataFrame(columns=["island_id", "overlap_weight"])
            weighted_scores = pair_scores.merge(
                pair_weights,
                on="island_id",
                how="inner",
                validate="many_to_one",
            )
            unweighted_data = _prepare(pair_scores, layer_covariates, pair_pattern)
            weighted_data = (
                _prepare(weighted_scores, layer_covariates, pair_pattern)
                if overlap_available
                else pd.DataFrame()
            )
            positivity = overlap_summary.loc[
                overlap_summary["context_layer"].eq(str(layer_name))
                & overlap_summary["context_a"].eq(context_a)
                & overlap_summary["context_b"].eq(context_b)
            ]
            positivity_pass = bool(
                not positivity.empty and positivity.iloc[0].get("positivity_pass", False)
            )

            for stratum in strata:
                for support_tier, threshold in support_tiers.items():
                    reference = _between_contexts(
                        unweighted_data,
                        stratum=stratum,
                        context_a=context_a,
                        context_b=context_b,
                        support_tier=support_tier,
                        threshold=threshold,
                        pattern_config=pair_pattern,
                    )
                    vector_rows.append(
                        {
                            "validation_mode": "reference_climate_main_effects",
                            "context_layer": str(layer_name),
                            "axis_set": primary_axis_set,
                            **reference,
                            "positivity_pass": np.nan,
                        }
                    )

                    if overlap_available:
                        overlap = _between_contexts(
                            weighted_data,
                            stratum=stratum,
                            context_a=context_a,
                            context_b=context_b,
                            support_tier=support_tier,
                            threshold=threshold,
                            pattern_config=pair_pattern,
                            weight_column="overlap_weight",
                        )
                    else:
                        overlap = {
                            "stratum": stratum,
                            "support_tier": support_tier,
                            "threshold": threshold,
                            "context_a": context_a,
                            "context_b": context_b,
                            "n_retained_syndromes": 0,
                            "retained_syndromes": "",
                            "status": "overlap_model_not_testable",
                        }
                    vector_rows.append(
                        {
                            "validation_mode": "outcome_blind_overlap_weighted",
                            "context_layer": str(layer_name),
                            "axis_set": primary_axis_set,
                            **overlap,
                            "positivity_pass": positivity_pass,
                        }
                    )

                    coefficients, climate = _climate_interaction_between_contexts(
                        unweighted_data,
                        stratum=stratum,
                        context_a=context_a,
                        context_b=context_b,
                        support_tier=support_tier,
                        threshold=threshold,
                        pattern_config=pair_pattern,
                        climate_predictors=climate_predictors,
                    )
                    vector_rows.append(
                        {
                            "validation_mode": "distance_x_climate_adjusted",
                            "context_layer": str(layer_name),
                            "axis_set": primary_axis_set,
                            **climate,
                            "positivity_pass": np.nan,
                        }
                    )
                    if not coefficients.empty:
                        coefficients.insert(0, "axis_set", primary_axis_set)
                        coefficients.insert(0, "context_layer", str(layer_name))
                        climate_coefficient_parts.append(coefficients)

    vector_tests = pd.DataFrame(vector_rows)
    if "p_value" not in vector_tests.columns:
        vector_tests["p_value"] = np.nan
    family = ["validation_mode", "axis_set", "stratum", "support_tier"]
    vector_tests["q_v1_contrast_family"] = vector_tests.groupby(
        family, group_keys=False
    )["p_value"].transform(_bh)
    alpha = float(spec["multiplicity"]["alpha"])
    vector_tests["v1_vector_difference_supported"] = (
        vector_tests["status"].eq("fit")
        & vector_tests["q_v1_contrast_family"].le(alpha)
    )
    climate_coefficients = (
        pd.concat(climate_coefficient_parts, ignore_index=True)
        if climate_coefficient_parts
        else pd.DataFrame()
    )
    return vector_tests, climate_coefficients


def run_climate_overlap_validation(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    explanation_config: dict[str, Any],
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
]:
    spec = _validation_spec(explanation_config)
    branch_scores = build_branch_scores(island_scores, branching_config)
    primary_axis_set = str(spec["primary_axis_set"])
    primary_axes = [
        str(x) for x in branching_config["axis_sets"][primary_axis_set]["axes"]
    ]
    primary_scores = branch_scores.loc[
        branch_scores["syndrome"].isin(primary_axes)
    ].copy()
    weights, overlap_coefficients, balance, overlap_summary = (
        build_context_overlap_weights(
            covariates,
            realm_assignment,
            pattern_config,
            branching_config,
            explanation_config,
        )
    )
    vector_tests, climate_coefficients = _run_pairwise_tests(
        primary_scores,
        covariates,
        realm_assignment,
        weights,
        overlap_summary,
        pattern_config,
        branching_config,
        explanation_config,
    )
    return (
        weights,
        overlap_coefficients,
        balance,
        overlap_summary,
        vector_tests,
        climate_coefficients,
    )


@app.command("run")
def run(
    island_scores_csv: Annotated[Path, typer.Option(exists=True)],
    covariates_csv: Annotated[Path, typer.Option(exists=True)],
    realm_assignment_csv: Annotated[Path, typer.Option(exists=True)],
    pattern_config_path: Annotated[Path, typer.Option(exists=True)],
    branching_config_path: Annotated[Path, typer.Option(exists=True)],
    explanation_config_path: Annotated[Path, typer.Option(exists=True)],
    evidence_scope: Annotated[str, typer.Option()],
    output_dir: Annotated[Path, typer.Option()],
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(
        branching_config_path.read_text(encoding="utf-8")
    )
    explanation_config = yaml.safe_load(
        explanation_config_path.read_text(encoding="utf-8")
    )
    outputs = run_climate_overlap_validation(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(realm_assignment_csv),
        pattern_config,
        branching_config,
        explanation_config,
    )
    names = [
        "v1_overlap_weights.csv.gz",
        "v1_overlap_model_coefficients.csv",
        "v1_overlap_balance_diagnostics.csv",
        "v1_overlap_summary.csv",
        "v1_vector_tests.csv",
        "v1_climate_interaction_coefficients.csv",
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame, name in zip(outputs, names, strict=True):
        frame.to_csv(output_dir / name, index=False)
    weights, _, balance, overlap_summary, vector_tests, _ = outputs
    finite_weight = pd.to_numeric(weights.get("overlap_weight"), errors="coerce")
    manifest = {
        "contract": "chapter1_V1_climate_overlap_validation_v1",
        "status": "post_baseline_prospective_falsification",
        "evidence_scope": evidence_scope,
        "primary_axis_set": _validation_spec(explanation_config)["primary_axis_set"],
        "n_overlap_weight_rows": len(weights),
        "n_balance_rows": len(balance),
        "n_vector_tests": len(vector_tests),
        "n_poor_overlap_contrasts": int(
            overlap_summary.get("poor_overlap", pd.Series(dtype=bool))
            .fillna(True)
            .astype(bool)
            .sum()
        ),
        "maximum_realized_overlap_weight": float(finite_weight.max()),
        "weight_model_uses_trait_availability": False,
        "weight_model_uses_trait_values": False,
        "weight_model_uses_effect_estimates": False,
        "weight_model_uses_p_values": False,
        "response_support_gates_lowered": False,
        "pollination_mechanism_identified": False,
        "baseline_PR141_rewritten": False,
    }
    (output_dir / "chapter1_v1_climate_overlap_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
