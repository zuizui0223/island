"""Within-Palearctic observation-selection sensitivity for restricted floral responses.

The global realm-interaction propensity model can fail to converge after exact SI / mating-
system restriction because the observed restricted score is sparse across all realms. The
scientific estimand here is narrower: the Palearctic isolation slope itself. This module
therefore fits the same outcome-blind availability model *within the Palearctic universe*,
using the same measured continuous predictors (distance, island area, climate), the same
propensity bounds and the same stabilized-IPW rules as PR138.

No response score, effect estimate or p-value enters the propensity model. The analysis is
run both on the full Palearctic universe and after removing the preidentified Aegean block
``lat12_lon20``. It is a missing-at-random sensitivity, not a proof of unbiased sampling.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design
from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.chapter1_pr138_outcrossing_selection_stress import (
    PREIDENTIFIED_AEGEAN_BLOCK,
    _bh,
    _standardize,
    build_restricted_attraction_scores,
)


def _expit(values: np.ndarray) -> np.ndarray:
    out = np.empty_like(values, dtype=float)
    positive = values >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-values[positive]))
    negative = ~positive
    exp_values = np.exp(values[negative])
    out[negative] = exp_values / (1.0 + exp_values)
    return out


def _palearctic_universe(
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    *,
    excluded_block: str | None,
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    universe = covariates.merge(
        realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    universe = universe.loc[universe["biogeographic_realm"].eq("Palearctic")].copy()
    if excluded_block is not None:
        universe = universe.loc[universe[cluster].astype(str).ne(excluded_block)].copy()
    for column in [geography, *baseline]:
        universe[column] = pd.to_numeric(universe[column], errors="coerce")
    universe[cluster] = universe[cluster].fillna("").astype(str)
    return universe.dropna(subset=[geography, *baseline]).loc[lambda x: x[cluster].ne("")].copy()


def _selection_design(
    universe: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> tuple[np.ndarray, list[str]]:
    geography = str(pattern_config["geography_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    predictors = [*baseline, geography]
    names = ["intercept"]
    columns = [np.ones(len(universe), dtype=float)]
    for predictor in predictors:
        values = pd.to_numeric(universe[predictor], errors="coerce").to_numpy(float)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=0))
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant selection predictor: {predictor}")
        names.append(f"z_{predictor}")
        columns.append((values - mean) / sd)
    return np.column_stack(columns), names


def build_within_palearctic_weights(
    observed_ids: set[str],
    universe: pd.DataFrame,
    pattern_config: dict[str, Any],
    selection_config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    cluster = str(pattern_config["cluster_column"])
    outcome = universe["island_id"].astype(str).isin(observed_ids).to_numpy(float)
    n_observed = int(outcome.sum())
    n_unobserved = int(len(outcome) - n_observed)
    minimum_observed = int(selection_config["selection_model"]["minimum_observed"])
    minimum_unobserved = int(selection_config["selection_model"]["minimum_unobserved"])
    unweighted = universe.loc[outcome.astype(bool), ["island_id"]].copy()
    unweighted["selection_mode"] = "unweighted"
    unweighted["analysis_weight"] = 1.0
    diagnostics = [
        {
            "selection_mode": "unweighted",
            "n_observed": n_observed,
            "n_unobserved": n_unobserved,
            "n_propensity_clipped": 0,
            "weight_min": 1.0,
            "weight_median": 1.0,
            "weight_max": 1.0,
            "effective_sample_size": float(n_observed),
            "effective_sample_fraction": 1.0 if n_observed else np.nan,
        }
    ]
    fit_base: dict[str, Any] = {"n_observed": n_observed, "n_unobserved": n_unobserved}
    if n_observed < minimum_observed or n_unobserved < minimum_unobserved:
        return unweighted, {**fit_base, "status": "below_selection_support"}, pd.DataFrame(diagnostics)

    design, names = _selection_design(universe, pattern_config)
    coefficients, fit, _ = _fit_grouped_binomial_design(
        outcome,
        np.ones(len(outcome), dtype=float),
        design,
        names,
        universe[cluster].astype(str).to_numpy(),
    )
    if coefficients.empty or fit.get("status") != "fit":
        return unweighted, {**fit_base, **fit}, pd.DataFrame(diagnostics)
    beta = coefficients.set_index("predictor").loc[names, "estimate_log_odds"].to_numpy(float)
    propensity = np.clip(_expit(design @ beta), 1e-8, 1 - 1e-8)
    prevalence = n_observed / len(universe)
    observed = universe.loc[outcome.astype(bool), ["island_id"]].copy()
    observed["propensity"] = propensity[outcome.astype(bool)]
    weight_parts = [unweighted]
    for mode, spec in selection_config["weight_modes"].items():
        if mode == "unweighted":
            continue
        raw = observed["propensity"].to_numpy(float)
        bounded = np.clip(raw, float(spec["minimum_propensity"]), float(spec["maximum_propensity"]))
        clipped = raw != bounded
        weights = prevalence / bounded
        if bool(selection_config["stabilization"]["normalize_within_axis_stratum_context"]):
            weights = weights / float(np.mean(weights))
        weights = np.minimum(weights, float(spec["maximum_weight"]))
        frame = observed[["island_id"]].copy()
        frame["selection_mode"] = str(mode)
        frame["analysis_weight"] = weights
        weight_parts.append(frame)
        sw = float(np.sum(weights))
        ess = sw**2 / float(np.sum(np.square(weights)))
        diagnostics.append(
            {
                "selection_mode": str(mode),
                "n_observed": n_observed,
                "n_unobserved": n_unobserved,
                "n_propensity_clipped": int(clipped.sum()),
                "weight_min": float(np.min(weights)),
                "weight_median": float(np.median(weights)),
                "weight_max": float(np.max(weights)),
                "effective_sample_size": ess,
                "effective_sample_fraction": ess / n_observed,
            }
        )
    fit_row = {
        **fit_base,
        **fit,
        "propensity_min": float(np.min(propensity)),
        "propensity_median": float(np.median(propensity)),
        "propensity_max": float(np.max(propensity)),
    }
    return pd.concat(weight_parts, ignore_index=True), fit_row, pd.DataFrame(diagnostics)


def _fit_response(
    scores: pd.DataFrame,
    weights: pd.DataFrame,
    universe: pd.DataFrame,
    pattern_config: dict[str, Any],
    selection_mode: str,
) -> dict[str, Any]:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    work = scores.merge(
        weights.loc[weights["selection_mode"].eq(selection_mode), ["island_id", "analysis_weight"]],
        on="island_id",
        how="inner",
        validate="one_to_one",
    )
    work = work.merge(
        universe[["island_id", geography, cluster, *baseline]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["syndrome_score", "analysis_weight", geography, *baseline]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cluster] = work[cluster].fillna("").astype(str)
    work = work.dropna(subset=["syndrome_score", "analysis_weight", geography, *baseline]).loc[
        lambda x: x[cluster].ne("")
    ].copy()
    threshold = int(pattern_config["support_tiers"]["confirmatory"])
    n_islands = int(work["island_id"].nunique())
    if n_islands < threshold:
        return {"status": "not_testable", "n_unique_islands": n_islands}
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    design = np.column_stack(
        [
            np.ones(len(work), dtype=float),
            _standardize(work[geography]),
            *[_standardize(work[x]) for x in baseline],
        ]
    )
    coef, _, fit = _fit_weighted_clustered_design(
        work["syndrome_score"].to_numpy(float),
        work["analysis_weight"].to_numpy(float),
        design,
        names,
        work[cluster].astype(str).to_numpy(),
    )
    if coef.empty:
        return {"status": str(fit.get("status", "fit_failed")), "n_unique_islands": n_islands}
    distance = coef.set_index("predictor").loc[f"z_{geography}"]
    w = work["analysis_weight"].to_numpy(float)
    sw = float(np.sum(w))
    ess = sw**2 / float(np.sum(np.square(w)))
    return {
        "status": "fit",
        "n_unique_islands": n_islands,
        "n_clusters": int(fit["n_clusters"]),
        "distance_estimate": float(distance["estimate"]),
        "distance_se": float(distance["cluster_robust_se"]),
        "distance_p": float(distance["p_value"]),
        "effective_sample_size": ess,
        "effective_sample_fraction": ess / len(work),
        "weight_max": float(np.max(w)),
    }


def run_palearctic_restricted_ipw(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    selection_config: dict[str, Any],
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    attraction = build_restricted_attraction_scores(adjusted_scores)
    restrictions = sorted(attraction["restriction"].dropna().astype(str).unique())
    source_modes = sorted(attraction["source_mode"].dropna().astype(str).unique())
    strata = [x for x in ["all_native", "native_nonendemic"] if x in set(attraction["stratum"])]
    universe_specs = {
        "full_palearctic": None,
        "palearctic_drop_aegean_lat12_lon20": PREIDENTIFIED_AEGEAN_BLOCK,
    }
    rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    diag_parts: list[pd.DataFrame] = []
    for universe_mode, excluded_block in universe_specs.items():
        universe = _palearctic_universe(
            covariates,
            realm_assignment,
            pattern_config,
            excluded_block=excluded_block,
        )
        valid_ids = set(universe["island_id"].astype(str))
        for restriction in restrictions:
            for source_mode in source_modes:
                scenario = attraction.loc[
                    attraction["restriction"].eq(restriction)
                    & attraction["source_mode"].eq(source_mode)
                    & attraction["island_id"].astype(str).isin(valid_ids)
                ].copy()
                for stratum in strata:
                    scores = scenario.loc[
                        scenario["stratum"].eq(stratum) & scenario["syndrome_score"].notna()
                    ].copy()
                    observed_ids = set(scores["island_id"].astype(str))
                    weights, fit, diagnostics = build_within_palearctic_weights(
                        observed_ids,
                        universe,
                        pattern_config,
                        selection_config,
                    )
                    fit_rows.append(
                        {
                            "universe_mode": universe_mode,
                            "restriction": restriction,
                            "source_mode": source_mode,
                            "stratum": stratum,
                            **fit,
                        }
                    )
                    diagnostics.insert(0, "stratum", stratum)
                    diagnostics.insert(0, "source_mode", source_mode)
                    diagnostics.insert(0, "restriction", restriction)
                    diagnostics.insert(0, "universe_mode", universe_mode)
                    diag_parts.append(diagnostics)
                    for selection_mode in selection_config["weight_modes"]:
                        response = _fit_response(
                            scores,
                            weights,
                            universe,
                            pattern_config,
                            str(selection_mode),
                        )
                        rows.append(
                            {
                                "universe_mode": universe_mode,
                                "excluded_block": excluded_block or "",
                                "restriction": restriction,
                                "source_mode": source_mode,
                                "stratum": stratum,
                                "selection_mode": str(selection_mode),
                                **response,
                            }
                        )
    result = pd.DataFrame(rows)
    result["distance_q_across_source_modes"] = np.nan
    fit_mask = result["status"].eq("fit")
    if fit_mask.any():
        family = ["universe_mode", "restriction", "stratum", "selection_mode"]
        result.loc[fit_mask, "distance_q_across_source_modes"] = (
            result.loc[fit_mask]
            .groupby(family, group_keys=False)["distance_p"]
            .apply(_bh)
        )
    manifest = {
        "contract": "chapter1_pr138_palearctic_restricted_ipw_v1",
        "target_population": "Palearctic islands only",
        "response": "source-adjusted attraction_shift within exact reproductive restriction",
        "selection_predictors": [
            *[str(x) for x in pattern_config["baseline_covariates"]],
            str(pattern_config["geography_column"]),
        ],
        "score_values_used_in_propensity": False,
        "effect_estimates_used_in_propensity": False,
        "p_values_used_in_propensity": False,
        "universe_modes": list(universe_specs),
        "preidentified_aegean_block": PREIDENTIFIED_AEGEAN_BLOCK,
        "claim_ceiling": "Within-Palearctic MAR sensitivity only; persistence does not identify pollinator loss or prove random observation.",
    }
    return {
        "results": result,
        "selection_fits": pd.DataFrame(fit_rows),
        "diagnostics": pd.concat(diag_parts, ignore_index=True) if diag_parts else pd.DataFrame(),
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adjusted-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--selection-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    outputs = run_palearctic_restricted_ipw(
        pd.read_csv(args.adjusted_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.selection_config_path.read_text(encoding="utf-8")),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in {
        "results": "palearctic_restricted_ipw_results.csv",
        "selection_fits": "palearctic_restricted_ipw_selection_fits.csv",
        "diagnostics": "palearctic_restricted_ipw_diagnostics.csv",
    }.items():
        frame = outputs[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(args.output_dir / filename, index=False)
    manifest = outputs["manifest"]
    assert isinstance(manifest, dict)
    (args.output_dir / "palearctic_restricted_ipw_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
