"""Species-list detection tipping-point sensitivity for Chapter 1.

This is deliberately separate from the frozen V5 trait-resolution MNAR analysis.
V5 asks whether, conditional on a recorded flora, trait resolution depends on the
unobserved trait state. V6 asks whether the island species list itself could be
state-selective and increasingly incomplete with geographic isolation.

The analysis does not estimate true island richness. It places an explicit sensitivity
surface over assumed list completeness and a state-dependent recording odds ratio, then
refits the frozen primary plant-response model while retaining the original trait-scored
species count as the information weight.
"""
from __future__ import annotations

import itertools
import json
import math
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_global_branching import build_branch_scores
from island_v2.chapter1_pr138_syndrome_analysis import _between_contexts, _bh, _prepare, _within_context
from island_v2.chapter1_trait_resolution_mnar import (
    build_total_species_counts,
    completed_prevalence_from_resolution_or,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_species_detection_tipping_v1"
SCOPE_PATHS = {
    "all_analysis_eligible": "all",
    "direct_only": "direct",
}


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected species-detection tipping contract")
    return config


def _expit(x: np.ndarray | float) -> np.ndarray | float:
    x_arr = np.asarray(x, dtype=float)
    out = np.empty_like(x_arr)
    positive = x_arr >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-x_arr[positive]))
    ex = np.exp(x_arr[~positive])
    out[~positive] = ex / (1.0 + ex)
    if np.ndim(x) == 0:
        return float(out)
    return out


def _logit(p: float) -> float:
    p = float(p)
    if not 0 < p < 1:
        raise ValueError("completeness must be strictly between zero and one")
    return math.log(p / (1.0 - p))


def distance_standardization(covariates: pd.DataFrame, geography: str) -> tuple[float, float]:
    if geography not in covariates.columns:
        raise typer.BadParameter(f"covariates missing distance column: {geography}")
    values = pd.to_numeric(covariates[geography], errors="coerce").dropna().to_numpy(float)
    if len(values) == 0:
        raise ValueError("no finite distance values")
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("distance standard deviation is invalid")
    return mean, sd


def completeness_from_distance(
    distance: np.ndarray,
    *,
    distance_mean: float,
    distance_sd: float,
    median_completeness: float,
    distance_completeness_odds_ratio: float,
) -> np.ndarray:
    """Return assumed species-list completeness for each island.

    ``distance_completeness_odds_ratio`` multiplies the odds of list completeness per
    one frozen all-island SD of log distance. Values below one encode lower completeness
    on more isolated islands.
    """
    if distance_completeness_odds_ratio <= 0:
        raise ValueError("distance completeness odds ratio must be positive")
    z = (np.asarray(distance, dtype=float) - float(distance_mean)) / float(distance_sd)
    eta = _logit(float(median_completeness)) + math.log(float(distance_completeness_odds_ratio)) * z
    result = np.asarray(_expit(eta), dtype=float)
    if not np.isfinite(result).all() or np.any(result <= 0) or np.any(result >= 1):
        raise ValueError("computed completeness must lie strictly between zero and one")
    return result


def build_recorded_stratum_counts(status_flora: pd.DataFrame, strata: list[str]) -> pd.DataFrame:
    counts = build_total_species_counts(status_flora, strata).rename(
        columns={"n_total_stratum_species": "n_recorded_stratum_species"}
    )
    counts["n_recorded_stratum_species"] = pd.to_numeric(
        counts["n_recorded_stratum_species"], errors="coerce"
    )
    return counts


def adjust_species_detection(
    island_scores: pd.DataFrame,
    recorded_counts: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    geography: str,
    distance_mean: float,
    distance_sd: float,
    median_completeness: float,
    distance_completeness_odds_ratio: float,
    state_recording_odds_ratio: float,
    affected_syndrome: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply one species-detection sensitivity scenario to one evidence ledger.

    The observed syndrome prevalence is first projected to the full *recorded* stratum
    under the V5 MAR point (OR_R=1). The species-list selection model then inverts from
    recorded to assumed true prevalence. The original ``n_species`` is retained exactly.
    """
    required = {"island_id", "stratum", "syndrome", "syndrome_score", "n_species"}
    if missing := required - set(island_scores.columns):
        raise typer.BadParameter(f"island scores missing columns: {sorted(missing)}")
    if geography not in covariates.columns:
        raise typer.BadParameter(f"covariates missing geography column: {geography}")

    scores = island_scores.copy()
    scores["island_id"] = scores["island_id"].astype(str)
    scores["stratum"] = scores["stratum"].astype(str)
    scores["syndrome"] = scores["syndrome"].astype(str)
    scores["syndrome_score"] = pd.to_numeric(scores["syndrome_score"], errors="coerce")
    scores["n_species"] = pd.to_numeric(scores["n_species"], errors="coerce")

    counts = recorded_counts.copy()
    counts["island_id"] = counts["island_id"].astype(str)
    counts["stratum"] = counts["stratum"].astype(str)
    cov = covariates[["island_id", geography]].drop_duplicates("island_id").copy()
    cov["island_id"] = cov["island_id"].astype(str)
    cov[geography] = pd.to_numeric(cov[geography], errors="coerce")

    work = scores.merge(counts, on=["island_id", "stratum"], how="left", validate="many_to_one")
    work = work.merge(cov, on="island_id", how="left", validate="many_to_one")
    affected = work["syndrome"].eq(str(affected_syndrome))
    focal = work.loc[affected].copy()
    if focal.empty:
        raise ValueError(f"affected syndrome absent: {affected_syndrome}")
    if focal[["n_recorded_stratum_species", geography]].isna().any().any():
        raise ValueError("affected scores lack recorded-richness or distance inputs")
    if (focal["n_species"] > focal["n_recorded_stratum_species"] + 1e-12).any():
        raise ValueError("trait-scored species exceed recorded stratum richness")

    completeness = completeness_from_distance(
        focal[geography].to_numpy(float),
        distance_mean=distance_mean,
        distance_sd=distance_sd,
        median_completeness=median_completeness,
        distance_completeness_odds_ratio=distance_completeness_odds_ratio,
    )
    recorded_total = focal["n_recorded_stratum_species"].to_numpy(float)
    true_total = recorded_total / completeness
    observed_membership = np.clip(
        (focal["syndrome_score"].to_numpy(float) + 1.0) / 2.0, 0.0, 1.0
    )
    recorded_focal_count = observed_membership * recorded_total

    completed = np.asarray(
        [
            completed_prevalence_from_resolution_or(
                total_species=n_true,
                resolved_species=n_recorded,
                resolved_focal_count=n_focal,
                resolution_odds_ratio=float(state_recording_odds_ratio),
            )
            for n_true, n_recorded, n_focal in zip(
                true_total, recorded_total, recorded_focal_count, strict=True
            )
        ],
        dtype=float,
    )
    adjusted = 2.0 * completed - 1.0
    work.loc[affected, "syndrome_score"] = adjusted

    diagnostics = focal[
        ["island_id", "stratum", "syndrome", "syndrome_score", "n_species", geography]
    ].copy()
    diagnostics["observed_soft_membership"] = observed_membership
    diagnostics["completed_soft_membership"] = completed
    diagnostics["adjusted_syndrome_score"] = adjusted
    diagnostics["n_recorded_stratum_species"] = recorded_total
    diagnostics["assumed_list_completeness"] = completeness
    diagnostics["assumed_true_stratum_species"] = true_total
    diagnostics["assumed_unrecorded_species"] = true_total - recorded_total
    diagnostics["median_distance_completeness"] = float(median_completeness)
    diagnostics["distance_completeness_odds_ratio"] = float(
        distance_completeness_odds_ratio
    )
    diagnostics["state_recording_odds_ratio"] = float(state_recording_odds_ratio)

    original_columns = list(island_scores.columns)
    return work[original_columns], diagnostics


def _reduced_branching_config(branching_config: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(branching_config)
    keep_axes = {"accessibility_generalization", "reproductive_assurance"}
    out["branch_axes"] = {
        key: value for key, value in out["branch_axes"].items() if key in keep_axes
    }
    out["axis_sets"] = {
        "universal_plant_response": deepcopy(
            branching_config["axis_sets"]["universal_plant_response"]
        )
    }
    return out


def fit_headline_family(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    *,
    strata: list[str],
    support_tier: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Refit only the frozen primary vector while preserving original FDR families.

    All within-context fits are retained for both context layers because their p-values
    define the original FDR family. All six analysis-regime context pairs are retained for
    the direct North--Tropical vector family. Realm between-context fits are unnecessary
    for the V6 headline targets and are not run.
    """
    branch_cfg = _reduced_branching_config(branching_config)
    branch_scores = build_branch_scores(island_scores, branch_cfg)
    alpha = float(branch_cfg["alpha"])
    threshold = int(pattern_config["support_tiers"][support_tier])

    slopes_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []

    for layer_name, layer_spec in branch_cfg["context_layers"].items():
        contexts = [str(x) for x in layer_spec["contexts"]]
        column = str(layer_spec["column"])
        layer_covariates = covariates.copy()
        if column not in layer_covariates.columns:
            assignment = realm_assignment[["island_id", column]].drop_duplicates("island_id")
            layer_covariates = layer_covariates.merge(
                assignment, on="island_id", how="left", validate="one_to_one"
            )
        layer_pattern = deepcopy(pattern_config)
        layer_pattern["context_column"] = column
        layer_pattern["contexts"] = contexts
        data = _prepare(branch_scores, layer_covariates, layer_pattern)

        for stratum in strata:
            for context in contexts:
                slopes, result = _within_context(
                    data,
                    stratum=stratum,
                    context_value=context,
                    support_tier=support_tier,
                    threshold=threshold,
                    pattern_config=layer_pattern,
                    syndrome_config={},
                )
                if not slopes.empty:
                    slopes.insert(0, "axis_set", "universal_plant_response")
                    slopes.insert(0, "context_layer", str(layer_name))
                    slopes_parts.append(slopes)
                result["axis_set"] = "universal_plant_response"
                result["context_layer"] = str(layer_name)
                within_rows.append(result)

            if str(layer_name) == "analysis_regime":
                for context_a, context_b in itertools.combinations(contexts, 2):
                    result = _between_contexts(
                        data,
                        stratum=stratum,
                        context_a=context_a,
                        context_b=context_b,
                        support_tier=support_tier,
                        threshold=threshold,
                        pattern_config=layer_pattern,
                    )
                    result["axis_set"] = "universal_plant_response"
                    result["context_layer"] = str(layer_name)
                    between_rows.append(result)

    slopes = pd.concat(slopes_parts, ignore_index=True) if slopes_parts else pd.DataFrame()
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    family = ["context_layer", "axis_set", "stratum", "support_tier"]
    if not slopes.empty:
        slopes["q_axis_family"] = slopes.groupby(family, group_keys=False)["p_value"].transform(_bh)
        slopes["axis_supported"] = slopes["q_axis_family"].le(alpha).fillna(False)
    if not within.empty:
        within["q_context_vector_family"] = within.groupby(family, group_keys=False)[
            "p_value"
        ].transform(_bh)
        within["context_vector_supported"] = within["q_context_vector_family"].le(alpha).fillna(False)
    if not between.empty:
        between["q_between_context_family"] = between.groupby(family, group_keys=False)[
            "p_value"
        ].transform(_bh)
        between["context_vector_difference_supported"] = between[
            "q_between_context_family"
        ].le(alpha).fillna(False)
    return slopes, within, between


def _target_rows(
    slopes: pd.DataFrame,
    between: pd.DataFrame,
    *,
    evidence_scope: str,
    median_completeness: float,
    distance_completeness_odds_ratio: float,
    state_recording_odds_ratio: float,
    strata: list[str],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for stratum in strata:
        for layer, context, expected in (
            ("biogeographic_realm", "Palearctic", "positive"),
            ("analysis_regime", "tropical", "negative"),
        ):
            subset = slopes.loc[
                slopes["context_layer"].astype(str).eq(layer)
                & slopes["stratum"].astype(str).eq(stratum)
                & slopes["support_tier"].astype(str).eq("confirmatory")
                & slopes["context"].astype(str).eq(context)
                & slopes["syndrome"].astype(str).eq("accessibility_generalization")
            ]
            if len(subset) != 1:
                rows.append(
                    {
                        "evidence_scope": evidence_scope,
                        "target": f"{context}_accessibility",
                        "stratum": stratum,
                        "status": "not_testable",
                    }
                )
                continue
            row = subset.iloc[0]
            estimate = float(row["distance_slope"])
            supported = bool(row["axis_supported"])
            sign_ok = estimate > 0 if expected == "positive" else estimate < 0
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "target": f"{context}_accessibility",
                    "stratum": stratum,
                    "status": "fit",
                    "expected_direction": expected,
                    "estimate": estimate,
                    "q_value": float(row["q_axis_family"]),
                    "supported": supported,
                    "sign_ok": bool(sign_ok),
                    "n_islands": int(row["n_islands"]),
                }
            )

        subset = between.loc[
            between["context_layer"].astype(str).eq("analysis_regime")
            & between["stratum"].astype(str).eq(stratum)
            & between["support_tier"].astype(str).eq("confirmatory")
            & between["context_a"].astype(str).eq("northern_midlatitude")
            & between["context_b"].astype(str).eq("tropical")
        ]
        if len(subset) != 1:
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "target": "north_tropical_vector_difference",
                    "stratum": stratum,
                    "status": "not_testable",
                }
            )
        else:
            row = subset.iloc[0]
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "target": "north_tropical_vector_difference",
                    "stratum": stratum,
                    "status": str(row.get("status", "fit")),
                    "expected_direction": "direct_difference",
                    "estimate": float("nan"),
                    "q_value": float(row["q_between_context_family"]),
                    "supported": bool(row["context_vector_difference_supported"]),
                    "sign_ok": True,
                    "n_islands": float(row.get("n_unique_islands", np.nan)),
                }
            )

    out = pd.DataFrame(rows)
    out["median_distance_completeness"] = float(median_completeness)
    out["distance_completeness_odds_ratio"] = float(distance_completeness_odds_ratio)
    out["state_recording_odds_ratio"] = float(state_recording_odds_ratio)
    return out


def build_tipping_surface(headline: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    group_cols = [
        "evidence_scope",
        "target",
        "stratum",
        "median_distance_completeness",
        "distance_completeness_odds_ratio",
    ]
    ordered = [
        float(x)
        for x in config["primary_bias_direction"][
            "positive_state_recording_odds_ratio_grid"
        ]
    ]
    for keys, group in headline.groupby(group_cols, dropna=False, sort=True):
        identity = dict(zip(group_cols, keys if isinstance(keys, tuple) else (keys,), strict=True))
        baseline = group.loc[np.isclose(group["state_recording_odds_ratio"], 1.0)]
        if len(baseline) != 1:
            raise ValueError(f"missing unique OR_D=1 baseline for {identity}")
        base = baseline.iloc[0]
        baseline_supported = bool(base.get("supported", False)) and bool(base.get("sign_ok", True))
        if not baseline_supported:
            rows.append(
                {
                    **identity,
                    "baseline_supported": False,
                    "tipping_state_recording_odds_ratio": np.nan,
                    "recording_advantage_nonfocal_over_focal": np.nan,
                    "tipping_event": "baseline_not_supported",
                    "verdict": "not_applicable",
                }
            )
            continue

        event_row: pd.Series | None = None
        event_label = ""
        for or_d in ordered:
            if np.isclose(or_d, 1.0):
                continue
            candidate = group.loc[np.isclose(group["state_recording_odds_ratio"], or_d)]
            if len(candidate) != 1:
                raise ValueError(f"missing OR_D={or_d} scenario for {identity}")
            row = candidate.iloc[0]
            lost = not bool(row.get("supported", False))
            sign_changed = not bool(row.get("sign_ok", True))
            if lost or sign_changed:
                event_row = row
                event_label = "|".join(
                    [name for name, flag in (("support_lost", lost), ("sign_changed", sign_changed)) if flag]
                )
                break
        if event_row is None:
            min_or = min(x for x in ordered if x < 1.0)
            rows.append(
                {
                    **identity,
                    "baseline_supported": True,
                    "tipping_state_recording_odds_ratio": np.nan,
                    "recording_advantage_nonfocal_over_focal": f">{1.0 / min_or:.6g}",
                    "tipping_event": "none_on_frozen_grid",
                    "verdict": "no_break_on_grid",
                }
            )
        else:
            or_d = float(event_row["state_recording_odds_ratio"])
            rows.append(
                {
                    **identity,
                    "baseline_supported": True,
                    "tipping_state_recording_odds_ratio": or_d,
                    "recording_advantage_nonfocal_over_focal": f"{1.0 / or_d:.6g}",
                    "tipping_event": event_label,
                    "verdict": "break_detected",
                }
            )
    return pd.DataFrame(rows)


def _observation_direction_diagnostic(coefficients: pd.DataFrame) -> dict[str, Any]:
    required = {"endpoint", "predictor", "estimate_log_odds", "p_value"}
    if missing := required - set(coefficients.columns):
        return {"status": "unavailable", "reason": f"missing columns {sorted(missing)}"}
    work = coefficients.loc[coefficients["endpoint"].astype(str).eq("flora_recorded")].copy()
    base_name = "z_log_distance_to_continent_km"
    base = work.loc[work["predictor"].astype(str).eq(base_name)]
    tropical = work.loc[
        work["predictor"].astype(str).eq(
            f"{base_name}:context[tropical]"
        )
    ]
    if len(base) != 1 or len(tropical) != 1:
        return {"status": "unavailable", "reason": "distance coefficients not uniquely found"}
    beta = float(base.iloc[0]["estimate_log_odds"])
    beta_tropical = beta + float(tropical.iloc[0]["estimate_log_odds"])
    return {
        "status": "descriptive_only",
        "endpoint": "any_flora_recorded_not_list_completeness",
        "reference_context": "northern_midlatitude",
        "northern_midlatitude_log_odds_per_sd": beta,
        "northern_midlatitude_odds_ratio_per_sd": math.exp(beta),
        "northern_midlatitude_base_p_value": float(base.iloc[0]["p_value"]),
        "tropical_combined_log_odds_per_sd": beta_tropical,
        "tropical_combined_odds_ratio_per_sd": math.exp(beta_tropical),
        "tropical_combined_p_value_available": False,
    }


def run_sensitivity(
    *,
    artifact_root: Path,
    config_path: Path,
    pattern_config_path: Path,
    branching_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    config = load_config(config_path)
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(branching_config_path.read_text(encoding="utf-8"))

    status_flora = pd.read_csv(
        artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz"
    )
    covariates = pd.read_csv(
        artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv"
    )
    realm_assignment = pd.read_csv(
        artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv"
    )
    observation_coefficients = pd.read_csv(
        artifact_root / "fixed/canonical/observation_bias/observation_selection_coefficients.csv"
    )

    geography = str(pattern_config["geography_column"])
    mean_distance, sd_distance = distance_standardization(covariates, geography)
    strata = [str(x) for x in config["strata"]]
    recorded_counts = build_recorded_stratum_counts(status_flora, strata)
    c0_grid = [
        float(x)
        for x in config["primary_bias_direction"]["median_distance_completeness_grid"]
    ]
    or_c_grid = [
        float(x)
        for x in config["primary_bias_direction"][
            "distance_completeness_odds_ratio_grid"
        ]
    ]
    or_d_grid = [
        float(x)
        for x in config["primary_bias_direction"][
            "positive_state_recording_odds_ratio_grid"
        ]
    ]
    affected = str(config["affected_response"]["syndrome"])
    support_tier = str(config["support_tier"])

    headline_parts: list[pd.DataFrame] = []
    diagnostic_parts: list[pd.DataFrame] = []
    maximum_or1_score_difference = 0.0

    for evidence_scope, short_scope in SCOPE_PATHS.items():
        scores = pd.read_csv(
            artifact_root / f"syndrome/{short_scope}/island_syndrome_scores.csv.gz"
        )
        baseline_fit: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame] | None = None
        baseline_scores = scores.copy()
        for c0 in c0_grid:
            for or_c in or_c_grid:
                for or_d in or_d_grid:
                    adjusted, diagnostics = adjust_species_detection(
                        scores,
                        recorded_counts,
                        covariates,
                        geography=geography,
                        distance_mean=mean_distance,
                        distance_sd=sd_distance,
                        median_completeness=c0,
                        distance_completeness_odds_ratio=or_c,
                        state_recording_odds_ratio=or_d,
                        affected_syndrome=affected,
                    )
                    diagnostics.insert(0, "evidence_scope", evidence_scope)
                    diagnostic_parts.append(diagnostics)
                    if np.isclose(or_d, 1.0):
                        original = baseline_scores.loc[
                            baseline_scores["syndrome"].astype(str).eq(affected),
                            ["island_id", "stratum", "syndrome_score"],
                        ].copy()
                        current = adjusted.loc[
                            adjusted["syndrome"].astype(str).eq(affected),
                            ["island_id", "stratum", "syndrome_score"],
                        ].copy()
                        check = original.merge(
                            current,
                            on=["island_id", "stratum"],
                            suffixes=("_observed", "_adjusted"),
                            validate="one_to_one",
                        )
                        maximum_or1_score_difference = max(
                            maximum_or1_score_difference,
                            float(
                                np.max(
                                    np.abs(
                                        check["syndrome_score_observed"].to_numpy(float)
                                        - check["syndrome_score_adjusted"].to_numpy(float)
                                    )
                                )
                            ),
                        )
                        if baseline_fit is None:
                            baseline_fit = fit_headline_family(
                                adjusted,
                                covariates,
                                realm_assignment,
                                pattern_config,
                                branching_config,
                                strata=strata,
                                support_tier=support_tier,
                            )
                        fit = baseline_fit
                    else:
                        fit = fit_headline_family(
                            adjusted,
                            covariates,
                            realm_assignment,
                            pattern_config,
                            branching_config,
                            strata=strata,
                            support_tier=support_tier,
                        )
                    headline_parts.append(
                        _target_rows(
                            fit[0],
                            fit[2],
                            evidence_scope=evidence_scope,
                            median_completeness=c0,
                            distance_completeness_odds_ratio=or_c,
                            state_recording_odds_ratio=or_d,
                            strata=strata,
                        )
                    )

    headline = pd.concat(headline_parts, ignore_index=True)
    diagnostics = pd.concat(diagnostic_parts, ignore_index=True)
    tipping = build_tipping_surface(headline, config)
    observation_direction = _observation_direction_diagnostic(observation_coefficients)

    output_dir.mkdir(parents=True, exist_ok=True)
    headline.to_csv(output_dir / "species_detection_headline_scenarios.csv.gz", index=False)
    tipping.to_csv(output_dir / "species_detection_tipping_surface.csv", index=False)
    diagnostics.to_csv(output_dir / "species_detection_score_diagnostics.csv.gz", index=False)
    (output_dir / "species_detection_observation_direction.json").write_text(
        json.dumps(observation_direction, indent=2) + "\n", encoding="utf-8"
    )

    manifest = {
        "contract": CONTRACT,
        "status": "species_detection_tipping_surface_complete",
        "pinned_workflow_run_id": int(config["pinned_input"]["workflow_run_id"]),
        "pinned_artifact_id": int(config["pinned_input"]["artifact_id"]),
        "pinned_artifact_digest": str(config["pinned_input"]["digest"]),
        "n_evidence_scopes": len(SCOPE_PATHS),
        "n_completeness_levels": len(c0_grid),
        "n_distance_completeness_levels": len(or_c_grid),
        "n_state_recording_levels": len(or_d_grid),
        "n_headline_rows": int(len(headline)),
        "n_tipping_rows": int(len(tipping)),
        "maximum_OR_D_1_score_difference": maximum_or1_score_difference,
        "trait_resolution_conditioning": "OR_R_equals_1_separate_from_V5",
        "hypothetical_species_increase_regression_precision": False,
        "observation_direction_diagnostic": observation_direction,
        "claim_ceiling": str(config["claim_ceiling"]),
    }
    (output_dir / "chapter1_species_detection_tipping_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_species_detection_tipping.yml"), exists=True, dir_okay=False
    ),
    pattern_config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True, dir_okay=False
    ),
    branching_config_path: Path = typer.Option(
        Path("config/chapter1_global_branching.yml"), exists=True, dir_okay=False
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            run_sensitivity(
                artifact_root=artifact_root,
                config_path=config_path,
                pattern_config_path=pattern_config_path,
                branching_config_path=branching_config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
