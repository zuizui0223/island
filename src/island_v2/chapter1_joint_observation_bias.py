"""Prospectively frozen joint observation-bias sensitivity for Chapter 1 P3.

The joint analysis combines two already-frozen observation processes without
reinterpreting either parent analysis:

* V5: trait-state-dependent resolution of ``selfing_core`` among recorded flora.
* V6: distance-dependent list completeness plus state-dependent recording of
  ``generalized_accessible`` species.

The processes act on distinct primary response components. Hypothetical or
imputed species never increase the original regression information weights.
The resulting grid is an assumption surface, not a posterior over missingness.
"""
from __future__ import annotations

import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_species_detection_tipping import (
    adjust_species_detection,
    build_recorded_stratum_counts,
    distance_standardization,
    fit_headline_family,
)
from island_v2.chapter1_trait_resolution_mnar import (
    _context_assignments,
    _scenario_rows,
    adjust_island_scores,
    build_total_species_counts,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
CONTRACT = "chapter1_joint_observation_bias_v1"
SCOPE_PATHS = {
    "all_analysis_eligible": "all",
    "direct_only": "direct",
}


def load_config(path: Path) -> dict[str, Any]:
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("contract") != CONTRACT:
        raise typer.BadParameter("unexpected joint observation-bias contract")
    return config


def _shared_v5_scenario(odds_ratio: float) -> dict[str, Any]:
    value = float(odds_ratio)
    if not math.isfinite(value) or value <= 0:
        raise ValueError("trait-resolution odds ratio must be positive")
    return {
        "scenario_id": f"joint_shared_or_r_{value:.10g}".replace(".", "p"),
        "scenario_family": "shared_state_dependent_resolution",
        "scenario_type": "selection_grid",
        "context_layer": "all_islands",
        "context_a": "",
        "context_b": "",
        "resolution_odds_ratio": value,
        "log_resolution_odds_ratio": math.log(value),
        "bound_mode": "",
    }


def _scope_scores(scores: pd.DataFrame, strata: list[str]) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score", "n_species"}
    if missing := required - set(scores.columns):
        raise ValueError(f"syndrome scores missing columns: {sorted(missing)}")
    scoped = scores.loc[scores["stratum"].astype(str).isin(set(strata))].copy()
    if scoped.empty:
        raise ValueError("no joint-contracted strata remain")
    observed = set(scoped["stratum"].astype(str))
    if observed != set(strata):
        raise ValueError(f"not all joint strata are present: {sorted(observed)}")
    return scoped


def _fit_targets(
    *,
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    branching_config: dict[str, Any],
    evidence_scope: str,
    strata: list[str],
    support_tier: str,
    or_r: float | None,
    c0: float,
    or_c: float,
    or_d: float,
    surface_type: str,
    v5_family: str,
    v5_bound_mode: str,
) -> pd.DataFrame:
    slopes, within, between = fit_headline_family(
        adjusted_scores,
        covariates,
        realm_assignment,
        pattern_config,
        branching_config,
        strata=strata,
        support_tier=support_tier,
    )
    rows: list[dict[str, Any]] = []

    def add_scalar(stratum: str, layer: str, context: str, target: str, expected: str) -> None:
        subset = slopes.loc[
            slopes["context_layer"].astype(str).eq(layer)
            & slopes["stratum"].astype(str).eq(stratum)
            & slopes["support_tier"].astype(str).eq(support_tier)
            & slopes["context"].astype(str).eq(context)
            & slopes["syndrome"].astype(str).eq("accessibility_generalization")
        ]
        if len(subset) != 1:
            rows.append({"target": target, "stratum": stratum, "status": "not_testable"})
            return
        row = subset.iloc[0]
        estimate = float(row["distance_slope"])
        supported = bool(row["axis_supported"])
        sign_ok = estimate > 0 if expected == "positive" else estimate < 0
        rows.append(
            {
                "target": target,
                "stratum": stratum,
                "status": "fit",
                "target_type": "scalar",
                "expected_direction": expected,
                "estimate": estimate,
                "q_value": float(row["q_axis_family"]),
                "supported": supported,
                "sign_ok": bool(sign_ok),
                "robust_cell": bool(supported and sign_ok),
                "n_islands": int(row["n_islands"]),
            }
        )

    def add_vector(stratum: str, layer: str, context: str, target: str) -> None:
        subset = within.loc[
            within["context_layer"].astype(str).eq(layer)
            & within["stratum"].astype(str).eq(stratum)
            & within["support_tier"].astype(str).eq(support_tier)
            & within["context"].astype(str).eq(context)
        ]
        if len(subset) != 1:
            rows.append({"target": target, "stratum": stratum, "status": "not_testable"})
            return
        row = subset.iloc[0]
        supported = bool(row["context_vector_supported"])
        rows.append(
            {
                "target": target,
                "stratum": stratum,
                "status": str(row.get("status", "fit")),
                "target_type": "vector",
                "expected_direction": "joint_support",
                "estimate": np.nan,
                "q_value": float(row["q_context_vector_family"]),
                "supported": supported,
                "sign_ok": True,
                "robust_cell": supported,
                "n_islands": float(row.get("n_unique_islands", np.nan)),
            }
        )

    def add_between(stratum: str) -> None:
        subset = between.loc[
            between["context_layer"].astype(str).eq("analysis_regime")
            & between["stratum"].astype(str).eq(stratum)
            & between["support_tier"].astype(str).eq(support_tier)
            & between["context_a"].astype(str).eq("northern_midlatitude")
            & between["context_b"].astype(str).eq("tropical")
        ]
        if len(subset) != 1:
            rows.append(
                {
                    "target": "north_tropical_vector_difference",
                    "stratum": stratum,
                    "status": "not_testable",
                }
            )
            return
        row = subset.iloc[0]
        supported = bool(row["context_vector_difference_supported"])
        rows.append(
            {
                "target": "north_tropical_vector_difference",
                "stratum": stratum,
                "status": str(row.get("status", "fit")),
                "target_type": "vector",
                "expected_direction": "joint_difference",
                "estimate": np.nan,
                "q_value": float(row["q_between_context_family"]),
                "supported": supported,
                "sign_ok": True,
                "robust_cell": supported,
                "n_islands": float(row.get("n_unique_islands", np.nan)),
            }
        )

    for stratum in strata:
        add_scalar(
            stratum,
            "biogeographic_realm",
            "Palearctic",
            "Palearctic_accessibility",
            "positive",
        )
        add_vector(stratum, "biogeographic_realm", "Palearctic", "Palearctic_vector")
        add_scalar(
            stratum,
            "analysis_regime",
            "tropical",
            "tropical_accessibility",
            "negative",
        )
        add_vector(stratum, "analysis_regime", "tropical", "tropical_vector")
        add_between(stratum)

    out = pd.DataFrame(rows)
    out.insert(0, "evidence_scope", evidence_scope)
    out["surface_type"] = surface_type
    out["trait_resolution_odds_ratio"] = np.nan if or_r is None else float(or_r)
    out["median_distance_completeness"] = float(c0)
    out["distance_completeness_odds_ratio"] = float(or_c)
    out["state_recording_odds_ratio"] = float(or_d)
    out["v5_scenario_family"] = v5_family
    out["v5_bound_mode"] = v5_bound_mode
    if "robust_cell" not in out:
        out["robust_cell"] = False
    out["robust_cell"] = out["robust_cell"].fillna(False).astype(bool)
    return out


def _summarize_primary(surface: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = ["evidence_scope", "target", "stratum", "target_type", "expected_direction"]
    for group_keys, group in surface.groupby(keys, dropna=False, sort=True):
        identity = dict(zip(keys, group_keys if isinstance(group_keys, tuple) else (group_keys,), strict=True))
        fit = group.loc[group["status"].astype(str).eq("fit")].copy()
        robust = fit["robust_cell"].astype(bool) if len(fit) else pd.Series(dtype=bool)
        estimates = pd.to_numeric(fit.get("estimate"), errors="coerce").dropna()
        rows.append(
            {
                **identity,
                "n_declared_cells": int(len(group)),
                "n_fit_cells": int(len(fit)),
                "n_robust_cells": int(robust.sum()) if len(robust) else 0,
                "robust_fraction_of_fit_grid": float(robust.mean()) if len(robust) else np.nan,
                "all_fit_cells_robust": bool(len(fit) > 0 and robust.all()),
                "any_fit_cell_fragile": bool(len(fit) > 0 and (~robust).any()),
                "estimate_min": float(estimates.min()) if len(estimates) else np.nan,
                "estimate_max": float(estimates.max()) if len(estimates) else np.nan,
                "grid_fraction_is_probability": False,
            }
        )
    return pd.DataFrame(rows)


def _summarize_envelope(bounds: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    keys = ["evidence_scope", "target", "stratum", "target_type", "expected_direction"]
    for group_keys, group in bounds.groupby(keys, dropna=False, sort=True):
        identity = dict(zip(keys, group_keys if isinstance(group_keys, tuple) else (group_keys,), strict=True))
        fit = group.loc[group["status"].astype(str).eq("fit")].copy()
        estimates = pd.to_numeric(fit.get("estimate"), errors="coerce").dropna()
        expected = str(identity["expected_direction"])
        if len(estimates):
            if expected == "positive":
                sign_identified = bool((estimates > 0).all())
            elif expected == "negative":
                sign_identified = bool((estimates < 0).all())
            else:
                sign_identified = np.nan
        else:
            sign_identified = np.nan
        support_identified = bool(len(fit) > 0 and fit["supported"].fillna(False).astype(bool).all())
        rows.append(
            {
                **identity,
                "n_bound_corner_cells": int(len(group)),
                "n_fit_cells": int(len(fit)),
                "estimate_lower": float(estimates.min()) if len(estimates) else np.nan,
                "estimate_upper": float(estimates.max()) if len(estimates) else np.nan,
                "expected_sign_identified": sign_identified,
                "support_identified_across_envelope": support_identified,
                "envelope_robust": bool(
                    support_identified and (True if pd.isna(sign_identified) else bool(sign_identified))
                ),
                "interpretation": "assumption_bounded_partial_identification_not_latent_truth",
            }
        )
    return pd.DataFrame(rows)


def run_joint_sensitivity(
    *,
    artifact_root: Path,
    config_path: Path,
    pattern_config_path: Path,
    branching_config_path: Path,
    explanation_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    config = load_config(config_path)
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    branching_config = yaml.safe_load(branching_config_path.read_text(encoding="utf-8"))
    explanation_config = yaml.safe_load(explanation_config_path.read_text(encoding="utf-8"))
    if explanation_config.get("contract") != "chapter1_explanation_gap_validation_v1":
        raise ValueError("unexpected V5 parent contract")

    status_flora = pd.read_csv(artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz")
    covariates = pd.read_csv(artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv")
    realm_assignment = pd.read_csv(
        artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv"
    )
    geography = str(pattern_config["geography_column"])
    mean_distance, sd_distance = distance_standardization(covariates, geography)

    joint_cfg = config["joint_primary_surface"]
    strata = [str(x) for x in joint_cfg["strata"]]
    support_tier = str(joint_cfg["support_tier"])
    recorded_counts = build_recorded_stratum_counts(status_flora, strata)
    total_counts = build_total_species_counts(status_flora, strata)
    assignments = _context_assignments(covariates, realm_assignment)

    v5_cfg = config["parent_contracts"]["V5_trait_resolution"]
    v6_cfg = config["parent_contracts"]["V6_species_detection"]
    or_r_grid = [float(x) for x in v5_cfg["resolution_odds_ratio_grid"]]
    c0_grid = [float(x) for x in v6_cfg["median_distance_completeness_grid"]]
    or_c_grid = [float(x) for x in v6_cfg["distance_completeness_odds_ratio_grid"]]
    or_d_grid = [float(x) for x in v6_cfg["state_recording_odds_ratio_grid"]]

    if len(or_r_grid) * len(c0_grid) * len(or_c_grid) * len(or_d_grid) != int(
        joint_cfg["n_surfaces_per_evidence_scope"]
    ):
        raise ValueError("joint primary surface size differs from frozen contract")

    primary_parts: list[pd.DataFrame] = []
    envelope_parts: list[pd.DataFrame] = []
    max_baseline_score_diff = 0.0

    v5_parent = explanation_config["validations"]["V5_trait_resolution_MNAR_tipping_point"]
    bound_scenarios = [
        scenario
        for scenario in _scenario_rows(v5_parent)
        if str(scenario["scenario_type"]) == "partial_identification_bound"
    ]
    frozen_bound_pairs = {
        (str(item["family"]), str(mode))
        for item in config["partial_identification_envelope"]["trait_bound_scenarios"]
        for mode in item["modes"]
    }
    bound_scenarios = [
        scenario
        for scenario in bound_scenarios
        if (str(scenario["scenario_family"]), str(scenario["bound_mode"])) in frozen_bound_pairs
    ]
    if len(bound_scenarios) != len(frozen_bound_pairs):
        raise ValueError("partial-identification V5 bounds do not match frozen contract")

    corner = config["partial_identification_envelope"]["v6_corner_grid"]
    corner_grid = list(
        itertools.product(
            [float(x) for x in corner["median_distance_completeness"]],
            [float(x) for x in corner["distance_completeness_odds_ratio"]],
            [float(x) for x in corner["state_recording_odds_ratio"]],
        )
    )
    expected_corner = int(
        config["partial_identification_envelope"]["n_corner_bound_surfaces_per_evidence_scope"]
    )
    if len(bound_scenarios) * len(corner_grid) != expected_corner:
        raise ValueError("partial-identification corner size differs from frozen contract")

    for evidence_scope, short_scope in SCOPE_PATHS.items():
        scores = _scope_scores(
            pd.read_csv(artifact_root / f"syndrome/{short_scope}/island_syndrome_scores.csv.gz"),
            strata,
        )
        base_index = scores.set_index(["island_id", "stratum", "syndrome"])["syndrome_score"]

        for or_r in or_r_grid:
            after_v5, _ = adjust_island_scores(
                scores,
                total_counts,
                assignments,
                _shared_v5_scenario(or_r),
                {str(v5_cfg["affected_syndrome"])},
            )
            for c0, or_c, or_d in itertools.product(c0_grid, or_c_grid, or_d_grid):
                adjusted, _ = adjust_species_detection(
                    after_v5,
                    recorded_counts,
                    covariates,
                    geography=geography,
                    distance_mean=mean_distance,
                    distance_sd=sd_distance,
                    median_completeness=c0,
                    distance_completeness_odds_ratio=or_c,
                    state_recording_odds_ratio=or_d,
                    affected_syndrome=str(v6_cfg["affected_syndrome"]),
                )
                if np.isclose(or_r, 1.0) and np.isclose(or_d, 1.0):
                    current = adjusted.set_index(["island_id", "stratum", "syndrome"])["syndrome_score"]
                    max_baseline_score_diff = max(
                        max_baseline_score_diff,
                        float((current - base_index).abs().max()),
                    )
                primary_parts.append(
                    _fit_targets(
                        adjusted_scores=adjusted,
                        covariates=covariates,
                        realm_assignment=realm_assignment,
                        pattern_config=pattern_config,
                        branching_config=branching_config,
                        evidence_scope=evidence_scope,
                        strata=strata,
                        support_tier=support_tier,
                        or_r=or_r,
                        c0=c0,
                        or_c=or_c,
                        or_d=or_d,
                        surface_type="joint_selection_grid",
                        v5_family="shared_state_dependent_resolution",
                        v5_bound_mode="",
                    )
                )

        for scenario in bound_scenarios:
            after_bound, _ = adjust_island_scores(
                scores,
                total_counts,
                assignments,
                scenario,
                {str(v5_cfg["affected_syndrome"])},
            )
            for c0, or_c, or_d in corner_grid:
                adjusted, _ = adjust_species_detection(
                    after_bound,
                    recorded_counts,
                    covariates,
                    geography=geography,
                    distance_mean=mean_distance,
                    distance_sd=sd_distance,
                    median_completeness=c0,
                    distance_completeness_odds_ratio=or_c,
                    state_recording_odds_ratio=or_d,
                    affected_syndrome=str(v6_cfg["affected_syndrome"]),
                )
                envelope_parts.append(
                    _fit_targets(
                        adjusted_scores=adjusted,
                        covariates=covariates,
                        realm_assignment=realm_assignment,
                        pattern_config=pattern_config,
                        branching_config=branching_config,
                        evidence_scope=evidence_scope,
                        strata=strata,
                        support_tier=support_tier,
                        or_r=None,
                        c0=c0,
                        or_c=or_c,
                        or_d=or_d,
                        surface_type="partial_identification_corner",
                        v5_family=str(scenario["scenario_family"]),
                        v5_bound_mode=str(scenario["bound_mode"]),
                    )
                )

    primary = pd.concat(primary_parts, ignore_index=True)
    envelope_cells = pd.concat(envelope_parts, ignore_index=True)
    robustness = _summarize_primary(primary)
    envelope = _summarize_envelope(envelope_cells)

    if not math.isfinite(max_baseline_score_diff) or max_baseline_score_diff > 1e-12:
        raise ValueError(f"joint no-differential baseline failed score reproduction: {max_baseline_score_diff}")

    output_dir.mkdir(parents=True, exist_ok=True)
    primary.to_csv(output_dir / "joint_surface_classification.csv.gz", index=False, compression="gzip")
    robustness.to_csv(output_dir / "joint_robustness_summary.csv", index=False)
    envelope_cells.to_csv(
        output_dir / "partial_identification_corner_cells.csv.gz", index=False, compression="gzip"
    )
    envelope.to_csv(output_dir / "partial_identification_envelope.csv", index=False)

    manifest = {
        "contract": CONTRACT,
        "status": "joint_surface_and_partial_identification_complete",
        "pinned_workflow_run_id": int(config["pinned_input"]["workflow_run_id"]),
        "pinned_artifact_id": int(config["pinned_input"]["artifact_id"]),
        "pinned_artifact_digest": str(config["pinned_input"]["digest"]),
        "n_primary_parameter_surfaces_per_scope": int(joint_cfg["n_surfaces_per_evidence_scope"]),
        "n_primary_target_rows": int(len(primary)),
        "n_partial_identification_surfaces_per_scope": expected_corner,
        "n_partial_identification_target_rows": int(len(envelope_cells)),
        "maximum_no_differential_score_difference": max_baseline_score_diff,
        "hypothetical_or_imputed_species_increase_regression_precision": False,
        "grid_fraction_is_probability": False,
        "pollinator_mechanism_promoted": False,
        "claim_ceiling": str(config["claim_ceiling"]),
    }
    (output_dir / "joint_observation_bias_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_joint_observation_bias.yml"), exists=True, dir_okay=False
    ),
    pattern_config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True, dir_okay=False
    ),
    branching_config_path: Path = typer.Option(
        Path("config/chapter1_global_branching.yml"), exists=True, dir_okay=False
    ),
    explanation_config_path: Path = typer.Option(
        Path("config/chapter1_explanation_gap_validation.yml"), exists=True, dir_okay=False
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            run_joint_sensitivity(
                artifact_root=artifact_root,
                config_path=config_path,
                pattern_config_path=pattern_config_path,
                branching_config_path=branching_config_path,
                explanation_config_path=explanation_config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
