"""Observed Chapter 1 response geometry after prospective V1/V2 audits.

This module opens the observed nonlinear response only after reading the V2 frozen
critical-D and cross-evidence shape-qualification tables. Geometry remains a plant-
response description; no mechanism is inferred from the selected shape.
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

from island_v2.chapter1_pr136_biogeographic_residual import (
    _fit_weighted_clustered_design,
)
from island_v2.chapter1_response_geometry_calibrated import (
    _score_shape_family,
)
from island_v2.chapter1_response_geometry_power import (
    _standardize,
    build_candidates,
    prepare_cell,
)


def _empirical_logit(frame: pd.DataFrame) -> np.ndarray:
    success = frame["successes"].to_numpy(float)
    trials = frame["trials"].to_numpy(float)
    failure = trials - success
    return np.log((success + 0.5) / (failure + 0.5))


def _base_design(
    frame: pd.DataFrame, config: dict[str, Any]
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    z = _standardize(frame[str(config["primary_exposure"])])
    columns = [np.ones(len(frame), dtype=float)]
    names = ["intercept"]
    for column in config["baseline_covariates"]:
        columns.append(_standardize(frame[str(column)]))
        names.append(f"z_{column}")
    return z, np.column_stack(columns), names


def _robust_effect(
    frame: pd.DataFrame,
    *,
    selected_shape: str,
    breakpoint_quantile: float,
    config: dict[str, Any],
) -> dict[str, Any]:
    y = _empirical_logit(frame)
    weights = frame["trials"].to_numpy(float)
    clusters = frame[str(config["cluster_column"])].astype(str).to_numpy()
    z, base, base_names = _base_design(frame, config)
    out: dict[str, Any] = {
        "effect_label": "not_applicable",
        "effect_estimate": np.nan,
        "effect_cluster_robust_se": np.nan,
        "effect_ci95_low": np.nan,
        "effect_ci95_high": np.nan,
        "turning_point_z": np.nan,
        "turning_point_z_ci95_low": np.nan,
        "turning_point_z_ci95_high": np.nan,
        "turning_point_quantile_ci95_low": np.nan,
        "turning_point_quantile_ci95_high": np.nan,
        "robust_fit_status": "not_requested",
    }
    if selected_shape == "G2_step":
        knot = float(np.quantile(z, breakpoint_quantile))
        step = (z > knot).astype(float)
        X = np.column_stack([base, step])
        names = [*base_names, "step_after_breakpoint"]
        coef, covariance, fit = _fit_weighted_clustered_design(
            y, weights, X, names, clusters
        )
        out["robust_fit_status"] = str(fit.get("status", "unknown"))
        if coef.empty:
            return out
        row = coef.set_index("predictor").loc["step_after_breakpoint"]
        estimate = float(row["estimate"])
        se = float(row["cluster_robust_se"])
        out.update(
            {
                "effect_label": "adjusted_step_empirical_logit",
                "effect_estimate": estimate,
                "effect_cluster_robust_se": se,
                "effect_ci95_low": estimate - 1.96 * se,
                "effect_ci95_high": estimate + 1.96 * se,
            }
        )
        return out

    if selected_shape == "G4_reversal":
        X = np.column_stack([base, z, z * z])
        names = [*base_names, "z_exposure", "z_exposure_squared"]
        coef, covariance, fit = _fit_weighted_clustered_design(
            y, weights, X, names, clusters
        )
        out["robust_fit_status"] = str(fit.get("status", "unknown"))
        if coef.empty:
            return out
        index = coef.set_index("predictor")
        b1 = float(index.loc["z_exposure", "estimate"])
        b2 = float(index.loc["z_exposure_squared", "estimate"])
        se2 = float(index.loc["z_exposure_squared", "cluster_robust_se"])
        out.update(
            {
                "effect_label": "adjusted_quadratic_curvature",
                "effect_estimate": b2,
                "effect_cluster_robust_se": se2,
                "effect_ci95_low": b2 - 1.96 * se2,
                "effect_ci95_high": b2 + 1.96 * se2,
            }
        )
        if abs(b2) <= 1e-12:
            return out
        turning = -b1 / (2.0 * b2)
        p1 = names.index("z_exposure")
        p2 = names.index("z_exposure_squared")
        sub = covariance[np.ix_([p1, p2], [p1, p2])]
        gradient = np.array(
            [-1.0 / (2.0 * b2), b1 / (2.0 * b2 * b2)], dtype=float
        )
        variance = float(gradient @ sub @ gradient)
        turn_se = math.sqrt(max(variance, 0.0))
        low = turning - 1.96 * turn_se
        high = turning + 1.96 * turn_se
        q_low = float(np.mean(z <= low))
        q_high = float(np.mean(z <= high))
        out.update(
            {
                "turning_point_z": turning,
                "turning_point_z_ci95_low": low,
                "turning_point_z_ci95_high": high,
                "turning_point_quantile_ci95_low": min(q_low, q_high),
                "turning_point_quantile_ci95_high": max(q_low, q_high),
            }
        )
    return out


def fit_scope_cell(
    frame: pd.DataFrame,
    *,
    evidence_scope: str,
    measurement_axis: str,
    outcome: str,
    context: str,
    stratum: str,
    critical_D: float,
    observed_config: dict[str, Any],
    v2_config: dict[str, Any],
) -> dict[str, Any]:
    z, weights, candidates = build_candidates(frame, v2_config)
    y = _empirical_logit(frame)[:, None]
    scored = _score_shape_family(y, z, weights, candidates, v2_config)
    D = float(scored["D"][0])
    best_monotonic = str(scored["best_monotonic"][0])
    best_nonlinear = str(scored["best_nonlinear"][0])
    breakpoint_q = float(scored["best_nonlinear_breakpoint"][0])
    exposure = str(observed_config["primary_exposure"])
    breakpoint_value = (
        float(np.quantile(frame[exposure].to_numpy(float), breakpoint_q))
        if np.isfinite(breakpoint_q)
        else np.nan
    )
    nonlinear_pass = bool(D > float(critical_D))
    effects = _robust_effect(
        frame,
        selected_shape=best_nonlinear,
        breakpoint_quantile=breakpoint_q,
        config=observed_config,
    )
    scores = {f"AICc_{shape}": float(scored["scores"][shape][0]) for shape in scored["scores"]}
    return {
        "evidence_scope": evidence_scope,
        "measurement_axis": measurement_axis,
        "outcome": outcome,
        "context": context,
        "stratum": stratum,
        "n_islands": int(len(frame)),
        "n_spatial_blocks": int(
            frame[str(observed_config["cluster_column"])].nunique()
        ),
        "critical_D": float(critical_D),
        "observed_D": D,
        "nonlinear_gate_pass": nonlinear_pass,
        "best_monotonic_shape": best_monotonic,
        "best_nonlinear_shape": best_nonlinear,
        "breakpoint_quantile": breakpoint_q,
        "breakpoint_exposure_value": breakpoint_value,
        **scores,
        **effects,
    }


def classify_cross_scope(
    scope: pd.DataFrame,
    cross_qualification: pd.DataFrame,
    observed_config: dict[str, Any],
) -> pd.DataFrame:
    keys = ["measurement_axis", "outcome", "context", "stratum"]
    tolerance = float(
        observed_config["uncertainty_reporting"]["breakpoint_consensus"][
            "maximum_scope_difference_quantile"
        ]
    )
    rows = []
    for values, part in scope.groupby(keys, sort=True):
        by_scope = part.set_index("evidence_scope")
        row = dict(zip(keys, values, strict=True))
        required = {"all_analysis_eligible", "direct_only"}
        if not required.issubset(set(by_scope.index.astype(str))):
            rows.append({**row, "classification": "monotonic_or_unresolved"})
            continue
        all_row = by_scope.loc["all_analysis_eligible"]
        direct_row = by_scope.loc["direct_only"]
        row.update(
            {
                "all_observed_D": float(all_row["observed_D"]),
                "all_critical_D": float(all_row["critical_D"]),
                "direct_observed_D": float(direct_row["observed_D"]),
                "direct_critical_D": float(direct_row["critical_D"]),
                "all_nonlinear_gate_pass": bool(all_row["nonlinear_gate_pass"]),
                "direct_nonlinear_gate_pass": bool(direct_row["nonlinear_gate_pass"]),
                "all_best_nonlinear_shape": str(all_row["best_nonlinear_shape"]),
                "direct_best_nonlinear_shape": str(direct_row["best_nonlinear_shape"]),
                "all_breakpoint_quantile": float(all_row["breakpoint_quantile"]),
                "direct_breakpoint_quantile": float(direct_row["breakpoint_quantile"]),
                "all_breakpoint_exposure_value": float(all_row["breakpoint_exposure_value"]),
                "direct_breakpoint_exposure_value": float(direct_row["breakpoint_exposure_value"]),
                "all_effect_estimate": float(all_row["effect_estimate"]),
                "all_effect_ci95_low": float(all_row["effect_ci95_low"]),
                "all_effect_ci95_high": float(all_row["effect_ci95_high"]),
                "direct_effect_estimate": float(direct_row["effect_estimate"]),
                "direct_effect_ci95_low": float(direct_row["effect_ci95_low"]),
                "direct_effect_ci95_high": float(direct_row["effect_ci95_high"]),
            }
        )
        if not (
            bool(all_row["nonlinear_gate_pass"])
            and bool(direct_row["nonlinear_gate_pass"])
        ):
            row["classification"] = "monotonic_or_unresolved"
            rows.append(row)
            continue
        all_shape = str(all_row["best_nonlinear_shape"])
        direct_shape = str(direct_row["best_nonlinear_shape"])
        if all_shape != direct_shape:
            row["classification"] = "nonlinear_shape_unresolved"
            rows.append(row)
            continue
        shape = all_shape
        qualified = cross_qualification.loc[
            cross_qualification["measurement_axis"].astype(str).eq(str(row["measurement_axis"]))
            & cross_qualification["outcome"].astype(str).eq(str(row["outcome"]))
            & cross_qualification["context"].astype(str).eq(str(row["context"]))
            & cross_qualification["stratum"].astype(str).eq(str(row["stratum"]))
            & cross_qualification["shape"].astype(str).eq(shape)
        ]
        shape_qualified = bool(
            len(qualified) == 1
            and bool(qualified.iloc[0]["headline_shape_qualified"])
        )
        row["matched_shape"] = shape
        row["shape_v2_qualified"] = shape_qualified
        if not shape_qualified or shape == "G3_hinge":
            row["classification"] = "nonlinear_shape_unresolved"
            rows.append(row)
            continue
        q_all = float(all_row["breakpoint_quantile"])
        q_direct = float(direct_row["breakpoint_quantile"])
        q_difference = abs(q_all - q_direct)
        row["breakpoint_scope_difference_quantile"] = q_difference
        if not np.isfinite(q_difference) or q_difference > tolerance:
            row["classification"] = "nonlinear_shape_unresolved"
            rows.append(row)
            continue
        row["consensus_breakpoint_quantile"] = (q_all + q_direct) / 2.0
        row["consensus_breakpoint_exposure_value"] = (
            float(all_row["breakpoint_exposure_value"])
            + float(direct_row["breakpoint_exposure_value"])
        ) / 2.0
        if shape == "G2_step":
            row["classification"] = "identified_step"
        elif shape == "G4_reversal":
            row["classification"] = "identified_reversal"
            row["all_turning_point_quantile_ci95_low"] = float(
                all_row["turning_point_quantile_ci95_low"]
            )
            row["all_turning_point_quantile_ci95_high"] = float(
                all_row["turning_point_quantile_ci95_high"]
            )
            row["direct_turning_point_quantile_ci95_low"] = float(
                direct_row["turning_point_quantile_ci95_low"]
            )
            row["direct_turning_point_quantile_ci95_high"] = float(
                direct_row["turning_point_quantile_ci95_high"]
            )
        else:
            row["classification"] = "nonlinear_shape_unresolved"
        rows.append(row)
    return pd.DataFrame(rows)


def run_observed(
    *,
    all_counts: Path,
    direct_counts: Path,
    covariates_csv: Path,
    v2_scope_qualification_csv: Path,
    v2_cross_qualification_csv: Path,
    observed_config_path: Path,
    v2_config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    observed_config = yaml.safe_load(
        observed_config_path.read_text(encoding="utf-8")
    )
    v2_config = yaml.safe_load(v2_config_path.read_text(encoding="utf-8"))
    scopes = {
        "all_analysis_eligible": pd.read_csv(all_counts),
        "direct_only": pd.read_csv(direct_counts),
    }
    covariates = pd.read_csv(covariates_csv)
    v2_scope = pd.read_csv(v2_scope_qualification_csv)
    v2_cross = pd.read_csv(v2_cross_qualification_csv)
    response_map = {
        str(outcome): str(axis)
        for axis, outcome in observed_config["responses"].items()
    }
    scope_rows = []
    for evidence_scope, counts in scopes.items():
        for outcome, axis in response_map.items():
            for context in observed_config["contexts"]:
                for stratum in observed_config["strata"]:
                    frame = prepare_cell(
                        counts,
                        covariates,
                        outcome=outcome,
                        context=str(context),
                        stratum=str(stratum),
                        config=v2_config,
                    )
                    critical_rows = v2_scope.loc[
                        v2_scope["evidence_scope"].astype(str).eq(evidence_scope)
                        & v2_scope["measurement_axis"].astype(str).eq(axis)
                        & v2_scope["outcome"].astype(str).eq(outcome)
                        & v2_scope["context"].astype(str).eq(str(context))
                        & v2_scope["stratum"].astype(str).eq(str(stratum))
                    ]
                    critical_values = critical_rows["critical_D"].dropna().unique()
                    if len(critical_values) != 1:
                        raise ValueError(
                            "V2 critical_D is missing or inconsistent for "
                            f"{evidence_scope}/{axis}/{context}/{stratum}"
                        )
                    scope_rows.append(
                        fit_scope_cell(
                            frame,
                            evidence_scope=evidence_scope,
                            measurement_axis=axis,
                            outcome=outcome,
                            context=str(context),
                            stratum=str(stratum),
                            critical_D=float(critical_values[0]),
                            observed_config=observed_config,
                            v2_config=v2_config,
                        )
                    )
    scope = pd.DataFrame(scope_rows)
    cross = classify_cross_scope(scope, v2_cross, observed_config)
    output_dir.mkdir(parents=True, exist_ok=True)
    scope.to_csv(output_dir / "observed_geometry_scope.csv", index=False)
    cross.to_csv(output_dir / "observed_geometry_cross_scope.csv", index=False)
    counts = cross["classification"].value_counts().to_dict()
    manifest = {
        "contract": observed_config["contract"],
        "status": "observed_geometry_opened_under_frozen_v2_gate",
        "n_scope_cells": int(len(scope)),
        "n_cross_scope_cells": int(len(cross)),
        "classification_counts": {str(k): int(v) for k, v in counts.items()},
        "observed_geometry_fitted": True,
        "mechanism_inferred": False,
        "physical_water_crossing_threshold_claimed": False,
        "claim_boundary": observed_config["claim_ceiling"].strip(),
    }
    (output_dir / "observed_geometry_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-counts", type=Path, required=True)
    parser.add_argument("--direct-counts", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--v2-scope-qualification-csv", type=Path, required=True)
    parser.add_argument("--v2-cross-qualification-csv", type=Path, required=True)
    parser.add_argument("--observed-config-path", type=Path, required=True)
    parser.add_argument("--v2-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = run_observed(
        all_counts=args.all_counts,
        direct_counts=args.direct_counts,
        covariates_csv=args.covariates_csv,
        v2_scope_qualification_csv=args.v2_scope_qualification_csv,
        v2_cross_qualification_csv=args.v2_cross_qualification_csv,
        observed_config_path=args.observed_config_path,
        v2_config_path=args.v2_config_path,
        output_dir=args.output_dir,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
