"""V2 held-out calibration for Chapter 1 nonlinear response geometry.

V1 showed that naive AICc selection can mistake spatial-block variation for a
breakpoint. V2 keeps the same candidate geometries but calibrates a cell-specific
nonlinearity statistic under simulated monotonic truths, then evaluates false
selection and nonlinear-shape recovery on independent simulations.

No observed nonlinear response is fit here.
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

from island_v2.chapter1_response_geometry_power import (
    GEOMETRIES,
    _aicc_vector,
    _shape_signal,
    _stable_seed,
    build_candidates,
    prepare_cell,
)

MONOTONIC = ("G0_flat", "G1_cline")
NONLINEAR = ("G2_step", "G3_hinge", "G4_reversal")


def _score_shape_family(
    y: np.ndarray,
    z: np.ndarray,
    weights: np.ndarray,
    candidates: list[Any],
    config: dict[str, Any],
) -> dict[str, Any]:
    """Return minimum AICc and breakpoint per shape for every replicate."""
    sqrtw = np.sqrt(weights)[:, None]
    yw = y * sqrtw
    total_ss = np.sum(yw * yw, axis=0)
    n, replicates = y.shape
    per_shape: dict[str, list[tuple[np.ndarray, np.ndarray]]] = {
        shape: [] for shape in GEOMETRIES
    }
    lo_q, hi_q = config["candidate_geometries"]["G4_reversal"][
        "turning_point_required_between_quantiles"
    ]
    z_lo, z_hi = np.quantile(z, [float(lo_q), float(hi_q)])
    z10, z90 = np.quantile(z, [0.10, 0.90])

    for candidate in candidates:
        projected = candidate.q_basis.T @ yw
        sse = (
            total_ss - np.sum(projected * projected, axis=0)
        ) / float(np.mean(weights))
        score = _aicc_vector(sse, n, candidate.n_parameters)
        breakpoint = np.full(replicates, np.nan)
        if candidate.shape == "G4_reversal":
            beta = candidate.beta_map @ y
            b1 = beta[-2]
            b2 = beta[-1]
            with np.errstate(divide="ignore", invalid="ignore"):
                turning = -b1 / (2.0 * b2)
            d_lo = b1 + 2.0 * b2 * z10
            d_hi = b1 + 2.0 * b2 * z90
            valid = (
                np.isfinite(turning)
                & (np.abs(b2) > 1e-12)
                & (turning >= z_lo)
                & (turning <= z_hi)
                & (d_lo * d_hi < 0)
            )
            score = np.where(valid, score, np.inf)
            if valid.any():
                breakpoint[valid] = np.mean(
                    z[:, None] <= turning[None, valid], axis=0
                )
        elif candidate.breakpoint_quantile is not None:
            breakpoint[:] = float(candidate.breakpoint_quantile)
        per_shape[candidate.shape].append((score, breakpoint))

    scores: dict[str, np.ndarray] = {}
    breakpoints: dict[str, np.ndarray] = {}
    for shape in GEOMETRIES:
        options = per_shape[shape]
        option_scores = np.vstack([item[0] for item in options])
        index = np.argmin(option_scores, axis=0)
        columns = np.arange(replicates)
        scores[shape] = option_scores[index, columns]
        option_breakpoints = np.vstack([item[1] for item in options])
        breakpoints[shape] = option_breakpoints[index, columns]

    monotonic_matrix = np.vstack([scores[name] for name in MONOTONIC])
    nonlinear_matrix = np.vstack([scores[name] for name in NONLINEAR])
    best_monotonic_index = np.argmin(monotonic_matrix, axis=0)
    best_nonlinear_index = np.argmin(nonlinear_matrix, axis=0)
    columns = np.arange(replicates)
    best_monotonic_score = monotonic_matrix[best_monotonic_index, columns]
    best_nonlinear_score = nonlinear_matrix[best_nonlinear_index, columns]
    best_monotonic = np.asarray(MONOTONIC, dtype=object)[best_monotonic_index]
    best_nonlinear = np.asarray(NONLINEAR, dtype=object)[best_nonlinear_index]
    nonlinear_breakpoint_matrix = np.vstack(
        [breakpoints[name] for name in NONLINEAR]
    )
    best_nonlinear_breakpoint = nonlinear_breakpoint_matrix[
        best_nonlinear_index, columns
    ]
    return {
        "scores": scores,
        "breakpoints": breakpoints,
        "best_monotonic": best_monotonic,
        "best_nonlinear": best_nonlinear,
        "best_nonlinear_breakpoint": best_nonlinear_breakpoint,
        "D": best_monotonic_score - best_nonlinear_score,
    }


def _simulate_empirical_logits(
    frame: pd.DataFrame,
    z: np.ndarray,
    *,
    shape: str,
    effect_logit_sd: float,
    breakpoint_quantile: float,
    cluster_sd: float,
    replicates: int,
    seed: int,
    config: dict[str, Any],
) -> np.ndarray:
    trials = frame["trials"].astype(int).to_numpy()
    pooled = float(
        (frame["successes"].sum() + 0.5) / (frame["trials"].sum() + 1.0)
    )
    pooled = min(max(pooled, 1e-4), 1 - 1e-4)
    intercept = math.log(pooled / (1.0 - pooled))
    blocks = frame[config["cluster_column"]].astype(str).to_numpy()
    unique_blocks, block_index = np.unique(blocks, return_inverse=True)
    rng = np.random.default_rng(seed)
    signal = _shape_signal(z, shape, breakpoint_quantile)
    signs = np.where(np.arange(replicates) % 2 == 0, 1.0, -1.0)
    random_effect = (
        rng.normal(0.0, cluster_sd, size=(len(unique_blocks), replicates))
        if cluster_sd > 0
        else np.zeros((len(unique_blocks), replicates))
    )
    eta = (
        intercept
        + signal[:, None] * float(effect_logit_sd) * signs[None, :]
        + random_effect[block_index]
    )
    probability = 1.0 / (1.0 + np.exp(-np.clip(eta, -30, 30)))
    successes = rng.binomial(trials[:, None], probability)
    failures = trials[:, None] - successes
    return np.log((successes + 0.5) / (failures + 0.5))


def _monotonic_scenarios(spec: dict[str, Any]) -> list[tuple[str, float]]:
    rows: list[tuple[str, float]] = []
    for shape, details in spec["monotonic_truths"].items():
        for effect in details["effect_logit_sd"]:
            rows.append((str(shape), float(effect)))
    return rows


def calibrate_cell(
    frame: pd.DataFrame,
    *,
    evidence_scope: str,
    measurement_axis: str,
    outcome: str,
    context: str,
    stratum: str,
    config: dict[str, Any],
) -> tuple[float, pd.DataFrame, np.ndarray, np.ndarray, list[Any]]:
    z, weights, candidates = build_candidates(frame, config)
    calibration = config["calibration"]
    replicates = int(calibration["replicates_per_scenario"])
    cluster_sd = float(calibration["cluster_random_intercept_sd"])
    quantile = float(calibration["nonlinear_gate_quantile"])
    rows = []
    critical_candidates = []
    for shape, effect in _monotonic_scenarios(calibration):
        y = _simulate_empirical_logits(
            frame,
            z,
            shape=shape,
            effect_logit_sd=effect,
            breakpoint_quantile=0.50,
            cluster_sd=cluster_sd,
            replicates=replicates,
            seed=_stable_seed(
                int(calibration["seed"]),
                "calibration",
                evidence_scope,
                measurement_axis,
                outcome,
                context,
                stratum,
                shape,
                str(effect),
            ),
            config=config,
        )
        scored = _score_shape_family(y, z, weights, candidates, config)
        critical = float(np.quantile(scored["D"], quantile))
        critical_candidates.append(critical)
        rows.append(
            {
                "evidence_scope": evidence_scope,
                "measurement_axis": measurement_axis,
                "outcome": outcome,
                "context": context,
                "stratum": stratum,
                "truth_shape": shape,
                "effect_logit_sd": effect,
                "cluster_random_intercept_sd": cluster_sd,
                "calibration_quantile": quantile,
                "scenario_critical_D": critical,
                "median_D": float(np.median(scored["D"])),
                "p95_D": float(np.quantile(scored["D"], 0.95)),
            }
        )
    critical_D = max(critical_candidates)
    table = pd.DataFrame(rows)
    table["cell_critical_D"] = critical_D
    return critical_D, table, z, weights, candidates


def validate_cell(
    frame: pd.DataFrame,
    *,
    evidence_scope: str,
    measurement_axis: str,
    outcome: str,
    context: str,
    stratum: str,
    critical_D: float,
    z: np.ndarray,
    weights: np.ndarray,
    candidates: list[Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    validation = config["validation"]
    replicates = int(validation["replicates_per_scenario"])
    cluster_sd = float(validation["cluster_random_intercept_sd"])
    seed0 = int(validation["seed"])
    rows: list[dict[str, Any]] = []

    for shape, effect in _monotonic_scenarios(validation):
        y = _simulate_empirical_logits(
            frame,
            z,
            shape=shape,
            effect_logit_sd=effect,
            breakpoint_quantile=0.50,
            cluster_sd=cluster_sd,
            replicates=replicates,
            seed=_stable_seed(
                seed0,
                "validation_monotonic",
                evidence_scope,
                measurement_axis,
                outcome,
                context,
                stratum,
                shape,
                str(effect),
            ),
            config=config,
        )
        scored = _score_shape_family(y, z, weights, candidates, config)
        nonlinear = scored["D"] > critical_D
        rows.append(
            {
                "evidence_scope": evidence_scope,
                "measurement_axis": measurement_axis,
                "outcome": outcome,
                "context": context,
                "stratum": stratum,
                "truth_family": "monotonic",
                "truth_shape": shape,
                "effect_logit_sd": effect,
                "truth_breakpoint_quantile": np.nan,
                "critical_D": critical_D,
                "false_nonlinear_rate": float(np.mean(nonlinear)),
                "shape_recovery_rate": float(
                    np.mean((~nonlinear) & (scored["best_monotonic"] == shape))
                ),
                "breakpoint_quantile_mae": np.nan,
            }
        )

    effect = float(validation["nonlinear_target_effect_logit_sd"])
    truth_breakpoints = [
        float(x) for x in validation["nonlinear_truth_breakpoint_quantiles"]
    ]
    for shape in NONLINEAR:
        shape_breakpoints = truth_breakpoints if shape != "G4_reversal" else [0.50]
        for breakpoint_quantile in shape_breakpoints:
            y = _simulate_empirical_logits(
                frame,
                z,
                shape=shape,
                effect_logit_sd=effect,
                breakpoint_quantile=breakpoint_quantile,
                cluster_sd=cluster_sd,
                replicates=replicates,
                seed=_stable_seed(
                    seed0,
                    "validation_nonlinear",
                    evidence_scope,
                    measurement_axis,
                    outcome,
                    context,
                    stratum,
                    shape,
                    str(breakpoint_quantile),
                ),
                config=config,
            )
            scored = _score_shape_family(y, z, weights, candidates, config)
            gate = scored["D"] > critical_D
            recovered = gate & (scored["best_nonlinear"] == shape)
            errors = np.abs(
                scored["best_nonlinear_breakpoint"][recovered]
                - breakpoint_quantile
            )
            rows.append(
                {
                    "evidence_scope": evidence_scope,
                    "measurement_axis": measurement_axis,
                    "outcome": outcome,
                    "context": context,
                    "stratum": stratum,
                    "truth_family": "nonlinear",
                    "truth_shape": shape,
                    "effect_logit_sd": effect,
                    "truth_breakpoint_quantile": breakpoint_quantile,
                    "critical_D": critical_D,
                    "false_nonlinear_rate": np.nan,
                    "shape_recovery_rate": float(np.mean(recovered)),
                    "breakpoint_quantile_mae": (
                        float(np.mean(errors)) if len(errors) else np.nan
                    ),
                }
            )

    detail = pd.DataFrame(rows)
    false_ceiling = float(
        config["qualification"]["maximum_validation_false_nonlinear_rate"]
    )
    recovery_floor = float(
        config["qualification"]["minimum_shape_recovery_rate"]
    )
    breakpoint_ceiling = float(
        config["qualification"]["maximum_breakpoint_quantile_mae"]
    )
    monotonic = detail.loc[detail["truth_family"].eq("monotonic")]
    max_false = float(monotonic["false_nonlinear_rate"].max())
    qualification_rows = []
    for shape in NONLINEAR:
        subset = detail.loc[detail["truth_shape"].eq(shape)].copy()
        min_recovery = float(subset["shape_recovery_rate"].min())
        max_mae = float(subset["breakpoint_quantile_mae"].max())
        passed = (
            max_false <= false_ceiling
            and min_recovery >= recovery_floor
            and np.isfinite(max_mae)
            and max_mae <= breakpoint_ceiling
        )
        qualification_rows.append(
            {
                "evidence_scope": evidence_scope,
                "measurement_axis": measurement_axis,
                "outcome": outcome,
                "context": context,
                "stratum": stratum,
                "n_islands": int(len(frame)),
                "n_spatial_blocks": int(
                    frame[config["cluster_column"]].nunique()
                ),
                "critical_D": critical_D,
                "shape": shape,
                "maximum_validation_false_nonlinear_rate": max_false,
                "minimum_shape_recovery_rate": min_recovery,
                "maximum_breakpoint_quantile_mae": max_mae,
                "scope_shape_qualified": bool(passed),
            }
        )
    return detail, pd.DataFrame(qualification_rows)


def cross_scope_qualification(scope: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "measurement_axis",
        "outcome",
        "context",
        "stratum",
        "shape",
    ]
    rows = []
    for values, part in scope.groupby(keys, sort=True):
        scopes = set(part["evidence_scope"].astype(str))
        all_ok = bool(
            part.loc[
                part["evidence_scope"].eq("all_analysis_eligible"),
                "scope_shape_qualified",
            ].all()
        ) if "all_analysis_eligible" in scopes else False
        direct_ok = bool(
            part.loc[
                part["evidence_scope"].eq("direct_only"),
                "scope_shape_qualified",
            ].all()
        ) if "direct_only" in scopes else False
        rows.append(
            {
                **dict(zip(keys, values, strict=True)),
                "n_evidence_scopes": len(scopes),
                "all_analysis_qualified": all_ok,
                "direct_only_qualified": direct_ok,
                "headline_shape_qualified": bool(
                    scopes == {"all_analysis_eligible", "direct_only"}
                    and all_ok
                    and direct_ok
                ),
            }
        )
    return pd.DataFrame(rows)


def run_audit(
    *,
    all_counts: Path,
    direct_counts: Path,
    covariates_csv: Path,
    config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    scopes = {
        "all_analysis_eligible": pd.read_csv(all_counts),
        "direct_only": pd.read_csv(direct_counts),
    }
    covariates = pd.read_csv(covariates_csv)
    response_to_axis = {
        str(spec["response"]): str(axis)
        for axis, spec in config["measurement_domain_representatives"].items()
    }
    minimum = int(config["minimum_islands_per_design_cell"])
    design_rows = []
    calibration_parts = []
    validation_parts = []
    qualification_parts = []

    for scope_name, counts in scopes.items():
        for outcome, axis in response_to_axis.items():
            for context in config["contexts"]:
                for stratum in config["floristic_strata"]:
                    frame = prepare_cell(
                        counts,
                        covariates,
                        outcome=outcome,
                        context=str(context),
                        stratum=str(stratum),
                        config=config,
                    )
                    supported = len(frame) >= minimum
                    design_rows.append(
                        {
                            "evidence_scope": scope_name,
                            "measurement_axis": axis,
                            "outcome": outcome,
                            "context": context,
                            "stratum": stratum,
                            "n_islands": len(frame),
                            "n_spatial_blocks": (
                                frame[config["cluster_column"]].nunique()
                                if len(frame)
                                else 0
                            ),
                            "support_qualified": supported,
                        }
                    )
                    if not supported:
                        continue
                    critical_D, calibration, z, weights, candidates = calibrate_cell(
                        frame,
                        evidence_scope=scope_name,
                        measurement_axis=axis,
                        outcome=outcome,
                        context=str(context),
                        stratum=str(stratum),
                        config=config,
                    )
                    detail, qualification = validate_cell(
                        frame,
                        evidence_scope=scope_name,
                        measurement_axis=axis,
                        outcome=outcome,
                        context=str(context),
                        stratum=str(stratum),
                        critical_D=critical_D,
                        z=z,
                        weights=weights,
                        candidates=candidates,
                        config=config,
                    )
                    calibration_parts.append(calibration)
                    validation_parts.append(detail)
                    qualification_parts.append(qualification)

    designs = pd.DataFrame(design_rows)
    if not qualification_parts:
        raise ValueError("no V2 geometry cells passed the support floor")
    calibration = pd.concat(calibration_parts, ignore_index=True)
    validation = pd.concat(validation_parts, ignore_index=True)
    scope_qualification = pd.concat(
        qualification_parts, ignore_index=True
    )
    cross = cross_scope_qualification(scope_qualification)

    output_dir.mkdir(parents=True, exist_ok=True)
    designs.to_csv(output_dir / "geometry_v2_design_support.csv", index=False)
    calibration.to_csv(output_dir / "geometry_v2_calibration.csv", index=False)
    validation.to_csv(output_dir / "geometry_v2_validation.csv", index=False)
    scope_qualification.to_csv(
        output_dir / "geometry_v2_scope_shape_qualification.csv", index=False
    )
    cross.to_csv(
        output_dir / "geometry_v2_cross_scope_shape_qualification.csv", index=False
    )
    manifest = {
        "contract": config["contract"],
        "status": "v2_recovery_audit_complete_observed_geometry_still_closed",
        "parent_v1_run": int(config["parent_failure"]["workflow_run_id"]),
        "n_design_cells": int(len(designs)),
        "n_support_qualified_cells": int(designs["support_qualified"].sum()),
        "n_scope_shape_tests": int(len(scope_qualification)),
        "n_cross_scope_shape_tests": int(len(cross)),
        "n_headline_shape_qualified": int(
            cross["headline_shape_qualified"].sum()
        ),
        "headline_qualified_by_shape": {
            shape: int(
                cross.loc[
                    cross["shape"].eq(shape), "headline_shape_qualified"
                ].sum()
            )
            for shape in NONLINEAR
        },
        "observed_geometry_fitted": False,
        "observed_breakpoint_reported": False,
        "calibration_and_validation_seeds_differ": (
            int(config["calibration"]["seed"])
            != int(config["validation"]["seed"])
        ),
        "claim_boundary": (
            "V2 qualifies only the ability to distinguish specified nonlinear shapes "
            "from monotonic responses under the frozen simulation assumptions. It does "
            "not establish an observed breakpoint, non-monotonicity or mechanism."
        ),
    }
    (output_dir / "geometry_v2_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all-counts", type=Path, required=True)
    parser.add_argument("--direct-counts", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    manifest = run_audit(
        all_counts=args.all_counts,
        direct_counts=args.direct_counts,
        covariates_csv=args.covariates_csv,
        config_path=args.config_path,
        output_dir=args.output_dir,
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
