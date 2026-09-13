"""Prospective response-geometry recovery audit for Chapter 1.

This module does not fit or report the observed nonlinear response. It asks whether
frozen island support can distinguish predeclared response shapes under simulation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

GEOMETRIES = ("G0_flat", "G1_cline", "G2_step", "G3_hinge", "G4_reversal")
NONMONOTONIC = {"G2_step", "G3_hinge", "G4_reversal"}
STEP_OR_HINGE = {"G2_step", "G3_hinge"}


def _standardize(values: pd.Series | np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _aicc_vector(sse: np.ndarray, n: int, k: int) -> np.ndarray:
    sse = np.clip(np.asarray(sse, dtype=float), 1e-12, None)
    aic = n * np.log(sse / n) + 2.0 * k
    denominator = n - k - 1
    if denominator <= 0:
        return np.full_like(sse, np.inf)
    return aic + (2.0 * k * (k + 1)) / denominator


@dataclass
class Candidate:
    shape: str
    breakpoint_quantile: float | None
    q_basis: np.ndarray
    beta_map: np.ndarray
    n_parameters: int
    valid: bool = True


def _candidate(
    shape: str,
    X: np.ndarray,
    weights: np.ndarray,
    *,
    breakpoint_quantile: float | None = None,
    extra_parameter_penalty: int = 0,
) -> Candidate:
    sqrtw = np.sqrt(weights)
    weighted = X * sqrtw[:, None]
    rank = int(np.linalg.matrix_rank(weighted))
    if rank != X.shape[1]:
        return Candidate(
            shape,
            breakpoint_quantile,
            np.empty((len(X), 0)),
            np.empty((0, len(X))),
            0,
            False,
        )
    q_basis, _ = np.linalg.qr(weighted, mode="reduced")
    xtwx = X.T @ (weights[:, None] * X)
    beta_map = np.linalg.solve(xtwx, X.T * weights[None, :])
    # Regression coefficients + residual variance + declared breakpoint-search penalty.
    n_parameters = X.shape[1] + 1 + int(extra_parameter_penalty)
    return Candidate(
        shape,
        breakpoint_quantile,
        q_basis,
        beta_map,
        n_parameters,
        True,
    )


def build_candidates(
    frame: pd.DataFrame, config: dict[str, Any]
) -> tuple[np.ndarray, np.ndarray, list[Candidate]]:
    exposure = str(config["primary_exposure"]["column"])
    z = _standardize(frame[exposure])
    base = [np.ones(len(frame), dtype=float)]
    for column in config["baseline_covariates"]:
        base.append(_standardize(frame[str(column)]))
    base_matrix = np.column_stack(base)
    weights = frame["trials"].to_numpy(float)
    candidates: list[Candidate] = []
    candidates.append(_candidate("G0_flat", base_matrix, weights))
    candidates.append(
        _candidate("G1_cline", np.column_stack([base_matrix, z]), weights)
    )

    for q in config["candidate_geometries"]["G2_step"][
        "breakpoint_search_quantiles"
    ]:
        knot = float(np.quantile(z, float(q)))
        step = (z > knot).astype(float)
        if step.min() == step.max():
            continue
        candidates.append(
            _candidate(
                "G2_step",
                np.column_stack([base_matrix, step]),
                weights,
                breakpoint_quantile=float(q),
                extra_parameter_penalty=int(
                    config["candidate_geometries"]["G2_step"][
                        "breakpoint_parameter_penalty"
                    ]
                ),
            )
        )

    for q in config["candidate_geometries"]["G3_hinge"][
        "breakpoint_search_quantiles"
    ]:
        knot = float(np.quantile(z, float(q)))
        hinge = np.maximum(0.0, z - knot)
        if float(np.std(hinge)) <= 0:
            continue
        candidates.append(
            _candidate(
                "G3_hinge",
                np.column_stack([base_matrix, z, hinge]),
                weights,
                breakpoint_quantile=float(q),
                extra_parameter_penalty=int(
                    config["candidate_geometries"]["G3_hinge"][
                        "breakpoint_parameter_penalty"
                    ]
                ),
            )
        )

    candidates.append(
        _candidate("G4_reversal", np.column_stack([base_matrix, z, z * z]), weights)
    )
    if not all(
        candidate.valid
        for candidate in candidates
        if candidate.shape in {"G0_flat", "G1_cline", "G4_reversal"}
    ):
        raise ValueError("rank-deficient mandatory geometry design")
    return z, weights, [candidate for candidate in candidates if candidate.valid]


def _shape_signal(
    z: np.ndarray, shape: str, breakpoint_quantile: float
) -> np.ndarray:
    knot = float(np.quantile(z, breakpoint_quantile))
    if shape == "G0_flat":
        raw = np.zeros_like(z)
    elif shape == "G1_cline":
        raw = z.copy()
    elif shape == "G2_step":
        raw = (z > knot).astype(float)
    elif shape == "G3_hinge":
        # Continuous erosion that becomes much shallower after the breakpoint.
        raw = np.where(z <= knot, z - knot, 0.2 * (z - knot))
    elif shape == "G4_reversal":
        raw = -((z - knot) ** 2)
    else:
        raise ValueError(f"unknown shape: {shape}")
    raw = raw - float(np.mean(raw))
    sd = float(np.std(raw, ddof=0))
    if shape != "G0_flat" and (not np.isfinite(sd) or sd <= 0):
        raise ValueError(f"invalid generated shape: {shape}")
    return raw if shape == "G0_flat" else raw / sd


def _score_batch(
    y: np.ndarray,
    z: np.ndarray,
    weights: np.ndarray,
    candidates: list[Candidate],
    config: dict[str, Any],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return selected shape, selected breakpoint quantile and delta-AICc per replicate."""
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

    shape_scores = []
    shape_breakpoints = []
    for shape in GEOMETRIES:
        options = per_shape[shape]
        scores = np.vstack([item[0] for item in options])
        index = np.argmin(scores, axis=0)
        selected_score = scores[index, np.arange(replicates)]
        breakpoints = np.vstack([item[1] for item in options])
        selected_breakpoint = breakpoints[index, np.arange(replicates)]
        shape_scores.append(selected_score)
        shape_breakpoints.append(selected_breakpoint)
    score_matrix = np.vstack(shape_scores)
    breakpoint_matrix = np.vstack(shape_breakpoints)
    order = np.argsort(score_matrix, axis=0)
    best_index = order[0]
    second_index = order[1]
    columns = np.arange(replicates)
    selected = np.asarray(GEOMETRIES, dtype=object)[best_index]
    selected_breakpoint = breakpoint_matrix[best_index, columns]
    delta = score_matrix[second_index, columns] - score_matrix[best_index, columns]
    return selected, selected_breakpoint, delta


def _stable_seed(base: int, *parts: str) -> int:
    token = "|".join(parts).encode("utf-8")
    digest = hashlib.sha256(token).hexdigest()[:8]
    return (int(base) + int(digest, 16)) % (2**32 - 1)


def prepare_cell(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    outcome: str,
    context: str,
    stratum: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    work = counts.loc[
        counts["outcome"].astype(str).eq(outcome)
        & counts["stratum"].astype(str).eq(stratum)
    ].copy()
    required_cov = [
        "island_id",
        "analysis_regime",
        config["cluster_column"],
        config["primary_exposure"]["column"],
        *config["baseline_covariates"],
    ]
    work = work.merge(
        covariates[required_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    work = work.loc[work["analysis_regime"].astype(str).eq(context)].copy()
    numeric = [
        "successes",
        "trials",
        config["primary_exposure"]["column"],
        *config["baseline_covariates"],
    ]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna(subset=[*numeric, config["cluster_column"]])
    work = work.loc[work["trials"].gt(0)].copy()
    return work.sort_values("island_id").reset_index(drop=True)


def simulate_cell(
    frame: pd.DataFrame,
    *,
    evidence_scope: str,
    measurement_axis: str,
    outcome: str,
    context: str,
    stratum: str,
    config: dict[str, Any],
) -> pd.DataFrame:
    z, weights, candidates = build_candidates(frame, config)
    trials = frame["trials"].astype(int).to_numpy()
    p0 = float(
        (frame["successes"].sum() + 0.5) / (frame["trials"].sum() + 1.0)
    )
    p0 = min(max(p0, 1e-4), 1 - 1e-4)
    intercept = math.log(p0 / (1.0 - p0))
    blocks = frame[config["cluster_column"]].astype(str).to_numpy()
    unique_blocks, block_index = np.unique(blocks, return_inverse=True)
    replicates = int(config["simulation"]["replicates"])
    rng = np.random.default_rng(
        _stable_seed(
            int(config["simulation"]["seed"]),
            evidence_scope,
            measurement_axis,
            outcome,
            context,
            stratum,
        )
    )
    rows: list[pd.DataFrame] = []
    for effect in config["simulation"]["effect_logit_sd_grid"]:
        effect = float(effect)
        for cluster_sd in config["simulation"]["cluster_random_intercept_sd_grid"]:
            cluster_sd = float(cluster_sd)
            for truth in GEOMETRIES:
                truth_quantiles = (
                    [
                        float(x)
                        for x in config["simulation"]["truth_breakpoint_quantiles"]
                    ]
                    if truth in {"G2_step", "G3_hinge"}
                    else [0.50]
                )
                for truth_quantile in truth_quantiles:
                    signal = _shape_signal(z, truth, truth_quantile)
                    # Alternate direction so recovery is not tied to one response sign.
                    signs = np.where(
                        np.arange(replicates) % 2 == 0, 1.0, -1.0
                    )
                    random_effects = (
                        rng.normal(
                            0.0,
                            cluster_sd,
                            size=(len(unique_blocks), replicates),
                        )
                        if cluster_sd > 0
                        else np.zeros((len(unique_blocks), replicates))
                    )
                    eta = (
                        intercept
                        + signal[:, None] * effect * signs[None, :]
                        + random_effects[block_index]
                    )
                    probability = 1.0 / (
                        1.0 + np.exp(-np.clip(eta, -30, 30))
                    )
                    successes = rng.binomial(trials[:, None], probability)
                    failures = trials[:, None] - successes
                    y = np.log((successes + 0.5) / (failures + 0.5))
                    selected, selected_q, delta = _score_batch(
                        y, z, weights, candidates, config
                    )
                    batch = pd.DataFrame(
                        {
                            "evidence_scope": evidence_scope,
                            "measurement_axis": measurement_axis,
                            "outcome": outcome,
                            "context": context,
                            "stratum": stratum,
                            "n_islands": len(frame),
                            "n_spatial_blocks": len(unique_blocks),
                            "effect_logit_sd": effect,
                            "cluster_random_intercept_sd": cluster_sd,
                            "truth_shape": truth,
                            "truth_breakpoint_quantile": truth_quantile,
                            "replicate": np.arange(replicates),
                            "selected_shape": selected,
                            "selected_breakpoint_quantile": selected_q,
                            "delta_AICc_to_second": delta,
                        }
                    )
                    batch["shape_recovered"] = batch["selected_shape"].eq(truth)
                    batch["strong_shape_recovered"] = batch[
                        "shape_recovered"
                    ] & batch["delta_AICc_to_second"].ge(
                        float(
                            config["selection_diagnostic"][
                                "strong_selection_delta_AICc"
                            ]
                        )
                    )
                    batch["breakpoint_abs_error"] = np.where(
                        batch["shape_recovered"]
                        & batch["truth_shape"].isin(["G2_step", "G3_hinge"]),
                        np.abs(
                            batch["selected_breakpoint_quantile"]
                            - truth_quantile
                        ),
                        np.nan,
                    )
                    rows.append(batch)
    return pd.concat(rows, ignore_index=True)


def summarize_simulations(
    simulations: pd.DataFrame, config: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    group = [
        "evidence_scope",
        "measurement_axis",
        "outcome",
        "context",
        "stratum",
        "n_islands",
        "n_spatial_blocks",
        "effect_logit_sd",
        "cluster_random_intercept_sd",
        "truth_shape",
        "truth_breakpoint_quantile",
    ]
    summary = (
        simulations.groupby(group, as_index=False)
        .agg(
            recovery_rate=("shape_recovered", "mean"),
            strong_recovery_rate=("strong_shape_recovered", "mean"),
            breakpoint_quantile_mae=("breakpoint_abs_error", "mean"),
            median_delta_AICc=("delta_AICc_to_second", "median"),
        )
    )
    target_effect = float(config["qualification"]["target_effect_logit_sd"])
    stress = float(config["qualification"]["clustered_stress_sd"])
    target = simulations.loc[
        np.isclose(simulations["effect_logit_sd"], target_effect)
        & np.isclose(simulations["cluster_random_intercept_sd"], stress)
    ].copy()
    keys = [
        "evidence_scope",
        "measurement_axis",
        "outcome",
        "context",
        "stratum",
        "n_islands",
        "n_spatial_blocks",
    ]
    rows: list[dict[str, Any]] = []
    for values, part in target.groupby(keys, sort=True):
        row = dict(zip(keys, values, strict=True))
        cline = part.loc[part["truth_shape"].eq("G1_cline")]
        false_nonmonotonic = float(
            cline["selected_shape"].isin(NONMONOTONIC).mean()
        )
        false_step_hinge = float(
            cline["selected_shape"].isin(STEP_OR_HINGE).mean()
        )
        recoveries: dict[str, float] = {}
        breakpoint_mae: dict[str, float] = {}
        for shape in ["G1_cline", "G2_step", "G3_hinge", "G4_reversal"]:
            subset = part.loc[part["truth_shape"].eq(shape)]
            recoveries[shape] = float(subset["shape_recovered"].mean())
            if shape in STEP_OR_HINGE:
                recovered = subset.loc[subset["shape_recovered"]]
                breakpoint_mae[shape] = (
                    float(recovered["breakpoint_abs_error"].mean())
                    if len(recovered)
                    else float("nan")
                )
        min_recovery = min(recoveries.values())
        finite_breakpoints = [
            value for value in breakpoint_mae.values() if np.isfinite(value)
        ]
        max_breakpoint_mae = (
            max(finite_breakpoints) if finite_breakpoints else float("nan")
        )
        passed = (
            min_recovery
            >= float(config["qualification"]["minimum_true_shape_recovery"])
            and false_nonmonotonic
            <= float(
                config["qualification"][
                    "maximum_false_nonmonotonic_selection_under_true_cline"
                ]
            )
            and false_step_hinge
            <= float(
                config["qualification"][
                    "maximum_false_step_or_hinge_selection_under_true_cline"
                ]
            )
            and np.isfinite(max_breakpoint_mae)
            and max_breakpoint_mae
            <= float(
                config["qualification"][
                    "maximum_breakpoint_quantile_mae_when_shape_recovered"
                ]
            )
        )
        row.update(
            {
                "target_effect_logit_sd": target_effect,
                "clustered_stress_sd": stress,
                **{f"recovery_{key}": value for key, value in recoveries.items()},
                "minimum_nonflat_recovery": min_recovery,
                "false_nonmonotonic_under_cline": false_nonmonotonic,
                "false_step_or_hinge_under_cline": false_step_hinge,
                "breakpoint_mae_G2_step": breakpoint_mae.get(
                    "G2_step", np.nan
                ),
                "breakpoint_mae_G3_hinge": breakpoint_mae.get(
                    "G3_hinge", np.nan
                ),
                "maximum_breakpoint_mae": max_breakpoint_mae,
                "scope_qualified": bool(passed),
            }
        )
        rows.append(row)
    scope_qualification = pd.DataFrame(rows)

    cross_keys = ["measurement_axis", "outcome", "context", "stratum"]
    cross_rows = []
    for values, part in scope_qualification.groupby(cross_keys, sort=True):
        scopes = set(part["evidence_scope"].astype(str))
        both = scopes == {"all_analysis_eligible", "direct_only"}
        all_analysis = bool(
            part.loc[
                part["evidence_scope"].eq("all_analysis_eligible"),
                "scope_qualified",
            ].all()
        ) if "all_analysis_eligible" in scopes else False
        direct = bool(
            part.loc[
                part["evidence_scope"].eq("direct_only"), "scope_qualified"
            ].all()
        ) if "direct_only" in scopes else False
        cross_rows.append(
            {
                **dict(zip(cross_keys, values, strict=True)),
                "n_evidence_scopes": len(scopes),
                "all_analysis_qualified": all_analysis,
                "direct_only_qualified": direct,
                "headline_geometry_qualified": bool(both and all_analysis and direct),
            }
        )
    return summary, pd.DataFrame(cross_rows), scope_qualification


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
    design_rows = []
    simulation_parts = []
    minimum = int(config["minimum_islands_per_design_cell"])
    for scope, counts in scopes.items():
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
                    design_rows.append(
                        {
                            "evidence_scope": scope,
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
                            "median_trials": (
                                float(frame["trials"].median())
                                if len(frame)
                                else np.nan
                            ),
                            "distance_min": (
                                float(
                                    frame[
                                        config["primary_exposure"]["column"]
                                    ].min()
                                )
                                if len(frame)
                                else np.nan
                            ),
                            "distance_max": (
                                float(
                                    frame[
                                        config["primary_exposure"]["column"]
                                    ].max()
                                )
                                if len(frame)
                                else np.nan
                            ),
                            "support_qualified": len(frame) >= minimum,
                        }
                    )
                    if len(frame) < minimum:
                        continue
                    simulation_parts.append(
                        simulate_cell(
                            frame,
                            evidence_scope=scope,
                            measurement_axis=axis,
                            outcome=outcome,
                            context=str(context),
                            stratum=str(stratum),
                            config=config,
                        )
                    )
    designs = pd.DataFrame(design_rows)
    if not simulation_parts:
        raise ValueError(
            "no response-geometry design cells meet the support threshold"
        )
    simulations = pd.concat(simulation_parts, ignore_index=True)
    summary, cross, scope_qualification = summarize_simulations(
        simulations, config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    designs.to_csv(output_dir / "geometry_design_support.csv", index=False)
    summary.to_csv(output_dir / "geometry_recovery_summary.csv", index=False)
    scope_qualification.to_csv(
        output_dir / "geometry_scope_qualification.csv", index=False
    )
    cross.to_csv(
        output_dir / "geometry_cross_scope_qualification.csv", index=False
    )
    # Replicate-level simulations are intentionally not persisted; gate summaries suffice.
    manifest = {
        "contract": config["contract"],
        "status": "recovery_audit_complete_observed_geometry_still_closed",
        "n_design_cells": int(len(designs)),
        "n_support_qualified_cells": int(designs["support_qualified"].sum()),
        "n_headline_geometry_qualified_cells": int(
            cross["headline_geometry_qualified"].sum()
        ),
        "n_cross_scope_cells": int(len(cross)),
        "observed_geometry_fitted": False,
        "observed_breakpoint_reported": False,
        "future_water_gap_status": config["future_physical_exposure"][
            "source_specific_water_gap"
        ]["status"],
        "claim_boundary": (
            "Recovery qualification is a design result only. A qualified cell may "
            "proceed to the separately tracked observed-geometry stage; it does not "
            "establish a threshold, reversal or mechanism in the biological data."
        ),
    }
    (output_dir / "geometry_recovery_manifest.json").write_text(
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
