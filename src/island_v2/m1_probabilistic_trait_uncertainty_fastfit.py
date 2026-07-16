"""Fast exact-equivalent fitting for probabilistic M1 trait imputations.

The species-level imputation logic is unchanged. This runner only removes repeated
DataFrame preparation from hundreds of grouped-binomial fits by freezing the
complete-case standardized design matrices and spatial-cluster codes once. IRLS,
AIC, and the one-way geographic cluster-robust sandwich covariance match the full
analysis equations.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.m1_probabilistic_trait_uncertainty import (
    _load_config,
    aggregate_species_states,
    build_probability_mapping,
    rubin_style_pool,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


@dataclass
class PreparedDesign:
    row_indices: np.ndarray
    X: np.ndarray
    trials: np.ndarray
    cluster_codes: np.ndarray
    n_clusters: int
    n_parameters: int
    bombus_parameter_index: int | None


def _expit(values: np.ndarray) -> np.ndarray:
    values = np.clip(values, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-values))


def prepare_design(
    island_table: pd.DataFrame,
    trials: np.ndarray,
    predictors: list[str],
    spatial_cluster: str,
    bombus_predictor: str | None,
) -> PreparedDesign:
    work = island_table[[spatial_cluster, *predictors]].copy()
    for predictor in predictors:
        work[predictor] = pd.to_numeric(work[predictor], errors="coerce")
    complete = work[predictors].notna().all(axis=1) & work[spatial_cluster].notna()
    complete &= pd.Series(trials > 0, index=work.index)
    row_indices = np.flatnonzero(complete.to_numpy())
    if len(row_indices) <= len(predictors) + 3:
        raise typer.BadParameter("too few complete islands for prepared design")
    subset = work.iloc[row_indices].copy()
    means = subset[predictors].mean(axis=0)
    sds = subset[predictors].std(axis=0, ddof=0)
    bad = sds.index[(~np.isfinite(sds)) | (sds <= 0)].tolist()
    if bad:
        raise typer.BadParameter(f"constant or invalid prepared predictors: {bad}")
    standardized = (subset[predictors] - means) / sds
    X = np.column_stack([np.ones(len(subset)), standardized.to_numpy(dtype=float)])
    cluster_codes, unique_clusters = pd.factorize(
        subset[spatial_cluster].astype(str), sort=False
    )
    bombus_parameter_index = (
        1 + predictors.index(bombus_predictor)
        if bombus_predictor is not None and bombus_predictor in predictors
        else None
    )
    return PreparedDesign(
        row_indices=row_indices,
        X=X,
        trials=trials[row_indices].astype(float),
        cluster_codes=cluster_codes.astype(int),
        n_clusters=int(len(unique_clusters)),
        n_parameters=int(X.shape[1]),
        bombus_parameter_index=bombus_parameter_index,
    )


def fit_prepared_design(
    design: PreparedDesign,
    all_successes: np.ndarray,
    need_cluster_se: bool,
    max_iter: int = 100,
    tolerance: float = 1e-9,
) -> dict[str, Any]:
    successes = all_successes[design.row_indices].astype(float)
    trials = design.trials
    X = design.X
    y = successes / trials
    beta = np.zeros(design.n_parameters, dtype=float)
    converged = False

    for iteration in range(1, max_iter + 1):
        probability = np.clip(_expit(X @ beta), 1e-8, 1.0 - 1e-8)
        variance = probability * (1.0 - probability)
        weights = trials * variance
        working_response = X @ beta + (y - probability) / variance
        xtwx = X.T @ (weights[:, None] * X)
        beta_new = np.linalg.pinv(xtwx) @ (X.T @ (weights * working_response))
        if float(np.max(np.abs(beta_new - beta))) < tolerance:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    probability = np.clip(_expit(X @ beta), 1e-10, 1.0 - 1e-10)
    log_likelihood = float(
        np.sum(
            successes * np.log(probability)
            + (trials - successes) * np.log(1.0 - probability)
        )
    )
    aic = float(-2.0 * log_likelihood + 2.0 * design.n_parameters)
    result: dict[str, Any] = {
        "status": "fitted" if converged else "max_iter_reached",
        "iterations": iteration,
        "n_islands": int(len(successes)),
        "n_spatial_blocks": design.n_clusters,
        "log_likelihood": log_likelihood,
        "aic": aic,
        "beta": beta,
    }
    if not need_cluster_se:
        return result
    if design.bombus_parameter_index is None:
        raise typer.BadParameter("cluster SE requested for design without Bombus predictor")

    hessian = X.T @ ((trials * probability * (1.0 - probability))[:, None] * X)
    bread = np.linalg.pinv(hessian)
    residual_counts = successes - trials * probability
    score_rows = residual_counts[:, None] * X
    cluster_scores = np.zeros((design.n_clusters, design.n_parameters), dtype=float)
    np.add.at(cluster_scores, design.cluster_codes, score_rows)
    meat = cluster_scores.T @ cluster_scores
    n_islands = len(successes)
    if design.n_clusters > 1 and n_islands > design.n_parameters:
        correction = (design.n_clusters / (design.n_clusters - 1.0)) * (
            (n_islands - 1.0) / (n_islands - design.n_parameters)
        )
        meat *= correction
    covariance = bread @ meat @ bread
    covariance = (covariance + covariance.T) / 2.0
    index = design.bombus_parameter_index
    variance = float(covariance[index, index])
    standard_error = float(math.sqrt(variance)) if variance > 0 else np.nan
    estimate = float(beta[index])
    z_value = float(estimate / standard_error) if standard_error > 0 else np.nan
    result.update(
        {
            "estimate": estimate,
            "within_se": standard_error,
            "z_value": z_value,
            "minimum_covariance_eigenvalue": float(np.linalg.eigvalsh(covariance).min()),
        }
    )
    return result


def run_fast_imputations(
    mapping: Any,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    n_imputations = int(config["imputation"]["n_imputations"])
    seed = int(config["imputation"]["random_seed"])
    z_threshold = float(config["nominal_z_threshold"])
    bombus = str(config["bombus_predictor"])
    spatial_cluster = str(config["spatial_cluster"])
    baselines = {
        str(key): [str(value) for value in values]
        for key, values in config["baseline_specs"].items()
    }

    designs: dict[tuple[str, str], PreparedDesign] = {}
    for baseline_name, baseline in baselines.items():
        designs[(baseline_name, "M0")] = prepare_design(
            mapping.island_table,
            mapping.trials,
            baseline,
            spatial_cluster,
            None,
        )
        designs[(baseline_name, "M1")] = prepare_design(
            mapping.island_table,
            mapping.trials,
            [*baseline, bombus],
            spatial_cluster,
            bombus,
        )

    probabilities = mapping.species_table["positive_probability"].to_numpy(dtype=float)
    hard_state = mapping.species_table["filled_value"].astype(str).eq(
        str(config["outcome"]["positive_value"])
    ).astype(float).to_numpy()
    expected_successes = aggregate_species_states(mapping, probabilities)
    hard_successes = aggregate_species_states(mapping, hard_state)

    def fit_pair(successes: np.ndarray, baseline_name: str) -> dict[str, Any]:
        m0 = fit_prepared_design(
            designs[(baseline_name, "M0")], successes, need_cluster_se=False
        )
        m1 = fit_prepared_design(
            designs[(baseline_name, "M1")], successes, need_cluster_se=True
        )
        return {
            "estimate": float(m1["estimate"]),
            "within_se": float(m1["within_se"]),
            "z_value": float(m1["z_value"]),
            "n_islands": int(m1["n_islands"]),
            "n_spatial_blocks": int(m1["n_spatial_blocks"]),
            "m0_aic": float(m0["aic"]),
            "m1_aic": float(m1["aic"]),
            "delta_aic_m1_minus_m0": float(m1["aic"] - m0["aic"]),
            "minimum_covariance_eigenvalue": float(
                m1["minimum_covariance_eigenvalue"]
            ),
        }

    comparator_rows: list[dict[str, Any]] = []
    for comparator, successes in (
        ("expected_probability_counts", expected_successes),
        ("hard_modal_grounded_taxonomic", hard_successes),
    ):
        for baseline_name in baselines:
            comparator_rows.append(
                {
                    "comparator": comparator,
                    "baseline_spec": baseline_name,
                    **fit_pair(successes, baseline_name),
                }
            )

    rng = np.random.default_rng(seed)
    draw_rows: list[dict[str, Any]] = []
    for imputation_index in range(1, n_imputations + 1):
        species_state = (rng.random(len(probabilities)) < probabilities).astype(float)
        successes = aggregate_species_states(mapping, species_state)
        for baseline_name in baselines:
            draw_rows.append(
                {
                    "imputation": imputation_index,
                    "baseline_spec": baseline_name,
                    **fit_pair(successes, baseline_name),
                }
            )
    draws = pd.DataFrame(draw_rows)
    pooled = pd.DataFrame(
        [
            rubin_style_pool(draws, baseline_name, z_threshold)
            for baseline_name in baselines
        ]
    )
    return draws, pooled, pd.DataFrame(comparator_rows)


@app.command("run")
def run(
    species_master_csv: Path = typer.Option(..., exists=True),
    reference_island_data_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/m1_probabilistic_trait_uncertainty.yml"), exists=True
    ),
) -> None:
    config = _load_config(config_path)
    master = pd.read_csv(
        species_master_csv,
        usecols=[
            "island_id",
            "accepted_species",
            "trait_name",
            "filled_value",
            "fill_tier",
            "value_distribution",
        ],
    )
    reference = pd.read_csv(reference_island_data_csv)
    mapping, metadata = build_probability_mapping(master, reference, config)
    draws, pooled, comparators = run_fast_imputations(mapping, config)

    minimum = int(config["reporting"].get("minimum_imputations", 100))
    if int(draws["imputation"].nunique()) < minimum:
        raise RuntimeError("fewer valid imputations than the declared minimum")

    output_dir.mkdir(parents=True, exist_ok=True)
    mapping.species_table.to_csv(
        output_dir / "probabilistic_species_trait_table.csv.gz",
        index=False,
        compression="gzip",
    )
    mapping.island_table.to_csv(
        output_dir / "probabilistic_analysis_islands.csv", index=False
    )
    draws.to_csv(output_dir / "imputation_m1_draws.csv", index=False)
    pooled.to_csv(output_dir / "imputation_m1_pooled.csv", index=False)
    comparators.to_csv(output_dir / "imputation_m1_comparators.csv", index=False)

    manifest = {
        "contract": str(config["analysis_contract"]),
        "runner": "precomputed_design_exact_equivalent",
        "analysis_role": str(config["analysis_role"]),
        **metadata,
        "n_imputations": int(draws["imputation"].nunique()),
        "pooled_results": pooled.to_dict("records"),
        "comparators": comparators.to_dict("records"),
        "interpretation": (
            "Species-direct cells are fixed; genus/family cells are sampled once per species per imputation. "
            "The fast runner reuses fixed complete-case standardized designs but retains the same grouped-binomial IRLS and geographic cluster sandwich equations."
        ),
    }
    (output_dir / "probabilistic_trait_uncertainty_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
