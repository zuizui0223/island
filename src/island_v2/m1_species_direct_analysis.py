"""Confirmatory M1 analysis on species-direct island x species observations.

The island-aggregated primary response is sparse (median six directly evidenced
species per island). This module therefore analyzes the directly evidenced
island-species rows themselves. The response is whether a directly evidenced
species belongs to the bumblebee pollination-functional guild.

Species traits are invariant across islands, so ordinary row-wise standard errors
would be pseudoreplicated. Inference uses Cameron-Gelbach-Miller style two-way
cluster-robust covariance over (species, geographic block), with a conservative
(family, geographic block) sensitivity. A second weighting scheme gives each
species equal total weight so geographically widespread species cannot dominate.
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
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


@dataclass
class WeightedLogitFit:
    beta: np.ndarray
    predictors: list[str]
    means: pd.Series
    sds: pd.Series
    probability: np.ndarray
    weights: np.ndarray
    design: np.ndarray
    score_rows: np.ndarray
    bread: np.ndarray
    log_likelihood: float
    converged: bool
    iterations: int


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("species-direct analysis config must be a mapping")
    return payload


def build_species_direct_data(
    species_master: pd.DataFrame,
    reference_island_data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    outcome_cfg = config["outcome"]
    required_master = {
        "island_id",
        "accepted_species",
        "family",
        "n_unique_gbif_ids",
        "trait_name",
        "filled_value",
        "fill_tier",
    }
    missing = required_master.difference(species_master.columns)
    if missing:
        raise typer.BadParameter(f"species master missing columns: {sorted(missing)}")

    # Observation-process proxy from the already materialized flora occurrences.
    # n_unique_gbif_ids is repeated across trait rows, so collapse to island x species first.
    occurrence = species_master[
        ["island_id", "accepted_species", "n_unique_gbif_ids"]
    ].copy()
    occurrence["n_unique_gbif_ids"] = pd.to_numeric(
        occurrence["n_unique_gbif_ids"], errors="coerce"
    ).fillna(0)
    occurrence = (
        occurrence.groupby(["island_id", "accepted_species"], as_index=False)[
            "n_unique_gbif_ids"
        ]
        .max()
    )
    effort = (
        occurrence.groupby("island_id", sort=True)["n_unique_gbif_ids"]
        .sum()
        .rename("gbif_flora_records")
        .reset_index()
    )
    effort["log_gbif_flora_records"] = np.log1p(effort["gbif_flora_records"])

    direct = species_master.loc[
        species_master["trait_name"].astype(str).eq(str(outcome_cfg["trait_name"]))
        & species_master["fill_tier"].astype(str).eq(
            str(outcome_cfg["required_fill_tier"])
        )
    ].copy()
    direct = direct.dropna(subset=["island_id", "accepted_species", "filled_value"])
    if direct.empty:
        raise typer.BadParameter("no species-direct outcome rows were found")
    if direct.duplicated(["island_id", "accepted_species"]).any():
        raise typer.BadParameter(
            "species-direct outcome must contain at most one row per island x species"
        )

    # A species-level direct trait must be internally stable across all island occurrences.
    inconsistent = (
        direct.groupby("accepted_species")["filled_value"].nunique().loc[lambda x: x > 1]
    )
    if len(inconsistent):
        raise typer.BadParameter(
            "species-direct outcome varies within species: "
            f"{inconsistent.head(20).to_dict()}"
        )
    direct["outcome_bumblebee_guild"] = direct["filled_value"].astype(str).eq(
        str(outcome_cfg["positive_value"])
    ).astype(int)
    direct["family"] = direct["family"].fillna("UNKNOWN_FAMILY").astype(str)

    baseline_full = [str(value) for value in config["baseline_full"]]
    baseline_no_latitude = [
        str(value) for value in config["baseline_sensitivity_no_abs_latitude"]
    ]
    reference_columns = list(
        dict.fromkeys(
            [
                "island_id",
                "spatial_block",
                str(config["bombus_predictor"]["column"]),
                *baseline_full,
                *baseline_no_latitude,
            ]
        )
    )
    missing_reference = set(reference_columns).difference(reference_island_data.columns)
    # log_gbif_flora_records is built above, not expected in the reference artifact.
    missing_reference.discard("log_gbif_flora_records")
    if missing_reference:
        raise typer.BadParameter(
            f"reference island data missing columns: {sorted(missing_reference)}"
        )
    reference_columns = [
        column for column in reference_columns if column != "log_gbif_flora_records"
    ]
    reference = reference_island_data[reference_columns].drop_duplicates("island_id")
    if reference["island_id"].duplicated().any():
        raise typer.BadParameter("reference island data must contain one row per island_id")

    result = direct.merge(reference, on="island_id", how="left", validate="many_to_one")
    result = result.merge(effort, on="island_id", how="left", validate="many_to_one")
    bombus_column = str(config["bombus_predictor"]["column"])
    if bombus_column != "bombus_channel_state":
        result["bombus_channel_state"] = pd.to_numeric(
            result[bombus_column], errors="coerce"
        )
    else:
        result["bombus_channel_state"] = pd.to_numeric(
            result["bombus_channel_state"], errors="coerce"
        )

    species_frequency = result.groupby("accepted_species").size()
    result["row_weighted"] = 1.0
    result["equal_species_weighted"] = result["accepted_species"].map(
        1.0 / species_frequency
    )
    # Normalize to mean one so likelihood magnitudes remain numerically comparable.
    result["equal_species_weighted"] /= float(result["equal_species_weighted"].mean())

    required_complete = list(
        dict.fromkeys(
            [
                "outcome_bumblebee_guild",
                "accepted_species",
                "family",
                "spatial_block",
                "bombus_channel_state",
                *baseline_full,
                *baseline_no_latitude,
            ]
        )
    )
    for column in required_complete:
        if column in {
            "accepted_species",
            "family",
            "spatial_block",
        }:
            continue
        result[column] = pd.to_numeric(result[column], errors="coerce")
    result = result.dropna(subset=required_complete).copy()
    result = result.loc[result["spatial_block"].astype(str).ne("missing")].copy()

    metadata = {
        "n_rows": int(len(result)),
        "n_islands": int(result["island_id"].nunique()),
        "n_species": int(result["accepted_species"].nunique()),
        "n_families": int(result["family"].nunique()),
        "n_spatial_blocks": int(result["spatial_block"].nunique()),
        "outcome_positive_rows": int(result["outcome_bumblebee_guild"].sum()),
        "outcome_positive_share": float(result["outcome_bumblebee_guild"].mean()),
        "median_islands_per_species": float(
            result.groupby("accepted_species").size().median()
        ),
        "max_islands_per_species": int(result.groupby("accepted_species").size().max()),
        "baseline_full": baseline_full,
        "baseline_no_abs_latitude": baseline_no_latitude,
        "observation_effort_covariate": "log_gbif_flora_records",
    }
    return result.reset_index(drop=True), metadata


def _expit(values: np.ndarray) -> np.ndarray:
    values = np.clip(values, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-values))


def _island_scaling(data: pd.DataFrame, predictors: list[str]) -> tuple[pd.Series, pd.Series]:
    unique_islands = data[["island_id", *predictors]].drop_duplicates("island_id")
    means = unique_islands[predictors].mean(axis=0)
    sds = unique_islands[predictors].std(axis=0, ddof=0)
    bad = sds.index[(~np.isfinite(sds)) | (sds <= 0)].tolist()
    if bad:
        raise typer.BadParameter(f"constant or invalid island predictors: {bad}")
    return means, sds


def fit_weighted_logit(
    data: pd.DataFrame,
    predictors: list[str],
    weight_column: str,
    means: pd.Series | None = None,
    sds: pd.Series | None = None,
    max_iter: int = 100,
    tolerance: float = 1e-10,
) -> WeightedLogitFit:
    columns = [
        "island_id",
        "outcome_bumblebee_guild",
        weight_column,
        *predictors,
    ]
    work = data[columns].copy()
    for column in ["outcome_bumblebee_guild", weight_column, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    if len(work) != len(data):
        raise typer.BadParameter("fit_weighted_logit received incomplete rows")
    if means is None or sds is None:
        means, sds = _island_scaling(data, predictors)

    standardized = (work[predictors] - means[predictors]) / sds[predictors]
    X = np.column_stack([np.ones(len(work)), standardized.to_numpy(dtype=float)])
    y = work["outcome_bumblebee_guild"].to_numpy(dtype=float)
    weights = work[weight_column].to_numpy(dtype=float)
    if (weights <= 0).any() or not np.isfinite(weights).all():
        raise typer.BadParameter(f"invalid analysis weights in {weight_column}")

    beta = np.zeros(X.shape[1], dtype=float)
    converged = False
    for iteration in range(1, max_iter + 1):
        probability = np.clip(_expit(X @ beta), 1e-8, 1.0 - 1e-8)
        variance = probability * (1.0 - probability)
        irls_weight = weights * variance
        working_response = X @ beta + (y - probability) / variance
        xtwx = X.T @ (irls_weight[:, None] * X)
        beta_new = np.linalg.pinv(xtwx) @ (X.T @ (irls_weight * working_response))
        if float(np.max(np.abs(beta_new - beta))) < tolerance:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    probability = np.clip(_expit(X @ beta), 1e-10, 1.0 - 1e-10)
    hessian = X.T @ ((weights * probability * (1.0 - probability))[:, None] * X)
    bread = np.linalg.pinv(hessian)
    score_rows = (weights * (y - probability))[:, None] * X
    log_likelihood = float(
        np.sum(
            weights
            * (y * np.log(probability) + (1.0 - y) * np.log(1.0 - probability))
        )
    )
    return WeightedLogitFit(
        beta=beta,
        predictors=predictors,
        means=means[predictors],
        sds=sds[predictors],
        probability=probability,
        weights=weights,
        design=X,
        score_rows=score_rows,
        bread=bread,
        log_likelihood=log_likelihood,
        converged=converged,
        iterations=iteration,
    )


def _cluster_meat(
    score_rows: np.ndarray,
    labels: pd.Series,
    n_parameters: int,
) -> tuple[np.ndarray, int]:
    codes, uniques = pd.factorize(labels.astype(str), sort=False)
    meat = np.zeros((n_parameters, n_parameters), dtype=float)
    for code in range(len(uniques)):
        cluster_score = score_rows[codes == code].sum(axis=0)
        meat += np.outer(cluster_score, cluster_score)
    n_clusters = int(len(uniques))
    return meat, n_clusters


def _finite_cluster_correction(
    meat: np.ndarray,
    n_clusters: int,
    n_rows: int,
    n_parameters: int,
) -> np.ndarray:
    if n_clusters <= 1 or n_rows <= n_parameters:
        return meat
    correction = (n_clusters / (n_clusters - 1.0)) * (
        (n_rows - 1.0) / (n_rows - n_parameters)
    )
    return meat * correction


def two_way_cluster_covariance(
    fit: WeightedLogitFit,
    cluster_a: pd.Series,
    cluster_b: pd.Series,
) -> tuple[np.ndarray, dict[str, Any]]:
    n_rows = len(fit.score_rows)
    n_parameters = fit.design.shape[1]
    if len(cluster_a) != n_rows or len(cluster_b) != n_rows:
        raise typer.BadParameter("cluster labels do not align with fitted rows")

    meat_a, n_a = _cluster_meat(fit.score_rows, cluster_a, n_parameters)
    meat_b, n_b = _cluster_meat(fit.score_rows, cluster_b, n_parameters)
    intersection = cluster_a.astype(str) + "||" + cluster_b.astype(str)
    meat_ab, n_ab = _cluster_meat(fit.score_rows, intersection, n_parameters)
    meat_a = _finite_cluster_correction(meat_a, n_a, n_rows, n_parameters)
    meat_b = _finite_cluster_correction(meat_b, n_b, n_rows, n_parameters)
    meat_ab = _finite_cluster_correction(meat_ab, n_ab, n_rows, n_parameters)

    covariance = fit.bread @ (meat_a + meat_b - meat_ab) @ fit.bread
    covariance = (covariance + covariance.T) / 2.0
    diagonal = np.diag(covariance)
    metadata = {
        "n_cluster_a": n_a,
        "n_cluster_b": n_b,
        "n_intersection_clusters": n_ab,
        "minimum_covariance_eigenvalue": float(np.linalg.eigvalsh(covariance).min()),
        "n_nonpositive_variance_diagonals": int((diagonal <= 0).sum()),
        "covariance_indefinite": bool(np.linalg.eigvalsh(covariance).min() < -1e-10),
    }
    return covariance, metadata


def _normal_two_sided_p(z_value: float) -> float:
    return float(math.erfc(abs(z_value) / math.sqrt(2.0)))


def coefficient_table(
    fit: WeightedLogitFit,
    covariance: np.ndarray,
    covariance_metadata: dict[str, Any],
    model: str,
    baseline_spec: str,
    weight_scheme: str,
    covariance_scheme: str,
    z_threshold: float,
) -> pd.DataFrame:
    names = ["intercept", *fit.predictors]
    variance = np.diag(covariance)
    rows: list[dict[str, Any]] = []
    for index, (name, estimate) in enumerate(zip(names, fit.beta, strict=True)):
        standard_error = float(math.sqrt(variance[index])) if variance[index] > 0 else np.nan
        z_value = float(estimate / standard_error) if standard_error > 0 else np.nan
        rows.append(
            {
                "model": model,
                "baseline_spec": baseline_spec,
                "weight_scheme": weight_scheme,
                "covariance_scheme": covariance_scheme,
                "predictor": name,
                "estimate_log_odds": float(estimate),
                "odds_ratio_per_predictor_sd": float(math.exp(estimate)),
                "two_way_cluster_robust_se": standard_error,
                "z_value": z_value,
                "two_sided_normal_p": (
                    _normal_two_sided_p(z_value) if math.isfinite(z_value) else np.nan
                ),
                "nominally_supported": (
                    bool(abs(z_value) >= z_threshold) if math.isfinite(z_value) else False
                ),
                "predictor_mean": 0.0 if index == 0 else float(fit.means[name]),
                "predictor_sd": 1.0 if index == 0 else float(fit.sds[name]),
                "n_rows": int(len(fit.probability)),
                **covariance_metadata,
            }
        )
    return pd.DataFrame(rows)


def fit_inference_grid(
    data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baselines = {
        "full": [str(value) for value in config["baseline_full"]],
        "no_abs_latitude": [
            str(value) for value in config["baseline_sensitivity_no_abs_latitude"]
        ],
    }
    weight_columns = {
        "row_weighted": "row_weighted",
        "equal_species_weighted": "equal_species_weighted",
    }
    covariance_schemes = {
        "species_x_spatial_block": ("accepted_species", "spatial_block"),
        "family_x_spatial_block": ("family", "spatial_block"),
    }
    z_threshold = float(config["reporting"].get("nominal_z_threshold", 1.96))

    coefficient_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    for baseline_name, baseline in baselines.items():
        means_m0, sds_m0 = _island_scaling(data, baseline)
        means_m1, sds_m1 = _island_scaling(data, [*baseline, "bombus_channel_state"])
        for weight_name, weight_column in weight_columns.items():
            for model, predictors, means, sds in (
                ("M0", baseline, means_m0, sds_m0),
                ("M1", [*baseline, "bombus_channel_state"], means_m1, sds_m1),
            ):
                fit = fit_weighted_logit(
                    data,
                    predictors,
                    weight_column,
                    means=means,
                    sds=sds,
                )
                fit_rows.append(
                    {
                        "model": model,
                        "baseline_spec": baseline_name,
                        "weight_scheme": weight_name,
                        "n_rows": int(len(data)),
                        "n_parameters": int(len(fit.beta)),
                        "weighted_log_likelihood": fit.log_likelihood,
                        "weighted_pseudo_aic": float(-2.0 * fit.log_likelihood + 2.0 * len(fit.beta)),
                        "converged": fit.converged,
                        "iterations": fit.iterations,
                    }
                )
                for covariance_name, (cluster_a, cluster_b) in covariance_schemes.items():
                    covariance, covariance_metadata = two_way_cluster_covariance(
                        fit,
                        data[cluster_a],
                        data[cluster_b],
                    )
                    coefficient_parts.append(
                        coefficient_table(
                            fit,
                            covariance,
                            covariance_metadata,
                            model,
                            baseline_name,
                            weight_name,
                            covariance_name,
                            z_threshold,
                        )
                    )
    coefficients = pd.concat(coefficient_parts, ignore_index=True)
    return coefficients, pd.DataFrame(fit_rows)


def assign_spatial_folds(
    island_blocks: pd.DataFrame,
    n_folds: int,
    seed: int,
) -> dict[str, int]:
    block_counts = (
        island_blocks.drop_duplicates(["island_id", "spatial_block"])
        .groupby("spatial_block", sort=True)
        .size()
        .rename("n_islands")
    )
    if len(block_counts) < n_folds:
        raise typer.BadParameter(
            f"fewer spatial blocks ({len(block_counts)}) than folds ({n_folds})"
        )
    rng = np.random.default_rng(seed)
    ordered = pd.DataFrame(
        {
            "spatial_block": block_counts.index.astype(str),
            "n_islands": block_counts.to_numpy(dtype=int),
            "tie_break": rng.random(len(block_counts)),
        }
    ).sort_values(["n_islands", "tie_break"], ascending=[False, True])
    loads = [0] * n_folds
    mapping: dict[str, int] = {}
    for row in ordered.itertuples(index=False):
        fold = min(range(n_folds), key=lambda index: (loads[index], index))
        mapping[str(row.spatial_block)] = fold
        loads[fold] += int(row.n_islands)
    return mapping


def _heldout_weighted_metrics(
    data: pd.DataFrame,
    probability: np.ndarray,
    weight_column: str,
) -> dict[str, float]:
    y = data["outcome_bumblebee_guild"].to_numpy(dtype=float)
    weights = data[weight_column].to_numpy(dtype=float)
    log_likelihood = float(
        np.sum(
            weights
            * (y * np.log(probability) + (1.0 - y) * np.log(1.0 - probability))
        )
    )
    brier = float(np.sum(weights * (y - probability) ** 2) / np.sum(weights))
    return {
        "weighted_log_likelihood": log_likelihood,
        "weighted_brier": brier,
        "sum_weights": float(np.sum(weights)),
    }


def spatial_cv(
    data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cv_cfg = config["spatial_cv"]
    n_folds = int(cv_cfg["folds"])
    repeats = int(cv_cfg["repeats"])
    seed = int(cv_cfg["seed"])
    baselines = {
        "full": [str(value) for value in config["baseline_full"]],
        "no_abs_latitude": [
            str(value) for value in config["baseline_sensitivity_no_abs_latitude"]
        ],
    }
    weight_columns = [str(value) for value in cv_cfg["weight_schemes"]]

    fold_rows: list[dict[str, Any]] = []
    island_blocks = data[["island_id", "spatial_block"]].drop_duplicates()
    for repeat in range(repeats):
        repeat_seed = seed + repeat * 1009
        mapping = assign_spatial_folds(island_blocks, n_folds, repeat_seed)
        fold_id = data["spatial_block"].astype(str).map(mapping)
        for fold in range(n_folds):
            test = data.loc[fold_id.eq(fold)].copy()
            train = data.loc[~fold_id.eq(fold)].copy()
            if test.empty or train.empty:
                continue
            train_species = set(train["accepted_species"].astype(str))
            test_species = set(test["accepted_species"].astype(str))
            unseen_test_species = test_species - train_species

            for baseline_name, baseline in baselines.items():
                for weight_column in weight_columns:
                    means_m0, sds_m0 = _island_scaling(train, baseline)
                    means_m1, sds_m1 = _island_scaling(
                        train, [*baseline, "bombus_channel_state"]
                    )
                    m0 = fit_weighted_logit(
                        train,
                        baseline,
                        weight_column,
                        means=means_m0,
                        sds=sds_m0,
                    )
                    m1 = fit_weighted_logit(
                        train,
                        [*baseline, "bombus_channel_state"],
                        weight_column,
                        means=means_m1,
                        sds=sds_m1,
                    )

                    def predict(fit: WeightedLogitFit) -> np.ndarray:
                        standardized = (
                            test[fit.predictors] - fit.means[fit.predictors]
                        ) / fit.sds[fit.predictors]
                        X = np.column_stack(
                            [np.ones(len(test)), standardized.to_numpy(dtype=float)]
                        )
                        return np.clip(_expit(X @ fit.beta), 1e-10, 1.0 - 1e-10)

                    m0_metrics = _heldout_weighted_metrics(
                        test, predict(m0), weight_column
                    )
                    m1_metrics = _heldout_weighted_metrics(
                        test, predict(m1), weight_column
                    )
                    fold_rows.append(
                        {
                            "repeat": repeat + 1,
                            "fold": fold + 1,
                            "seed": repeat_seed,
                            "baseline_spec": baseline_name,
                            "weight_scheme": weight_column,
                            "n_train_rows": int(len(train)),
                            "n_test_rows": int(len(test)),
                            "n_train_islands": int(train["island_id"].nunique()),
                            "n_test_islands": int(test["island_id"].nunique()),
                            "n_test_spatial_blocks": int(test["spatial_block"].nunique()),
                            "n_test_species": int(len(test_species)),
                            "n_test_species_unseen_in_train": int(len(unseen_test_species)),
                            "share_test_species_unseen_in_train": float(
                                len(unseen_test_species) / len(test_species)
                                if test_species
                                else np.nan
                            ),
                            "test_weight_sum": m1_metrics["sum_weights"],
                            "delta_log_likelihood_m1_minus_m0": float(
                                m1_metrics["weighted_log_likelihood"]
                                - m0_metrics["weighted_log_likelihood"]
                            ),
                            "delta_log_likelihood_per_1000_weight": float(
                                1000.0
                                * (
                                    m1_metrics["weighted_log_likelihood"]
                                    - m0_metrics["weighted_log_likelihood"]
                                )
                                / m1_metrics["sum_weights"]
                            ),
                            "delta_brier_m1_minus_m0": float(
                                m1_metrics["weighted_brier"]
                                - m0_metrics["weighted_brier"]
                            ),
                        }
                    )

    folds = pd.DataFrame(fold_rows)
    summary_rows: list[dict[str, Any]] = []
    for (baseline_spec, weight_scheme), group in folds.groupby(
        ["baseline_spec", "weight_scheme"], sort=True
    ):
        repeat_summary = (
            group.groupby("repeat", sort=True)
            .apply(
                lambda g: pd.Series(
                    {
                        "delta_log_likelihood_per_1000_weight": float(
                            1000.0
                            * g["delta_log_likelihood_m1_minus_m0"].sum()
                            / g["test_weight_sum"].sum()
                        ),
                        "weighted_delta_brier": float(
                            np.average(
                                g["delta_brier_m1_minus_m0"],
                                weights=g["test_weight_sum"],
                            )
                        ),
                    }
                ),
                include_groups=False,
            )
            .reset_index()
        )
        ll = repeat_summary["delta_log_likelihood_per_1000_weight"].to_numpy(
            dtype=float
        )
        brier = repeat_summary["weighted_delta_brier"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "baseline_spec": baseline_spec,
                "weight_scheme": weight_scheme,
                "n_repeats": int(len(repeat_summary)),
                "median_delta_log_likelihood_per_1000_weight": float(np.median(ll)),
                "min_delta_log_likelihood_per_1000_weight": float(np.min(ll)),
                "max_delta_log_likelihood_per_1000_weight": float(np.max(ll)),
                "proportion_repeats_m1_better_loglik": float(np.mean(ll > 0)),
                "median_delta_brier_m1_minus_m0": float(np.median(brier)),
                "proportion_repeats_m1_better_brier": float(np.mean(brier < 0)),
                "predictive_gain_consistent": bool(np.all(ll > 0) and np.all(brier < 0)),
            }
        )
    return folds, pd.DataFrame(summary_rows)


@app.command("run")
def run(
    species_master_csv: Path = typer.Option(..., exists=True),
    reference_island_data_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/m1_species_direct_analysis.yml"), exists=True
    ),
) -> None:
    config = _load_config(config_path)
    species_master = pd.read_csv(species_master_csv)
    reference = pd.read_csv(reference_island_data_csv)
    data, metadata = build_species_direct_data(species_master, reference, config)

    reporting = config["reporting"]
    if len(data) < int(reporting["minimum_rows"]):
        raise typer.BadParameter("too few species-direct rows")
    if data["accepted_species"].nunique() < int(reporting["minimum_species_clusters"]):
        raise typer.BadParameter("too few species clusters")
    if data["family"].nunique() < int(reporting["minimum_family_clusters"]):
        raise typer.BadParameter("too few family clusters")
    if data["spatial_block"].nunique() < int(reporting["minimum_spatial_blocks"]):
        raise typer.BadParameter("too few spatial blocks")

    coefficients, fit = fit_inference_grid(data, config)
    cv_folds, cv_summary = spatial_cv(data, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(output_dir / "species_direct_analysis_rows.csv.gz", index=False, compression="gzip")
    coefficients.to_csv(output_dir / "species_direct_coefficients.csv", index=False)
    fit.to_csv(output_dir / "species_direct_model_fit.csv", index=False)
    cv_folds.to_csv(output_dir / "species_direct_spatial_cv_folds.csv", index=False)
    cv_summary.to_csv(output_dir / "species_direct_spatial_cv_summary.csv", index=False)

    bombus_coefficients = coefficients.loc[
        coefficients["predictor"].eq("bombus_channel_state")
    ].copy()
    manifest = {
        "contract": str(config["analysis_contract"]),
        **metadata,
        "bombus_coefficients": bombus_coefficients.to_dict("records"),
        "spatial_cv": cv_summary.to_dict("records"),
        "primary_inference": {
            "baseline_spec": "full",
            "weight_scheme": "row_weighted",
            "covariance_scheme": "species_x_spatial_block",
        },
        "sensitivity_inference": [
            "equal_species_weighted",
            "family_x_spatial_block",
            "no_abs_latitude",
        ],
        "interpretation": (
            "Direct-evidence island-species composition model. Repeated species and geographic blocks are explicitly clustered; "
            "Bombus is environmental compatibility, not observed absence."
        ),
    }
    (output_dir / "species_direct_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
