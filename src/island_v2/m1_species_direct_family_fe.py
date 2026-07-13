"""Family-fixed-effect sensitivity for the species-direct M1 analysis.

The direct island-species outcome is taxonomically structured. This sensitivity
keeps only families containing both directly evidenced bumblebee-guild and
non-bumblebee-guild species, then adds family fixed effects to M0/M1. Inference
still uses two-way species x geographic-block cluster-robust covariance.

The coefficient therefore asks whether Bombus environmental compatibility is
associated with island compositional filtering *within informative families*, not
merely with turnover among families.
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

from island_v2.m1_species_direct_analysis import two_way_cluster_covariance

app = typer.Typer(add_completion=False, no_args_is_help=True)


@dataclass
class FamilyFEFit:
    beta: np.ndarray
    names: list[str]
    numeric_predictors: list[str]
    means: pd.Series
    sds: pd.Series
    family_levels: list[str]
    reference_family: str
    design: np.ndarray
    score_rows: np.ndarray
    bread: np.ndarray
    log_likelihood: float
    converged: bool
    iterations: int


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("species-direct config must be a mapping")
    return payload


def informative_family_subset(data: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {"accepted_species", "family", "outcome_bumblebee_guild"}
    missing = required.difference(data.columns)
    if missing:
        raise typer.BadParameter(f"species-direct rows missing columns: {sorted(missing)}")
    species = data[
        ["accepted_species", "family", "outcome_bumblebee_guild"]
    ].drop_duplicates("accepted_species")
    if species["accepted_species"].duplicated().any():
        raise typer.BadParameter("species outcome table is not unique by species")
    variation = species.groupby("family")["outcome_bumblebee_guild"].nunique()
    informative = sorted(variation.loc[variation > 1].index.astype(str))
    subset = data.loc[data["family"].astype(str).isin(informative)].copy()
    metadata = {
        "n_families_total": int(species["family"].nunique()),
        "n_families_informative": int(len(informative)),
        "n_species_total": int(species["accepted_species"].nunique()),
        "n_species_informative_families": int(
            subset["accepted_species"].nunique()
        ),
        "n_rows_informative_families": int(len(subset)),
        "n_islands_informative_families": int(subset["island_id"].nunique()),
        "n_spatial_blocks_informative_families": int(
            subset["spatial_block"].nunique()
        ),
    }
    return subset.reset_index(drop=True), metadata


def _expit(values: np.ndarray) -> np.ndarray:
    values = np.clip(values, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-values))


def _build_design(
    data: pd.DataFrame,
    numeric_predictors: list[str],
    means: pd.Series | None = None,
    sds: pd.Series | None = None,
    family_levels: list[str] | None = None,
    reference_family: str | None = None,
) -> tuple[np.ndarray, list[str], pd.Series, pd.Series, list[str], str]:
    numeric = data[numeric_predictors].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        raise typer.BadParameter("family-FE model received incomplete numeric predictors")
    if means is None or sds is None:
        island_numeric = (
            data[["island_id", *numeric_predictors]]
            .drop_duplicates("island_id")[numeric_predictors]
            .apply(pd.to_numeric, errors="coerce")
        )
        means = island_numeric.mean(axis=0)
        sds = island_numeric.std(axis=0, ddof=0)
    bad = sds.index[(~np.isfinite(sds)) | (sds <= 0)].tolist()
    if bad:
        raise typer.BadParameter(f"constant or invalid numeric predictors: {bad}")
    numeric_z = (numeric - means[numeric_predictors]) / sds[numeric_predictors]

    observed_families = sorted(data["family"].astype(str).unique())
    if family_levels is None:
        family_levels = observed_families
    if not family_levels:
        raise typer.BadParameter("no family levels available")
    if reference_family is None:
        reference_family = family_levels[0]
    dummy_levels = [family for family in family_levels if family != reference_family]
    family_values = data["family"].astype(str).to_numpy()
    dummies = np.column_stack(
        [(family_values == family).astype(float) for family in dummy_levels]
    ) if dummy_levels else np.empty((len(data), 0), dtype=float)
    X = np.column_stack(
        [np.ones(len(data)), numeric_z.to_numpy(dtype=float), dummies]
    )
    names = [
        "intercept",
        *numeric_predictors,
        *[f"family_FE::{family}" for family in dummy_levels],
    ]
    return X, names, means, sds, family_levels, reference_family


def fit_family_fe_logit(
    data: pd.DataFrame,
    numeric_predictors: list[str],
    weight_column: str,
    means: pd.Series | None = None,
    sds: pd.Series | None = None,
    family_levels: list[str] | None = None,
    reference_family: str | None = None,
    max_iter: int = 150,
    tolerance: float = 1e-10,
) -> FamilyFEFit:
    X, names, means, sds, family_levels, reference_family = _build_design(
        data,
        numeric_predictors,
        means=means,
        sds=sds,
        family_levels=family_levels,
        reference_family=reference_family,
    )
    y = pd.to_numeric(data["outcome_bumblebee_guild"], errors="coerce").to_numpy(
        dtype=float
    )
    weights = pd.to_numeric(data[weight_column], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(y).all() or not np.isfinite(weights).all() or (weights <= 0).any():
        raise typer.BadParameter("invalid outcome or weights in family-FE model")

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
    return FamilyFEFit(
        beta=beta,
        names=names,
        numeric_predictors=numeric_predictors,
        means=means,
        sds=sds,
        family_levels=family_levels,
        reference_family=reference_family,
        design=X,
        score_rows=score_rows,
        bread=bread,
        log_likelihood=log_likelihood,
        converged=converged,
        iterations=iteration,
    )


def family_fe_inference(
    data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    baselines = {
        "full": [str(value) for value in config["baseline_full"]],
        "no_abs_latitude": [
            str(value) for value in config["baseline_sensitivity_no_abs_latitude"]
        ],
    }
    z_threshold = float(config["reporting"].get("nominal_z_threshold", 1.96))
    coefficient_rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []

    for baseline_name, baseline in baselines.items():
        for weight_column in ("row_weighted", "equal_species_weighted"):
            m0 = fit_family_fe_logit(data, baseline, weight_column)
            m1 = fit_family_fe_logit(
                data,
                [*baseline, "bombus_channel_state"],
                weight_column,
                family_levels=m0.family_levels,
                reference_family=m0.reference_family,
            )
            for model_name, fit in (("M0_family_FE", m0), ("M1_family_FE", m1)):
                fit_rows.append(
                    {
                        "model": model_name,
                        "baseline_spec": baseline_name,
                        "weight_scheme": weight_column,
                        "n_rows": int(len(data)),
                        "n_parameters": int(len(fit.beta)),
                        "n_family_fixed_effect_levels": int(len(fit.family_levels)),
                        "reference_family": fit.reference_family,
                        "weighted_log_likelihood": fit.log_likelihood,
                        "weighted_pseudo_aic": float(
                            -2.0 * fit.log_likelihood + 2.0 * len(fit.beta)
                        ),
                        "converged": fit.converged,
                        "iterations": fit.iterations,
                    }
                )
            fit_proxy = type(
                "FitProxy",
                (),
                {
                    "score_rows": m1.score_rows,
                    "design": m1.design,
                    "bread": m1.bread,
                },
            )()
            covariance, covariance_metadata = two_way_cluster_covariance(
                fit_proxy,
                data["accepted_species"],
                data["spatial_block"],
            )
            bombus_index = m1.names.index("bombus_channel_state")
            variance = float(covariance[bombus_index, bombus_index])
            se = float(math.sqrt(variance)) if variance > 0 else np.nan
            estimate = float(m1.beta[bombus_index])
            z_value = float(estimate / se) if se > 0 else np.nan
            coefficient_rows.append(
                {
                    "model": "M1_family_FE",
                    "baseline_spec": baseline_name,
                    "weight_scheme": weight_column,
                    "covariance_scheme": "species_x_spatial_block",
                    "predictor": "bombus_channel_state",
                    "estimate_log_odds": estimate,
                    "odds_ratio_per_predictor_sd": float(math.exp(estimate)),
                    "two_way_cluster_robust_se": se,
                    "z_value": z_value,
                    "two_sided_normal_p": (
                        float(math.erfc(abs(z_value) / math.sqrt(2.0)))
                        if math.isfinite(z_value)
                        else np.nan
                    ),
                    "nominally_supported": (
                        bool(abs(z_value) >= z_threshold)
                        if math.isfinite(z_value)
                        else False
                    ),
                    "n_rows": int(len(data)),
                    "n_family_fixed_effect_levels": int(len(m1.family_levels)),
                    **covariance_metadata,
                }
            )
    return pd.DataFrame(coefficient_rows), pd.DataFrame(fit_rows)


def _predict(fit: FamilyFEFit, test: pd.DataFrame) -> np.ndarray:
    X, _, _, _, _, _ = _build_design(
        test,
        fit.numeric_predictors,
        means=fit.means,
        sds=fit.sds,
        family_levels=fit.family_levels,
        reference_family=fit.reference_family,
    )
    return np.clip(_expit(X @ fit.beta), 1e-10, 1.0 - 1e-10)


def spatial_cv_family_fe(
    data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cv = config["spatial_cv"]
    n_folds = int(cv["folds"])
    repeats = int(cv["repeats"])
    seed = int(cv["seed"])
    baselines = {
        "full": [str(value) for value in config["baseline_full"]],
        "no_abs_latitude": [
            str(value) for value in config["baseline_sensitivity_no_abs_latitude"]
        ],
    }
    island_blocks = data[["island_id", "spatial_block"]].drop_duplicates()
    block_counts = island_blocks.groupby("spatial_block").size()
    if len(block_counts) < n_folds:
        raise typer.BadParameter("too few blocks for family-FE spatial CV")

    fold_rows: list[dict[str, Any]] = []
    for repeat in range(repeats):
        rng = np.random.default_rng(seed + repeat * 1009)
        order = pd.DataFrame(
            {
                "spatial_block": block_counts.index.astype(str),
                "n_islands": block_counts.to_numpy(dtype=int),
                "tie": rng.random(len(block_counts)),
            }
        ).sort_values(["n_islands", "tie"], ascending=[False, True])
        loads = [0] * n_folds
        mapping: dict[str, int] = {}
        for row in order.itertuples(index=False):
            fold = min(range(n_folds), key=lambda index: (loads[index], index))
            mapping[str(row.spatial_block)] = fold
            loads[fold] += int(row.n_islands)
        fold_id = data["spatial_block"].astype(str).map(mapping)

        for fold in range(n_folds):
            train = data.loc[~fold_id.eq(fold)].copy()
            test = data.loc[fold_id.eq(fold)].copy()
            for baseline_name, baseline in baselines.items():
                for weight_column in ("row_weighted", "equal_species_weighted"):
                    m0 = fit_family_fe_logit(train, baseline, weight_column)
                    m1 = fit_family_fe_logit(
                        train,
                        [*baseline, "bombus_channel_state"],
                        weight_column,
                        family_levels=m0.family_levels,
                        reference_family=m0.reference_family,
                    )
                    y = test["outcome_bumblebee_guild"].to_numpy(dtype=float)
                    weights = test[weight_column].to_numpy(dtype=float)
                    p0 = _predict(m0, test)
                    p1 = _predict(m1, test)
                    ll0 = float(
                        np.sum(
                            weights
                            * (y * np.log(p0) + (1.0 - y) * np.log(1.0 - p0))
                        )
                    )
                    ll1 = float(
                        np.sum(
                            weights
                            * (y * np.log(p1) + (1.0 - y) * np.log(1.0 - p1))
                        )
                    )
                    brier0 = float(np.sum(weights * (y - p0) ** 2) / np.sum(weights))
                    brier1 = float(np.sum(weights * (y - p1) ** 2) / np.sum(weights))
                    fold_rows.append(
                        {
                            "repeat": repeat + 1,
                            "fold": fold + 1,
                            "baseline_spec": baseline_name,
                            "weight_scheme": weight_column,
                            "n_train_rows": int(len(train)),
                            "n_test_rows": int(len(test)),
                            "test_weight_sum": float(np.sum(weights)),
                            "delta_log_likelihood_m1_minus_m0": ll1 - ll0,
                            "delta_brier_m1_minus_m0": brier1 - brier0,
                        }
                    )

    folds = pd.DataFrame(fold_rows)
    summary_rows: list[dict[str, Any]] = []
    for (baseline_spec, weight_scheme), group in folds.groupby(
        ["baseline_spec", "weight_scheme"], sort=True
    ):
        repeat_rows: list[dict[str, float]] = []
        for _, repeat_group in group.groupby("repeat", sort=True):
            total_weight = float(repeat_group["test_weight_sum"].sum())
            repeat_rows.append(
                {
                    "delta_log_likelihood_per_1000_weight": float(
                        1000.0
                        * repeat_group["delta_log_likelihood_m1_minus_m0"].sum()
                        / total_weight
                    ),
                    "weighted_delta_brier": float(
                        np.average(
                            repeat_group["delta_brier_m1_minus_m0"],
                            weights=repeat_group["test_weight_sum"],
                        )
                    ),
                }
            )
        repeat_table = pd.DataFrame(repeat_rows)
        ll = repeat_table["delta_log_likelihood_per_1000_weight"].to_numpy(dtype=float)
        brier = repeat_table["weighted_delta_brier"].to_numpy(dtype=float)
        summary_rows.append(
            {
                "baseline_spec": baseline_spec,
                "weight_scheme": weight_scheme,
                "n_repeats": int(len(repeat_table)),
                "median_delta_log_likelihood_per_1000_weight": float(np.median(ll)),
                "proportion_repeats_m1_better_loglik": float(np.mean(ll > 0)),
                "median_delta_brier_m1_minus_m0": float(np.median(brier)),
                "proportion_repeats_m1_better_brier": float(np.mean(brier < 0)),
                "predictive_gain_consistent": bool(np.all(ll > 0) and np.all(brier < 0)),
            }
        )
    return folds, pd.DataFrame(summary_rows)


@app.command("run")
def run(
    species_direct_rows_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/m1_species_direct_analysis.yml"), exists=True
    ),
) -> None:
    config = _load_config(config_path)
    data = pd.read_csv(species_direct_rows_csv)
    subset, subset_metadata = informative_family_subset(data)
    coefficients, fit = family_fe_inference(subset, config)
    cv_folds, cv_summary = spatial_cv_family_fe(subset, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    subset.to_csv(
        output_dir / "species_direct_informative_family_rows.csv.gz",
        index=False,
        compression="gzip",
    )
    coefficients.to_csv(output_dir / "family_fe_bombus_coefficients.csv", index=False)
    fit.to_csv(output_dir / "family_fe_model_fit.csv", index=False)
    cv_folds.to_csv(output_dir / "family_fe_spatial_cv_folds.csv", index=False)
    cv_summary.to_csv(output_dir / "family_fe_spatial_cv_summary.csv", index=False)

    manifest = {
        "contract": "m1_species_direct_family_fixed_effect_sensitivity_v1",
        **subset_metadata,
        "coefficients": coefficients.to_dict("records"),
        "spatial_cv": cv_summary.to_dict("records"),
        "interpretation": (
            "Family fixed effects restrict inference to within-family compositional sorting among families with both outcome classes. "
            "Species and geographic blocks remain the two cluster dimensions."
        ),
    }
    (output_dir / "family_fe_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
