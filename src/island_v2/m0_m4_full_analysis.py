"""Fit the locked-data M0-M4 full path-analysis ladder.

This is the first inferential model layer beyond the exploratory OLS snapshot. It
uses grouped-binomial responses (species counts, not naked proportions), a frozen
baseline covariate set, and cluster-robust sandwich uncertainty over geographic
blocks. The Bombus predictor remains environmental compatibility and must not be
reported as observed biological absence.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("full-analysis config must be a mapping")
    return payload


def _proxy_counts(
    composition: pd.DataFrame,
    trait_name: str,
    positive_values: list[str],
    prefix: str,
) -> pd.DataFrame:
    required = {"island_id", "trait_name", "filled_value", "n_species"}
    missing = required.difference(composition.columns)
    if missing:
        raise typer.BadParameter(f"trait composition missing columns: {sorted(missing)}")
    work = composition.loc[composition["trait_name"].astype(str).eq(str(trait_name))].copy()
    if work.empty:
        raise typer.BadParameter(f"trait composition has no rows for {trait_name}")
    work["n_species"] = pd.to_numeric(work["n_species"], errors="coerce")
    work = work.dropna(subset=["n_species"])
    work = work.loc[work["n_species"] >= 0].copy()
    positives = {str(value) for value in positive_values}
    work["is_positive"] = work["filled_value"].astype(str).isin(positives)

    total = work.groupby("island_id", sort=True)["n_species"].sum().rename(f"{prefix}_trials")
    success = (
        work.loc[work["is_positive"]]
        .groupby("island_id", sort=True)["n_species"]
        .sum()
        .rename(f"{prefix}_successes")
    )
    result = (
        total.to_frame()
        .join(success, how="left")
        .fillna({f"{prefix}_successes": 0.0})
        .reset_index()
    )
    trials = pd.to_numeric(result[f"{prefix}_trials"], errors="coerce")
    successes = pd.to_numeric(result[f"{prefix}_successes"], errors="coerce")
    result[f"{prefix}_share"] = successes / trials.replace(0, np.nan)
    return result


def build_full_analysis_data(
    composition: pd.DataFrame,
    island_metrics: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, dict[str, Any]]:
    for label, table in {
        "island_metrics": island_metrics,
        "covariates": covariates,
    }.items():
        if "island_id" not in table.columns:
            raise typer.BadParameter(f"{label} missing island_id")
        if table["island_id"].duplicated().any():
            raise typer.BadParameter(f"{label} must contain one row per island_id")

    proxy_cfg = config["proxies"]
    proxy_tables = []
    for prefix, key in (
        ("floral_phenotype", "floral_phenotype"),
        ("mixed_mating", "mixed_mating"),
        ("self_compatibility", "self_compatibility"),
        ("alternative_pollen_vector", "alternative_pollen_vector"),
    ):
        spec = proxy_cfg[key]
        proxy_tables.append(
            _proxy_counts(
                composition,
                str(spec["trait_name"]),
                [str(value) for value in spec["positive_values"]],
                prefix,
            )
        )

    result = island_metrics.drop_duplicates("island_id").merge(
        covariates.drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    for proxy_table in proxy_tables:
        result = result.merge(proxy_table, on="island_id", how="left", validate="one_to_one")

    predictor = str(config["bombus_channel"]["predictor"])
    if predictor not in result.columns:
        raise typer.BadParameter(f"island metrics missing Bombus predictor {predictor}")
    result["bombus_channel_state"] = pd.to_numeric(result[predictor], errors="coerce")
    result["log_trait_species_richness"] = np.log1p(
        pd.to_numeric(result.get("n_trait_species"), errors="coerce")
    )
    result["alternative_pollen_vector"] = pd.to_numeric(
        result["alternative_pollen_vector_share"], errors="coerce"
    )
    result["bombus_x_alternative"] = (
        pd.to_numeric(result["bombus_channel_state"], errors="coerce")
        * result["alternative_pollen_vector"]
    )

    baseline = [str(value) for value in config["baseline_covariates"]]
    missing_baseline = [column for column in baseline if column not in result.columns]
    metadata = {
        "baseline_covariates": baseline,
        "missing_baseline_covariates": missing_baseline,
        "baseline_complete": not missing_baseline,
        "bombus_predictor": predictor,
        "mediators": ["mixed_mating", "self_compatibility"],
        "non_informative_traits": config.get("non_informative_traits", {}),
    }
    return result, metadata


def _expit(values: np.ndarray) -> np.ndarray:
    values = np.clip(values, -35.0, 35.0)
    return 1.0 / (1.0 + np.exp(-values))


def _prepare_predictors(
    data: pd.DataFrame,
    predictors: list[str],
) -> tuple[np.ndarray, list[str], dict[str, dict[str, float]]]:
    columns = []
    scaling: dict[str, dict[str, float]] = {}
    for predictor in predictors:
        values = pd.to_numeric(data[predictor], errors="coerce").to_numpy(dtype=float)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=0))
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant or invalid predictor in complete-case data: {predictor}")
        columns.append((values - mean) / sd)
        scaling[predictor] = {"mean": mean, "sd": sd}
    X = np.column_stack([np.ones(len(data)), *columns])
    return X, ["intercept", *predictors], scaling


def fit_grouped_binomial_clustered(
    data: pd.DataFrame,
    successes_column: str,
    trials_column: str,
    predictors: list[str],
    cluster_column: str,
    z_threshold: float,
    max_iter: int = 100,
    tolerance: float = 1e-9,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    columns = [successes_column, trials_column, cluster_column, *predictors]
    work = data[columns].copy()
    for column in [successes_column, trials_column, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    work = work.loc[
        (work[trials_column] > 0)
        & (work[successes_column] >= 0)
        & (work[successes_column] <= work[trials_column])
    ].copy()

    n_islands = int(len(work))
    if n_islands < max(10, len(predictors) + 3):
        return pd.DataFrame(), {
            "status": "insufficient_complete_islands",
            "n_islands": n_islands,
            "n_clusters": int(work[cluster_column].nunique()),
        }

    successes = work[successes_column].to_numpy(dtype=float)
    trials = work[trials_column].to_numpy(dtype=float)
    y = successes / trials
    try:
        X, names, scaling = _prepare_predictors(work, predictors)
    except ValueError as exc:
        return pd.DataFrame(), {
            "status": "constant_or_invalid_predictor",
            "n_islands": n_islands,
            "n_clusters": int(work[cluster_column].nunique()),
            "error": str(exc),
        }
    beta = np.zeros(X.shape[1], dtype=float)
    converged = False

    for iteration in range(1, max_iter + 1):
        eta = X @ beta
        probability = np.clip(_expit(eta), 1e-8, 1.0 - 1e-8)
        variance = probability * (1.0 - probability)
        weights = trials * variance
        working_response = eta + (y - probability) / variance
        xtwx = X.T @ (weights[:, None] * X)
        beta_new = np.linalg.pinv(xtwx) @ (X.T @ (weights * working_response))
        if float(np.max(np.abs(beta_new - beta))) < tolerance:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    probability = np.clip(_expit(X @ beta), 1e-10, 1.0 - 1e-10)
    weights = trials * probability * (1.0 - probability)
    bread = np.linalg.pinv(X.T @ (weights[:, None] * X))

    residual_counts = successes - trials * probability
    clusters = work[cluster_column].fillna("missing").astype(str).to_numpy()
    unique_clusters = np.unique(clusters)
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for cluster in unique_clusters:
        mask = clusters == cluster
        score = X[mask].T @ residual_counts[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n_clusters = int(len(unique_clusters))
    p = int(X.shape[1])
    if n_clusters > 1 and n_islands > p:
        correction = (n_clusters / (n_clusters - 1.0)) * (
            (n_islands - 1.0) / (n_islands - p)
        )
        covariance *= correction
    standard_errors = np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    log_likelihood = float(
        np.sum(
            successes * np.log(probability)
            + (trials - successes) * np.log(1.0 - probability)
        )
    )
    mean_probability = float(
        np.clip(successes.sum() / trials.sum(), 1e-10, 1.0 - 1e-10)
    )
    null_log_likelihood = float(
        np.sum(
            successes * np.log(mean_probability)
            + (trials - successes) * np.log(1.0 - mean_probability)
        )
    )
    aic = float(-2.0 * log_likelihood + 2.0 * p)
    pseudo_r2 = (
        float(1.0 - log_likelihood / null_log_likelihood)
        if null_log_likelihood != 0
        else np.nan
    )
    pearson = float(
        np.sum(
            (successes - trials * probability) ** 2
            / np.maximum(trials * probability * (1.0 - probability), 1e-12)
        )
    )
    dispersion = float(pearson / max(n_islands - p, 1))

    rows: list[dict[str, Any]] = []
    for index, (name, estimate, standard_error) in enumerate(
        zip(names, beta, standard_errors, strict=True)
    ):
        z_value = float(estimate / standard_error) if standard_error > 0 else np.nan
        rows.append(
            {
                "predictor": name,
                "estimate_log_odds": float(estimate),
                "cluster_robust_se": float(standard_error),
                "z_value": z_value,
                "nominally_supported": (
                    bool(abs(z_value) >= z_threshold) if math.isfinite(z_value) else False
                ),
                "predictor_mean": 0.0 if index == 0 else scaling[name]["mean"],
                "predictor_sd": 1.0 if index == 0 else scaling[name]["sd"],
                "n_islands": n_islands,
                "n_clusters": n_clusters,
            }
        )
    fit = {
        "status": "fitted" if converged else "max_iter_reached",
        "n_islands": n_islands,
        "n_clusters": n_clusters,
        "iterations": iteration,
        "log_likelihood": log_likelihood,
        "aic": aic,
        "mcfadden_pseudo_r2": pseudo_r2,
        "pearson_dispersion": dispersion,
    }
    return pd.DataFrame(rows), fit


def _equations(baseline: list[str]) -> dict[str, list[tuple[str, str, str, list[str]]]]:
    return {
        "M0": [
            (
                "floral_phenotype",
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                baseline,
            ),
        ],
        "M1": [
            (
                "floral_phenotype",
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [*baseline, "bombus_channel_state"],
            ),
        ],
        "M2": [
            (
                "mixed_mating",
                "mixed_mating_successes",
                "mixed_mating_trials",
                [*baseline, "bombus_channel_state"],
            ),
            (
                "self_compatibility",
                "self_compatibility_successes",
                "self_compatibility_trials",
                [*baseline, "bombus_channel_state"],
            ),
            (
                "floral_phenotype",
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [*baseline, "mixed_mating_share", "self_compatibility_share"],
            ),
        ],
        "M3": [
            (
                "mixed_mating",
                "mixed_mating_successes",
                "mixed_mating_trials",
                [*baseline, "bombus_channel_state"],
            ),
            (
                "self_compatibility",
                "self_compatibility_successes",
                "self_compatibility_trials",
                [*baseline, "bombus_channel_state"],
            ),
            (
                "floral_phenotype",
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [
                    *baseline,
                    "bombus_channel_state",
                    "mixed_mating_share",
                    "self_compatibility_share",
                ],
            ),
        ],
        "M4": [
            (
                "floral_phenotype",
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [
                    *baseline,
                    "bombus_channel_state",
                    "alternative_pollen_vector",
                    "bombus_x_alternative",
                ],
            ),
        ],
    }


def fit_model_ladder(
    data: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    baseline = [str(value) for value in config["baseline_covariates"]]
    missing = [column for column in baseline if column not in data.columns]
    if missing:
        raise typer.BadParameter(f"full-analysis data missing baseline covariates: {missing}")
    cluster_column = "spatial_block"
    if cluster_column not in data.columns:
        raise typer.BadParameter("full-analysis data missing spatial_block")
    z_threshold = float(config["reporting"].get("nominal_z_threshold", 1.96))

    coefficient_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    lookup: dict[tuple[str, str, str], tuple[float, float]] = {}

    for model, equations in _equations(baseline).items():
        for equation_number, (outcome, successes, trials, predictors) in enumerate(
            equations, start=1
        ):
            coefficients, fit = fit_grouped_binomial_clustered(
                data,
                successes,
                trials,
                predictors,
                cluster_column,
                z_threshold,
            )
            if not coefficients.empty:
                coefficients.insert(0, "outcome", outcome)
                coefficients.insert(0, "equation", equation_number)
                coefficients.insert(0, "model", model)
                coefficient_parts.append(coefficients)
                for row in coefficients.itertuples(index=False):
                    lookup[(model, outcome, row.predictor)] = (
                        float(row.estimate_log_odds),
                        float(row.cluster_robust_se),
                    )
            fit_rows.append(
                {
                    "model": model,
                    "equation": equation_number,
                    "outcome": outcome,
                    "predictors": " + ".join(predictors),
                    **fit,
                }
            )

    indirect_rows: list[dict[str, Any]] = []
    for model in ("M2", "M3"):
        products: list[tuple[str, float, float]] = []
        for mediator, mediator_predictor in (
            ("mixed_mating", "mixed_mating_share"),
            ("self_compatibility", "self_compatibility_share"),
        ):
            a = lookup.get((model, mediator, "bombus_channel_state"))
            b = lookup.get((model, "floral_phenotype", mediator_predictor))
            if a is None or b is None:
                estimate = np.nan
                standard_error = np.nan
            else:
                estimate = float(a[0] * b[0])
                standard_error = float(
                    math.sqrt((b[0] ** 2) * (a[1] ** 2) + (a[0] ** 2) * (b[1] ** 2))
                )
            products.append((mediator, estimate, standard_error))
            indirect_rows.append(
                {
                    "model": model,
                    "mediator": mediator,
                    "effect_scale": "product_of_standardized_logit_coefficients",
                    "indirect_effect": estimate,
                    "delta_method_se_approx": standard_error,
                }
            )
        finite = [value for _, value, _ in products if math.isfinite(value)]
        indirect_rows.append(
            {
                "model": model,
                "mediator": "parallel_total",
                "effect_scale": "sum_of_product_coefficients_logit_scale",
                "indirect_effect": float(sum(finite)) if finite else np.nan,
                "delta_method_se_approx": np.nan,
            }
        )

    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True) if coefficient_parts else pd.DataFrame()
    )
    fits = pd.DataFrame(fit_rows)
    indirect = pd.DataFrame(indirect_rows)
    model_summary = (
        fits.groupby("model", sort=True)
        .agg(
            n_equations=("equation", "count"),
            total_aic=("aic", "sum"),
            total_log_likelihood=("log_likelihood", "sum"),
            minimum_equation_islands=("n_islands", "min"),
            minimum_spatial_blocks=("n_clusters", "min"),
        )
        .reset_index()
    )
    return coefficients, fits, indirect, model_summary


def summarize_tier_robustness(
    coefficient_tables: dict[str, pd.DataFrame],
    z_threshold: float = 1.96,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    key_columns = ["model", "equation", "outcome", "predictor"]
    indexed = {
        tier: table.set_index(key_columns) for tier, table in coefficient_tables.items()
    }
    all_keys = sorted(set().union(*(set(table.index) for table in indexed.values())))
    for key in all_keys:
        estimates: dict[str, float] = {}
        z_values: dict[str, float] = {}
        for tier, table in indexed.items():
            if key in table.index:
                row = table.loc[key]
                if isinstance(row, pd.DataFrame):
                    row = row.iloc[0]
                estimates[tier] = float(row["estimate_log_odds"])
                z_values[tier] = float(row["z_value"])
        finite_estimates = [value for value in estimates.values() if math.isfinite(value)]
        nonzero_signs = {int(np.sign(value)) for value in finite_estimates if value != 0}
        rows.append(
            {
                "model": key[0],
                "equation": key[1],
                "outcome": key[2],
                "predictor": key[3],
                "n_tiers_available": len(finite_estimates),
                "direction_consistent_across_tiers": (
                    len(finite_estimates) >= 2 and len(nonzero_signs) <= 1
                ),
                "all_tiers_nominally_supported": (
                    len(z_values) >= 2
                    and all(abs(value) >= z_threshold for value in z_values.values())
                ),
                "max_absolute_coefficient_difference": (
                    float(max(finite_estimates) - min(finite_estimates))
                    if finite_estimates
                    else np.nan
                ),
                **{
                    f"{tier}_estimate": estimates.get(tier, np.nan)
                    for tier in coefficient_tables
                },
                **{f"{tier}_z": z_values.get(tier, np.nan) for tier in coefficient_tables},
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    trait_composition_csv: Path = typer.Option(..., exists=True),
    island_metrics_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    analysis_tier: str = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_full_analysis.yml"), exists=True),
) -> None:
    config = _load_config(config_path)
    composition = pd.read_csv(trait_composition_csv)
    metrics = pd.read_csv(island_metrics_csv)
    covariates = pd.read_csv(covariates_csv)
    data, metadata = build_full_analysis_data(composition, metrics, covariates, config)

    minimum = int(config["reporting"].get("minimum_complete_islands", 50))
    if int(data["island_id"].nunique()) < minimum:
        raise typer.BadParameter(f"fewer than {minimum} islands available before model fitting")

    coefficients, fits, indirect, model_summary = fit_model_ladder(data, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    data.to_csv(output_dir / "full_analysis_island_data.csv", index=False)
    coefficients.to_csv(output_dir / "full_analysis_path_coefficients.csv", index=False)
    fits.to_csv(output_dir / "full_analysis_equation_fit.csv", index=False)
    indirect.to_csv(output_dir / "full_analysis_indirect_effects.csv", index=False)
    model_summary.to_csv(output_dir / "full_analysis_model_summary.csv", index=False)

    summary_lookup = model_summary.set_index("model")["total_aic"].to_dict()
    manifest = {
        "contract": str(config["analysis_contract"]),
        "analysis_tier": analysis_tier,
        "n_islands_input": int(data["island_id"].nunique()),
        "n_spatial_blocks_input": int(data["spatial_block"].nunique()),
        **metadata,
        "m1_minus_m0_total_aic": (
            float(summary_lookup["M1"] - summary_lookup["M0"])
            if "M1" in summary_lookup and "M0" in summary_lookup
            else None
        ),
        "m3_minus_m2_total_aic": (
            float(summary_lookup["M3"] - summary_lookup["M2"])
            if "M3" in summary_lookup and "M2" in summary_lookup
            else None
        ),
        "m4_alternative_proxy": "pollen_vector_mode: abiotic_wind",
        "m4_previous_complement_proxy_removed": True,
        "interpretation": (
            "Grouped-binomial island-level path models with geographic-cluster robust uncertainty. "
            "Bombus is environmental compatibility, not observed absence."
        ),
    }
    (output_dir / "full_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
