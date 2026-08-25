"""Canonical Chapter 1 M0-M2 context-dependent grouped-binomial analysis.

M0: baseline covariates only.
M1: M0 + isolation (universal-isolation baseline).
M2: M1 + biogeographic context + isolation x context (primary Chapter 1 test).

The fitted outcomes are retained trait categories. No Bombus, pollinator occurrence,
or pollination-syndrome variable is required or admitted as a primary predictor.
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


def _expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -35.0, 35.0)))


def _two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _fit_grouped_binomial_design(
    successes: np.ndarray,
    trials: np.ndarray,
    X: np.ndarray,
    names: list[str],
    clusters: np.ndarray,
    *,
    max_iter: int = 100,
    tolerance: float = 1e-9,
) -> tuple[pd.DataFrame, dict[str, Any], np.ndarray]:
    beta = np.zeros(X.shape[1], dtype=float)
    y = successes / trials
    converged = False
    iteration = 0
    for iteration in range(1, max_iter + 1):
        eta = X @ beta
        p = np.clip(_expit(eta), 1e-8, 1.0 - 1e-8)
        variance = p * (1.0 - p)
        weights = trials * variance
        z = eta + (y - p) / variance
        xtwx = X.T @ (weights[:, None] * X)
        new_beta = np.linalg.pinv(xtwx) @ (X.T @ (weights * z))
        if float(np.max(np.abs(new_beta - beta))) < tolerance:
            beta = new_beta
            converged = True
            break
        beta = new_beta

    p = np.clip(_expit(X @ beta), 1e-10, 1.0 - 1e-10)
    weights = trials * p * (1.0 - p)
    bread = np.linalg.pinv(X.T @ (weights[:, None] * X))
    residual_counts = successes - trials * p
    unique_clusters = np.unique(clusters.astype(str))
    meat = np.zeros((X.shape[1], X.shape[1]), dtype=float)
    for cluster in unique_clusters:
        mask = clusters.astype(str) == cluster
        score = X[mask].T @ residual_counts[mask]
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    n = len(successes)
    k = X.shape[1]
    g = len(unique_clusters)
    if g > 1 and n > k:
        covariance *= (g / (g - 1.0)) * ((n - 1.0) / (n - k))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    rows = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z_value = float(estimate / stderr) if stderr > 0 else float("nan")
        rows.append({
            "predictor": name,
            "estimate_log_odds": float(estimate),
            "cluster_robust_se": float(stderr),
            "z_value": z_value,
            "p_value": _two_sided_p(z_value) if math.isfinite(z_value) else float("nan"),
        })

    log_likelihood = float(np.sum(
        successes * np.log(p) + (trials - successes) * np.log(1.0 - p)
    ))
    aic = float(-2.0 * log_likelihood + 2.0 * k)
    fit = {
        "status": "fit" if converged else "max_iter_reached",
        "n_islands": int(n),
        "n_clusters": int(g),
        "n_parameters": int(k),
        "iterations": int(iteration),
        "log_likelihood": log_likelihood,
        "aic": aic,
    }
    return pd.DataFrame(rows), fit, covariance


def _z(values: pd.Series) -> tuple[np.ndarray, dict[str, float]]:
    x = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not math.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid continuous predictor")
    return (x - mean) / sd, {"mean": mean, "sd": sd}


def _fit_one_category(
    data: pd.DataFrame,
    *,
    baseline: list[str],
    isolation_column: str,
    context_column: str,
    cluster_column: str,
    reference_context: str,
    min_islands: int,
    min_context_islands: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    required = {"successes", "trials", isolation_column, context_column, cluster_column, *baseline}
    missing = required - set(data.columns)
    if missing:
        raise typer.BadParameter(f"Chapter 1 model input missing columns: {sorted(missing)}")

    work = data.copy()
    numeric = ["successes", "trials", isolation_column, *baseline]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[context_column] = work[context_column].fillna("").astype(str)
    work[cluster_column] = work[cluster_column].fillna("").astype(str)
    work = work.dropna(subset=numeric)
    work = work.loc[
        (work["trials"] > 0)
        & (work["successes"] >= 0)
        & (work["successes"] <= work["trials"])
        & work[context_column].ne("")
        & work[cluster_column].ne("")
    ].copy()
    if len(work) < min_islands:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {
            "status": "below_min_islands", "n_islands": int(len(work))
        }

    context_counts = work[context_column].value_counts()
    eligible_contexts = sorted(context_counts.loc[context_counts >= min_context_islands].index.astype(str))
    work = work.loc[work[context_column].isin(eligible_contexts)].copy()
    if len(work) < min_islands or len(eligible_contexts) < 2:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {
            "status": "insufficient_context_replication",
            "n_islands": int(len(work)),
            "eligible_contexts": eligible_contexts,
        }
    if reference_context not in eligible_contexts:
        reference_context = eligible_contexts[0]

    scaling: dict[str, dict[str, float]] = {}
    columns: dict[str, np.ndarray] = {}
    for predictor in baseline + [isolation_column]:
        columns[f"z_{predictor}"], scaling[predictor] = _z(work[predictor])

    contexts = [c for c in eligible_contexts if c != reference_context]
    for context in contexts:
        dummy = work[context_column].eq(context).astype(float).to_numpy()
        columns[f"context[{context}]"] = dummy
        columns[f"z_isolation:context[{context}]"] = columns[f"z_{isolation_column}"] * dummy

    intercept = np.ones(len(work))
    m0_names = ["intercept", *[f"z_{p}" for p in baseline]]
    m1_names = [*m0_names, f"z_{isolation_column}"]
    m2_names = [
        *m1_names,
        *[f"context[{c}]" for c in contexts],
        *[f"z_isolation:context[{c}]" for c in contexts],
    ]

    def matrix(names: list[str]) -> np.ndarray:
        arrays = []
        for name in names:
            arrays.append(intercept if name == "intercept" else columns[name])
        return np.column_stack(arrays)

    successes = work["successes"].to_numpy(dtype=float)
    trials = work["trials"].to_numpy(dtype=float)
    clusters = work[cluster_column].to_numpy(dtype=str)

    coefficient_parts = []
    fit_rows = []
    fits: dict[str, tuple[pd.DataFrame, dict[str, Any], np.ndarray, list[str]]] = {}
    for model, names in (("M0", m0_names), ("M1", m1_names), ("M2", m2_names)):
        coef, fit, covariance = _fit_grouped_binomial_design(
            successes, trials, matrix(names), names, clusters
        )
        coef.insert(0, "model", model)
        coefficient_parts.append(coef)
        fit_rows.append({"model": model, **fit})
        fits[model] = (coef, fit, covariance, names)

    fit_table = pd.DataFrame(fit_rows)
    m1_aic = float(fit_table.loc[fit_table["model"].eq("M1"), "aic"].iloc[0])
    m2_aic = float(fit_table.loc[fit_table["model"].eq("M2"), "aic"].iloc[0])
    fit_table["delta_aic_vs_m1"] = fit_table["aic"] - m1_aic
    fit_table["m2_improvement_over_m1_aic"] = m1_aic - m2_aic

    coef_m2, _, cov_m2, names_m2 = fits["M2"]
    beta = coef_m2.set_index("predictor")["estimate_log_odds"]
    iso_name = f"z_{isolation_column}"
    iso_idx = names_m2.index(iso_name)
    slope_rows = []
    for context in eligible_contexts:
        estimate = float(beta[iso_name])
        variance = float(cov_m2[iso_idx, iso_idx])
        if context != reference_context:
            interaction = f"z_isolation:context[{context}]"
            j = names_m2.index(interaction)
            estimate += float(beta[interaction])
            variance += float(cov_m2[j, j] + 2.0 * cov_m2[iso_idx, j])
        stderr = math.sqrt(max(variance, 0.0))
        z_value = estimate / stderr if stderr > 0 else float("nan")
        slope_rows.append({
            "context": context,
            "reference_context": reference_context,
            "isolation_slope_log_odds_per_sd": estimate,
            "cluster_robust_se": stderr,
            "z_value": z_value,
            "p_value": _two_sided_p(z_value) if math.isfinite(z_value) else float("nan"),
            "n_context_islands": int(context_counts.get(context, 0)),
        })

    metadata = {
        "status": "fit",
        "reference_context": reference_context,
        "eligible_contexts": eligible_contexts,
        "scaling": scaling,
        "n_islands": int(len(work)),
    }
    return pd.concat(coefficient_parts, ignore_index=True), fit_table, pd.DataFrame(slope_rows), metadata


def fit_chapter1_context_models(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    required_comp = {"island_id", "stratum", "trait_name", "trait_state", "successes", "trials"}
    missing = required_comp - set(composition.columns)
    if missing:
        raise typer.BadParameter(f"composition missing columns: {sorted(missing)}")
    if "island_id" not in covariates.columns or covariates["island_id"].duplicated().any():
        raise typer.BadParameter("covariates must contain one unique row per island_id")

    baseline = [str(x) for x in config["baseline_covariates"]]
    isolation = str(config["isolation_column"])
    context = str(config["context_column"])
    cluster = str(config["cluster_column"])
    reference = str(config["reference_context"])
    min_islands = int(config.get("min_islands_per_fit", 30))
    min_context = int(config.get("min_islands_per_context", 8))

    joined = composition.merge(covariates, on="island_id", how="left", validate="many_to_one")
    coefficient_parts = []
    fit_parts = []
    slope_parts = []
    support_rows = []
    for (stratum, trait, state), subset in joined.groupby(
        ["stratum", "trait_name", "trait_state"], sort=True
    ):
        coef, fits, slopes, meta = _fit_one_category(
            subset,
            baseline=baseline,
            isolation_column=isolation,
            context_column=context,
            cluster_column=cluster,
            reference_context=reference,
            min_islands=min_islands,
            min_context_islands=min_context,
        )
        support_rows.append({
            "stratum": stratum,
            "trait_name": trait,
            "trait_state": state,
            "status": meta.get("status"),
            "n_islands": int(meta.get("n_islands", len(subset))),
            "eligible_contexts": "|".join(meta.get("eligible_contexts", [])),
            "reference_context": meta.get("reference_context", ""),
        })
        if coef.empty:
            continue
        for table in (coef, fits, slopes):
            table.insert(0, "trait_state", state)
            table.insert(0, "trait_name", trait)
            table.insert(0, "stratum", stratum)
        coefficient_parts.append(coef)
        fit_parts.append(fits)
        slope_parts.append(slopes)

    coefficients = pd.concat(coefficient_parts, ignore_index=True) if coefficient_parts else pd.DataFrame()
    fits = pd.concat(fit_parts, ignore_index=True) if fit_parts else pd.DataFrame()
    slopes = pd.concat(slope_parts, ignore_index=True) if slope_parts else pd.DataFrame()
    support = pd.DataFrame(support_rows)
    return coefficients, fits, slopes, support


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    required = {
        "baseline_covariates", "isolation_column", "context_column",
        "cluster_column", "reference_context"
    }
    if not isinstance(payload, dict) or not required.issubset(payload):
        raise typer.BadParameter(f"config must contain {sorted(required)}")
    forbidden = {"bombus_deficit", "bombus_channel_state", "bombus_environmental_compatibility"}
    declared = set(map(str, payload.get("baseline_covariates", []))) | {
        str(payload.get("isolation_column", "")), str(payload.get("context_column", ""))
    }
    if declared & forbidden:
        raise typer.BadParameter("Bombus predictors are prohibited from canonical Chapter 1 config")
    return payload


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/chapter1_context_analysis.yml"), exists=True),
) -> None:
    config = _load_config(config_path)
    coefficients, fits, slopes, support = fit_chapter1_context_models(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output_dir / "chapter1_m0_m2_coefficients.csv", index=False)
    fits.to_csv(output_dir / "chapter1_m0_m2_fit_stats.csv", index=False)
    slopes.to_csv(output_dir / "chapter1_context_isolation_slopes.csv", index=False)
    support.to_csv(output_dir / "chapter1_model_support.csv", index=False)
    # M4 is the category-preserving view of the same M2 simple slopes.
    slopes.to_csv(output_dir / "chapter1_m4_category_decomposition.csv", index=False)
    manifest = {
        "contract": "chapter1_context_dependent_models_v1",
        "model_ladder": ["M0_baseline", "M1_universal_isolation", "M2_isolation_x_context"],
        "m3": "separate genus-fixed status/lineage residual layer",
        "m4": "category-preserving decomposition of M2 context slopes",
        "pollinator_primary_predictors": False,
        "config": config,
        "n_fitted_categories": int(slopes[["stratum", "trait_name", "trait_state"]].drop_duplicates().shape[0]) if not slopes.empty else 0,
    }
    (output_dir / "chapter1_context_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
