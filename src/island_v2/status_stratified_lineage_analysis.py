"""Fit isolation slopes after a genus-composition-preserving trait null.

This module is a lineage sensitivity layer. It consumes island-level output from
``genus_fixed_trait_null`` and asks whether the observed-minus-genus-null trait
composition still covaries with isolation after area and climate are represented.

It does not impute missing trait states and does not treat genus-null residuals as
causal evolutionary effects. Native non-endemic and corroborated-endemic strata
are fitted separately so lineage turnover and endemicity can be separated before
pollination-channel terms are interpreted.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

DEFAULT_BASELINE = [
    "log_distance_to_continent_km",
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
]
DEFAULT_STRATA = ("all_native", "native_nonendemic", "endemic")
DEFAULT_REGIMES = (
    "northern_midlatitude",
    "tropical",
    "southern_extratropical",
)


def _normal_two_sided_p(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0))


def _prepare_design(
    data: pd.DataFrame, predictors: list[str]
) -> tuple[np.ndarray, list[str], dict[str, dict[str, float]]]:
    columns: list[np.ndarray] = []
    scaling: dict[str, dict[str, float]] = {}
    for predictor in predictors:
        values = pd.to_numeric(data[predictor], errors="coerce").to_numpy(dtype=float)
        mean = float(np.mean(values))
        sd = float(np.std(values, ddof=0))
        if not math.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant or invalid predictor: {predictor}")
        columns.append((values - mean) / sd)
        scaling[predictor] = {"mean": mean, "sd": sd}
    X = np.column_stack([np.ones(len(data)), *columns])
    return X, ["intercept", *predictors], scaling


def fit_weighted_linear_clustered(
    data: pd.DataFrame,
    *,
    response_column: str,
    weight_column: str,
    predictors: list[str],
    cluster_column: str,
    z_threshold: float = 1.96,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Weighted least squares with a cluster-robust sandwich covariance."""
    required = {response_column, weight_column, cluster_column, *predictors}
    missing = required - set(data.columns)
    if missing:
        raise typer.BadParameter(f"lineage analysis missing columns: {sorted(missing)}")

    work = data[list(required)].copy()
    for column in [response_column, weight_column, *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    work = work.loc[work[weight_column] > 0].copy()
    n = int(len(work))
    p = len(predictors) + 1
    n_clusters = int(work[cluster_column].nunique())
    if n < max(10, p + 3) or n_clusters < 2:
        return pd.DataFrame(), {
            "status": "insufficient_complete_islands",
            "n_islands": n,
            "n_clusters": n_clusters,
        }

    try:
        X, names, scaling = _prepare_design(work, predictors)
    except ValueError as exc:
        return pd.DataFrame(), {
            "status": "constant_or_invalid_predictor",
            "n_islands": n,
            "n_clusters": n_clusters,
            "error": str(exc),
        }

    y = work[response_column].to_numpy(dtype=float)
    w = work[weight_column].to_numpy(dtype=float)
    xtwx = X.T @ (w[:, None] * X)
    bread = np.linalg.pinv(xtwx)
    beta = bread @ (X.T @ (w * y))
    residual = y - X @ beta

    meat = np.zeros((p, p), dtype=float)
    clusters = work[cluster_column].astype(str).to_numpy()
    for cluster in np.unique(clusters):
        mask = clusters == cluster
        score = X[mask].T @ (w[mask] * residual[mask])
        meat += np.outer(score, score)
    covariance = bread @ meat @ bread
    if n_clusters > 1 and n > p:
        covariance *= (n_clusters / (n_clusters - 1.0)) * ((n - 1.0) / (n - p))
    se = np.sqrt(np.clip(np.diag(covariance), 0.0, None))

    rows: list[dict[str, Any]] = []
    for name, estimate, stderr in zip(names, beta, se, strict=True):
        z = float(estimate / stderr) if stderr > 0 else float("nan")
        p_value = _normal_two_sided_p(z) if math.isfinite(z) else float("nan")
        rows.append(
            {
                "predictor": name,
                "estimate": float(estimate),
                "cluster_robust_se": float(stderr),
                "z_value": z,
                "p_value": p_value,
                "nominally_supported": bool(math.isfinite(z) and abs(z) >= z_threshold),
                "n_islands": n,
                "n_clusters": n_clusters,
            }
        )
    return pd.DataFrame(rows), {
        "status": "fit",
        "n_islands": n,
        "n_clusters": n_clusters,
        "scaling": scaling,
        "weight_column": weight_column,
    }


def fit_status_stratified_lineage_models(
    genus_null: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    strata: tuple[str, ...] = DEFAULT_STRATA,
    regimes: tuple[str, ...] = DEFAULT_REGIMES,
    baseline: list[str] | None = None,
    cluster_column: str = "spatial_block",
    pilot_min_islands: int = 30,
    confirmatory_min_islands: int = 50,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_null = {
        "island_id",
        "outcome",
        "stratum",
        "observed_n_species",
        "deviation_observed_minus_null",
    }
    missing = required_null - set(genus_null.columns)
    if missing:
        raise typer.BadParameter(f"genus-null table missing columns: {sorted(missing)}")
    predictors = list(baseline or DEFAULT_BASELINE)
    required_cov = {"island_id", "analysis_regime", cluster_column, *predictors}
    missing_cov = required_cov - set(covariates.columns)
    if missing_cov:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing_cov)}")

    data = genus_null.merge(
        covariates[list(required_cov)].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    coefficient_tables: list[pd.DataFrame] = []
    support_rows: list[dict[str, Any]] = []
    for stratum in strata:
        for regime in regimes:
            outcomes = sorted(
                data.loc[data["stratum"].eq(stratum), "outcome"].dropna().astype(str).unique()
            )
            for outcome in outcomes:
                subset = data.loc[
                    data["stratum"].eq(stratum)
                    & data["analysis_regime"].eq(regime)
                    & data["outcome"].eq(outcome)
                ].copy()
                n_islands = int(len(subset))
                if n_islands >= confirmatory_min_islands:
                    support_class = "confirmatory_count_met"
                elif n_islands >= pilot_min_islands:
                    support_class = "targeted_acquisition_zone"
                else:
                    support_class = "below_pilot_support"
                support_rows.append(
                    {
                        "stratum": stratum,
                        "regime": regime,
                        "outcome": outcome,
                        "n_islands": n_islands,
                        "n_spatial_blocks": int(subset[cluster_column].dropna().nunique()),
                        "median_direct_trials": (
                            float(subset["observed_n_species"].median()) if n_islands else np.nan
                        ),
                        "support_class": support_class,
                    }
                )
                coefficients, fit = fit_weighted_linear_clustered(
                    subset,
                    response_column="deviation_observed_minus_null",
                    weight_column="observed_n_species",
                    predictors=predictors,
                    cluster_column=cluster_column,
                )
                if coefficients.empty:
                    continue
                coefficients.insert(0, "outcome", outcome)
                coefficients.insert(0, "regime", regime)
                coefficients.insert(0, "stratum", stratum)
                coefficients["support_class"] = support_class
                coefficients["fit_status"] = str(fit.get("status", ""))
                coefficient_tables.append(coefficients)

    coefficients = (
        pd.concat(coefficient_tables, ignore_index=True)
        if coefficient_tables
        else pd.DataFrame()
    )
    return coefficients, pd.DataFrame(support_rows)


@app.command("run")
def run(
    genus_null_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    null = pd.read_csv(genus_null_csv)
    covariates = pd.read_csv(covariates_csv)
    coefficients, support = fit_status_stratified_lineage_models(null, covariates)
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output_dir / "status_stratified_lineage_coefficients.csv", index=False)
    support.to_csv(output_dir / "status_stratified_lineage_support.csv", index=False)
    isolation = (
        coefficients.loc[coefficients["predictor"].eq("log_distance_to_continent_km")].copy()
        if not coefficients.empty
        else coefficients.copy()
    )
    isolation.to_csv(output_dir / "status_stratified_lineage_isolation.csv", index=False)
    manifest = {
        "contract": "status_stratified_genus_fixed_lineage_v1",
        "response": "observed trait share minus genus-fixed null mean",
        "weights": "direct species trials per island",
        "strata": list(DEFAULT_STRATA),
        "regimes": list(DEFAULT_REGIMES),
        "baseline": DEFAULT_BASELINE,
        "interpretation": (
            "Lineage sensitivity only. A residual slope indicates information beyond the "
            "observed genus composition under this null; it does not identify a causal "
            "pollination mechanism."
        ),
    }
    (output_dir / "status_stratified_lineage_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
