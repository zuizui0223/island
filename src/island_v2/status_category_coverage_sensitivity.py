"""Stress-test supported category slopes against direct-trait coverage imbalance.

Only category contrasts that already survive BH-FDR with >=50 islands enter this
sensitivity. The analysis does not fill missing traits. It measures each island's
direct-trait fraction within the relevant floristic stratum, then checks whether
the isolation coefficient survives (1) continuous coverage adjustment and (2)
prespecified higher-coverage subsets.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.status_stratified_lineage_analysis import (
    DEFAULT_BASELINE,
    fit_weighted_linear_clustered,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)
THRESHOLDS = (0.25, 0.30, 0.40, 0.50)


def _text(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip()


def _stratum_mask(status: pd.DataFrame, stratum: str) -> pd.Series:
    origin = _text(status["origin_status"]).str.lower()
    endemic = _text(status["endemic_status"]).str.lower()
    if stratum == "all_native":
        return origin.eq("native")
    if stratum == "native_nonendemic":
        return origin.eq("native") & endemic.eq("nonendemic")
    if stratum == "endemic":
        return origin.eq("native") & endemic.eq("endemic")
    raise typer.BadParameter(f"unsupported stratum: {stratum}")


def _distance_row(
    data: pd.DataFrame,
    predictors: list[str],
    *,
    model: str,
    cluster: str,
) -> dict[str, Any]:
    coefficients, fit = fit_weighted_linear_clustered(
        data,
        response_column="deviation_observed_minus_genus",
        weight_column="trials",
        predictors=predictors,
        cluster_column=cluster,
    )
    if coefficients.empty:
        return {
            "model": model,
            "fit_status": str(fit.get("status", "")),
            "n_islands": int(fit.get("n_islands", 0)),
            "n_clusters": int(fit.get("n_clusters", 0)),
        }
    row = coefficients.loc[
        coefficients["predictor"].eq("log_distance_to_continent_km")
    ].iloc[0]
    result = {
        "model": model,
        "fit_status": str(fit.get("status", "fit")),
        "n_islands": int(row["n_islands"]),
        "n_clusters": int(row["n_clusters"]),
        "distance_beta": float(row["estimate"]),
        "distance_se": float(row["cluster_robust_se"]),
        "distance_p": float(row["p_value"]),
    }
    iso = pd.to_numeric(data["log_distance_to_continent_km"], errors="coerce").dropna()
    coverage = pd.to_numeric(data["direct_fraction"], errors="coerce").dropna()
    if len(iso):
        q05, q50, q95 = np.quantile(iso, [0.05, 0.5, 0.95])
        result.update(
            isolation_q05=float(q05), isolation_median=float(q50), isolation_q95=float(q95)
        )
    if len(coverage):
        result["median_direct_fraction"] = float(coverage.median())
    return result


def run_coverage_sensitivity(
    residuals: pd.DataFrame,
    category_fdr: pd.DataFrame,
    status_ledger: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    cluster: str = "spatial_block",
    thresholds: tuple[float, ...] = THRESHOLDS,
) -> pd.DataFrame:
    required_fdr = {
        "stratum",
        "regime",
        "domain",
        "category",
        "fdr_supported",
        "n_islands",
    }
    missing = required_fdr - set(category_fdr.columns)
    if missing:
        raise typer.BadParameter(f"category FDR table missing columns: {sorted(missing)}")
    focus = category_fdr.loc[
        category_fdr["fdr_supported"].astype(bool)
        & pd.to_numeric(category_fdr["n_islands"], errors="coerce").ge(50)
    ].drop_duplicates(["stratum", "regime", "domain", "category"])
    if focus.empty:
        return pd.DataFrame()

    required_status = {"island_id", "origin_status", "endemic_status"}
    missing = required_status - set(status_ledger.columns)
    if missing:
        raise typer.BadParameter(f"status ledger missing columns: {sorted(missing)}")
    status = status_ledger.copy()
    status["island_id"] = _text(status["island_id"])
    required_cov = {"island_id", "analysis_regime", cluster, *DEFAULT_BASELINE}
    missing = required_cov - set(covariates.columns)
    if missing:
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")

    rows: list[dict[str, Any]] = []
    for key in focus.itertuples(index=False):
        totals = (
            status.loc[_stratum_mask(status, str(key.stratum))]
            .groupby("island_id")
            .size()
            .rename("n_stratum_species")
            .reset_index()
        )
        data = residuals.loc[
            residuals["stratum"].eq(key.stratum)
            & residuals["domain"].eq(key.domain)
            & residuals["category"].eq(key.category)
        ].copy()
        data = data.merge(totals, on="island_id", how="left", validate="one_to_one")
        data = data.merge(
            covariates[list(required_cov)].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="one_to_one",
        )
        data = data.loc[data["analysis_regime"].eq(key.regime)].copy()
        data["direct_fraction"] = (
            pd.to_numeric(data["trials"], errors="coerce")
            / pd.to_numeric(data["n_stratum_species"], errors="coerce").replace(0, np.nan)
        )

        base = _distance_row(data, list(DEFAULT_BASELINE), model="baseline", cluster=cluster)
        adjusted = _distance_row(
            data,
            [*DEFAULT_BASELINE, "direct_fraction"],
            model="continuous_coverage_adjusted",
            cluster=cluster,
        )
        for result in (base, adjusted):
            rows.append(
                {
                    "stratum": key.stratum,
                    "regime": key.regime,
                    "domain": key.domain,
                    "category": key.category,
                    **result,
                }
            )
        for threshold in thresholds:
            subset = data.loc[data["direct_fraction"].ge(float(threshold))].copy()
            if len(subset) < 30:
                continue
            result = _distance_row(
                subset,
                list(DEFAULT_BASELINE),
                model=f"direct_fraction_ge_{threshold:.2f}",
                cluster=cluster,
            )
            rows.append(
                {
                    "stratum": key.stratum,
                    "regime": key.regime,
                    "domain": key.domain,
                    "category": key.category,
                    **result,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["same_direction_as_baseline"] = False
        for keys, idx in out.groupby(["stratum", "regime", "domain", "category"]).groups.items():
            group = out.loc[idx]
            baseline = group.loc[group["model"].eq("baseline"), "distance_beta"]
            if baseline.empty or not math.isfinite(float(baseline.iloc[0])):
                continue
            sign = np.sign(float(baseline.iloc[0]))
            out.loc[idx, "same_direction_as_baseline"] = (
                np.sign(pd.to_numeric(group["distance_beta"], errors="coerce")) == sign
            ).to_numpy()
    return out


@app.command("run")
def run(
    residuals_csv: Path = typer.Option(..., exists=True),
    category_fdr_csv: Path = typer.Option(..., exists=True),
    status_ledger_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_csv: Path = typer.Option(...),
) -> None:
    result = run_coverage_sensitivity(
        pd.read_csv(residuals_csv),
        pd.read_csv(category_fdr_csv),
        pd.read_csv(status_ledger_csv),
        pd.read_csv(covariates_csv),
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)


if __name__ == "__main__":
    app()
