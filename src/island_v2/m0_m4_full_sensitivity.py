"""Sensitivity analyses for the locked full M0-M4 island path models.

This module answers two reviewer-facing questions without acquiring new data:
1. Is the Bombus environmental-compatibility association only a linear restatement
   of area, latitude, climate PCs, and trait richness?
2. Are the key path estimates stable when geographic 10-degree blocks, rather than
   individual islands, are resampled?

The Bombus variable remains environmental compatibility, never observed absence.
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

from island_v2.m0_m4_full_analysis import fit_grouped_binomial_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("full-analysis config must be a mapping")
    return payload


def _numeric_complete(data: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    work = data[columns].copy()
    for column in columns:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    return work.dropna()


def residualize_bombus(
    data: pd.DataFrame,
    baseline: list[str],
    bombus_column: str = "bombus_channel_state",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Regress Bombus compatibility on the frozen baseline and retain raw residuals."""
    columns = [bombus_column, *baseline]
    complete = _numeric_complete(data, columns)
    if len(complete) <= len(baseline) + 2:
        raise typer.BadParameter("too few complete islands to residualize Bombus compatibility")

    y = complete[bombus_column].to_numpy(dtype=float)
    X = np.column_stack(
        [np.ones(len(complete)), *[complete[column].to_numpy(dtype=float) for column in baseline]]
    )
    beta = np.linalg.pinv(X.T @ X) @ (X.T @ y)
    fitted = X @ beta
    residuals = y - fitted
    sst = float(np.sum((y - y.mean()) ** 2))
    sse = float(np.sum(residuals**2))
    r2 = float(1.0 - sse / sst) if sst > 0 else np.nan

    result = data.copy()
    result["bombus_residualized"] = np.nan
    result.loc[complete.index, "bombus_residualized"] = residuals
    result["bombus_residualized_x_alternative"] = (
        pd.to_numeric(result["bombus_residualized"], errors="coerce")
        * pd.to_numeric(result["alternative_pollen_vector"], errors="coerce")
    )

    metadata = {
        "n_complete": int(len(complete)),
        "baseline_predictors": baseline,
        "bombus_column": bombus_column,
        "baseline_r2_for_bombus": r2,
        "residual_mean": float(np.mean(residuals)),
        "residual_sd": float(np.std(residuals, ddof=0)),
        "coefficients": {
            "intercept": float(beta[0]),
            **{column: float(value) for column, value in zip(baseline, beta[1:], strict=True)},
        },
    }
    return result, metadata


def collinearity_diagnostics(
    data: pd.DataFrame,
    predictors: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Return complete-case correlations, VIFs, and standardized design condition number."""
    complete = _numeric_complete(data, predictors)
    if len(complete) <= len(predictors) + 1:
        raise typer.BadParameter("too few complete islands for collinearity diagnostics")

    means = complete.mean(axis=0)
    sds = complete.std(axis=0, ddof=0)
    nonconstant = sds[(sds > 0) & np.isfinite(sds)].index.tolist()
    dropped = [column for column in predictors if column not in nonconstant]
    if len(nonconstant) < 2:
        raise typer.BadParameter("fewer than two nonconstant predictors for diagnostics")

    z = (complete[nonconstant] - means[nonconstant]) / sds[nonconstant]
    correlation = z.corr()
    corr_matrix = correlation.to_numpy(dtype=float)
    vif_values = np.diag(np.linalg.pinv(corr_matrix))
    singular_values = np.linalg.svd(z.to_numpy(dtype=float), compute_uv=False)
    positive = singular_values[singular_values > np.finfo(float).eps]
    condition_number = float(positive.max() / positive.min()) if len(positive) else np.nan

    correlation_long = (
        correlation.rename_axis(index="predictor_a", columns="predictor_b")
        .stack()
        .rename("correlation")
        .reset_index()
    )
    vif = pd.DataFrame(
        {
            "predictor": nonconstant,
            "vif": [float(value) for value in vif_values],
        }
    ).sort_values("vif", ascending=False, ignore_index=True)
    summary = {
        "n_complete": int(len(complete)),
        "predictors_requested": predictors,
        "predictors_used": nonconstant,
        "constant_predictors_dropped": dropped,
        "condition_number_standardized_design": condition_number,
        "max_vif": float(vif["vif"].max()),
        "max_absolute_pairwise_correlation": float(
            correlation.where(~np.eye(len(correlation), dtype=bool)).abs().max().max()
        ),
    }
    return correlation_long, vif, summary


def _coefficient(
    coefficients: pd.DataFrame,
    predictor: str,
) -> float:
    row = coefficients.loc[coefficients["predictor"].eq(predictor)]
    return float(row.iloc[0]["estimate_log_odds"]) if len(row) else np.nan


def _fit_estimate(
    data: pd.DataFrame,
    successes: str,
    trials: str,
    predictors: list[str],
    target_predictor: str,
    z_threshold: float,
) -> float:
    coefficients, fit = fit_grouped_binomial_clustered(
        data,
        successes,
        trials,
        predictors,
        "spatial_block",
        z_threshold,
    )
    if fit.get("status") not in {"fitted", "max_iter_reached"} or coefficients.empty:
        return np.nan
    return _coefficient(coefficients, target_predictor)


def fit_sensitivity_snapshot(
    data: pd.DataFrame,
    baseline: list[str],
    z_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit raw and baseline-residualized versions of the key M1/M3/M4 paths."""
    specifications = [
        (
            "M1_raw_bombus",
            "floral_phenotype_successes",
            "floral_phenotype_trials",
            [*baseline, "bombus_channel_state"],
            "bombus_channel_state",
        ),
        (
            "M1_residualized_bombus",
            "floral_phenotype_successes",
            "floral_phenotype_trials",
            [*baseline, "bombus_residualized"],
            "bombus_residualized",
        ),
        (
            "M3_raw_direct_bombus",
            "floral_phenotype_successes",
            "floral_phenotype_trials",
            [*baseline, "bombus_channel_state", "mixed_mating_share", "self_compatibility_share"],
            "bombus_channel_state",
        ),
        (
            "M3_residualized_direct_bombus",
            "floral_phenotype_successes",
            "floral_phenotype_trials",
            [*baseline, "bombus_residualized", "mixed_mating_share", "self_compatibility_share"],
            "bombus_residualized",
        ),
        (
            "M4_raw_interaction",
            "floral_phenotype_successes",
            "floral_phenotype_trials",
            [*baseline, "bombus_channel_state", "alternative_pollen_vector", "bombus_x_alternative"],
            "bombus_x_alternative",
        ),
        (
            "M4_residualized_interaction",
            "floral_phenotype_successes",
            "floral_phenotype_trials",
            [
                *baseline,
                "bombus_residualized",
                "alternative_pollen_vector",
                "bombus_residualized_x_alternative",
            ],
            "bombus_residualized_x_alternative",
        ),
    ]

    coefficient_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    for analysis, successes, trials, predictors, target in specifications:
        coefficients, fit = fit_grouped_binomial_clustered(
            data, successes, trials, predictors, "spatial_block", z_threshold
        )
        if not coefficients.empty:
            coefficients.insert(0, "analysis", analysis)
            coefficients["target_predictor"] = target
            coefficient_parts.append(coefficients)
        fit_rows.append(
            {
                "analysis": analysis,
                "target_predictor": target,
                "predictors": " + ".join(predictors),
                **fit,
            }
        )
    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True) if coefficient_parts else pd.DataFrame()
    )
    return coefficients, pd.DataFrame(fit_rows)


def _indirect_effects(
    data: pd.DataFrame,
    baseline: list[str],
    bombus_predictor: str,
    z_threshold: float,
) -> dict[str, float]:
    a_mixed = _fit_estimate(
        data,
        "mixed_mating_successes",
        "mixed_mating_trials",
        [*baseline, bombus_predictor],
        bombus_predictor,
        z_threshold,
    )
    a_sc = _fit_estimate(
        data,
        "self_compatibility_successes",
        "self_compatibility_trials",
        [*baseline, bombus_predictor],
        bombus_predictor,
        z_threshold,
    )
    b_predictors = [
        *baseline,
        bombus_predictor,
        "mixed_mating_share",
        "self_compatibility_share",
    ]
    b_mixed = _fit_estimate(
        data,
        "floral_phenotype_successes",
        "floral_phenotype_trials",
        b_predictors,
        "mixed_mating_share",
        z_threshold,
    )
    b_sc = _fit_estimate(
        data,
        "floral_phenotype_successes",
        "floral_phenotype_trials",
        b_predictors,
        "self_compatibility_share",
        z_threshold,
    )
    mixed = float(a_mixed * b_mixed) if math.isfinite(a_mixed) and math.isfinite(b_mixed) else np.nan
    sc = float(a_sc * b_sc) if math.isfinite(a_sc) and math.isfinite(b_sc) else np.nan
    total = float(mixed + sc) if math.isfinite(mixed) and math.isfinite(sc) else np.nan
    return {
        "mixed_mating_indirect": mixed,
        "self_compatibility_indirect": sc,
        "parallel_total_indirect": total,
    }


def _sample_spatial_blocks(
    data: pd.DataFrame,
    rng: np.random.Generator,
) -> pd.DataFrame:
    blocks = sorted(data["spatial_block"].dropna().astype(str).unique())
    if len(blocks) < 2:
        raise typer.BadParameter("fewer than two spatial blocks for bootstrap")
    draws = rng.choice(blocks, size=len(blocks), replace=True)
    parts: list[pd.DataFrame] = []
    for draw_index, block in enumerate(draws):
        part = data.loc[data["spatial_block"].astype(str).eq(str(block))].copy()
        if part.empty:
            continue
        # Duplicate block draws must remain distinct bootstrap clusters.
        part["spatial_block"] = f"bootstrap_draw_{draw_index:04d}"
        part["island_id"] = part["island_id"].astype(str) + f"__boot{draw_index:04d}"
        parts.append(part)
    return pd.concat(parts, ignore_index=True) if parts else data.head(0).copy()


def spatial_block_bootstrap(
    data: pd.DataFrame,
    baseline: list[str],
    replicates: int,
    seed: int,
    z_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resample geographic blocks and refit key raw/residualized paths."""
    if replicates < 20:
        raise typer.BadParameter("bootstrap replicates must be at least 20")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []

    for replicate in range(1, replicates + 1):
        sample = _sample_spatial_blocks(data, rng)
        metrics = {
            "M1_raw_bombus": _fit_estimate(
                sample,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [*baseline, "bombus_channel_state"],
                "bombus_channel_state",
                z_threshold,
            ),
            "M1_residualized_bombus": _fit_estimate(
                sample,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [*baseline, "bombus_residualized"],
                "bombus_residualized",
                z_threshold,
            ),
            "M3_raw_direct_bombus": _fit_estimate(
                sample,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [
                    *baseline,
                    "bombus_channel_state",
                    "mixed_mating_share",
                    "self_compatibility_share",
                ],
                "bombus_channel_state",
                z_threshold,
            ),
            "M3_residualized_direct_bombus": _fit_estimate(
                sample,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [
                    *baseline,
                    "bombus_residualized",
                    "mixed_mating_share",
                    "self_compatibility_share",
                ],
                "bombus_residualized",
                z_threshold,
            ),
            "M4_raw_interaction": _fit_estimate(
                sample,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [
                    *baseline,
                    "bombus_channel_state",
                    "alternative_pollen_vector",
                    "bombus_x_alternative",
                ],
                "bombus_x_alternative",
                z_threshold,
            ),
            "M4_residualized_interaction": _fit_estimate(
                sample,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                [
                    *baseline,
                    "bombus_residualized",
                    "alternative_pollen_vector",
                    "bombus_residualized_x_alternative",
                ],
                "bombus_residualized_x_alternative",
                z_threshold,
            ),
        }
        metrics.update(
            {
                f"M3_raw_{key}": value
                for key, value in _indirect_effects(
                    sample, baseline, "bombus_channel_state", z_threshold
                ).items()
            }
        )
        metrics.update(
            {
                f"M3_residualized_{key}": value
                for key, value in _indirect_effects(
                    sample, baseline, "bombus_residualized", z_threshold
                ).items()
            }
        )
        for metric, estimate in metrics.items():
            rows.append(
                {
                    "replicate": replicate,
                    "metric": metric,
                    "estimate": estimate,
                    "n_rows": int(len(sample)),
                    "n_bootstrap_blocks": int(sample["spatial_block"].nunique()),
                }
            )

    draws = pd.DataFrame(rows)
    summaries: list[dict[str, Any]] = []
    for metric, group in draws.groupby("metric", sort=True):
        values = pd.to_numeric(group["estimate"], errors="coerce").dropna().to_numpy(dtype=float)
        if not len(values):
            summaries.append(
                {
                    "metric": metric,
                    "n_valid_replicates": 0,
                    "bootstrap_median": np.nan,
                    "ci_2_5": np.nan,
                    "ci_97_5": np.nan,
                    "proportion_positive": np.nan,
                    "ci_excludes_zero": False,
                }
            )
            continue
        lower, median, upper = np.quantile(values, [0.025, 0.5, 0.975])
        summaries.append(
            {
                "metric": metric,
                "n_valid_replicates": int(len(values)),
                "bootstrap_median": float(median),
                "ci_2_5": float(lower),
                "ci_97_5": float(upper),
                "proportion_positive": float(np.mean(values > 0)),
                "ci_excludes_zero": bool(lower > 0 or upper < 0),
            }
        )
    return draws, pd.DataFrame(summaries)


@app.command("run")
def run(
    island_data_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    analysis_tier: str = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_full_analysis.yml"), exists=True),
    bootstrap_replicates: int = typer.Option(200, min=20),
    bootstrap_seed: int = typer.Option(20260713),
) -> None:
    config = _load_config(config_path)
    data = pd.read_csv(island_data_csv)
    baseline = [str(value) for value in config["baseline_covariates"]]
    missing = [column for column in [*baseline, "bombus_channel_state", "spatial_block"] if column not in data.columns]
    if missing:
        raise typer.BadParameter(f"full-analysis island data missing columns: {missing}")

    residualized, residualization = residualize_bombus(data, baseline)
    diagnostic_predictors = [*baseline, "bombus_channel_state"]
    correlations, vif, collinearity = collinearity_diagnostics(
        residualized, diagnostic_predictors
    )
    z_threshold = float(config["reporting"].get("nominal_z_threshold", 1.96))
    snapshot_coefficients, snapshot_fit = fit_sensitivity_snapshot(
        residualized, baseline, z_threshold
    )
    bootstrap_draws, bootstrap_summary = spatial_block_bootstrap(
        residualized,
        baseline,
        replicates=bootstrap_replicates,
        seed=bootstrap_seed,
        z_threshold=z_threshold,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    residualized.to_csv(output_dir / "sensitivity_island_data.csv", index=False)
    correlations.to_csv(output_dir / "predictor_correlations.csv", index=False)
    vif.to_csv(output_dir / "predictor_vif.csv", index=False)
    snapshot_coefficients.to_csv(output_dir / "sensitivity_path_coefficients.csv", index=False)
    snapshot_fit.to_csv(output_dir / "sensitivity_path_fit.csv", index=False)
    bootstrap_draws.to_csv(output_dir / "spatial_block_bootstrap_draws.csv.gz", index=False, compression="gzip")
    bootstrap_summary.to_csv(output_dir / "spatial_block_bootstrap_summary.csv", index=False)

    manifest = {
        "contract": "m0_m4_full_sensitivity_v1",
        "analysis_tier": analysis_tier,
        "n_islands_input": int(data["island_id"].nunique()),
        "n_spatial_blocks_input": int(data["spatial_block"].nunique()),
        "residualization": residualization,
        "collinearity": collinearity,
        "bootstrap": {
            "unit": "10_degree_spatial_block",
            "replicates": bootstrap_replicates,
            "seed": bootstrap_seed,
            "n_metrics": int(bootstrap_summary["metric"].nunique()),
        },
        "interpretation": (
            "Residualized Bombus sensitivity tests incremental association beyond the frozen linear baseline. "
            "Spatial block bootstrap resamples geographic blocks, not individual islands."
        ),
    }
    (output_dir / "full_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
