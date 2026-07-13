"""Reviewer-grade sensitivity analysis with within-bootstrap Bombus residualization.

Each geographic block-bootstrap replicate re-fits the linear baseline model for
Bombus environmental compatibility before the residualized path models are fit.
This propagates uncertainty from the residualization step instead of treating a
single set of residuals as fixed.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer

from island_v2.m0_m4_full_sensitivity import (
    _fit_estimate,
    _indirect_effects,
    _load_config,
    _sample_spatial_blocks,
    collinearity_diagnostics,
    fit_sensitivity_snapshot,
    residualize_bombus,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def spatial_block_bootstrap_refit_residualization(
    data: pd.DataFrame,
    baseline: list[str],
    replicates: int,
    seed: int,
    z_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Resample spatial blocks and re-fit residualization inside every replicate."""
    if replicates < 20:
        raise typer.BadParameter("bootstrap replicates must be at least 20")
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []

    for replicate in range(1, replicates + 1):
        sampled = _sample_spatial_blocks(data, rng)
        sample, residualization = residualize_bombus(sampled, baseline)
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
                    "bootstrap_bombus_baseline_r2": float(
                        residualization["baseline_r2_for_bombus"]
                    ),
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


def _validate_bootstrap(
    summary: pd.DataFrame,
    expected_replicates: int,
) -> None:
    if summary.empty:
        raise RuntimeError("bootstrap summary is empty")
    minimum_valid = int(math.floor(expected_replicates * 0.90))
    if int(summary["n_valid_replicates"].min()) < minimum_valid:
        failed = summary.loc[summary["n_valid_replicates"] < minimum_valid]
        raise RuntimeError(
            "fewer than 90% valid bootstrap replicates for one or more metrics: "
            f"{failed[['metric', 'n_valid_replicates']].to_dict('records')}"
        )


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
    required = [*baseline, "bombus_channel_state", "spatial_block"]
    missing = [column for column in required if column not in data.columns]
    if missing:
        raise typer.BadParameter(f"full-analysis island data missing columns: {missing}")

    residualized, residualization = residualize_bombus(data, baseline)
    correlations, vif, collinearity = collinearity_diagnostics(
        residualized, [*baseline, "bombus_channel_state"]
    )
    z_threshold = float(config["reporting"].get("nominal_z_threshold", 1.96))
    snapshot_coefficients, snapshot_fit = fit_sensitivity_snapshot(
        residualized, baseline, z_threshold
    )
    bootstrap_draws, bootstrap_summary = spatial_block_bootstrap_refit_residualization(
        data,
        baseline,
        replicates=bootstrap_replicates,
        seed=bootstrap_seed,
        z_threshold=z_threshold,
    )
    _validate_bootstrap(bootstrap_summary, bootstrap_replicates)

    output_dir.mkdir(parents=True, exist_ok=True)
    residualized.to_csv(output_dir / "sensitivity_island_data.csv", index=False)
    correlations.to_csv(output_dir / "predictor_correlations.csv", index=False)
    vif.to_csv(output_dir / "predictor_vif.csv", index=False)
    snapshot_coefficients.to_csv(output_dir / "sensitivity_path_coefficients.csv", index=False)
    snapshot_fit.to_csv(output_dir / "sensitivity_path_fit.csv", index=False)
    bootstrap_draws.to_csv(
        output_dir / "spatial_block_bootstrap_draws.csv.gz",
        index=False,
        compression="gzip",
    )
    bootstrap_summary.to_csv(
        output_dir / "spatial_block_bootstrap_summary.csv", index=False
    )

    bootstrap_r2 = (
        bootstrap_draws[["replicate", "bootstrap_bombus_baseline_r2"]]
        .drop_duplicates("replicate")["bootstrap_bombus_baseline_r2"]
        .dropna()
    )
    manifest = {
        "contract": "m0_m4_full_sensitivity_v2_refit_residualization",
        "analysis_tier": analysis_tier,
        "n_islands_input": int(data["island_id"].nunique()),
        "n_spatial_blocks_input": int(data["spatial_block"].nunique()),
        "residualization": residualization,
        "collinearity": collinearity,
        "bootstrap": {
            "unit": "10_degree_spatial_block",
            "replicates": bootstrap_replicates,
            "seed": bootstrap_seed,
            "residualization_refit_within_each_replicate": True,
            "minimum_valid_fraction_required": 0.90,
            "n_metrics": int(bootstrap_summary["metric"].nunique()),
            "bombus_baseline_r2_bootstrap_median": float(bootstrap_r2.median()),
            "bombus_baseline_r2_bootstrap_2_5": float(bootstrap_r2.quantile(0.025)),
            "bombus_baseline_r2_bootstrap_97_5": float(bootstrap_r2.quantile(0.975)),
        },
        "interpretation": (
            "Residualized Bombus sensitivity tests incremental association beyond the frozen linear baseline. "
            "Every spatial block-bootstrap replicate re-fits the residualization model before path fitting."
        ),
    }
    (output_dir / "full_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
