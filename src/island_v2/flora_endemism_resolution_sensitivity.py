"""Sensitivity of the isolation-endemism model to status-resolution quality."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.flora_endemism_analysis import prepare_endemism_response
from island_v2.m0_m4_full_analysis import fit_grouped_binomial_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)

DEFAULT_THRESHOLDS = (0.0, 0.70, 0.80, 0.90, 0.95)


def run_resolution_sensitivity(
    raw: pd.DataFrame,
    config: dict[str, Any],
    thresholds: tuple[float, ...] = DEFAULT_THRESHOLDS,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    prepared = prepare_endemism_response(raw, config)
    regime_column = str(config["regime_column"])
    cluster = str(config["spatial_cluster"])
    baseline = [str(value) for value in config["baseline_covariates"]]
    z_threshold = float(config.get("nominal_z_threshold", 1.96))

    coefficient_tables: list[pd.DataFrame] = []
    support_rows: list[dict[str, Any]] = []
    for threshold in thresholds:
        subset_threshold = prepared.loc[
            prepared["endemism_support_eligible"]
            & prepared["endemism_resolution_fraction"].ge(float(threshold))
        ].copy()
        for regime in [str(value) for value in config["regimes"]]:
            subset = subset_threshold.loc[
                subset_threshold[regime_column].astype(str).eq(regime)
            ].copy()
            support_rows.append(
                {
                    "minimum_endemism_resolution_fraction": float(threshold),
                    "regime": regime,
                    "n_islands": int(len(subset)),
                    "n_spatial_blocks": int(subset[cluster].dropna().nunique()),
                    "median_endemism_resolution_fraction": (
                        float(subset["endemism_resolution_fraction"].median())
                        if len(subset)
                        else float("nan")
                    ),
                }
            )
            coefficients, fit = fit_grouped_binomial_clustered(
                subset,
                successes_column="endemism_successes",
                trials_column="endemism_trials",
                predictors=baseline,
                cluster_column=cluster,
                z_threshold=z_threshold,
            )
            if coefficients.empty:
                continue
            coefficients.insert(0, "regime", regime)
            coefficients.insert(0, "minimum_endemism_resolution_fraction", float(threshold))
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
    endemism_input_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/flora_endemism_analysis.yml"), exists=True),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    raw = pd.read_csv(endemism_input_csv)
    coefficients, support = run_resolution_sensitivity(raw, config)
    output_dir.mkdir(parents=True, exist_ok=True)
    coefficients.to_csv(output_dir / "endemism_resolution_sensitivity_coefficients.csv", index=False)
    support.to_csv(output_dir / "endemism_resolution_sensitivity_support.csv", index=False)
    distance = coefficients.loc[
        coefficients["predictor"].eq("log_distance_to_continent_km")
    ].copy() if not coefficients.empty else coefficients.copy()
    distance.to_csv(output_dir / "endemism_resolution_isolation_coefficients.csv", index=False)
    manifest = {
        "contract": "flora_endemism_resolution_sensitivity_v1",
        "thresholds": list(DEFAULT_THRESHOLDS),
        "policy": (
            "Status resolution declines with isolation in the real checkpoint, so the "
            "isolation-endemism association must be stable across stricter resolved-fraction "
            "subsets before it is called robust."
        ),
        "support": support.to_dict("records"),
    }
    (output_dir / "endemism_resolution_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
