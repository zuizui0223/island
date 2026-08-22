"""Isolation-endemism analysis on source-backed native floristic status.

Endemism is modelled among *endemism-resolved native species only*:

    endemic / (endemic + confirmed native non-endemic)

Native species whose endemic status is unresolved are never counted as
non-endemics. Their fraction is reported as a measurement diagnostic. The model
is an island-level floristic association, not an estimate of speciation rate.
"""

from __future__ import annotations

import json
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
    required = {"regime_column", "regimes", "baseline_covariates", "spatial_cluster"}
    if not isinstance(payload, dict) or not required.issubset(payload):
        raise typer.BadParameter(f"config must contain {sorted(required)}")
    return payload


def prepare_endemism_response(counts: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    required = {
        "island_id",
        "n_native_species",
        "n_native_nonendemic_species",
        "n_endemic_species",
        str(config["regime_column"]),
        str(config["spatial_cluster"]),
    }
    required |= set(str(value) for value in config["baseline_covariates"])
    missing = required - set(counts.columns)
    if missing:
        raise typer.BadParameter(f"endemism input missing columns: {sorted(missing)}")

    frame = counts.copy()
    for column in (
        "n_native_species",
        "n_native_nonendemic_species",
        "n_endemic_species",
    ):
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["endemism_successes"] = frame["n_endemic_species"]
    frame["endemism_trials"] = (
        frame["n_endemic_species"] + frame["n_native_nonendemic_species"]
    )
    frame["endemism_resolution_fraction"] = (
        frame["endemism_trials"] / frame["n_native_species"].replace(0, np.nan)
    )
    frame["endemic_share_resolved_native"] = (
        frame["endemism_successes"] / frame["endemism_trials"].replace(0, np.nan)
    )
    min_resolved = int(config.get("min_resolved_native_species", 1))
    frame["endemism_support_eligible"] = frame["endemism_trials"].ge(min_resolved)
    regimes = {str(value) for value in config["regimes"]}
    frame = frame.loc[frame[str(config["regime_column"])].astype(str).isin(regimes)].copy()
    return frame


def _fit_one(
    data: pd.DataFrame,
    predictors: list[str],
    cluster_column: str,
    z_threshold: float,
    label: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    coef, fit = fit_grouped_binomial_clustered(
        data,
        successes_column="endemism_successes",
        trials_column="endemism_trials",
        predictors=predictors,
        cluster_column=cluster_column,
        z_threshold=z_threshold,
    )
    if not coef.empty:
        coef.insert(0, "model", label)
    return coef, {"model": label, **fit}


def fit_endemism_models(
    prepared: pd.DataFrame, config: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    regime_column = str(config["regime_column"])
    cluster = str(config["spatial_cluster"])
    baseline = [str(value) for value in config["baseline_covariates"]]
    z_threshold = float(config.get("nominal_z_threshold", 1.96))
    eligible = prepared.loc[prepared["endemism_support_eligible"]].copy()

    coefficient_tables: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    support_rows: list[dict[str, Any]] = []
    for regime in [str(value) for value in config["regimes"]]:
        subset = eligible.loc[eligible[regime_column].astype(str).eq(regime)].copy()
        support_rows.append(
            {
                "regime": regime,
                "n_islands": int(len(subset)),
                "n_spatial_blocks": int(subset[cluster].dropna().nunique()),
                "total_resolved_native_species": float(subset["endemism_trials"].sum()),
                "median_endemism_resolution_fraction": float(
                    subset["endemism_resolution_fraction"].median()
                )
                if len(subset)
                else np.nan,
            }
        )
        coef, fit = _fit_one(
            subset,
            baseline,
            cluster,
            z_threshold,
            f"endemism_{regime}",
        )
        if not coef.empty:
            coef.insert(1, "regime", regime)
            coefficient_tables.append(coef)
        fit_rows.append({"regime": regime, **fit})

    coefficients = (
        pd.concat(coefficient_tables, ignore_index=True)
        if coefficient_tables
        else pd.DataFrame()
    )
    return coefficients, pd.DataFrame(fit_rows), pd.DataFrame(support_rows)


@app.command("run")
def run(
    endemism_input_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/flora_endemism_analysis.yml"), exists=True),
) -> None:
    config = _load_config(config_path)
    raw = pd.read_csv(endemism_input_csv)
    prepared = prepare_endemism_response(raw, config)
    coefficients, fits, support = fit_endemism_models(prepared, config)

    output_dir.mkdir(parents=True, exist_ok=True)
    prepared.to_csv(output_dir / "endemism_analysis_input.csv", index=False)
    coefficients.to_csv(output_dir / "endemism_coefficients.csv", index=False)
    fits.to_csv(output_dir / "endemism_model_fits.csv", index=False)
    support.to_csv(output_dir / "endemism_support.csv", index=False)
    manifest = {
        "contract": "flora_endemism_analysis_v1",
        "response": "endemic / (endemic + confirmed native non-endemic)",
        "unresolved_endemism_policy": "excluded from response denominator; resolution fraction reported",
        "regimes": [str(value) for value in config["regimes"]],
        "baseline_covariates": [str(value) for value in config["baseline_covariates"]],
        "interpretation": (
            "Island-level floristic association only; Isolation-endemism association does not by "
            "itself distinguish differential colonization, extinction, or in-situ speciation."
        ),
        "support": support.to_dict("records"),
    }
    (output_dir / "endemism_analysis_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
