"""Diagnose whether M0-M4 conclusions depend on inferred trait evidence layers.

The trait matrix is rectangular but not uniform in evidence quality. This module
separates an outcome into species-direct, genus-inferred, and family-inferred layers,
then quantifies both
layer-specific Bombus paths and primary-tier predictive performance across a full
predeclared grid of direct-evidence coverage thresholds.
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
from island_v2.m0_m4_spatial_cv import (
    assign_balanced_block_folds,
    fit_grouped_binomial,
    heldout_metrics,
    predict_probability,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

FILL_LAYERS = ("species_direct", "genus_inference", "family_inference")
PRIMARY_COVERAGE_THRESHOLDS = (1, 3, 5, 10, 20, 30, 50)
MEDIATOR_COVERAGE_THRESHOLDS = (1, 2, 3, 5, 10)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("full-analysis config must be a mapping")
    return payload


def aggregate_outcome_by_fill_layer(
    species_master: pd.DataFrame,
    trait_name: str = "pollination_functional_guild",
    positive_value: str = "bumblebees",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required = {"island_id", "accepted_species", "trait_name", "filled_value", "fill_tier"}
    missing = required.difference(species_master.columns)
    if missing:
        raise typer.BadParameter(f"species master missing columns: {sorted(missing)}")
    work = species_master.loc[species_master["trait_name"].astype(str).eq(trait_name)].copy()
    if work.empty:
        raise typer.BadParameter(f"species master has no rows for {trait_name}")
    work = work.drop_duplicates(["island_id", "accepted_species", "trait_name"])
    work["is_positive"] = work["filled_value"].astype(str).eq(positive_value)

    rows: list[dict[str, Any]] = []
    island_parts: list[pd.DataFrame] = []
    for fill_layer in FILL_LAYERS:
        layer = work.loc[work["fill_tier"].astype(str).eq(fill_layer)].copy()
        grouped = (
            layer.groupby("island_id", sort=True)
            .agg(
                trials=("accepted_species", "nunique"),
                successes=("is_positive", "sum"),
            )
            .reset_index()
        )
        grouped["fill_layer"] = fill_layer
        grouped["success_share"] = grouped["successes"] / grouped["trials"].replace(0, np.nan)
        island_parts.append(grouped)
        rows.append(
            {
                "fill_layer": fill_layer,
                "n_island_species_rows": int(len(layer)),
                "n_unique_species": int(layer["accepted_species"].nunique()),
                "n_islands": int(grouped["island_id"].nunique()),
                "median_species_per_island": float(grouped["trials"].median()) if len(grouped) else np.nan,
                "q25_species_per_island": float(grouped["trials"].quantile(0.25)) if len(grouped) else np.nan,
                "q75_species_per_island": float(grouped["trials"].quantile(0.75)) if len(grouped) else np.nan,
                "pooled_bumblebee_share": (
                    float(grouped["successes"].sum() / grouped["trials"].sum())
                    if grouped["trials"].sum() > 0
                    else np.nan
                ),
            }
        )
    island_counts = pd.concat(island_parts, ignore_index=True) if island_parts else pd.DataFrame()
    return island_counts, pd.DataFrame(rows)


def fit_fill_layer_bombus_paths(
    island_counts: pd.DataFrame,
    reference_island_data: pd.DataFrame,
    baseline: list[str],
    z_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    required_reference = {"island_id", "spatial_block", "bombus_channel_state", *baseline}
    missing = required_reference.difference(reference_island_data.columns)
    if missing:
        raise typer.BadParameter(f"reference island data missing columns: {sorted(missing)}")

    coefficient_parts: list[pd.DataFrame] = []
    fit_rows: list[dict[str, Any]] = []
    reference = reference_island_data[
        ["island_id", "spatial_block", "bombus_channel_state", *baseline]
    ].drop_duplicates("island_id")
    for fill_layer in FILL_LAYERS:
        counts = island_counts.loc[island_counts["fill_layer"].eq(fill_layer)].copy()
        data = reference.merge(counts, on="island_id", how="inner", validate="one_to_one")
        coefficients, fit = fit_grouped_binomial_clustered(
            data,
            "successes",
            "trials",
            [*baseline, "bombus_channel_state"],
            "spatial_block",
            z_threshold,
        )
        if not coefficients.empty:
            coefficients.insert(0, "fill_layer", fill_layer)
            coefficient_parts.append(coefficients)
        fit_rows.append({"fill_layer": fill_layer, **fit})
    coefficients = (
        pd.concat(coefficient_parts, ignore_index=True) if coefficient_parts else pd.DataFrame()
    )
    return coefficients, pd.DataFrame(fit_rows)


def _m1_spatial_cv_for_threshold(
    data: pd.DataFrame,
    baseline: list[str],
    threshold: int,
    n_folds: int,
    repeats: int,
    seed: int,
) -> list[dict[str, Any]]:
    subset = data.loc[
        pd.to_numeric(data["floral_phenotype_trials"], errors="coerce").ge(threshold)
    ].copy()
    predictors = [*baseline, "bombus_channel_state"]
    columns = [
        "island_id",
        "spatial_block",
        "floral_phenotype_successes",
        "floral_phenotype_trials",
        *predictors,
    ]
    work = subset[columns].copy()
    for column in ["floral_phenotype_successes", "floral_phenotype_trials", *predictors]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work.dropna()
    work = work.loc[work["floral_phenotype_trials"] > 0].copy()
    if len(work) < 50 or work["spatial_block"].nunique() < n_folds:
        return [
            {
                "threshold": threshold,
                "status": "insufficient_support",
                "n_islands": int(len(work)),
                "n_spatial_blocks": int(work["spatial_block"].nunique()),
            }
        ]

    rows: list[dict[str, Any]] = []
    for repeat in range(repeats):
        repeat_seed = seed + repeat * 1009
        mapping = assign_balanced_block_folds(work, n_folds, repeat_seed)
        fold_id = work["spatial_block"].astype(str).map(mapping)
        delta_ll = 0.0
        delta_brier_weighted = 0.0
        total_trials = 0.0
        for fold in range(n_folds):
            test = work.loc[fold_id.eq(fold)].copy()
            train = work.loc[~fold_id.eq(fold)].copy()
            reduced_model = fit_grouped_binomial(
                train,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                baseline,
            )
            full_model = fit_grouped_binomial(
                train,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                predictors,
            )
            reduced = heldout_metrics(
                test,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                predict_probability(reduced_model, test),
            )
            full = heldout_metrics(
                test,
                "floral_phenotype_successes",
                "floral_phenotype_trials",
                predict_probability(full_model, test),
            )
            delta_ll += full["log_likelihood"] - reduced["log_likelihood"]
            delta_brier_weighted += (
                full["weighted_brier"] - reduced["weighted_brier"]
            ) * full["n_trials"]
            total_trials += full["n_trials"]
        rows.append(
            {
                "threshold": threshold,
                "status": "fitted",
                "repeat": repeat + 1,
                "n_islands": int(len(work)),
                "n_spatial_blocks": int(work["spatial_block"].nunique()),
                "median_direct_species_per_island": float(
                    work["floral_phenotype_trials"].median()
                ),
                "delta_log_likelihood_per_1000_trials": float(
                    1000.0 * delta_ll / total_trials
                ),
                "delta_brier_full_minus_reduced": float(
                    delta_brier_weighted / total_trials
                ),
            }
        )
    return rows


def primary_coverage_threshold_cv(
    primary_island_data: pd.DataFrame,
    baseline: list[str],
    n_folds: int = 10,
    repeats: int = 5,
    seed: int = 20260713,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for threshold in PRIMARY_COVERAGE_THRESHOLDS:
        rows.extend(
            _m1_spatial_cv_for_threshold(
                primary_island_data,
                baseline,
                threshold,
                n_folds,
                repeats,
                seed,
            )
        )
    repeat_table = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    for threshold, group in repeat_table.groupby("threshold", sort=True):
        fitted = group.loc[group["status"].eq("fitted")]
        if fitted.empty:
            first = group.iloc[0]
            summary_rows.append(
                {
                    "threshold": int(threshold),
                    "status": "insufficient_support",
                    "n_islands": int(first["n_islands"]),
                    "n_spatial_blocks": int(first["n_spatial_blocks"]),
                    "n_repeats": 0,
                    "median_delta_log_likelihood_per_1000_trials": np.nan,
                    "proportion_repeats_positive_loglik_gain": np.nan,
                    "median_delta_brier_full_minus_reduced": np.nan,
                    "proportion_repeats_brier_improved": np.nan,
                }
            )
            continue
        summary_rows.append(
            {
                "threshold": int(threshold),
                "status": "fitted",
                "n_islands": int(fitted["n_islands"].iloc[0]),
                "n_spatial_blocks": int(fitted["n_spatial_blocks"].iloc[0]),
                "n_repeats": int(len(fitted)),
                "median_delta_log_likelihood_per_1000_trials": float(
                    fitted["delta_log_likelihood_per_1000_trials"].median()
                ),
                "proportion_repeats_positive_loglik_gain": float(
                    (fitted["delta_log_likelihood_per_1000_trials"] > 0).mean()
                ),
                "median_delta_brier_full_minus_reduced": float(
                    fitted["delta_brier_full_minus_reduced"].median()
                ),
                "proportion_repeats_brier_improved": float(
                    (fitted["delta_brier_full_minus_reduced"] < 0).mean()
                ),
            }
        )
    return repeat_table, pd.DataFrame(summary_rows)


def mediator_direct_evidence_coverage(primary_island_data: pd.DataFrame) -> pd.DataFrame:
    required = {
        "floral_phenotype_trials",
        "mixed_mating_trials",
        "self_compatibility_trials",
        "spatial_block",
    }
    missing = required.difference(primary_island_data.columns)
    if missing:
        raise typer.BadParameter(f"primary island data missing columns: {sorted(missing)}")
    rows: list[dict[str, Any]] = []
    for threshold in MEDIATOR_COVERAGE_THRESHOLDS:
        mask = (
            pd.to_numeric(primary_island_data["floral_phenotype_trials"], errors="coerce").ge(threshold)
            & pd.to_numeric(primary_island_data["mixed_mating_trials"], errors="coerce").ge(threshold)
            & pd.to_numeric(primary_island_data["self_compatibility_trials"], errors="coerce").ge(threshold)
        )
        subset = primary_island_data.loc[mask]
        rows.append(
            {
                "minimum_direct_species_per_required_trait": threshold,
                "n_islands": int(len(subset)),
                "n_spatial_blocks": int(subset["spatial_block"].nunique()),
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    species_master_csv: Path = typer.Option(..., exists=True),
    primary_island_data_csv: Path = typer.Option(..., exists=True),
    reference_island_data_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_full_analysis.yml"), exists=True),
) -> None:
    config = _load_config(config_path)
    master = pd.read_csv(
        species_master_csv,
        usecols=["island_id", "accepted_species", "trait_name", "filled_value", "fill_tier"],
    )
    primary = pd.read_csv(primary_island_data_csv)
    reference = pd.read_csv(reference_island_data_csv)
    baseline = [str(value) for value in config["baseline_covariates"]]
    z_threshold = float(config["reporting"].get("nominal_z_threshold", 1.96))

    island_counts, layer_summary = aggregate_outcome_by_fill_layer(master)
    layer_coefficients, layer_fit = fit_fill_layer_bombus_paths(
        island_counts, reference, baseline, z_threshold
    )
    threshold_repeats, threshold_summary = primary_coverage_threshold_cv(
        primary, baseline
    )
    mediator_coverage = mediator_direct_evidence_coverage(primary)

    output_dir.mkdir(parents=True, exist_ok=True)
    island_counts.to_csv(output_dir / "fill_layer_island_counts.csv.gz", index=False, compression="gzip")
    layer_summary.to_csv(output_dir / "fill_layer_summary.csv", index=False)
    layer_coefficients.to_csv(output_dir / "fill_layer_bombus_path_coefficients.csv", index=False)
    layer_fit.to_csv(output_dir / "fill_layer_bombus_path_fit.csv", index=False)
    threshold_repeats.to_csv(output_dir / "primary_coverage_threshold_cv_repeats.csv", index=False)
    threshold_summary.to_csv(output_dir / "primary_coverage_threshold_cv_summary.csv", index=False)
    mediator_coverage.to_csv(output_dir / "primary_mediator_direct_evidence_coverage.csv", index=False)

    bombus_rows = layer_coefficients.loc[
        layer_coefficients["predictor"].eq("bombus_channel_state")
    ].copy()
    manifest = {
        "contract": "m0_m4_trait_evidence_layer_diagnostics_v1",
        "outcome_trait": "pollination_functional_guild",
        "positive_outcome": "bumblebees",
        "fill_layers": layer_summary.to_dict("records"),
        "layer_specific_bombus_paths": bombus_rows.to_dict("records"),
        "primary_coverage_thresholds": threshold_summary.to_dict("records"),
        "primary_mediator_direct_evidence_coverage": mediator_coverage.to_dict("records"),
        "interpretation": (
            "This diagnostic separates evidence quantity from evidence tier. Threshold results are reported as a full grid and must not be cherry-picked."
        ),
    }
    (output_dir / "evidence_layer_diagnostics_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
