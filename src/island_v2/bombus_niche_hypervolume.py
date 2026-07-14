"""Species-level Bombus environmental niche hypervolumes for island compatibility.

This restores the production estimator used by the successful PR #112 real-data run:
winsorized and standardized species environments, ridge-regularized covariance,
Mahalanobis ellipsoidal hypervolumes, empirical tail support, and extrapolation
 diagnostics. The score is environmental compatibility, not realized occurrence or
Bombus absence.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)
ID_COLUMNS = {"island_id", "bombus_species"}


@app.callback()
def main() -> None:
    """Fit and score species-specific Bombus environmental niche envelopes."""


def _environment_columns(
    occurrences: pd.DataFrame,
    island_species: pd.DataFrame,
    requested: list[str] | None = None,
) -> list[str]:
    if requested:
        missing = [
            column
            for column in requested
            if column not in occurrences.columns or column not in island_species.columns
        ]
        if missing:
            raise typer.BadParameter(
                f"environment columns missing from one or both tables: {missing}"
            )
        return requested
    shared = [
        column
        for column in occurrences.columns
        if column in island_species.columns and column not in ID_COLUMNS
    ]
    numeric: list[str] = []
    for column in shared:
        left = pd.to_numeric(occurrences[column], errors="coerce")
        right = pd.to_numeric(island_species[column], errors="coerce")
        if left.notna().any() and right.notna().any():
            numeric.append(column)
    if not numeric:
        raise typer.BadParameter("no shared numeric environmental columns were found")
    return numeric


def _regularized_inverse_covariance(
    matrix: np.ndarray,
    ridge_fraction: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    covariance = np.atleast_2d(np.cov(matrix, rowvar=False)).astype(float)
    dimensions = covariance.shape[0]
    average_variance = float(np.trace(covariance) / max(dimensions, 1))
    ridge = max(average_variance * ridge_fraction, 1e-9)
    regularized = covariance + np.eye(dimensions) * ridge
    return regularized, np.linalg.pinv(regularized), float(np.linalg.cond(regularized))


def _mahalanobis_d2(
    values: np.ndarray,
    center: np.ndarray,
    inverse_covariance: np.ndarray,
) -> np.ndarray:
    delta = values - center
    return np.einsum("ij,jk,ik->i", delta, inverse_covariance, delta)


def _standardize_training_environment(
    training: pd.DataFrame,
    environment_columns: list[str],
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    list[str],
    list[str],
]:
    lower = training.quantile(0.01).to_numpy(dtype=float)
    upper = training.quantile(0.99).to_numpy(dtype=float)
    matrix = training.clip(lower=lower, upper=upper, axis="columns").to_numpy(dtype=float)
    center = np.nanmedian(matrix, axis=0)
    scale = np.nanstd(matrix, axis=0, ddof=1)
    usable = np.isfinite(center) & np.isfinite(scale) & (scale > 0)
    active = [c for c, keep in zip(environment_columns, usable, strict=True) if keep]
    dropped = [c for c, keep in zip(environment_columns, usable, strict=True) if not keep]
    if not active:
        empty = np.empty(0, dtype=float)
        return np.empty((len(training), 0), dtype=float), empty, empty, empty, empty, active, dropped
    active_center = center[usable]
    active_scale = scale[usable]
    standardized = (matrix[:, usable] - active_center) / active_scale
    return standardized, active_center, active_scale, lower[usable], upper[usable], active, dropped


def _unresolved_row(
    target: pd.Series,
    species: str,
    n_records: int,
    n_dimensions: int,
    n_requested_dimensions: int,
    dropped_text: str,
    status: str,
    threshold: float | object = pd.NA,
) -> dict[str, object]:
    return {
        "island_id": target["island_id"],
        "bombus_species": species,
        "environmental_compatibility": pd.NA,
        "inside_hypervolume": pd.NA,
        "inside_ellipsoidal_envelope": pd.NA,
        "mahalanobis_d2": pd.NA,
        "hypervolume_d2_threshold": threshold,
        "ellipsoidal_d2_threshold": threshold,
        "hypervolume_distance_ratio": pd.NA,
        "ellipsoidal_distance_ratio": pd.NA,
        "empirical_training_tail_support": pd.NA,
        "univariate_extrapolation_count": pd.NA,
        "univariate_extrapolation_fraction": pd.NA,
        "environmental_extrapolation": pd.NA,
        "covariance_condition_number": pd.NA,
        "n_occurrence_records": n_records,
        "n_environmental_dimensions": n_dimensions,
        "n_environmental_dimensions_requested": n_requested_dimensions,
        "n_environmental_dimensions_active": n_dimensions,
        "dropped_environmental_dimensions": dropped_text,
        "model_status": status,
        "model_type": "ridge_regularized_mahalanobis_ellipsoidal_hypervolume",
    }


def score_niche_hypervolumes(
    occurrences: pd.DataFrame,
    island_species: pd.DataFrame,
    environment_columns: list[str] | None = None,
    min_occurrences: int = 10,
    hypervolume_quantile: float = 0.95,
    ridge_fraction: float = 0.05,
    **_: object,
) -> pd.DataFrame:
    required_occurrence = {"bombus_species"}
    required_island = {"island_id", "bombus_species"}
    if missing := required_occurrence.difference(occurrences.columns):
        raise typer.BadParameter(f"occurrence table missing columns: {sorted(missing)}")
    if missing := required_island.difference(island_species.columns):
        raise typer.BadParameter(f"island source-pool table missing columns: {sorted(missing)}")
    if min_occurrences < 3:
        raise typer.BadParameter("min_occurrences must be at least 3")
    if not 0.5 < hypervolume_quantile < 1.0:
        raise typer.BadParameter("hypervolume_quantile must be between 0.5 and 1")
    if ridge_fraction <= 0:
        raise typer.BadParameter("ridge_fraction must be positive")

    env = _environment_columns(occurrences, island_species, environment_columns)
    occurrence_work = occurrences.copy()
    island_work = island_species.copy()
    occurrence_work["bombus_species"] = occurrence_work["bombus_species"].fillna("").astype(str).str.strip()
    occurrence_work = occurrence_work.loc[occurrence_work["bombus_species"].ne("")].copy()
    for column in ("island_id", "bombus_species"):
        island_work[column] = island_work[column].fillna("").astype(str).str.strip()
    if island_work[["island_id", "bombus_species"]].eq("").any(axis=None):
        raise typer.BadParameter("island source-pool table contains blank island_id or bombus_species")
    if island_work.duplicated(["island_id", "bombus_species"]).any():
        raise typer.BadParameter("island source-pool table must have one row per island_id x bombus_species")
    for column in env:
        occurrence_work[column] = pd.to_numeric(occurrence_work[column], errors="coerce")
        island_work[column] = pd.to_numeric(island_work[column], errors="coerce")

    rows: list[dict[str, object]] = []
    for species, targets in island_work.groupby("bombus_species", sort=True):
        training = occurrence_work.loc[occurrence_work["bombus_species"].eq(species), env].dropna()
        n_records = int(len(training))
        n_requested = len(env)
        if n_records < min_occurrences:
            rows.extend(
                _unresolved_row(t, species, n_records, n_requested, n_requested, "", "insufficient_occurrences")
                for _, t in targets.iterrows()
            )
            continue

        matrix, center, scale, lower, upper, active_env, dropped_env = _standardize_training_environment(training, env)
        n_dimensions = len(active_env)
        dropped_text = "|".join(dropped_env)
        if n_dimensions == 0:
            rows.extend(
                _unresolved_row(t, species, n_records, 0, n_requested, dropped_text, "insufficient_environmental_variation")
                for _, t in targets.iterrows()
            )
            continue
        if n_records < n_dimensions + 2:
            rows.extend(
                _unresolved_row(t, species, n_records, n_dimensions, n_requested, dropped_text, "insufficient_occurrences")
                for _, t in targets.iterrows()
            )
            continue

        _, inverse_covariance, condition_number = _regularized_inverse_covariance(matrix, ridge_fraction)
        standardized_center = np.zeros(n_dimensions, dtype=float)
        training_d2 = _mahalanobis_d2(matrix, standardized_center, inverse_covariance)
        threshold = max(float(np.quantile(training_d2, hypervolume_quantile)), 1e-9)

        for _, target in targets.iterrows():
            target_values = target[active_env]
            if target_values.isna().any():
                rows.append(
                    _unresolved_row(t=target, species=species, n_records=n_records, n_dimensions=n_dimensions, n_requested_dimensions=n_requested, dropped_text=dropped_text, status="missing_island_environment", threshold=threshold)
                )
                continue
            raw_vector = target_values.to_numpy(dtype=float)
            vector = ((raw_vector - center) / scale).reshape(1, -1)
            d2 = float(_mahalanobis_d2(vector, standardized_center, inverse_covariance)[0])
            ratio = d2 / threshold
            compatibility = math.exp(-math.log(2.0) * ratio)
            empirical_tail = (float(np.count_nonzero(training_d2 >= d2)) + 0.5) / (len(training_d2) + 1.0)
            outside_axes = (raw_vector < lower) | (raw_vector > upper)
            extrapolation_count = int(np.count_nonzero(outside_axes))
            inside = bool(d2 <= threshold)
            rows.append({
                "island_id": target["island_id"],
                "bombus_species": species,
                "environmental_compatibility": max(0.0, min(1.0, compatibility)),
                "inside_hypervolume": inside,
                "inside_ellipsoidal_envelope": inside,
                "mahalanobis_d2": d2,
                "hypervolume_d2_threshold": threshold,
                "ellipsoidal_d2_threshold": threshold,
                "hypervolume_distance_ratio": ratio,
                "ellipsoidal_distance_ratio": ratio,
                "empirical_training_tail_support": empirical_tail,
                "univariate_extrapolation_count": extrapolation_count,
                "univariate_extrapolation_fraction": extrapolation_count / n_dimensions,
                "environmental_extrapolation": bool(extrapolation_count > 0),
                "covariance_condition_number": condition_number,
                "n_occurrence_records": n_records,
                "n_environmental_dimensions": n_dimensions,
                "n_environmental_dimensions_requested": n_requested,
                "n_environmental_dimensions_active": n_dimensions,
                "dropped_environmental_dimensions": dropped_text,
                "model_status": "scored",
                "model_type": "ridge_regularized_mahalanobis_ellipsoidal_hypervolume",
            })
    return pd.DataFrame(rows)


@app.command("run")
def run(
    occurrence_environment_csv: Path = typer.Option(..., exists=True),
    island_source_pool_environment_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    environment_columns: str = typer.Option(""),
    min_occurrences: int = typer.Option(10, min=3),
    hypervolume_quantile: float = typer.Option(0.95, min=0.500001, max=0.999999),
    ridge_fraction: float = typer.Option(0.05, min=1e-9),
    pca_variance_retained: float = typer.Option(0.95, min=0.5, max=1.0),
) -> None:
    del pca_variance_retained  # retained only for temporary CLI compatibility with #120
    occurrences = pd.read_csv(occurrence_environment_csv, dtype=str).fillna("")
    island_species = pd.read_csv(island_source_pool_environment_csv, dtype=str).fillna("")
    requested = [value.strip() for value in environment_columns.split(",") if value.strip()]
    result = score_niche_hypervolumes(
        occurrences,
        island_species,
        environment_columns=requested or None,
        min_occurrences=min_occurrences,
        hypervolume_quantile=hypervolume_quantile,
        ridge_fraction=ridge_fraction,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    canonical = output_dir / "bombus_species_environmental_compatibility.csv"
    result.to_csv(canonical, index=False)
    # Compatibility aliases for the short-lived #120 naming change.
    result.to_csv(output_dir / "bombus_regularized_ellipsoidal_niche_scores.csv", index=False)
    result.to_csv(output_dir / "bombus_niche_hypervolume_scores.csv", index=False)
    status_counts = result["model_status"].value_counts(dropna=False).to_dict()
    scored = result.loc[result["model_status"].eq("scored")]
    summary = {
        "n_rows": int(len(result)),
        "n_islands": int(result["island_id"].nunique()),
        "n_bombus_species": int(result["bombus_species"].nunique()),
        "environment_columns": requested or "inferred_shared_numeric",
        "min_occurrences": min_occurrences,
        "hypervolume_quantile": hypervolume_quantile,
        "ridge_fraction": ridge_fraction,
        "model_status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "n_scored_extrapolated": int(scored["environmental_extrapolation"].fillna(False).astype(bool).sum()),
        "compatibility_calibration": "0.5 at the fitted hypervolume boundary",
        "method": "species-specific winsorized and standardized ridge-regularized Mahalanobis ellipsoidal hypervolume with empirical-tail and extrapolation diagnostics",
        "interpretation": "environmental compatibility only; not realized occurrence probability or Bombus absence",
        "reference_implementation": "PR #112 successful real-data estimator",
    }
    (output_dir / "bombus_niche_hypervolume_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
