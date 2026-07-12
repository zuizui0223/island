"""Species-level Bombus environmental niche hypervolumes for island compatibility.

The estimator is deliberately dependency-light and auditable. For each Bombus species,
it fits a regularized multivariate ellipsoidal niche from occurrence-site environmental
variables, then scores island environments against that species-specific hypervolume.
Species scores remain explicit and can later be aggregated within a predeclared source pool.
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
) -> tuple[np.ndarray, np.ndarray]:
    covariance = np.cov(matrix, rowvar=False)
    covariance = np.atleast_2d(covariance).astype(float)
    dimensions = covariance.shape[0]
    average_variance = float(np.trace(covariance) / max(dimensions, 1))
    ridge = max(average_variance * ridge_fraction, 1e-9)
    regularized = covariance + np.eye(dimensions) * ridge
    return regularized, np.linalg.pinv(regularized)


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[str]]:
    """Winsorize and standardize a species' niche variables before ridge fitting.

    Ridge regularization on raw variables is unit-dependent: a precipitation
    column measured in millimetres can dominate a temperature column measured in
    degrees. Species-specific standardization makes the fitted distance invariant
    to harmless unit rescaling. Constant variables carry no fitted niche
    information and are reported as dropped rather than regularized arbitrarily.
    """
    lower = training.quantile(0.01)
    upper = training.quantile(0.99)
    clipped = training.clip(lower=lower, upper=upper, axis="columns")
    matrix = clipped.to_numpy(dtype=float)
    center = np.nanmedian(matrix, axis=0)
    scale = np.nanstd(matrix, axis=0, ddof=1)
    usable = np.isfinite(center) & np.isfinite(scale) & (scale > 0)
    active_columns = [
        column for column, keep in zip(environment_columns, usable, strict=True) if keep
    ]
    dropped_columns = [
        column for column, keep in zip(environment_columns, usable, strict=True) if not keep
    ]
    if not active_columns:
        return (
            np.empty((len(training), 0), dtype=float),
            np.empty(0, dtype=float),
            np.empty(0, dtype=float),
            active_columns,
            dropped_columns,
        )
    active_center = center[usable]
    active_scale = scale[usable]
    standardized = (matrix[:, usable] - active_center) / active_scale
    return standardized, active_center, active_scale, active_columns, dropped_columns


def score_niche_hypervolumes(
    occurrences: pd.DataFrame,
    island_species: pd.DataFrame,
    environment_columns: list[str] | None = None,
    min_occurrences: int = 10,
    hypervolume_quantile: float = 0.95,
    ridge_fraction: float = 0.05,
) -> pd.DataFrame:
    """Score island x source-pool Bombus rows against species-specific niche volumes."""
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
    occurrence_work["bombus_species"] = (
        occurrence_work["bombus_species"].fillna("").astype(str).str.strip()
    )
    occurrence_work = occurrence_work.loc[occurrence_work["bombus_species"].ne("")].copy()
    for column in ("island_id", "bombus_species"):
        island_work[column] = island_work[column].fillna("").astype(str).str.strip()
    if island_work[["island_id", "bombus_species"]].eq("").any(axis=None):
        raise typer.BadParameter(
            "island source-pool table contains blank island_id or bombus_species"
        )
    if island_work.duplicated(["island_id", "bombus_species"]).any():
        raise typer.BadParameter(
            "island source-pool table must have one row per island_id x bombus_species"
        )
    for column in env:
        occurrence_work[column] = pd.to_numeric(occurrence_work[column], errors="coerce")
        island_work[column] = pd.to_numeric(island_work[column], errors="coerce")

    rows: list[dict[str, object]] = []
    for species, targets in island_work.groupby("bombus_species", sort=True):
        training = occurrence_work.loc[
            occurrence_work["bombus_species"].eq(species), env
        ].dropna()
        n_records = int(len(training))
        n_requested_dimensions = len(env)

        if n_records < min_occurrences:
            for _, target in targets.iterrows():
                rows.append(
                    {
                        "island_id": target["island_id"],
                        "bombus_species": species,
                        "environmental_compatibility": pd.NA,
                        "inside_hypervolume": pd.NA,
                        "mahalanobis_d2": pd.NA,
                        "hypervolume_d2_threshold": pd.NA,
                        "n_occurrence_records": n_records,
                        "n_environmental_dimensions": n_requested_dimensions,
                        "n_environmental_dimensions_requested": n_requested_dimensions,
                        "dropped_environmental_dimensions": "",
                        "model_status": "insufficient_occurrences",
                    }
                )
            continue

        matrix, center, scale, active_env, dropped_env = _standardize_training_environment(
            training,
            env,
        )
        n_dimensions = len(active_env)
        dropped_text = "|".join(dropped_env)
        if n_dimensions == 0:
            for _, target in targets.iterrows():
                rows.append(
                    {
                        "island_id": target["island_id"],
                        "bombus_species": species,
                        "environmental_compatibility": pd.NA,
                        "inside_hypervolume": pd.NA,
                        "mahalanobis_d2": pd.NA,
                        "hypervolume_d2_threshold": pd.NA,
                        "n_occurrence_records": n_records,
                        "n_environmental_dimensions": 0,
                        "n_environmental_dimensions_requested": n_requested_dimensions,
                        "dropped_environmental_dimensions": dropped_text,
                        "model_status": "insufficient_environmental_variation",
                    }
                )
            continue
        if n_records < n_dimensions + 2:
            for _, target in targets.iterrows():
                rows.append(
                    {
                        "island_id": target["island_id"],
                        "bombus_species": species,
                        "environmental_compatibility": pd.NA,
                        "inside_hypervolume": pd.NA,
                        "mahalanobis_d2": pd.NA,
                        "hypervolume_d2_threshold": pd.NA,
                        "n_occurrence_records": n_records,
                        "n_environmental_dimensions": n_dimensions,
                        "n_environmental_dimensions_requested": n_requested_dimensions,
                        "dropped_environmental_dimensions": dropped_text,
                        "model_status": "insufficient_occurrences",
                    }
                )
            continue

        _, inverse_covariance = _regularized_inverse_covariance(
            matrix,
            ridge_fraction=ridge_fraction,
        )
        standardized_center = np.zeros(n_dimensions, dtype=float)
        training_d2 = _mahalanobis_d2(
            matrix,
            standardized_center,
            inverse_covariance,
        )
        threshold = float(np.quantile(training_d2, hypervolume_quantile))
        threshold = max(threshold, 1e-9)

        for _, target in targets.iterrows():
            target_values = target[active_env]
            if target_values.isna().any():
                rows.append(
                    {
                        "island_id": target["island_id"],
                        "bombus_species": species,
                        "environmental_compatibility": pd.NA,
                        "inside_hypervolume": pd.NA,
                        "mahalanobis_d2": pd.NA,
                        "hypervolume_d2_threshold": threshold,
                        "n_occurrence_records": n_records,
                        "n_environmental_dimensions": n_dimensions,
                        "n_environmental_dimensions_requested": n_requested_dimensions,
                        "dropped_environmental_dimensions": dropped_text,
                        "model_status": "missing_island_environment",
                    }
                )
                continue

            vector = (
                (target_values.to_numpy(dtype=float) - center) / scale
            ).reshape(1, -1)
            d2 = float(
                _mahalanobis_d2(
                    vector,
                    standardized_center,
                    inverse_covariance,
                )[0]
            )
            scaled_distance = d2 / threshold
            compatibility = math.exp(-0.5 * scaled_distance)
            rows.append(
                {
                    "island_id": target["island_id"],
                    "bombus_species": species,
                    "environmental_compatibility": max(
                        0.0,
                        min(1.0, compatibility),
                    ),
                    "inside_hypervolume": bool(d2 <= threshold),
                    "mahalanobis_d2": d2,
                    "hypervolume_d2_threshold": threshold,
                    "n_occurrence_records": n_records,
                    "n_environmental_dimensions": n_dimensions,
                    "n_environmental_dimensions_requested": n_requested_dimensions,
                    "dropped_environmental_dimensions": dropped_text,
                    "model_status": "scored",
                }
            )

    return pd.DataFrame(rows)


@app.command("run")
def run(
    occurrence_environment_csv: Path = typer.Option(..., exists=True),
    island_source_pool_environment_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    environment_columns: str = typer.Option(
        "",
        help=(
            "Comma-separated environmental columns. Empty infers shared numeric columns "
            "for local exploration only; production workflows require an explicit list."
        ),
    ),
    min_occurrences: int = typer.Option(10, min=3),
    hypervolume_quantile: float = typer.Option(0.95, min=0.500001, max=0.999999),
    ridge_fraction: float = typer.Option(0.05, min=1e-9),
) -> None:
    occurrences = pd.read_csv(occurrence_environment_csv, dtype=str).fillna("")
    island_species = pd.read_csv(
        island_source_pool_environment_csv,
        dtype=str,
    ).fillna("")
    requested = [
        value.strip()
        for value in environment_columns.split(",")
        if value.strip()
    ]
    result = score_niche_hypervolumes(
        occurrences,
        island_species,
        environment_columns=requested or None,
        min_occurrences=min_occurrences,
        hypervolume_quantile=hypervolume_quantile,
        ridge_fraction=ridge_fraction,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "bombus_species_environmental_compatibility.csv"
    result.to_csv(output_path, index=False)
    status_counts = result["model_status"].value_counts(dropna=False).to_dict()
    summary = {
        "n_rows": int(len(result)),
        "n_islands": int(result["island_id"].nunique()),
        "n_bombus_species": int(result["bombus_species"].nunique()),
        "environment_columns": requested or "inferred_shared_numeric",
        "min_occurrences": min_occurrences,
        "hypervolume_quantile": hypervolume_quantile,
        "ridge_fraction": ridge_fraction,
        "model_status_counts": {str(key): int(value) for key, value in status_counts.items()},
        "method": (
            "species-specific winsorized and standardized regularized Mahalanobis "
            "ellipsoidal hypervolume"
        ),
    }
    (output_dir / "bombus_niche_hypervolume_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
