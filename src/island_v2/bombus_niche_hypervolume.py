"""Species-level regularized ellipsoidal environmental-niche compatibility.

The historical module name is retained for CLI compatibility. The fitted object is NOT
claimed to be a generic ecological hypervolume, realized occurrence probability,
source-pool membership, or Bombus absence. Correlated predictors are compressed to a
species-specific PCA basis before the regularized Mahalanobis envelope is fitted.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)
ID_COLUMNS = {"island_id", "bombus_species", "projection_scope", "environment_point_estimate"}


def _environment_columns(occurrences: pd.DataFrame, island_species: pd.DataFrame, requested: list[str] | None = None) -> list[str]:
    if requested:
        missing = [c for c in requested if c not in occurrences.columns or c not in island_species.columns]
        if missing:
            raise typer.BadParameter(f"environment columns missing from one or both tables: {missing}")
        return requested
    shared = [c for c in occurrences.columns if c in island_species.columns and c not in ID_COLUMNS]
    numeric = []
    for c in shared:
        if pd.to_numeric(occurrences[c], errors="coerce").notna().any() and pd.to_numeric(island_species[c], errors="coerce").notna().any():
            numeric.append(c)
    if not numeric:
        raise typer.BadParameter("no shared numeric environmental columns were found")
    return numeric


def _standardize_training_environment(training: pd.DataFrame, environment_columns: list[str]):
    lower = training.quantile(0.01).to_numpy(dtype=float)
    upper = training.quantile(0.99).to_numpy(dtype=float)
    matrix = training.clip(lower=lower, upper=upper, axis="columns").to_numpy(dtype=float)
    center = np.nanmedian(matrix, axis=0)
    scale = np.nanstd(matrix, axis=0, ddof=1)
    usable = np.isfinite(center) & np.isfinite(scale) & (scale > 0)
    active = [c for c, keep in zip(environment_columns, usable, strict=True) if keep]
    dropped = [c for c, keep in zip(environment_columns, usable, strict=True) if not keep]
    if not active:
        return np.empty((len(training), 0)), np.empty(0), np.empty(0), np.empty(0), np.empty(0), active, dropped
    return (matrix[:, usable] - center[usable]) / scale[usable], center[usable], scale[usable], lower[usable], upper[usable], active, dropped


def _pca_basis(matrix: np.ndarray, variance_retained: float = 0.95) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if matrix.ndim != 2 or matrix.shape[1] == 0:
        return matrix, np.empty((0, 0)), np.empty(0)
    _, singular, vt = np.linalg.svd(matrix, full_matrices=False)
    eigen = (singular ** 2) / max(matrix.shape[0] - 1, 1)
    fractions = eigen / eigen.sum() if eigen.sum() > 0 else np.ones_like(eigen) / len(eigen)
    k = int(np.searchsorted(np.cumsum(fractions), variance_retained) + 1)
    k = max(1, min(k, matrix.shape[1], matrix.shape[0] - 2))
    loadings = vt[:k].T
    return matrix @ loadings, loadings, fractions[:k]


def _regularized_inverse_covariance(matrix: np.ndarray, ridge_fraction: float):
    covariance = np.atleast_2d(np.cov(matrix, rowvar=False)).astype(float)
    dimensions = covariance.shape[0]
    average_variance = float(np.trace(covariance) / max(dimensions, 1))
    ridge = max(average_variance * ridge_fraction, 1e-9)
    regularized = covariance + np.eye(dimensions) * ridge
    return np.linalg.pinv(regularized), float(np.linalg.cond(regularized))


def _mahalanobis_d2(values: np.ndarray, inverse_covariance: np.ndarray) -> np.ndarray:
    return np.einsum("ij,jk,ik->i", values, inverse_covariance, values)


def _unresolved_row(target: pd.Series, species: str, n_records: int, status: str, n_requested: int, n_active: int = 0, n_pca: int = 0, threshold: float | object = pd.NA):
    return {
        "island_id": target["island_id"], "bombus_species": species,
        "environmental_compatibility": pd.NA, "inside_ellipsoidal_envelope": pd.NA,
        "mahalanobis_d2": pd.NA, "ellipsoidal_d2_threshold": threshold,
        "ellipsoidal_distance_ratio": pd.NA, "empirical_training_tail_support": pd.NA,
        "univariate_extrapolation_count": pd.NA, "univariate_extrapolation_fraction": pd.NA,
        "environmental_extrapolation": pd.NA, "covariance_condition_number": pd.NA,
        "n_occurrence_records": n_records, "n_environmental_dimensions_requested": n_requested,
        "n_environmental_dimensions_active": n_active, "n_pca_dimensions": n_pca,
        "model_status": status, "model_type": "regularized_ellipsoidal_environmental_niche",
    }


def score_niche_hypervolumes(occurrences: pd.DataFrame, island_species: pd.DataFrame, environment_columns: list[str] | None = None, min_occurrences: int = 50, envelope_quantile: float = 0.95, ridge_fraction: float = 0.05, pca_variance_retained: float = 0.95) -> pd.DataFrame:
    if missing := {"bombus_species"}.difference(occurrences.columns):
        raise typer.BadParameter(f"occurrence table missing columns: {sorted(missing)}")
    if missing := {"island_id", "bombus_species"}.difference(island_species.columns):
        raise typer.BadParameter(f"island target table missing columns: {sorted(missing)}")
    if min_occurrences < 20:
        raise typer.BadParameter("min_occurrences must be at least 20")
    if not 0.5 < envelope_quantile < 1.0:
        raise typer.BadParameter("envelope_quantile must be between 0.5 and 1")
    if not 0.5 <= pca_variance_retained <= 1.0:
        raise typer.BadParameter("pca_variance_retained must be between 0.5 and 1")
    if ridge_fraction <= 0:
        raise typer.BadParameter("ridge_fraction must be positive")

    env = _environment_columns(occurrences, island_species, environment_columns)
    occurrence_work = occurrences.copy()
    island_work = island_species.copy()
    occurrence_work["bombus_species"] = occurrence_work["bombus_species"].fillna("").astype(str).str.strip()
    occurrence_work = occurrence_work[occurrence_work["bombus_species"].ne("")]
    for c in ("island_id", "bombus_species"):
        island_work[c] = island_work[c].fillna("").astype(str).str.strip()
    if island_work[["island_id", "bombus_species"]].eq("").any(axis=None):
        raise typer.BadParameter("island target table contains blank island_id or bombus_species")
    if island_work.duplicated(["island_id", "bombus_species"]).any():
        raise typer.BadParameter("island target table must have one row per island_id x bombus_species")
    for c in env:
        occurrence_work[c] = pd.to_numeric(occurrence_work[c], errors="coerce")
        island_work[c] = pd.to_numeric(island_work[c], errors="coerce")

    rows = []
    for species, targets in island_work.groupby("bombus_species", sort=True):
        training = occurrence_work.loc[occurrence_work["bombus_species"].eq(species), env].dropna()
        n_records = len(training)
        if n_records < min_occurrences:
            rows.extend(_unresolved_row(t, species, n_records, "insufficient_occurrences_after_quality_filtering", len(env)) for _, t in targets.iterrows())
            continue
        standardized, center, scale, lower, upper, active_env, dropped_env = _standardize_training_environment(training, env)
        if not active_env:
            rows.extend(_unresolved_row(t, species, n_records, "insufficient_environmental_variation", len(env)) for _, t in targets.iterrows())
            continue
        pca_training, loadings, explained = _pca_basis(standardized, pca_variance_retained)
        n_pca = pca_training.shape[1]
        if n_records < max(min_occurrences, 5 * n_pca):
            rows.extend(_unresolved_row(t, species, n_records, "insufficient_records_for_effective_dimension", len(env), len(active_env), n_pca) for _, t in targets.iterrows())
            continue
        inverse_covariance, condition_number = _regularized_inverse_covariance(pca_training, ridge_fraction)
        training_d2 = _mahalanobis_d2(pca_training, inverse_covariance)
        threshold = max(float(np.quantile(training_d2, envelope_quantile)), 1e-9)
        for _, target in targets.iterrows():
            values = target[active_env]
            if values.isna().any():
                rows.append(_unresolved_row(target, species, n_records, "missing_island_environment", len(env), len(active_env), n_pca, threshold))
                continue
            raw = values.to_numpy(dtype=float)
            pca_target = (((raw - center) / scale).reshape(1, -1)) @ loadings
            d2 = float(_mahalanobis_d2(pca_target, inverse_covariance)[0])
            ratio = d2 / threshold
            compatibility = math.exp(-math.log(2.0) * ratio)
            empirical_tail = (float(np.count_nonzero(training_d2 >= d2)) + 0.5) / (len(training_d2) + 1.0)
            outside = (raw < lower) | (raw > upper)
            rows.append({
                "island_id": target["island_id"], "bombus_species": species,
                "environmental_compatibility": max(0.0, min(1.0, compatibility)),
                "inside_ellipsoidal_envelope": bool(d2 <= threshold),
                "mahalanobis_d2": d2, "ellipsoidal_d2_threshold": threshold,
                "ellipsoidal_distance_ratio": ratio, "empirical_training_tail_support": empirical_tail,
                "univariate_extrapolation_count": int(np.count_nonzero(outside)),
                "univariate_extrapolation_fraction": float(np.count_nonzero(outside) / len(active_env)),
                "environmental_extrapolation": bool(np.count_nonzero(outside) > 0),
                "covariance_condition_number": condition_number,
                "n_occurrence_records": n_records, "n_environmental_dimensions_requested": len(env),
                "n_environmental_dimensions_active": len(active_env), "n_pca_dimensions": n_pca,
                "pca_variance_retained_actual": float(explained.sum()),
                "dropped_environmental_dimensions": "|".join(dropped_env),
                "model_status": "scored", "model_type": "regularized_ellipsoidal_environmental_niche",
            })
    return pd.DataFrame(rows)


@app.command("run")
def run(
    occurrence_environment_csv: Path = typer.Option(..., exists=True),
    island_source_pool_environment_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    environment_columns: str = typer.Option(""),
    min_occurrences: int = typer.Option(50, min=20),
    hypervolume_quantile: float = typer.Option(0.95, min=0.500001, max=0.999999),
    ridge_fraction: float = typer.Option(0.05, min=1e-9),
    pca_variance_retained: float = typer.Option(0.95, min=0.5, max=1.0),
) -> None:
    occurrences = pd.read_csv(occurrence_environment_csv, dtype=str).fillna("")
    targets = pd.read_csv(island_source_pool_environment_csv, dtype=str).fillna("")
    requested = [x.strip() for x in environment_columns.split(",") if x.strip()]
    result = score_niche_hypervolumes(occurrences, targets, requested or None, min_occurrences, hypervolume_quantile, ridge_fraction, pca_variance_retained)
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "bombus_regularized_ellipsoidal_niche_scores.csv", index=False)
    result.to_csv(output_dir / "bombus_niche_hypervolume_scores.csv", index=False)
    manifest = {
        "model_type": "regularized_ellipsoidal_environmental_niche",
        "interpretation": "relative climatic-environmental compatibility only; not occurrence probability, source-pool membership, absence, abundance, or pollination service",
        "n_rows": int(len(result)), "n_scored": int(result["model_status"].eq("scored").sum()),
        "min_occurrences_after_quality_filtering": min_occurrences,
        "envelope_quantile": hypervolume_quantile, "ridge_fraction": ridge_fraction,
        "pca_variance_retained": pca_variance_retained, "environment_columns": requested or _environment_columns(occurrences, targets),
    }
    (output_dir / "bombus_ellipsoidal_niche_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest))


if __name__ == "__main__":
    app()
