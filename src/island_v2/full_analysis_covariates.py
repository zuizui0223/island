"""Build deterministic baseline covariates for the locked M0-M4 full analysis.

No new biological occurrence data are acquired here. Island area comes from the
frozen GSHHG island manifest, geography and climate come from the already locked
PR #112 island-environment products, and climate PCs are fitted once over the
frozen global island universe.
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

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _load_config(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise typer.BadParameter("full-analysis config must be a mapping")
    return payload


def _circular_mean_longitude(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if not len(numeric):
        return float("nan")
    radians = np.deg2rad(numeric)
    angle = math.atan2(float(np.sin(radians).mean()), float(np.cos(radians).mean()))
    return float(np.rad2deg(angle))


def _island_geography(samples: pd.DataFrame, block_degrees: float) -> pd.DataFrame:
    required = {"island_id", "decimal_longitude", "decimal_latitude"}
    missing = required.difference(samples.columns)
    if missing:
        raise typer.BadParameter(f"environment samples missing columns: {sorted(missing)}")
    if block_degrees <= 0 or block_degrees > 90:
        raise typer.BadParameter("spatial block_degrees must be in (0, 90]")

    rows: list[dict[str, Any]] = []
    for island_id, group in samples.groupby("island_id", sort=True):
        latitude = pd.to_numeric(group["decimal_latitude"], errors="coerce").median()
        longitude = _circular_mean_longitude(group["decimal_longitude"])
        if pd.isna(latitude) or not math.isfinite(longitude):
            lat_block = pd.NA
            lon_block = pd.NA
            spatial_block = "missing"
        else:
            lat_block = int(math.floor((float(latitude) + 90.0) / block_degrees))
            lon_block = int(math.floor((float(longitude) + 180.0) / block_degrees))
            spatial_block = f"lat{lat_block:02d}_lon{lon_block:02d}"
        rows.append(
            {
                "island_id": str(island_id),
                "island_latitude": float(latitude) if pd.notna(latitude) else np.nan,
                "island_longitude": longitude,
                "abs_latitude": abs(float(latitude)) if pd.notna(latitude) else np.nan,
                "spatial_block": spatial_block,
            }
        )
    return pd.DataFrame(rows)


def _fit_climate_pca(
    environment: pd.DataFrame,
    climate_columns: list[str],
    n_pcs: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    missing = {"island_id", *climate_columns}.difference(environment.columns)
    if missing:
        raise typer.BadParameter(f"island environment missing columns: {sorted(missing)}")
    if n_pcs < 1 or n_pcs > len(climate_columns):
        raise typer.BadParameter("n_pcs must be between 1 and the number of climate columns")

    work = environment[["island_id", *climate_columns]].copy()
    values = work[climate_columns].apply(pd.to_numeric, errors="coerce")
    complete = values.notna().all(axis=1)
    if int(complete.sum()) <= n_pcs:
        raise typer.BadParameter("too few complete islands to fit climate PCA")

    fitted = values.loc[complete].copy()
    means = fitted.mean(axis=0)
    sds = fitted.std(axis=0, ddof=0)
    if (sds <= 0).any() or (~np.isfinite(sds)).any():
        bad = sds.index[(sds <= 0) | (~np.isfinite(sds))].tolist()
        raise typer.BadParameter(f"constant or invalid climate columns: {bad}")

    z = (fitted - means) / sds
    _, singular_values, vt = np.linalg.svd(z.to_numpy(dtype=float), full_matrices=False)
    loadings = vt[:n_pcs].T.copy()
    scores = z.to_numpy(dtype=float) @ loadings

    # SVD signs are arbitrary. Orient each PC deterministically so the loading with
    # the largest absolute magnitude is positive.
    for pc_index in range(n_pcs):
        anchor = int(np.argmax(np.abs(loadings[:, pc_index])))
        if loadings[anchor, pc_index] < 0:
            loadings[:, pc_index] *= -1.0
            scores[:, pc_index] *= -1.0

    eigenvalues = singular_values**2
    variance_ratio = eigenvalues / eigenvalues.sum()

    score_table = work[["island_id"]].copy()
    for pc_index in range(n_pcs):
        column = f"climate_pc{pc_index + 1}"
        score_table[column] = np.nan
        score_table.loc[complete, column] = scores[:, pc_index]

    loading_rows: list[dict[str, Any]] = []
    for variable_index, variable in enumerate(climate_columns):
        row: dict[str, Any] = {
            "variable": variable,
            "mean": float(means[variable]),
            "sd": float(sds[variable]),
        }
        for pc_index in range(n_pcs):
            row[f"climate_pc{pc_index + 1}_loading"] = float(loadings[variable_index, pc_index])
        loading_rows.append(row)

    metadata = {
        "n_islands": int(len(work)),
        "n_complete_for_pca": int(complete.sum()),
        "climate_columns": climate_columns,
        "n_pcs": n_pcs,
        "explained_variance_ratio": {
            f"climate_pc{i + 1}": float(variance_ratio[i]) for i in range(n_pcs)
        },
        "cumulative_explained_variance": float(variance_ratio[:n_pcs].sum()),
        "orientation_rule": "largest-absolute-loading-positive",
    }
    return score_table, pd.DataFrame(loading_rows), metadata


def build_covariates(
    island_environment: pd.DataFrame,
    environment_samples: pd.DataFrame,
    island_manifest: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    if "island_id" not in island_manifest.columns or "area_km2" not in island_manifest.columns:
        raise typer.BadParameter("island manifest must contain island_id and area_km2")
    if island_manifest["island_id"].duplicated().any():
        raise typer.BadParameter("island manifest must contain one row per island_id")

    climate_cfg = config["climate"]
    climate_columns = [str(value) for value in climate_cfg["columns"]]
    n_pcs = int(climate_cfg["n_pcs"])
    block_degrees = float(config["spatial_inference"]["block_degrees"])

    geography = _island_geography(environment_samples, block_degrees)
    pca_scores, pca_loadings, pca_metadata = _fit_climate_pca(
        island_environment, climate_columns, n_pcs
    )
    area = island_manifest[["island_id", "area_km2"]].copy()
    area["area_km2"] = pd.to_numeric(area["area_km2"], errors="coerce")
    area["log_island_area_km2"] = np.log1p(area["area_km2"])

    covariates = area.merge(geography, on="island_id", how="outer", validate="one_to_one")
    covariates = covariates.merge(pca_scores, on="island_id", how="outer", validate="one_to_one")
    covariates = covariates.sort_values("island_id").reset_index(drop=True)

    expected_min = float(climate_cfg.get("expected_min_cumulative_variance", 0.0))
    if pca_metadata["cumulative_explained_variance"] < expected_min:
        raise RuntimeError(
            "climate PCA retained less variance than declared minimum: "
            f"{pca_metadata['cumulative_explained_variance']:.3f} < {expected_min:.3f}"
        )

    manifest = {
        "contract": "full_analysis_covariates_v1",
        "n_islands": int(covariates["island_id"].nunique()),
        "n_area_nonmissing": int(covariates["area_km2"].notna().sum()),
        "n_geography_nonmissing": int(covariates["abs_latitude"].notna().sum()),
        "n_spatial_blocks": int(covariates["spatial_block"].replace("missing", pd.NA).nunique()),
        "spatial_block_degrees": block_degrees,
        "pca": pca_metadata,
        "isolation_status": "not_connected_no_ad_hoc_proxy_used",
    }
    return covariates, pca_loadings, manifest


@app.command("run")
def run(
    island_environment_csv: Path = typer.Option(..., exists=True),
    environment_samples_csv: Path = typer.Option(..., exists=True),
    island_manifest_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(Path("config/m0_m4_full_analysis.yml"), exists=True),
) -> None:
    config = _load_config(config_path)
    environment = pd.read_csv(island_environment_csv)
    samples = pd.read_csv(environment_samples_csv)
    island_manifest = pd.read_csv(island_manifest_csv)
    covariates, loadings, manifest = build_covariates(
        environment, samples, island_manifest, config
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    covariates.to_csv(output_dir / "full_analysis_covariates.csv", index=False)
    loadings.to_csv(output_dir / "climate_pca_loadings.csv", index=False)
    (output_dir / "full_analysis_covariate_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
