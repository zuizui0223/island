"""Single prequalified observed H5c pollination-mode specificity test."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _standardize(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def build_observed_design(
    *,
    status_flora: pd.DataFrame,
    species_scores: pd.DataFrame,
    mode_species: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
) -> pd.DataFrame:
    scores = species_scores.loc[
        species_scores["syndrome"].astype(str).eq("generalized_accessible")
    ][["accepted_species", "syndrome_concordance"]].copy()
    scores["accepted_species"] = scores["accepted_species"].astype(str)
    scores["syndrome_concordance"] = pd.to_numeric(scores["syndrome_concordance"], errors="coerce")
    scores = scores.dropna().drop_duplicates("accepted_species")

    modes = mode_species[["accepted_species", "pollen_vector_mode"]].copy()
    modes["accepted_species"] = modes["accepted_species"].astype(str)
    modes = modes.loc[modes["pollen_vector_mode"].astype(str).isin(["biotic", "abiotic_wind"])]
    modes = modes.drop_duplicates("accepted_species")

    flora = status_flora.loc[stratum_mask(status_flora, "native_nonendemic")].copy()
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    merged = flora.merge(scores, on="accepted_species", how="inner", validate="many_to_one")
    merged = merged.merge(modes, on="accepted_species", how="inner", validate="many_to_one")
    merged = merged.drop_duplicates(["island_id", "accepted_species", "pollen_vector_mode"])
    design = (
        merged.groupby(["island_id", "pollen_vector_mode"], as_index=False)
        .agg(response=("syndrome_concordance", "mean"), n_species=("accepted_species", "nunique"))
    )

    baselines = ["log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    required = {"island_id", "log_distance_to_continent_km", "spatial_block", *baselines}
    if missing := required - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    cov = covariates[list(required)].drop_duplicates("island_id")
    design = design.merge(cov, on="island_id", how="left", validate="many_to_one")
    design = design.merge(
        realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    design = design.loc[design["biogeographic_realm"].astype(str).eq("Palearctic")].copy()
    numeric = ["response", "n_species", "log_distance_to_continent_km", *baselines]
    for column in numeric:
        design[column] = pd.to_numeric(design[column], errors="coerce")
    design["spatial_block"] = design["spatial_block"].fillna("").astype(str)
    design = design.dropna(subset=numeric)
    design = design.loc[design["spatial_block"].ne("")]
    return design.reset_index(drop=True)


def fit_observed(design: pd.DataFrame) -> dict[str, Any]:
    geo = _standardize(design["log_distance_to_continent_km"])
    biotic = design["pollen_vector_mode"].astype(str).eq("biotic").to_numpy(float)
    names = ["intercept", "z_distance", "biotic", "distance_x_biotic"]
    columns = [np.ones(len(design)), geo, biotic, geo * biotic]
    for predictor in ("log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"):
        names.append(f"z_{predictor}")
        columns.append(_standardize(design[predictor]))
    X = np.column_stack(columns)
    weights = design["n_species"].to_numpy(float)
    clusters = design["spatial_block"].astype(str).to_numpy()
    coef, covariance, fit = _fit_weighted_clustered_design(
        design["response"].to_numpy(float), weights, X, names, clusters
    )
    if coef.empty or fit.get("status") != "fit":
        raise RuntimeError(f"H5c observed fit failed: {fit}")
    indexed = coef.set_index("predictor")
    b_d = float(indexed.loc["z_distance", "estimate"])
    b_i = float(indexed.loc["distance_x_biotic", "estimate"])
    p_i = float(indexed.loc["distance_x_biotic", "p_value"])
    se_i = float(indexed.loc["distance_x_biotic", "cluster_robust_se"])
    i_d = names.index("z_distance")
    i_i = names.index("distance_x_biotic")
    wind_slope = b_d
    wind_var = float(covariance[i_d, i_d])
    biotic_slope = b_d + b_i
    biotic_var = float(covariance[i_d, i_d] + covariance[i_i, i_i] + 2 * covariance[i_d, i_i])
    wind_se = math.sqrt(max(wind_var, 0.0))
    biotic_se = math.sqrt(max(biotic_var, 0.0))
    if p_i <= 0.05 and b_i > 0:
        classification = "expected_biotic_specificity_supported"
    elif p_i <= 0.05 and b_i < 0:
        classification = "opposite_pollination_mode_specificity_supported"
    else:
        classification = "no_pollination_mode_specificity_support"
    return {
        "classification": classification,
        "interaction_estimate": b_i,
        "interaction_se": se_i,
        "interaction_ci_low": b_i - 1.96 * se_i,
        "interaction_ci_high": b_i + 1.96 * se_i,
        "interaction_p_value": p_i,
        "wind_distance_slope": wind_slope,
        "wind_distance_se": wind_se,
        "wind_distance_ci_low": wind_slope - 1.96 * wind_se,
        "wind_distance_ci_high": wind_slope + 1.96 * wind_se,
        "biotic_distance_slope": biotic_slope,
        "biotic_distance_se": biotic_se,
        "biotic_distance_ci_low": biotic_slope - 1.96 * biotic_se,
        "biotic_distance_ci_high": biotic_slope + 1.96 * biotic_se,
        "n_design_rows": int(len(design)),
        "n_unique_islands": int(design["island_id"].nunique()),
        "n_clusters": int(fit["n_clusters"]),
        "n_islands_biotic": int(design.loc[biotic == 1, "island_id"].nunique()),
        "n_islands_wind": int(design.loc[biotic == 0, "island_id"].nunique()),
        "n_species_sum_biotic": int(design.loc[biotic == 1, "n_species"].sum()),
        "n_species_sum_wind": int(design.loc[biotic == 0, "n_species"].sum()),
    }


def run_observed(
    *, artifact_root: Path, mode_species_csv: Path, config_path: Path, output_dir: Path
) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_h5c_observed_pollination_mode_specificity_v1":
        raise ValueError("unexpected H5c observed contract")
    status = pd.read_csv(artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz")
    scores = pd.read_csv(artifact_root / "syndrome/direct/species_syndrome_concordance.csv.gz")
    cov = pd.read_csv(artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv")
    realm = pd.read_csv(artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv")
    modes = pd.read_csv(mode_species_csv)
    design = build_observed_design(
        status_flora=status,
        species_scores=scores,
        mode_species=modes,
        covariates=cov,
        realm_assignment=realm,
    )
    result = fit_observed(design)
    output_dir.mkdir(parents=True, exist_ok=True)
    design.to_csv(output_dir / "h5c_observed_design.csv.gz", index=False, compression="gzip")
    payload = {
        "contract": config["contract"],
        "qualified_cell": config["only_cell"],
        **result,
        "claim_ceiling": config["claim_ceiling"],
        "N1_N2_reopened": False,
    }
    (output_dir / "h5c_observed_result.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([payload]).to_csv(output_dir / "h5c_observed_result.csv", index=False)
    return payload


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    mode_species_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(Path("config/chapter1_h5c_observed_specificity.yml"), exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(run_observed(artifact_root=artifact_root, mode_species_csv=mode_species_csv, config_path=config_path, output_dir=output_dir), indent=2))


if __name__ == "__main__":
    app()
