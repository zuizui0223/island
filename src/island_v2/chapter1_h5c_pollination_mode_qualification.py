"""Outcome-blind qualification for Chapter 1 H5c pollination-mode specificity.

This module never fits the observed biotic-vs-wind interaction. It uses independently
reported ``pollen_vector_mode`` only to define coverage/support, then simulates response
values on the realized island, mode, covariate, species-support and spatial-block design.
Observed H5c effects may be opened only if the frozen qualification gate passes.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _standardize(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.nanmean(x))
    sd = float(np.nanstd(x))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid simulated-design predictor")
    return (x - mean) / sd


def collapse_pollination_mode(gift: pd.DataFrame) -> pd.DataFrame:
    required = {"scientific_name", "trait_name", "trait_value"}
    missing = required - set(gift.columns)
    if missing:
        raise typer.BadParameter(f"GIFT prepared table missing columns: {sorted(missing)}")
    x = gift.loc[gift["trait_name"].astype(str).eq("pollen_vector_mode")].copy()
    x["scientific_name"] = x["scientific_name"].fillna("").astype(str).str.strip()
    x["trait_value"] = x["trait_value"].fillna("").astype(str).str.strip()
    x = x.loc[x["scientific_name"].ne("")]
    rows: list[dict[str, str]] = []
    for species, group in x.groupby("scientific_name", sort=True):
        states = set(group["trait_value"])
        if "mixed" in states:
            continue
        allowed = states & {"biotic", "abiotic_wind"}
        if len(allowed) != 1 or states - {"biotic", "abiotic_wind"}:
            continue
        rows.append({"accepted_species": str(species), "pollen_vector_mode": next(iter(allowed))})
    return pd.DataFrame(rows, columns=["accepted_species", "pollen_vector_mode"])


def build_mode_design(
    *,
    status_flora: pd.DataFrame,
    eligible_species: pd.DataFrame,
    mode_species: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    stratum: str,
    context_layer: str,
    context: str,
    minimum_species_per_island_mode: int,
) -> pd.DataFrame:
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
        "floristic_status",
    }
    if missing := required_flora - set(status_flora.columns):
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    eligible = set(
        eligible_species.loc[
            eligible_species["syndrome"].astype(str).eq("generalized_accessible"),
            "accepted_species",
        ].astype(str)
    )
    flora = status_flora.loc[stratum_mask(status_flora, stratum)].copy()
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.loc[flora["accepted_species"].isin(eligible)]
    flora = flora.merge(mode_species, on="accepted_species", how="inner", validate="many_to_one")
    flora = flora.drop_duplicates(["island_id", "accepted_species", "pollen_vector_mode"])
    counts = (
        flora.groupby(["island_id", "pollen_vector_mode"], as_index=False)["accepted_species"]
        .nunique()
        .rename(columns={"accepted_species": "n_species"})
    )
    counts = counts.loc[counts["n_species"].ge(int(minimum_species_per_island_mode))].copy()

    geography = "log_distance_to_continent_km"
    baselines = ["log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"]
    cov_required = {"island_id", geography, "analysis_regime", "spatial_block", *baselines}
    if missing := cov_required - set(covariates.columns):
        raise typer.BadParameter(f"covariates missing columns: {sorted(missing)}")
    cov = covariates[list(cov_required)].drop_duplicates("island_id").copy()
    if context_layer == "biogeographic_realm":
        if not {"island_id", "biogeographic_realm"} <= set(realm_assignment.columns):
            raise typer.BadParameter("realm assignment missing required columns")
        cov = cov.merge(
            realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="one_to_one",
        )
    elif context_layer != "analysis_regime":
        raise ValueError(f"unsupported context layer: {context_layer}")
    out = counts.merge(cov, on="island_id", how="left", validate="many_to_one")
    out = out.loc[out[context_layer].astype(str).eq(str(context))].copy()
    numeric = ["n_species", geography, *baselines]
    for column in numeric:
        out[column] = pd.to_numeric(out[column], errors="coerce")
    out["spatial_block"] = out["spatial_block"].fillna("").astype(str)
    out = out.dropna(subset=numeric)
    out = out.loc[out["spatial_block"].ne("")]
    return out


def _support_summary(design: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for mode in ("biotic", "abiotic_wind"):
        x = design.loc[design["pollen_vector_mode"].astype(str).eq(mode)]
        summary[f"n_islands_{mode}"] = int(x["island_id"].nunique())
        summary[f"n_blocks_{mode}"] = int(x["spatial_block"].nunique())
        summary[f"n_species_sum_{mode}"] = int(pd.to_numeric(x["n_species"], errors="coerce").sum())
    summary["n_design_rows"] = int(len(design))
    summary["n_unique_islands"] = int(design["island_id"].nunique())
    return summary


def _simulation_design(design: pd.DataFrame) -> tuple[np.ndarray, list[str], np.ndarray, np.ndarray, np.ndarray]:
    geo = _standardize(design["log_distance_to_continent_km"])
    biotic = design["pollen_vector_mode"].astype(str).eq("biotic").to_numpy(float)
    columns = [np.ones(len(design)), geo, biotic, geo * biotic]
    names = ["intercept", "z_distance", "biotic", "distance_x_biotic"]
    for predictor in ("log_island_area_km2", "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"):
        columns.append(_standardize(design[predictor]))
        names.append(f"z_{predictor}")
    X = np.column_stack(columns)
    weights = pd.to_numeric(design["n_species"], errors="coerce").to_numpy(float)
    clusters = design["spatial_block"].astype(str).to_numpy()
    return X, names, weights, clusters, geo * biotic


def simulate_cell(
    design: pd.DataFrame,
    *,
    effects: list[float],
    replicates: int,
    seed: int,
    cluster_sd: float,
    residual_sd: float,
    alpha: float,
) -> pd.DataFrame:
    X, names, weights, clusters, interaction = _simulation_design(design)
    unique_clusters = np.unique(clusters)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    interaction_index = names.index("distance_x_biotic")
    for effect in effects:
        detected = 0
        fit_count = 0
        for _ in range(int(replicates)):
            cluster_values = {
                cluster: rng.normal(0.0, float(cluster_sd)) for cluster in unique_clusters
            }
            cluster_noise = np.asarray([cluster_values[c] for c in clusters], dtype=float)
            residual = rng.normal(0.0, float(residual_sd) / np.sqrt(np.maximum(weights, 1.0)))
            y = float(effect) * interaction + cluster_noise + residual
            coef, _, fit = _fit_weighted_clustered_design(y, weights, X, names, clusters)
            if coef.empty or fit.get("status") != "fit":
                continue
            fit_count += 1
            row = coef.iloc[interaction_index]
            p_value = float(row["p_value"])
            estimate = float(row["estimate"])
            if np.isfinite(p_value) and p_value <= float(alpha) and (effect == 0 or estimate > 0):
                detected += 1
        rows.append(
            {
                "interaction_effect_sd": float(effect),
                "n_replicates": int(replicates),
                "n_fit": int(fit_count),
                "detection_rate": float(detected / fit_count) if fit_count else float("nan"),
            }
        )
    return pd.DataFrame(rows)


def run_qualification(
    *,
    artifact_root: Path,
    gift_prepared_csv: Path,
    config_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if config.get("contract") != "chapter1_h5c_pollination_mode_specificity_qualification_v1":
        raise ValueError("unexpected H5c contract")
    status = pd.read_csv(artifact_root / "fixed/canonical/input/chapter1_status_flora.csv.gz")
    covariates = pd.read_csv(artifact_root / "fixed/isolation/results/purpose_shortest_island_data.csv")
    realm = pd.read_csv(artifact_root / "fixed/realm/realm/island_biogeographic_realm_assignment.csv")
    gift = pd.read_csv(gift_prepared_csv)
    modes = collapse_pollination_mode(gift)

    scope_paths = {
        "all_analysis_eligible": artifact_root / "syndrome/all/species_syndrome_concordance.csv.gz",
        "direct_only": artifact_root / "syndrome/direct/species_syndrome_concordance.csv.gz",
    }
    support_rows: list[dict[str, Any]] = []
    power_rows: list[pd.DataFrame] = []
    qual_rows: list[dict[str, Any]] = []
    sim = config["simulation"]
    gate = config["support_gate"]
    qgate = config["qualification"]
    cell_index = 0
    for scope in config["evidence_scopes"]:
        eligible = pd.read_csv(scope_paths[str(scope)])
        for stratum in config["strata"]:
            for context_spec in config["contexts"]:
                layer = str(context_spec["context_layer"])
                context = str(context_spec["context"])
                design = build_mode_design(
                    status_flora=status,
                    eligible_species=eligible,
                    mode_species=modes,
                    covariates=covariates,
                    realm_assignment=realm,
                    stratum=str(stratum),
                    context_layer=layer,
                    context=context,
                    minimum_species_per_island_mode=int(gate["minimum_species_per_island_mode"]),
                )
                support = _support_summary(design)
                support_ok = (
                    support["n_islands_biotic"] >= int(gate["minimum_islands_per_mode"])
                    and support["n_islands_abiotic_wind"] >= int(gate["minimum_islands_per_mode"])
                    and support["n_blocks_biotic"] >= int(gate["minimum_spatial_blocks_per_mode"])
                    and support["n_blocks_abiotic_wind"] >= int(gate["minimum_spatial_blocks_per_mode"])
                )
                identity = {
                    "evidence_scope": str(scope),
                    "stratum": str(stratum),
                    "context_layer": layer,
                    "context": context,
                }
                support_rows.append({**identity, **support, "support_gate_pass": bool(support_ok)})
                if support_ok:
                    power = simulate_cell(
                        design,
                        effects=[float(x) for x in sim["interaction_effect_sd_grid"]],
                        replicates=int(sim["replicates"]),
                        seed=int(sim["seed"]) + cell_index * 1000,
                        cluster_sd=float(sim["cluster_random_effect_sd"]),
                        residual_sd=float(sim["residual_sd"]),
                        alpha=float(sim["alpha"]),
                    )
                    for key, value in reversed(list(identity.items())):
                        power.insert(0, key, value)
                    power_rows.append(power)
                    null = power.loc[np.isclose(power["interaction_effect_sd"], 0.0)].iloc[0]
                    target = power.loc[
                        np.isclose(power["interaction_effect_sd"], float(qgate["target_effect_sd"]))
                    ].iloc[0]
                    type1 = float(null["detection_rate"])
                    target_power = float(target["detection_rate"])
                    qualified = (
                        type1 <= float(qgate["maximum_type1_error"])
                        and target_power >= float(qgate["minimum_target_recovery"])
                    )
                else:
                    type1 = float("nan")
                    target_power = float("nan")
                    qualified = False
                qual_rows.append(
                    {
                        **identity,
                        "support_gate_pass": bool(support_ok),
                        "type1_error": type1,
                        "target_effect_sd": float(qgate["target_effect_sd"]),
                        "target_recovery": target_power,
                        "qualified": bool(qualified),
                    }
                )
                cell_index += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    support_df = pd.DataFrame(support_rows)
    power_df = pd.concat(power_rows, ignore_index=True) if power_rows else pd.DataFrame()
    qualification_df = pd.DataFrame(qual_rows)
    support_df.to_csv(output_dir / "h5c_pollination_mode_support.csv", index=False)
    power_df.to_csv(output_dir / "h5c_pollination_mode_power.csv", index=False)
    qualification_df.to_csv(output_dir / "h5c_pollination_mode_qualification.csv", index=False)
    modes.to_csv(output_dir / "h5c_independent_pollination_mode_species.csv.gz", index=False, compression="gzip")
    manifest = {
        "contract": config["contract"],
        "observed_effect_opened": False,
        "n_independent_mode_species": int(len(modes)),
        "n_biotic_species": int((modes["pollen_vector_mode"] == "biotic").sum()),
        "n_wind_species": int((modes["pollen_vector_mode"] == "abiotic_wind").sum()),
        "n_design_cells": int(len(qualification_df)),
        "n_support_qualified_cells": int(qualification_df["support_gate_pass"].sum()),
        "n_fully_qualified_cells": int(qualification_df["qualified"].sum()),
        "any_observed_test_may_open": bool(qualification_df["qualified"].any()),
        "claim_ceiling": config["claim_ceiling"],
    }
    (output_dir / "h5c_pollination_mode_qualification_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    gift_prepared_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    config_path: Path = typer.Option(
        Path("config/chapter1_h5c_pollination_mode_qualification.yml"), exists=True
    ),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(
        json.dumps(
            run_qualification(
                artifact_root=artifact_root,
                gift_prepared_csv=gift_prepared_csv,
                config_path=config_path,
                output_dir=output_dir,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
