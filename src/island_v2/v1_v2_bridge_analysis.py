"""Minimal v1-to-v2 bridge analysis on the locked v2 species-trait master.

This is not a sensitivity-analysis grid. It asks one narrow question: does the v1
geographic sign reappear when the v2 data are restricted and recoded to match the
v1 analysis as closely as the v2 ontology permits?
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from island_v2.m0_m4_full_analysis import fit_grouped_binomial_clustered
from island_v2.purpose_shortest_analysis import (
    TROPIC,
    compute_distance_metrics,
    replace_trait_richness,
)
from island_v2.trait_bombus_analysis import (
    GROUNDED_FILL_TIERS,
    load_trait_rows_for_grounded_analysis,
    require_grounded_trait_rows,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

ANIMAL_GUILDS = {
    "bumblebees",
    "other_bees",
    "flies",
    "butterflies",
    "moths",
    "birds",
    "bats",
    "beetles_wasps_ants",
    "mixed_or_generalist",
}
V1_LIKE_PLAIN = {"white", "green_brown_inconspicuous", "other_described"}
V1_LIKE_CONSPICUOUS = {"yellow_orange", "red_pink", "blue_purple"}
V1_LIKE_GENERALIZED = {
    "open_radial",
    "bell_campanulate",
    "brush_puff",
    "composite_head",
    "reduced_wind",
}
V1_LIKE_SPECIALIZED = {
    "tubular",
    "salverform",
    "funnel_trumpet",
    "urn_urceolate",
    "zygomorphic",
    "papilionaceous",
    "spurred",
}


def build_v1_like_animal_counts(species_master: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "trait_name", "filled_value"}
    missing = required.difference(species_master.columns)
    if missing:
        raise typer.BadParameter(f"species master missing columns: {sorted(missing)}")

    species_master = require_grounded_trait_rows(species_master)
    keep_traits = {"pollination_functional_guild", "flower_primary_color", "floral_form"}
    work = species_master.loc[
        species_master["trait_name"].astype(str).isin(keep_traits),
        ["island_id", "accepted_species", "trait_name", "filled_value"],
    ].copy()
    work = work.drop_duplicates(["island_id", "accepted_species", "trait_name"])
    wide = work.pivot(
        index=["island_id", "accepted_species"],
        columns="trait_name",
        values="filled_value",
    ).reset_index()

    animal = wide.loc[
        wide["pollination_functional_guild"].astype(str).isin(ANIMAL_GUILDS)
    ].copy()

    color_allowed = V1_LIKE_PLAIN | V1_LIKE_CONSPICUOUS
    form_allowed = V1_LIKE_GENERALIZED | V1_LIKE_SPECIALIZED
    animal["color_usable"] = animal["flower_primary_color"].astype(str).isin(color_allowed)
    animal["plain"] = animal["flower_primary_color"].astype(str).isin(V1_LIKE_PLAIN)
    animal["form_usable"] = animal["floral_form"].astype(str).isin(form_allowed)
    animal["generalized"] = animal["floral_form"].astype(str).isin(V1_LIKE_GENERALIZED)

    counts = animal.groupby("island_id", sort=True).agg(
        n_animal_species=("accepted_species", "nunique"),
        plain_color_trials=("color_usable", "sum"),
        plain_color_successes=("plain", lambda x: int(x[animal.loc[x.index, "color_usable"]].sum())),
        generalized_form_trials=("form_usable", "sum"),
        generalized_form_successes=("generalized", lambda x: int(x[animal.loc[x.index, "form_usable"]].sum())),
    ).reset_index()
    return counts


def fit_bridge(
    data: pd.DataFrame,
    *,
    outcome: str,
    successes: str,
    trials: str,
    distance_predictor: str,
    model: str,
) -> tuple[pd.DataFrame, dict[str, object]]:
    coefficients, fit = fit_grouped_binomial_clustered(
        data,
        successes_column=successes,
        trials_column=trials,
        predictors=[distance_predictor, "log_island_area_km2"],
        cluster_column="spatial_block",
        z_threshold=1.96,
    )
    if not coefficients.empty:
        coefficients.insert(0, "outcome", outcome)
        coefficients.insert(0, "model", model)
    return coefficients, {"model": model, "outcome": outcome, **fit}


@app.command("run")
def run(
    species_master_csv: Path = typer.Option(..., exists=True),
    island_metrics_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    natural_earth_land_zip: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    master = load_trait_rows_for_grounded_analysis(species_master_csv)
    grounded = require_grounded_trait_rows(master)
    grounded_richness = (
        grounded.groupby("island_id", sort=True)["accepted_species"]
        .nunique()
        .rename("n_trait_species")
        .reset_index()
    )
    metrics = replace_trait_richness(pd.read_csv(island_metrics_csv), grounded_richness)
    covariates = pd.read_csv(covariates_csv)
    distance = compute_distance_metrics(islands_gpkg, natural_earth_land_zip)
    counts = build_v1_like_animal_counts(grounded)

    data = (
        metrics.drop_duplicates("island_id")
        .merge(covariates.drop_duplicates("island_id"), on="island_id", how="left", validate="one_to_one")
        .merge(distance.drop_duplicates("island_id"), on="island_id", how="left", validate="one_to_one")
        .merge(counts, on="island_id", how="left", validate="one_to_one")
    )

    # Match the v1 geographic analysis set as closely as possible:
    # northern non-tropical islands, >=20 km2, and legacy distance > 0.
    bridge = data.loc[
        (pd.to_numeric(data["island_latitude"], errors="coerce") >= TROPIC)
        & (pd.to_numeric(data["area_km2"], errors="coerce") >= 20.0)
        & (pd.to_numeric(data["legacy_v1_all_land_distance_km"], errors="coerce") > 0.0)
    ].copy()

    coefficient_tables: list[pd.DataFrame] = []
    fits: list[dict[str, object]] = []
    for outcome, successes, trials in (
        ("plain_color", "plain_color_successes", "plain_color_trials"),
        ("generalized_form", "generalized_form_successes", "generalized_form_trials"),
    ):
        for model, distance_predictor in (
            ("v1_like_legacy_distance", "legacy_v1_log_all_land_distance_km"),
            ("same_scope_corrected_continent_distance", "log_distance_to_continent_km"),
        ):
            coef, fit = fit_bridge(
                bridge,
                outcome=outcome,
                successes=successes,
                trials=trials,
                distance_predictor=distance_predictor,
                model=model,
            )
            coefficient_tables.append(coef)
            fits.append(fit)

    coefficients = pd.concat(
        [table for table in coefficient_tables if not table.empty], ignore_index=True
    )
    fit_table = pd.DataFrame(fits)

    output_dir.mkdir(parents=True, exist_ok=True)
    bridge.to_csv(output_dir / "v1_v2_bridge_island_data.csv", index=False)
    coefficients.to_csv(output_dir / "v1_v2_bridge_coefficients.csv", index=False)
    fit_table.to_csv(output_dir / "v1_v2_bridge_model_fit.csv", index=False)

    key = coefficients.loc[
        coefficients["predictor"].isin(
            ["legacy_v1_log_all_land_distance_km", "log_distance_to_continent_km"]
        )
    ].to_dict("records")
    summary = {
        "contract": "v1_v2_minimal_bridge_v1",
        "analysis_tier": "sensitivity_all_grounded_available_values",
        "grounded_input_filter": {
            "required_fill_tiers": sorted(GROUNDED_FILL_TIERS),
            "n_species_master_rows": int(len(master)),
            "n_grounded_rows": int(len(grounded)),
            "n_excluded_rows": int(len(master) - len(grounded)),
            "forbidden_rows_after_filter": int(
                (~grounded["fill_tier"].astype(str).isin(GROUNDED_FILL_TIERS)).sum()
            ),
        },
        "scope": "northern_non_tropical_area_ge_20km2_legacy_distance_gt_0",
        "plant_filter": "animal_pollinated_functional_guilds_only",
        "model": "grouped binomial outcome ~ standardized distance + standardized island area",
        "n_scope_islands": int(bridge["island_id"].nunique()),
        "distance_effects": key,
        "interpretation": (
            "This bridge is for sign comparability only. It approximates v1 trait categories using the v2 ontology; "
            "it does not recreate the original v1 trait table or PDI predictor."
        ),
    }
    (output_dir / "v1_v2_bridge_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
