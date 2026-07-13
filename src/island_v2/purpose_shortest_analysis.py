"""Run the shortest purpose-aligned global analysis on the locked all-data tier.

This is intentionally simple. It restores mainland/continent distance as the v1
continuity variable, keeps the Northern mid-latitude Bombus test separate from
tropical and Southern extratropical regime comparisons, and avoids mediation or
sensitivity layers.

Interpretation guardrails:
- the legacy v1 all-land distance is retained for audit only;
- the fitted models use distance to seeded major continental landmasses, not distance
  to every Natural Earth land polygon;
- Bombus is environmental compatibility, not observed biological absence;
- tropical/Southern models use plant pollination-guild composition as exploratory
  functional-replacement correlates, not direct evidence of local pollinator abundance.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import geopandas as gpd
import numpy as np
import pandas as pd
import typer
from shapely.geometry import Point

from island_v2.m0_m4_full_analysis import fit_grouped_binomial_clustered

app = typer.Typer(add_completion=False, no_args_is_help=True)

TROPIC = 23.4366
ARCTIC = 66.5634

# Seed points identify the six major continental landmasses in Natural Earth 110m.
# This avoids the v1 implementation problem where distance to *all* land polygons can
# become zero when the focal island itself is represented as a Natural Earth polygon.
CONTINENT_SEEDS = {
    "africa": (20.0, 0.0),
    "eurasia": (80.0, 50.0),
    "north_america": (-100.0, 40.0),
    "south_america": (-60.0, -15.0),
    "australia": (135.0, -25.0),
    "antarctica": (0.0, -80.0),
}

PLAIN_COLORS = {"white", "green_brown_inconspicuous"}
CONSPICUOUS_COLORS = {"yellow_orange", "red_pink", "blue_purple"}
GENERALIZED_FORMS = {"open_radial", "brush_puff", "composite_head"}
SPECIALIZED_FORMS = {
    "tubular",
    "salverform",
    "funnel_trumpet",
    "urn_urceolate",
    "zygomorphic",
    "papilionaceous",
    "spurred",
}
POLLEN_VECTOR_VALUES = {"biotic", "abiotic_wind", "mixed"}
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
SHOWY_ALTERNATIVE_GUILDS = {"butterflies", "moths", "birds", "bats"}
OTHER_BEE_GUILDS = {"other_bees"}
GENERALIST_INSECT_GUILDS = {"flies", "beetles_wasps_ants", "mixed_or_generalist"}


def _contrast_counts(
    composition: pd.DataFrame,
    *,
    trait_name: str,
    positive_values: Iterable[str],
    denominator_values: Iterable[str],
    prefix: str,
) -> pd.DataFrame:
    work = composition.loc[
        composition["trait_name"].astype(str).eq(trait_name),
        ["island_id", "filled_value", "n_species"],
    ].copy()
    allowed = {str(value) for value in denominator_values}
    positive = {str(value) for value in positive_values}
    work = work.loc[work["filled_value"].astype(str).isin(allowed)].copy()
    work["n_species"] = pd.to_numeric(work["n_species"], errors="coerce")
    work = work.dropna(subset=["n_species"])
    work["is_positive"] = work["filled_value"].astype(str).isin(positive)
    total = work.groupby("island_id")["n_species"].sum().rename(f"{prefix}_trials")
    success = (
        work.loc[work["is_positive"]]
        .groupby("island_id")["n_species"]
        .sum()
        .rename(f"{prefix}_successes")
    )
    out = total.to_frame().join(success, how="left").fillna({f"{prefix}_successes": 0.0})
    out = out.reset_index()
    out[f"{prefix}_share"] = out[f"{prefix}_successes"] / out[f"{prefix}_trials"].replace(0, np.nan)
    return out


def compute_distance_metrics(
    islands_gpkg: Path,
    natural_earth_land_zip: Path,
) -> pd.DataFrame:
    """Return true continent distance plus the legacy v1 all-land comparator."""
    islands = gpd.read_file(islands_gpkg)
    if "island_id" not in islands.columns:
        raise typer.BadParameter("islands GPKG must contain island_id")
    land = gpd.read_file(f"zip://{natural_earth_land_zip.resolve()}")
    if islands.crs is None or land.crs is None:
        raise typer.BadParameter("island and Natural Earth layers must have CRS metadata")

    islands_m = islands[["island_id", "geometry"]].to_crs(3857)
    land_wgs = land.to_crs(4326)
    land_m = land_wgs.to_crs(3857)

    # Legacy comparator: minimum distance to every Natural Earth land polygon.
    legacy_land_union = land_m.geometry.union_all()
    legacy_distance_km = (
        islands_m.geometry.distance(legacy_land_union).to_numpy(dtype=float) / 1000.0
    )

    # Select only the major continental land polygons using fixed geographic seeds.
    selected_indices: set[object] = set()
    for name, (longitude, latitude) in CONTINENT_SEEDS.items():
        point = Point(longitude, latitude)
        covered = land_wgs.geometry.covers(point)
        if not bool(covered.any()):
            raise typer.BadParameter(
                f"Natural Earth land layer did not contain the {name} continent seed"
            )
        selected_indices.add(land_wgs.index[covered][0])

    continent_m = land_wgs.loc[sorted(selected_indices)].to_crs(3857)
    continent_union = continent_m.geometry.union_all()
    continent_distance_km = (
        islands_m.geometry.distance(continent_union).to_numpy(dtype=float) / 1000.0
    )

    return pd.DataFrame(
        {
            "island_id": islands_m["island_id"].astype(str).to_numpy(),
            "distance_to_continent_km": continent_distance_km,
            "log_distance_to_continent_km": np.log1p(continent_distance_km),
            "legacy_v1_all_land_distance_km": legacy_distance_km,
            "legacy_v1_log_all_land_distance_km": np.log1p(legacy_distance_km),
        }
    )


def _assign_regime(latitude: pd.Series) -> pd.Series:
    lat = pd.to_numeric(latitude, errors="coerce")
    regime = pd.Series("unresolved", index=lat.index, dtype="object")
    regime.loc[lat.abs() < TROPIC] = "tropical"
    regime.loc[(lat >= TROPIC) & (lat < ARCTIC)] = "northern_midlatitude"
    regime.loc[lat <= -TROPIC] = "southern_extratropical"
    regime.loc[lat >= ARCTIC] = "northern_high_latitude"
    return regime


def _fit(
    data: pd.DataFrame,
    *,
    model: str,
    regime: str,
    outcome: str,
    successes: str,
    trials: str,
    predictors: list[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    subset = data.loc[data["analysis_regime"].eq(regime)].copy()
    coefficients, fit = fit_grouped_binomial_clustered(
        subset,
        successes_column=successes,
        trials_column=trials,
        predictors=predictors,
        cluster_column="spatial_block",
        z_threshold=1.96,
    )
    if not coefficients.empty:
        coefficients.insert(0, "outcome", outcome)
        coefficients.insert(0, "regime", regime)
        coefficients.insert(0, "model", model)
    fit = {"model": model, "regime": regime, "outcome": outcome, **fit}
    return coefficients, fit


def run_analysis(
    composition: pd.DataFrame,
    island_metrics: pd.DataFrame,
    covariates: pd.DataFrame,
    distance: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = (
        island_metrics.drop_duplicates("island_id")
        .merge(covariates.drop_duplicates("island_id"), on="island_id", how="left", validate="one_to_one")
        .merge(distance.drop_duplicates("island_id"), on="island_id", how="left", validate="one_to_one")
    )
    base["analysis_regime"] = _assign_regime(base["island_latitude"])
    base["bombus_channel_state"] = pd.to_numeric(
        base["environmental_compatibility_max"], errors="coerce"
    )

    contrasts = [
        _contrast_counts(
            composition,
            trait_name="flower_primary_color",
            positive_values=PLAIN_COLORS,
            denominator_values=PLAIN_COLORS | CONSPICUOUS_COLORS,
            prefix="plain_color",
        ),
        _contrast_counts(
            composition,
            trait_name="flower_primary_color",
            positive_values=CONSPICUOUS_COLORS,
            denominator_values=PLAIN_COLORS | CONSPICUOUS_COLORS,
            prefix="conspicuous_color",
        ),
        _contrast_counts(
            composition,
            trait_name="floral_form",
            positive_values=GENERALIZED_FORMS,
            denominator_values=GENERALIZED_FORMS | SPECIALIZED_FORMS,
            prefix="generalized_form",
        ),
        _contrast_counts(
            composition,
            trait_name="pollen_vector_mode",
            positive_values={"abiotic_wind", "mixed"},
            denominator_values=POLLEN_VECTOR_VALUES,
            prefix="wind_mixed",
        ),
        _contrast_counts(
            composition,
            trait_name="pollination_functional_guild",
            positive_values=SHOWY_ALTERNATIVE_GUILDS,
            denominator_values=ANIMAL_GUILDS,
            prefix="showy_alternative_guild",
        ),
        _contrast_counts(
            composition,
            trait_name="pollination_functional_guild",
            positive_values=OTHER_BEE_GUILDS,
            denominator_values=ANIMAL_GUILDS,
            prefix="other_bee_guild",
        ),
        _contrast_counts(
            composition,
            trait_name="pollination_functional_guild",
            positive_values=GENERALIST_INSECT_GUILDS,
            denominator_values=ANIMAL_GUILDS,
            prefix="generalist_insect_guild",
        ),
    ]
    data = base.copy()
    for table in contrasts:
        data = data.merge(table, on="island_id", how="left", validate="one_to_one")

    baseline = [
        "log_distance_to_continent_km",
        "log_island_area_km2",
        "climate_pc1",
        "climate_pc2",
        "climate_pc3",
        "climate_pc4",
    ]

    coefficient_tables: list[pd.DataFrame] = []
    fits: list[dict[str, object]] = []
    for outcome, successes, trials in (
        ("plain_color", "plain_color_successes", "plain_color_trials"),
        ("generalized_form", "generalized_form_successes", "generalized_form_trials"),
    ):
        for model, predictors in (
            ("M0_geography", baseline),
            ("M1_add_bombus", [*baseline, "bombus_channel_state"]),
        ):
            coef, fit = _fit(
                data,
                model=model,
                regime="northern_midlatitude",
                outcome=outcome,
                successes=successes,
                trials=trials,
                predictors=predictors,
            )
            coefficient_tables.append(coef)
            fits.append(fit)

    for regime in ("tropical", "southern_extratropical"):
        coef, fit = _fit(
            data,
            model="M4_replacement_shortest",
            regime=regime,
            outcome="conspicuous_color",
            successes="conspicuous_color_successes",
            trials="conspicuous_color_trials",
            predictors=[
                *baseline,
                "showy_alternative_guild_share",
                "other_bee_guild_share",
                "generalist_insect_guild_share",
                "wind_mixed_share",
            ],
        )
        coefficient_tables.append(coef)
        fits.append(fit)

    coefficients = pd.concat(
        [table for table in coefficient_tables if not table.empty], ignore_index=True
    )
    fit_table = pd.DataFrame(fits)

    region_summary = (
        data.groupby("analysis_regime", dropna=False)
        .agg(
            n_islands=("island_id", "nunique"),
            mean_distance_to_continent_km=("distance_to_continent_km", "mean"),
            mean_legacy_v1_all_land_distance_km=("legacy_v1_all_land_distance_km", "mean"),
            mean_bombus_compatibility=("bombus_channel_state", "mean"),
            mean_plain_color_share=("plain_color_share", "mean"),
            mean_conspicuous_color_share=("conspicuous_color_share", "mean"),
            mean_generalized_form_share=("generalized_form_share", "mean"),
            mean_wind_mixed_share=("wind_mixed_share", "mean"),
            mean_showy_alternative_guild_share=("showy_alternative_guild_share", "mean"),
            mean_other_bee_guild_share=("other_bee_guild_share", "mean"),
            mean_generalist_insect_guild_share=("generalist_insect_guild_share", "mean"),
        )
        .reset_index()
    )
    return data, coefficients, fit_table, region_summary


@app.command("run")
def run(
    trait_composition_csv: Path = typer.Option(..., exists=True),
    island_metrics_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    islands_gpkg: Path = typer.Option(..., exists=True),
    natural_earth_land_zip: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    composition = pd.read_csv(trait_composition_csv)
    island_metrics = pd.read_csv(island_metrics_csv)
    covariates = pd.read_csv(covariates_csv)
    distance = compute_distance_metrics(islands_gpkg, natural_earth_land_zip)
    data, coefficients, fits, region_summary = run_analysis(
        composition, island_metrics, covariates, distance
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    distance.to_csv(output_dir / "continent_and_legacy_distance.csv", index=False)
    data.to_csv(output_dir / "purpose_shortest_island_data.csv", index=False)
    coefficients.to_csv(output_dir / "purpose_shortest_coefficients.csv", index=False)
    fits.to_csv(output_dir / "purpose_shortest_model_fit.csv", index=False)
    region_summary.to_csv(output_dir / "purpose_shortest_region_summary.csv", index=False)

    summary = {
        "contract": "purpose_shortest_distance_regime_analysis_v3",
        "analysis_tier": "sensitivity_all_all_available_filled_data",
        "n_islands": int(data["island_id"].nunique()),
        "distance_definition": {
            "model_predictor": "distance_to_continent_km",
            "continents": sorted(CONTINENT_SEEDS),
            "legacy_all_land_distance_retained_for_audit": True,
            "legacy_zero_distance_share": float(
                (distance["legacy_v1_all_land_distance_km"] == 0).mean()
            ),
            "continent_zero_distance_share": float(
                (distance["distance_to_continent_km"] == 0).mean()
            ),
        },
        "regime_counts": {
            str(key): int(value)
            for key, value in data["analysis_regime"].value_counts(dropna=False).items()
        },
        "models": [
            "northern_midlatitude: M0 geography -> plain color / generalized form",
            "northern_midlatitude: M1 geography + Bombus environmental compatibility",
            "tropical: conspicuous color ~ geography + showy alternatives + other bees + generalist insects + wind/mixed",
            "southern_extratropical: same replacement model",
        ],
        "guardrails": [
            "fitted distance is to seeded major continental landmasses; the v1 all-land metric is audit-only",
            "northern_midlatitude is a geographic proxy, not a frozen global Bombus applicability registry",
            "Bombus predictor is environmental compatibility, not observed biological absence",
            "alternative guild shares are plant pollination-guild composition, not direct island pollinator abundance",
            "all fill tiers including global fallback are used because this is the requested first all-data exploratory run",
        ],
    }
    (output_dir / "purpose_shortest_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(summary))


if __name__ == "__main__":
    app()
