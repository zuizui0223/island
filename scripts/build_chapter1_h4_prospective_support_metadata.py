"""Build outcome-blind prospective H4 support metadata from Methods candidates."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_h5_glopl_global_distance import build_preflight_table

app = typer.Typer(add_completion=False, no_args_is_help=True)

PARENT_GLOPL_CONFIG = Path("config/chapter1_h5_glopl_global_distance_v1.yml")
FROZEN_DISTANCE_MEAN = 1.1074920173877423
FROZEN_DISTANCE_SD = 2.4207834708622444


def _key(*values: object) -> str:
    payload = "|".join(str(value) for value in values)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def _parent_config() -> dict:
    value = yaml.safe_load(PARENT_GLOPL_CONFIG.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise typer.BadParameter("invalid parent GloPL config")
    return value


def build_support_metadata(
    methods: pd.DataFrame,
    traits: pd.DataFrame,
    land: gpd.GeoDataFrame,
) -> pd.DataFrame:
    needed = {
        "doi",
        "pmcid",
        "publication_date",
        "title",
        "matched_species",
        "coordinate_pairs_json",
        "candidate_status",
    }
    if missing := needed - set(methods.columns):
        raise typer.BadParameter(f"Methods screen missing columns: {sorted(missing)}")
    if "accepted_species" not in traits.columns:
        raise typer.BadParameter("trait table missing accepted_species")

    frozen_species = set(traits["accepted_species"].dropna().astype(str))
    expanded: list[dict[str, object]] = []
    admitted = methods.loc[
        methods["candidate_status"].astype(str).eq("automatic_methods_candidate")
    ].copy()
    for row in admitted.itertuples(index=False):
        species_values = [
            value.strip()
            for value in str(row.matched_species or "").split("|")
            if value.strip()
        ]
        species_values = [value for value in species_values if value in frozen_species]
        try:
            coordinates = json.loads(str(row.coordinate_pairs_json or "[]"))
        except json.JSONDecodeError as exc:
            raise typer.BadParameter("invalid coordinate_pairs_json") from exc
        publication_id = str(row.doi or "").strip().casefold() or str(row.pmcid or "").strip()
        if not publication_id:
            continue
        date = str(row.publication_date or "").strip()
        year = date[:4]
        for species in sorted(set(species_values)):
            for coordinate in coordinates:
                lat = float(coordinate["lat"])
                lon = float(coordinate["lon"])
                expanded.append(
                    {
                        "publication_id": publication_id,
                        "doi": str(row.doi or "").strip().casefold(),
                        "pmcid": str(row.pmcid or "").strip(),
                        "publication_date": date,
                        "title": str(row.title or ""),
                        "accepted_species": species,
                        "Latitude": lat,
                        "Longitude": lon,
                        "Author": "",
                        "Year": year,
                    }
                )

    base_columns = [
        "publication_id",
        "doi",
        "pmcid",
        "publication_date",
        "title",
        "accepted_species",
        "Latitude",
        "Longitude",
        "Author",
        "Year",
    ]
    raw = pd.DataFrame(expanded, columns=base_columns)
    output_columns = [
        "experiment_key",
        "publication_id",
        "doi",
        "publication_date",
        "study_key",
        "site_key",
        "accepted_species",
        "analysis_regime",
        "z_distance",
        "title",
        "source_database",
        "source_record_id",
        "latitude",
        "longitude",
        "distance_to_major_continent_km",
        "log1p_distance_to_major_continent_km",
    ]
    if raw.empty:
        return pd.DataFrame(columns=output_columns)

    raw = raw.drop_duplicates(
        ["publication_id", "accepted_species", "Latitude", "Longitude"]
    ).reset_index(drop=True)
    geo = build_preflight_table(
        raw[["Latitude", "Longitude", "doi", "Author", "Year"]]
        .rename(columns={"doi": "DOI"}),
        land,
        _parent_config(),
    )
    raw["study_key"] = geo["study_key"].astype(str)
    raw["site_key"] = geo["site_key"].astype(str)
    raw["analysis_regime"] = geo["analysis_regime"].astype(str)
    raw["distance_to_major_continent_km"] = pd.to_numeric(
        geo["distance_to_major_continent_km"], errors="coerce"
    )
    raw["log1p_distance_to_major_continent_km"] = pd.to_numeric(
        geo["log1p_distance_to_major_continent_km"], errors="coerce"
    )
    raw["z_distance"] = (
        raw["log1p_distance_to_major_continent_km"] - FROZEN_DISTANCE_MEAN
    ) / FROZEN_DISTANCE_SD
    raw["experiment_key"] = [
        _key(pub, species, lat, lon)
        for pub, species, lat, lon in zip(
            raw["publication_id"],
            raw["accepted_species"],
            raw["Latitude"],
            raw["Longitude"],
            strict=True,
        )
    ]
    raw["source_database"] = "Europe_PMC_methods_only"
    raw["source_record_id"] = raw["pmcid"]
    raw["latitude"] = raw["Latitude"]
    raw["longitude"] = raw["Longitude"]

    valid = (
        np.isfinite(raw["z_distance"].to_numpy(float))
        & raw["study_key"].ne("")
        & raw["site_key"].ne("")
        & raw["analysis_regime"].ne("unresolved")
    )
    return (
        raw.loc[valid, output_columns]
        .sort_values(["publication_id", "accepted_species", "site_key"])
        .reset_index(drop=True)
    )


@app.command("run")
def run(
    methods_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    trait_states_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    land_geojson: Path = typer.Option(..., exists=True, dir_okay=False),
    output_csv: Path = typer.Option(...),
) -> None:
    methods = pd.read_csv(methods_csv).fillna("")
    traits = pd.read_csv(trait_states_csv)
    land = gpd.read_file(land_geojson)
    result = build_support_metadata(methods, traits, land)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)
    typer.echo(
        json.dumps(
            {
                "n_experiment_keys": int(result["experiment_key"].nunique()),
                "n_publications": int(result["study_key"].nunique()),
                "n_species": int(result["accepted_species"].nunique()),
                "n_sites": int(result["site_key"].nunique()),
                "by_context": result["analysis_regime"].value_counts().to_dict(),
                "distance_standardization": {
                    "mean": FROZEN_DISTANCE_MEAN,
                    "sd": FROZEN_DISTANCE_SD,
                    "source": "frozen_parent_GloPL_global_distance",
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    app()
