"""Acquire quality-filtered Bombus occurrence environments and island targets.

This module builds a GLOBAL ENVIRONMENTAL SCREEN only. It does not define island
source pools, realized Bombus presence, or pollination service.
"""
from __future__ import annotations

import io
import json
import math
import zipfile
from pathlib import Path

import geopandas as gpd
import httpx
import pandas as pd
import rasterio
import typer

from island_v2.island_environment_summary import sample_island_environment

app = typer.Typer(add_completion=False, no_args_is_help=True)
GBIF = "https://api.gbif.org/v1"
WORLDCLIM_RESOLUTION = "5m"
WORLDCLIM_BIO = f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_{WORLDCLIM_RESOLUTION}_bio.zip"
WORLDCLIM_ELEV = f"https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_{WORLDCLIM_RESOLUTION}_elev.zip"
ENVIRONMENT_COLUMNS = ["bio1", "bio4", "bio5", "bio6", "bio12", "bio15", "bio18", "bio19", "elevation"]
BIO_NUMBERS = {"bio1": 1, "bio4": 4, "bio5": 5, "bio6": 6, "bio12": 12, "bio15": 15, "bio18": 18, "bio19": 19}
EXCLUDED_BASIS = {"FOSSIL_SPECIMEN", "MACHINE_OBSERVATION"}


def _download_zip_member(url: str, destination: Path, suffix: str) -> None:
    if destination.exists():
        return
    response = httpx.get(url, timeout=240.0, follow_redirects=True)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        member = next(name for name in archive.namelist() if name.endswith(suffix))
        destination.write_bytes(archive.read(member))


def _download_worldclim(cache_dir: Path) -> dict[str, Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    rasters: dict[str, Path] = {}
    for key, number in BIO_NUMBERS.items():
        path = cache_dir / f"wc2.1_{WORLDCLIM_RESOLUTION}_bio_{number}.tif"
        _download_zip_member(WORLDCLIM_BIO, path, path.name)
        rasters[key] = path
    elevation = cache_dir / f"wc2.1_{WORLDCLIM_RESOLUTION}_elev.tif"
    _download_zip_member(WORLDCLIM_ELEV, elevation, elevation.name)
    rasters["elevation"] = elevation
    return rasters


def _gbif_genus_key() -> int:
    response = httpx.get(f"{GBIF}/species/match", params={"name": "Bombus", "rank": "GENUS"}, timeout=60.0)
    response.raise_for_status()
    payload = response.json()
    if payload.get("matchType") == "NONE" or "usageKey" not in payload:
        raise typer.BadParameter("GBIF could not match genus Bombus")
    return int(payload["usageKey"])


def discover_bombus_species(min_occurrence_records: int, max_species: int) -> pd.DataFrame:
    response = httpx.get(
        f"{GBIF}/occurrence/search",
        params={"taxon_key": _gbif_genus_key(), "has_coordinate": "true", "occurrence_status": "PRESENT", "limit": 0, "facet": "speciesKey", "facet_limit": max(max_species * 3, 500)},
        timeout=120.0,
    )
    response.raise_for_status()
    facets = response.json().get("facets", [])
    counts = facets[0].get("counts", []) if facets else []
    rows = []
    for item in counts:
        count = int(item.get("count", 0))
        if count < min_occurrence_records:
            continue
        species_key = int(item["name"])
        r = httpx.get(f"{GBIF}/species/{species_key}", timeout=60.0)
        r.raise_for_status()
        p = r.json()
        name = str(p.get("canonicalName") or p.get("scientificName") or "").strip()
        if name.startswith("Bombus "):
            rows.append({"species_key": species_key, "bombus_species": name, "gbif_coordinate_record_count": count})
    return pd.DataFrame(rows).sort_values(["gbif_coordinate_record_count", "bombus_species"], ascending=[False, True]).head(max_species).reset_index(drop=True)


def _gbif_occurrences(species_key: int, name: str, limit: int) -> pd.DataFrame:
    rows = []
    offset = 0
    while len(rows) < limit:
        page_size = min(300, limit - len(rows))
        response = httpx.get(
            f"{GBIF}/occurrence/search",
            params={"taxon_key": species_key, "has_coordinate": "true", "occurrence_status": "PRESENT", "limit": page_size, "offset": offset},
            timeout=90.0,
        )
        response.raise_for_status()
        payload = response.json()
        batch = payload.get("results", [])
        if not batch:
            break
        for item in batch:
            lat, lon = item.get("decimalLatitude"), item.get("decimalLongitude")
            if lat is None or lon is None:
                continue
            rows.append({
                "bombus_species": name,
                "species_key": species_key,
                "gbif_id": item.get("key", ""),
                "decimal_latitude": lat,
                "decimal_longitude": lon,
                "year": item.get("year", ""),
                "dataset_key": item.get("datasetKey", ""),
                "basis_of_record": item.get("basisOfRecord", ""),
                "coordinate_uncertainty_m": item.get("coordinateUncertaintyInMeters", ""),
                "country_code": item.get("countryCode", ""),
                "issues": "|".join(item.get("issues", []) or []),
            })
        if payload.get("endOfRecords", False):
            break
        offset += page_size
    return pd.DataFrame(rows)


def quality_filter_and_thin_occurrences(frame: pd.DataFrame, *, grid_degrees: float = 0.25, max_coordinate_uncertainty_m: float = 20_000, min_year: int = 1950) -> pd.DataFrame:
    """Apply transparent record-quality rules and one deterministic record per species x grid cell."""
    d = frame.copy()
    d["decimal_latitude"] = pd.to_numeric(d["decimal_latitude"], errors="coerce")
    d["decimal_longitude"] = pd.to_numeric(d["decimal_longitude"], errors="coerce")
    d["year_numeric"] = pd.to_numeric(d.get("year"), errors="coerce")
    d["uncertainty_numeric"] = pd.to_numeric(d.get("coordinate_uncertainty_m"), errors="coerce")
    d = d[d["decimal_latitude"].between(-90, 90) & d["decimal_longitude"].between(-180, 180)]
    d = d[~d["basis_of_record"].fillna("").isin(EXCLUDED_BASIS)]
    d = d[d["year_numeric"].isna() | (d["year_numeric"] >= min_year)]
    d = d[d["uncertainty_numeric"].isna() | (d["uncertainty_numeric"] <= max_coordinate_uncertainty_m)]
    d = d.drop_duplicates(["bombus_species", "decimal_latitude", "decimal_longitude"])
    d["grid_x"] = (d["decimal_longitude"] / grid_degrees).apply(math.floor)
    d["grid_y"] = (d["decimal_latitude"] / grid_degrees).apply(math.floor)
    d = d.sort_values(["bombus_species", "grid_x", "grid_y", "uncertainty_numeric", "year_numeric"], na_position="last")
    return d.drop_duplicates(["bombus_species", "grid_x", "grid_y"], keep="first").drop(columns=["grid_x", "grid_y", "year_numeric", "uncertainty_numeric"])


def _sample(paths: dict[str, Path], frame: pd.DataFrame, lon: str, lat: str) -> pd.DataFrame:
    result = frame.copy()
    coordinates = list(zip(result[lon].astype(float), result[lat].astype(float), strict=True))
    for key, path in paths.items():
        with rasterio.open(path) as dataset:
            values = [float(v[0]) for v in dataset.sample(coordinates)]
            nodata = dataset.nodata
        s = pd.Series(values, dtype="float64")
        if nodata is not None:
            s = s.mask(s.eq(float(nodata)))
        result[key] = s
    return result


def build_global_targets(islands: gpd.GeoDataFrame, species: pd.DataFrame, island_id_column: str) -> pd.DataFrame:
    if island_id_column not in islands.columns:
        raise typer.BadParameter(f"island polygons missing {island_id_column}")
    ids = islands[island_id_column].astype(str).drop_duplicates().sort_values().reset_index(drop=True)
    target = ids.to_frame(name="island_id").merge(species[["bombus_species"]], how="cross")
    target["projection_scope"] = "global_environmental_screen_not_source_pool"
    return target


@app.command("acquire-global")
def acquire_global(
    islands_gpkg: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    island_id_column: str = typer.Option("island_id"),
    min_occurrence_records: int = typer.Option(50, min=20),
    max_species: int = typer.Option(300, min=1),
    max_occurrences_per_species: int = typer.Option(10000, min=20),
    max_environment_samples_per_island: int = typer.Option(128, min=1, max=512),
    environment_grid_size: int = typer.Option(16, min=2, max=80),
    thinning_grid_degrees: float = typer.Option(0.25, min=0.05, max=2.0),
    max_coordinate_uncertainty_m: float = typer.Option(20_000, min=0),
    min_year: int = typer.Option(1950, min=1800, max=2100),
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rasters = _download_worldclim(output_dir / "worldclim_cache")
    species = discover_bombus_species(min_occurrence_records, max_species)
    if species.empty:
        raise typer.BadParameter("no Bombus species passed the global occurrence threshold")
    species.to_csv(output_dir / "bombus_global_species_registry.csv", index=False)

    raw = pd.concat([_gbif_occurrences(int(r.species_key), str(r.bombus_species), max_occurrences_per_species) for r in species.itertuples(index=False)], ignore_index=True)
    raw.to_csv(output_dir / "bombus_occurrences_raw.csv.gz", index=False, compression="gzip")
    occurrences = quality_filter_and_thin_occurrences(raw, grid_degrees=thinning_grid_degrees, max_coordinate_uncertainty_m=max_coordinate_uncertainty_m, min_year=min_year)
    occurrences = _sample(rasters, occurrences, "decimal_longitude", "decimal_latitude").dropna(subset=ENVIRONMENT_COLUMNS)
    counts = occurrences.groupby("bombus_species").size()
    eligible = counts[counts >= min_occurrence_records].index
    occurrences = occurrences[occurrences["bombus_species"].isin(eligible)].copy()
    occurrences.to_csv(output_dir / "bombus_occurrence_environment.csv", index=False)

    islands = gpd.read_file(islands_gpkg).to_crs(4326)
    island_env, island_samples = sample_island_environment(islands, rasters, island_id_column=island_id_column, max_points_per_island=max_environment_samples_per_island, grid_size=environment_grid_size)
    island_env.to_csv(output_dir / "global_island_environment.csv", index=False)
    island_samples.to_csv(output_dir / "global_island_environment_samples.csv.gz", index=False, compression="gzip")

    scored_species = species[species["bombus_species"].isin(eligible)]
    targets = build_global_targets(islands, scored_species, island_id_column)
    env_cols = ["island_id", "n_environment_samples", *ENVIRONMENT_COLUMNS, *[f"{c}_p10" for c in ENVIRONMENT_COLUMNS], *[f"{c}_p90" for c in ENVIRONMENT_COLUMNS]]
    scored_targets = targets.merge(island_env[env_cols], on="island_id", how="left", validate="many_to_one")
    scored_targets["environment_point_estimate"] = "polygon_sample_median_with_p10_p90_diagnostics"
    scored_targets.to_csv(output_dir / "island_global_environmental_screen.csv", index=False)

    manifest = {
        "scope": "global_environmental_screen_not_source_pool",
        "model_name": "regularized_ellipsoidal_environmental_niche",
        "sources": {"occurrences": "GBIF occurrence search API", "climate": f"WorldClim 2.1 {WORLDCLIM_RESOLUTION} bioclimatic and elevation layers"},
        "n_species_discovered": int(len(species)),
        "n_species_after_quality_filter_and_thinning": int(len(scored_species)),
        "n_occurrence_rows_raw": int(len(raw)),
        "n_occurrence_rows_model": int(len(occurrences)),
        "n_islands": int(island_env["island_id"].nunique()),
        "n_island_species_rows": int(len(scored_targets)),
        "environment_columns": ENVIRONMENT_COLUMNS,
        "occurrence_quality": {"spatial_thinning_degrees": thinning_grid_degrees, "max_coordinate_uncertainty_m": max_coordinate_uncertainty_m, "min_year": min_year, "excluded_basis_of_record": sorted(EXCLUDED_BASIS)},
        "island_environment_summary": "deterministic polygon samples; median projection with p10/p90 retained; not local microclimate",
        "max_environment_samples_per_island": max_environment_samples_per_island,
        "environment_grid_size": environment_grid_size,
        "min_occurrence_records_after_filtering": min_occurrence_records,
    }
    (output_dir / "bombus_niche_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest))


if __name__ == "__main__":
    app()
