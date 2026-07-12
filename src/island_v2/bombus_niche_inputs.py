"""Acquire global Bombus occurrence environments and global island targets."""

from __future__ import annotations

import io
import json
import zipfile
from pathlib import Path

import geopandas as gpd
import httpx
import pandas as pd
import rasterio
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)
GBIF = "https://api.gbif.org/v1"
WORLDCLIM_BIO = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_bio.zip"
WORLDCLIM_ELEV = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_10m_elev.zip"
ENVIRONMENT_COLUMNS = ["bio1", "bio4", "bio5", "bio6", "bio12", "bio15", "bio18", "bio19", "elevation"]
BIO_NUMBERS = {"bio1": 1, "bio4": 4, "bio5": 5, "bio6": 6, "bio12": 12, "bio15": 15, "bio18": 18, "bio19": 19}


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
        path = cache_dir / f"wc2.1_10m_bio_{number}.tif"
        _download_zip_member(WORLDCLIM_BIO, path, path.name)
        rasters[key] = path
    elevation = cache_dir / "wc2.1_10m_elev.tif"
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
    """Discover globally recorded Bombus species from a genus-level GBIF facet."""
    response = httpx.get(
        f"{GBIF}/occurrence/search",
        params={
            "taxon_key": _gbif_genus_key(),
            "has_coordinate": "true",
            "occurrence_status": "PRESENT",
            "limit": 0,
            "facet": "speciesKey",
            "facet_limit": max(max_species * 3, 500),
        },
        timeout=120.0,
    )
    response.raise_for_status()
    facets = response.json().get("facets", [])
    counts = facets[0].get("counts", []) if facets else []
    rows: list[dict[str, object]] = []
    for item in counts:
        count = int(item.get("count", 0))
        if count < min_occurrence_records:
            continue
        species_key = int(item["name"])
        species_response = httpx.get(f"{GBIF}/species/{species_key}", timeout=60.0)
        species_response.raise_for_status()
        species_payload = species_response.json()
        scientific_name = str(species_payload.get("canonicalName") or species_payload.get("scientificName") or "").strip()
        if not scientific_name.startswith("Bombus "):
            continue
        rows.append({"species_key": species_key, "bombus_species": scientific_name, "gbif_coordinate_record_count": count})
    return pd.DataFrame(rows).sort_values(["gbif_coordinate_record_count", "bombus_species"], ascending=[False, True]).head(max_species).reset_index(drop=True)


def _gbif_occurrences(species_key: int, name: str, limit: int) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    offset = 0
    while len(rows) < limit:
        page_size = min(300, limit - len(rows))
        response = httpx.get(
            f"{GBIF}/occurrence/search",
            params={
                "taxon_key": species_key,
                "has_coordinate": "true",
                "occurrence_status": "PRESENT",
                "limit": page_size,
                "offset": offset,
            },
            timeout=90.0,
        )
        response.raise_for_status()
        payload = response.json()
        batch = payload.get("results", [])
        if not batch:
            break
        for item in batch:
            lat = item.get("decimalLatitude")
            lon = item.get("decimalLongitude")
            if lat is None or lon is None:
                continue
            rows.append(
                {
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
                }
            )
        if payload.get("endOfRecords", False):
            break
        offset += page_size
    return pd.DataFrame(rows).drop_duplicates(["bombus_species", "decimal_latitude", "decimal_longitude"])


def _sample(paths: dict[str, Path], frame: pd.DataFrame, lon: str, lat: str) -> pd.DataFrame:
    result = frame.copy()
    coordinates = list(zip(result[lon].astype(float), result[lat].astype(float), strict=True))
    for key, path in paths.items():
        with rasterio.open(path) as dataset:
            values = [float(value[0]) for value in dataset.sample(coordinates)]
            nodata = dataset.nodata
        series = pd.Series(values, dtype="float64")
        if nodata is not None:
            series = series.mask(series.eq(float(nodata)))
        result[key] = series
    return result


def build_global_targets(islands: gpd.GeoDataFrame, species: pd.DataFrame, island_id_column: str) -> pd.DataFrame:
    if island_id_column not in islands.columns:
        raise typer.BadParameter(f"island polygons missing {island_id_column}")
    island_ids = islands[island_id_column].astype(str).drop_duplicates().sort_values().reset_index(drop=True)
    target = island_ids.to_frame(name="island_id").merge(species[["bombus_species"]], how="cross")
    target["projection_scope"] = "global_environmental_screen_not_source_pool"
    return target


@app.command("acquire-global")
def acquire_global(
    islands_gpkg: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    island_id_column: str = typer.Option("island_id"),
    min_occurrence_records: int = typer.Option(50, min=20),
    max_species: int = typer.Option(300, min=1),
    max_occurrences_per_species: int = typer.Option(5000, min=20),
) -> None:
    """Create global Bombus hypervolume inputs from public occurrence and climate data."""
    output_dir.mkdir(parents=True, exist_ok=True)
    rasters = _download_worldclim(output_dir / "worldclim_cache")
    species = discover_bombus_species(min_occurrence_records, max_species)
    if species.empty:
        raise typer.BadParameter("no Bombus species passed the global occurrence threshold")
    species.to_csv(output_dir / "bombus_global_species_registry.csv", index=False)

    occurrence_parts = [
        _gbif_occurrences(int(row.species_key), str(row.bombus_species), max_occurrences_per_species)
        for row in species.itertuples(index=False)
    ]
    occurrences = pd.concat(occurrence_parts, ignore_index=True)
    occurrences = _sample(rasters, occurrences, "decimal_longitude", "decimal_latitude")
    occurrences = occurrences.dropna(subset=ENVIRONMENT_COLUMNS)
    occurrences.to_csv(output_dir / "bombus_occurrence_environment.csv", index=False)

    islands = gpd.read_file(islands_gpkg).to_crs(4326)
    points = islands.geometry.representative_point()
    island_env = pd.DataFrame(
        {
            "island_id": islands[island_id_column].astype(str),
            "decimal_longitude": points.x,
            "decimal_latitude": points.y,
        }
    ).drop_duplicates("island_id")
    island_env = _sample(rasters, island_env, "decimal_longitude", "decimal_latitude")
    island_env.to_csv(output_dir / "global_island_environment.csv", index=False)

    targets = build_global_targets(islands, species, island_id_column)
    scored_targets = targets.merge(island_env[["island_id", *ENVIRONMENT_COLUMNS]], on="island_id", how="left", validate="many_to_one")
    scored_targets.to_csv(output_dir / "island_source_pool_environment.csv", index=False)

    manifest = {
        "scope": "global_environmental_screen_not_source_pool",
        "sources": {
            "occurrences": "GBIF occurrence search API",
            "climate": "WorldClim 2.1 10-minute bioclimatic and elevation layers",
        },
        "n_species": int(len(species)),
        "n_occurrence_rows": int(len(occurrences)),
        "n_islands": int(island_env["island_id"].nunique()),
        "n_island_species_rows": int(len(scored_targets)),
        "environment_columns": ENVIRONMENT_COLUMNS,
        "min_occurrence_records": min_occurrence_records,
        "max_species": max_species,
        "max_occurrences_per_species": max_occurrences_per_species,
    }
    (output_dir / "bombus_niche_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest))


if __name__ == "__main__":
    app()
