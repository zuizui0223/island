"""Acquire real Bombus occurrence environments and island target environments."""

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


def _download_worldclim(cache_dir: Path) -> dict[str, Path]:
    cache_dir.mkdir(parents=True, exist_ok=True)
    wanted = {"bio1": cache_dir / "wc2.1_10m_bio_1.tif", "bio12": cache_dir / "wc2.1_10m_bio_12.tif"}
    if all(path.exists() for path in wanted.values()):
        return wanted
    response = httpx.get(WORLDCLIM_BIO, timeout=180.0, follow_redirects=True)
    response.raise_for_status()
    with zipfile.ZipFile(io.BytesIO(response.content)) as archive:
        for key, path in wanted.items():
            member = next(name for name in archive.namelist() if name.endswith(path.name))
            path.write_bytes(archive.read(member))
    return wanted


def _gbif_taxon_key(name: str) -> int:
    response = httpx.get(f"{GBIF}/species/match", params={"name": name}, timeout=60.0)
    response.raise_for_status()
    payload = response.json()
    if payload.get("matchType") == "NONE" or "usageKey" not in payload:
        raise typer.BadParameter(f"GBIF could not match species: {name}")
    return int(payload["usageKey"])


def _gbif_occurrences(name: str, limit: int) -> pd.DataFrame:
    key = _gbif_taxon_key(name)
    rows: list[dict[str, object]] = []
    offset = 0
    while len(rows) < limit:
        page_size = min(300, limit - len(rows))
        response = httpx.get(
            f"{GBIF}/occurrence/search",
            params={
                "taxon_key": key,
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
                    "gbif_id": item.get("key", ""),
                    "decimal_latitude": lat,
                    "decimal_longitude": lon,
                    "year": item.get("year", ""),
                    "dataset_key": item.get("datasetKey", ""),
                    "basis_of_record": item.get("basisOfRecord", ""),
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


@app.command("acquire")
def acquire(
    islands_gpkg: Path = typer.Option(..., exists=True),
    island_species_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    island_id_column: str = typer.Option("island_id"),
    max_occurrences_per_species: int = typer.Option(5000, min=20),
) -> None:
    """Create the two production hypervolume inputs from public real data."""
    targets = pd.read_csv(island_species_csv, dtype=str).fillna("")
    required = {"island_id", "bombus_species"}
    if missing := required.difference(targets.columns):
        raise typer.BadParameter(f"island-species table missing columns: {sorted(missing)}")
    species = sorted(value for value in targets["bombus_species"].unique() if value)
    if not species:
        raise typer.BadParameter("no Bombus candidate species supplied")

    output_dir.mkdir(parents=True, exist_ok=True)
    rasters = _download_worldclim(output_dir / "worldclim_cache")

    occurrence_parts = [_gbif_occurrences(name, max_occurrences_per_species) for name in species]
    occurrences = pd.concat(occurrence_parts, ignore_index=True)
    occurrences = _sample(rasters, occurrences, "decimal_longitude", "decimal_latitude")
    occurrences = occurrences.dropna(subset=["bio1", "bio12"])
    occurrences.to_csv(output_dir / "bombus_occurrence_environment.csv", index=False)

    islands = gpd.read_file(islands_gpkg)
    if island_id_column not in islands.columns:
        raise typer.BadParameter(f"island polygons missing {island_id_column}")
    islands = islands[[island_id_column, "geometry"]].rename(columns={island_id_column: "island_id"})
    islands = islands.to_crs(4326)
    points = islands.geometry.representative_point()
    island_env = pd.DataFrame(
        {
            "island_id": islands["island_id"].astype(str),
            "decimal_longitude": points.x,
            "decimal_latitude": points.y,
        }
    )
    island_env = _sample(rasters, island_env, "decimal_longitude", "decimal_latitude")
    scored_targets = targets.merge(island_env[["island_id", "bio1", "bio12"]], on="island_id", how="left", validate="many_to_one")
    scored_targets.to_csv(output_dir / "island_source_pool_environment.csv", index=False)

    manifest = {
        "source": {"occurrences": "GBIF occurrence search API", "climate": "WorldClim 2.1 10-minute bioclim"},
        "species": species,
        "n_occurrence_rows": int(len(occurrences)),
        "n_island_species_rows": int(len(scored_targets)),
        "environment_columns": ["bio1", "bio12"],
        "max_occurrences_per_species": max_occurrences_per_species,
    }
    (output_dir / "bombus_niche_input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    typer.echo(json.dumps(manifest))


if __name__ == "__main__":
    app()
