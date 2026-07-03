"""Download public GSHHG shorelines and derive documented island polygons."""

from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import UTC, datetime
from pathlib import Path

import geopandas as gpd
import httpx
import typer
import yaml
from shapely.ops import unary_union

from island_v2.gbif_flora import stable_island_table

app = typer.Typer(add_completion=False, no_args_is_help=True)


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    with httpx.stream("GET", url, timeout=600, follow_redirects=True) as response:
        response.raise_for_status()
        with temporary.open("wb") as handle:
            for block in response.iter_bytes():
                handle.write(block)
    temporary.replace(destination)


def extract_l1(archive: Path, directory: Path, resolution: str) -> Path:
    stem = f"GSHHS_{resolution}_L1"
    directory.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as zipped:
        members = [name for name in zipped.namelist() if Path(name).name.startswith(stem + ".")]
        if not members:
            raise RuntimeError(f"Missing {stem} in archive")
        for member in members:
            target = directory / Path(member).name
            if not target.exists():
                target.write_bytes(zipped.read(member))
    return directory / f"{stem}.shp"


def dissolve_ids(frame: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    id_name = next((name for name in ("id", "ID") if name in frame.columns), None)
    if id_name is None:
        return frame[["geometry"]].copy()
    rows = []
    for identifier, group in frame.groupby(id_name, dropna=False, sort=True):
        rows.append({"gshhg_polygon_id": str(identifier), "geometry": unary_union(group.geometry.tolist())})
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=frame.crs)


@app.command()
def build(
    config_path: Path = typer.Option(Path("config/island_source_gshhg.yml")),
    output_dir: Path = typer.Option(Path("data/v2/external/islands/gshhg")),
) -> None:
    """Acquire GSHHG and write a v2-compatible exact-polygon GeoPackage."""
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    resolution = str(config["resolution"])
    if resolution not in {"f", "h", "i", "l", "c"}:
        raise typer.BadParameter("GSHHG resolution must be f, h, i, l, or c")
    archive = output_dir / "source" / f"gshhg-shp-{config['source_version']}.zip"
    download(str(config["source_url"]), archive)
    shapefile = extract_l1(archive, output_dir / "source" / f"L1_{resolution}", resolution)
    coast = gpd.read_file(shapefile)
    if coast.crs is None:
        coast = coast.set_crs(4326)
    dissolved = dissolve_ids(coast)
    islands = stable_island_table(
        dissolved,
        source_label=f"gshhg_{config['source_version']}_{resolution}",
        min_area_km2=float(config["min_area_km2"]),
        name_field=None,
        parent_id_field="gshhg_polygon_id" if "gshhg_polygon_id" in dissolved.columns else None,
    )
    cutoff = float(config["mainland_area_threshold_km2"])
    islands = islands.loc[islands["area_km2"] <= cutoff].copy()
    islands["landmass_rule"] = f"L1 area <= {cutoff} km2"
    prepared = output_dir / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    islands.to_file(prepared / "islands_v2.gpkg", layer="islands", driver="GPKG")
    manifest = islands.drop(columns="geometry").copy()
    manifest["source_url"] = config["source_url"]
    manifest["source_version"] = config["source_version"]
    manifest["source_resolution"] = resolution
    manifest["source_sha256"] = checksum(archive)
    manifest["prepared_at_utc"] = datetime.now(UTC).isoformat()
    manifest.to_csv(prepared / "island_manifest.csv", index=False)
    (prepared / "source_policy.json").write_text(
        json.dumps({**config, "archive_sha256": checksum(archive), "n_islands": len(islands)}, indent=2),
        encoding="utf-8",
    )
    typer.echo(f"Prepared {len(islands)} exact GSHHG island polygons")


if __name__ == "__main__":
    app()
