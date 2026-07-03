"""Download public GSHHG shorelines and derive documented island polygons."""

from __future__ import annotations

import hashlib
import json
import zipfile
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import geopandas as gpd
import httpx
import typer
import yaml
from shapely import make_valid
from shapely.geometry import MultiPolygon, Polygon
from shapely.ops import unary_union

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.callback()
def main() -> None:
    """Acquire and prepare public GSHHG island polygons for island-floral v2."""


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def geometry_hash(geometry: Any) -> str:
    return hashlib.sha256(geometry.wkb).hexdigest()[:20]


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    try:
        with httpx.stream("GET", url, timeout=600, follow_redirects=True) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for block in response.iter_bytes():
                    handle.write(block)
        temporary.replace(destination)
    except httpx.HTTPError as exc:
        temporary.unlink(missing_ok=True)
        raise typer.BadParameter(f"GSHHG archive download failed for {url}: {exc}") from exc


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
    """Recombine components of each GSHHG landmass before any area filtering."""
    id_name = next((name for name in ("id", "ID") if name in frame.columns), None)
    if id_name is None:
        frame = frame.copy()
        frame["gshhg_polygon_id"] = [str(index) for index in frame.index]
        return frame[["gshhg_polygon_id", "geometry"]]
    rows = []
    for identifier, group in frame.groupby(id_name, dropna=False, sort=True):
        geometry = make_valid(unary_union(group.geometry.tolist()))
        if geometry.is_empty:
            continue
        rows.append({"gshhg_polygon_id": str(identifier), "geometry": geometry})
    return gpd.GeoDataFrame(rows, geometry="geometry", crs=frame.crs)


def polygonal_geometry(geometry: Any) -> Polygon | MultiPolygon | None:
    repaired = make_valid(geometry)
    if isinstance(repaired, (Polygon, MultiPolygon)):
        return repaired
    parts = [part for part in getattr(repaired, "geoms", []) if isinstance(part, Polygon)]
    if not parts:
        return None
    return MultiPolygon(parts) if len(parts) > 1 else parts[0]


def make_island_units(
    dissolved: gpd.GeoDataFrame,
    source_label: str,
    min_area_km2: float,
    mainland_area_threshold_km2: float,
) -> gpd.GeoDataFrame:
    """Keep one GSHHG L1 landmass as one analysis unit, including MultiPolygons."""
    if dissolved.crs is None:
        raise typer.BadParameter("GSHHG input has no CRS")
    units = dissolved.to_crs(4326).copy()
    units["geometry"] = units.geometry.map(polygonal_geometry)
    units = units.loc[units.geometry.notna() & ~units.geometry.is_empty].copy()
    units["area_km2"] = units.to_crs(6933).area / 1_000_000
    units = units.loc[
        (units["area_km2"] >= min_area_km2)
        & (units["area_km2"] <= mainland_area_threshold_km2)
    ].copy()
    units["island_id"] = [f"{source_label}_{geometry_hash(geometry)}" for geometry in units.geometry]
    if units["island_id"].duplicated().any():
        raise RuntimeError("GSHHG geometry-derived island IDs are not unique")
    units["source_label"] = source_label
    units["parent_feature_id"] = units["gshhg_polygon_id"].astype(str)
    units["part_index"] = 1
    units["island_name"] = ""
    units["geometry_sha256"] = [geometry_hash(geometry) for geometry in units.geometry]
    units["landmass_rule"] = (
        f"L1 landmass {min_area_km2} <= area_km2 <= {mainland_area_threshold_km2}"
    )
    columns = [
        "island_id",
        "source_label",
        "parent_feature_id",
        "part_index",
        "island_name",
        "area_km2",
        "geometry_sha256",
        "landmass_rule",
        "geometry",
    ]
    return units[columns].sort_values("island_id").reset_index(drop=True)


@app.command()
def build(
    config_path: Path = typer.Option(Path("config/island_source_gshhg.yml")),
    output_dir: Path = typer.Option(Path("data/v2/external/islands/gshhg")),
    source_url: str | None = typer.Option(
        None,
        help="Optional archive URL override; records the effective URL in the output manifest.",
    ),
) -> None:
    """Acquire GSHHG and write a v2-compatible exact-polygon GeoPackage."""
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    resolution = str(config["resolution"])
    if resolution not in {"f", "h", "i", "l", "c"}:
        raise typer.BadParameter("GSHHG resolution must be f, h, i, l, or c")
    effective_source_url = source_url or str(config["source_url"])
    archive = output_dir / "source" / f"gshhg-shp-{config['source_version']}.zip"
    download(effective_source_url, archive)
    shapefile = extract_l1(archive, output_dir / "source" / f"L1_{resolution}", resolution)
    coast = gpd.read_file(shapefile)
    if coast.crs is None:
        coast = coast.set_crs(4326)
    islands = make_island_units(
        dissolved=dissolve_ids(coast),
        source_label=f"gshhg_{config['source_version']}_{resolution}",
        min_area_km2=float(config["min_area_km2"]),
        mainland_area_threshold_km2=float(config["mainland_area_threshold_km2"]),
    )
    prepared = output_dir / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    islands.to_file(prepared / "islands_v2.gpkg", layer="islands", driver="GPKG")
    manifest = islands.drop(columns="geometry").copy()
    manifest["source_url"] = effective_source_url
    manifest["source_version"] = config["source_version"]
    manifest["source_resolution"] = resolution
    manifest["source_sha256"] = checksum(archive)
    manifest["prepared_at_utc"] = datetime.now(UTC).isoformat()
    manifest.to_csv(prepared / "island_manifest.csv", index=False)
    (prepared / "source_policy.json").write_text(
        json.dumps(
            {**config, "source_url_used": effective_source_url, "archive_sha256": checksum(archive), "n_islands": len(islands)},
            indent=2,
        ),
        encoding="utf-8",
    )
    typer.echo(f"Prepared {len(islands)} exact GSHHG island polygons")


if __name__ == "__main__":
    app()
