"""Acquire public island polygons with a recorded GSHHG-to-Natural-Earth fallback.

GSHHG is the preferred source because it explicitly encodes hierarchical
shoreline polygons.  Some runners cannot reach its legacy distribution hosts;
when that happens this module falls back to a pinned Natural Earth 10m land
GeoJSON stored on GitHub.  The fallback is never silent: backend, URL, checksum,
and the GSHHG failure message are written to source_policy.json.
"""

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

# Pinned upstream commit, rather than a moving default branch.
NATURAL_EARTH_COMMIT = "ca96624a56bd078437bca8184e78163e5039ad19"
NATURAL_EARTH_LAND_URL = (
    "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
    f"{NATURAL_EARTH_COMMIT}/geojson/ne_10m_land.geojson"
)


@app.callback()
def main() -> None:
    """Acquire and prepare public island polygons for island-floral v2."""


def checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def geometry_hash(geometry: Any) -> str:
    return hashlib.sha256(geometry.wkb).hexdigest()[:20]


def download(url: str, destination: Path) -> None:
    """Download a public source with bounded connect waits and clean failures."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 0:
        return
    temporary = destination.with_suffix(destination.suffix + ".part")
    timeout = httpx.Timeout(connect=20.0, read=180.0, write=60.0, pool=20.0)
    try:
        with httpx.stream(
            "GET",
            url,
            timeout=timeout,
            follow_redirects=True,
            headers={"User-Agent": "island-floral-v2/0.1"},
        ) as response:
            response.raise_for_status()
            with temporary.open("wb") as handle:
                for block in response.iter_bytes():
                    handle.write(block)
        temporary.replace(destination)
    except (httpx.HTTPError, OSError) as exc:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"download failed for {url}: {exc}") from exc


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
        frame["source_feature_id"] = [str(index) for index in frame.index]
        return frame[["source_feature_id", "geometry"]]
    rows = []
    for identifier, group in frame.groupby(id_name, dropna=False, sort=True):
        geometry = make_valid(unary_union(group.geometry.tolist()))
        if not geometry.is_empty:
            rows.append({"source_feature_id": str(identifier), "geometry": geometry})
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
    features: gpd.GeoDataFrame,
    source_label: str,
    min_area_km2: float,
    mainland_area_threshold_km2: float,
) -> gpd.GeoDataFrame:
    """Create one analysis unit per input landmass geometry."""
    if features.crs is None:
        raise typer.BadParameter("Public island source has no CRS")
    units = features.to_crs(4326).copy()
    units["geometry"] = units.geometry.map(polygonal_geometry)
    units = units.loc[units.geometry.notna() & ~units.geometry.is_empty].copy()
    units["area_km2"] = units.to_crs(6933).area / 1_000_000
    units = units.loc[
        (units["area_km2"] >= min_area_km2)
        & (units["area_km2"] <= mainland_area_threshold_km2)
    ].copy()
    units["island_id"] = [f"{source_label}_{geometry_hash(geometry)}" for geometry in units.geometry]
    if units["island_id"].duplicated().any():
        raise RuntimeError("Geometry-derived island IDs are not unique")
    units["source_label"] = source_label
    units["parent_feature_id"] = units["source_feature_id"].astype(str)
    units["part_index"] = 1
    units["island_name"] = ""
    units["geometry_sha256"] = [geometry_hash(geometry) for geometry in units.geometry]
    units["landmass_rule"] = (
        f"{min_area_km2} <= area_km2 <= {mainland_area_threshold_km2}"
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


def natural_earth_parts(frame: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Split Natural Earth land MultiPolygons into individual candidate islands."""
    if frame.crs is None:
        frame = frame.set_crs(4326)
    polygons: list[dict[str, Any]] = []
    for feature_index, geometry in enumerate(frame.geometry):
        repaired = polygonal_geometry(geometry)
        if repaired is None:
            continue
        parts = list(repaired.geoms) if isinstance(repaired, MultiPolygon) else [repaired]
        for part_index, part in enumerate(parts, start=1):
            polygons.append(
                {
                    "source_feature_id": f"ne_{feature_index}_{part_index}",
                    "geometry": part,
                }
            )
    return gpd.GeoDataFrame(polygons, geometry="geometry", crs=frame.crs)


def build_from_gshhg(
    source_url: str,
    config: dict[str, Any],
    output_dir: Path,
) -> tuple[gpd.GeoDataFrame, Path]:
    resolution = str(config["resolution"])
    if resolution not in {"f", "h", "i", "l", "c"}:
        raise typer.BadParameter("GSHHG resolution must be f, h, i, l, or c")
    archive = output_dir / "source" / f"gshhg-shp-{config['source_version']}.zip"
    download(source_url, archive)
    shapefile = extract_l1(archive, output_dir / "source" / f"L1_{resolution}", resolution)
    coast = gpd.read_file(shapefile)
    if coast.crs is None:
        coast = coast.set_crs(4326)
    units = make_island_units(
        features=dissolve_ids(coast),
        source_label=f"gshhg_{config['source_version']}_{resolution}",
        min_area_km2=float(config["min_area_km2"]),
        mainland_area_threshold_km2=float(config["mainland_area_threshold_km2"]),
    )
    return units, archive


def build_from_natural_earth(
    config: dict[str, Any],
    output_dir: Path,
) -> tuple[gpd.GeoDataFrame, Path]:
    geojson_path = output_dir / "source" / f"natural_earth_10m_land_{NATURAL_EARTH_COMMIT}.geojson"
    download(NATURAL_EARTH_LAND_URL, geojson_path)
    land = gpd.read_file(geojson_path)
    units = make_island_units(
        features=natural_earth_parts(land),
        source_label=f"natural_earth_10m_{NATURAL_EARTH_COMMIT[:12]}",
        min_area_km2=float(config["min_area_km2"]),
        mainland_area_threshold_km2=float(config["mainland_area_threshold_km2"]),
    )
    return units, geojson_path


@app.command()
def build(
    config_path: Path = typer.Option(Path("config/island_source_gshhg.yml")),
    output_dir: Path = typer.Option(Path("data/v2/external/islands/gshhg")),
    source_url: str | None = typer.Option(
        None,
        help="Optional GSHHG archive URL override; retained in source_policy.json.",
    ),
    allow_natural_earth_fallback: bool = typer.Option(
        True,
        help="Use pinned Natural Earth 10m land only when GSHHG acquisition fails.",
    ),
) -> None:
    """Acquire GSHHG, or a recorded public fallback, and write island polygons."""
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    gshhg_url = source_url or str(config["source_url"])
    gshhg_error: str | None = None
    backend = "gshhg"
    source_url_used = gshhg_url
    try:
        islands, source_file = build_from_gshhg(gshhg_url, config, output_dir)
    except Exception as exc:
        gshhg_error = str(exc)
        if not allow_natural_earth_fallback:
            raise typer.BadParameter(gshhg_error) from exc
        backend = "natural_earth_10m_fallback"
        source_url_used = NATURAL_EARTH_LAND_URL
        islands, source_file = build_from_natural_earth(config, output_dir)

    prepared = output_dir / "prepared"
    prepared.mkdir(parents=True, exist_ok=True)
    islands.to_file(prepared / "islands_v2.gpkg", layer="islands", driver="GPKG")
    manifest = islands.drop(columns="geometry").copy()
    manifest["source_backend"] = backend
    manifest["source_url"] = source_url_used
    manifest["source_file_sha256"] = checksum(source_file)
    manifest["prepared_at_utc"] = datetime.now(UTC).isoformat()
    manifest.to_csv(prepared / "island_manifest.csv", index=False)
    policy = {
        **config,
        "source_backend": backend,
        "source_url_used": source_url_used,
        "source_file_sha256": checksum(source_file),
        "gshhg_attempt_error": gshhg_error,
        "natural_earth_fallback_url": NATURAL_EARTH_LAND_URL,
        "n_islands": len(islands),
    }
    (prepared / "source_policy.json").write_text(
        json.dumps(policy, indent=2), encoding="utf-8"
    )
    typer.echo(f"Prepared {len(islands)} exact island polygons using {backend}")


if __name__ == "__main__":
    app()
