"""Summarize environmental rasters across island polygons.

Representative-point values are retained only as a fallback for very small polygons.
The production summary uses deterministic interior sampling and records median and
10th/90th percentiles for every environmental axis.
"""

from __future__ import annotations

import math
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from shapely.geometry import Point


def deterministic_polygon_points(
    geometry: object,
    *,
    max_points: int = 64,
    grid_size: int = 12,
) -> list[Point]:
    """Return deterministic interior sample points, with representative-point fallback."""
    if geometry is None or getattr(geometry, "is_empty", True):
        return []
    minx, miny, maxx, maxy = geometry.bounds
    if not all(math.isfinite(value) for value in (minx, miny, maxx, maxy)):
        return []
    if minx == maxx or miny == maxy:
        return [geometry.representative_point()]

    xs = np.linspace(minx, maxx, grid_size + 2)[1:-1]
    ys = np.linspace(miny, maxy, grid_size + 2)[1:-1]
    candidates = [Point(float(x), float(y)) for y in ys for x in xs]
    inside = [point for point in candidates if geometry.covers(point)]
    if not inside:
        return [geometry.representative_point()]
    if len(inside) <= max_points:
        return inside
    indices = np.linspace(0, len(inside) - 1, max_points).round().astype(int)
    return [inside[int(index)] for index in indices]


def sample_island_environment(
    islands: gpd.GeoDataFrame,
    rasters: dict[str, Path],
    *,
    island_id_column: str = "island_id",
    max_points_per_island: int = 64,
    grid_size: int = 12,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sample raster values within each island and return summaries plus raw samples."""
    if island_id_column not in islands.columns:
        raise ValueError(f"island polygons missing {island_id_column}")
    if islands.crs is None:
        raise ValueError("island polygons must have a declared CRS")
    frame = islands.to_crs(4326)

    rows: list[dict[str, object]] = []
    for record in frame[[island_id_column, "geometry"]].itertuples(index=False):
        island_id = str(getattr(record, island_id_column))
        points = deterministic_polygon_points(
            record.geometry,
            max_points=max_points_per_island,
            grid_size=grid_size,
        )
        for sample_index, point in enumerate(points):
            rows.append(
                {
                    "island_id": island_id,
                    "sample_index": sample_index,
                    "decimal_longitude": point.x,
                    "decimal_latitude": point.y,
                }
            )
    samples = pd.DataFrame(rows)
    if samples.empty:
        return pd.DataFrame(columns=["island_id"]), samples

    coordinates = list(
        zip(
            samples["decimal_longitude"].astype(float),
            samples["decimal_latitude"].astype(float),
            strict=True,
        )
    )
    for key, path in rasters.items():
        with rasterio.open(path) as dataset:
            values = [float(value[0]) for value in dataset.sample(coordinates)]
            nodata = dataset.nodata
        series = pd.Series(values, dtype="float64")
        if nodata is not None:
            series = series.mask(series.eq(float(nodata)))
        samples[key] = series

    summaries: list[dict[str, object]] = []
    for island_id, group in samples.groupby("island_id", sort=True):
        row: dict[str, object] = {
            "island_id": island_id,
            "n_environment_samples": int(len(group)),
        }
        for key in rasters:
            values = pd.to_numeric(group[key], errors="coerce").dropna()
            row[f"{key}_median"] = float(values.median()) if not values.empty else np.nan
            row[f"{key}_p10"] = float(values.quantile(0.10)) if not values.empty else np.nan
            row[f"{key}_p90"] = float(values.quantile(0.90)) if not values.empty else np.nan
            row[key] = row[f"{key}_median"]
        summaries.append(row)
    return pd.DataFrame(summaries), samples
