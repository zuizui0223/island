from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin
from shapely.geometry import Polygon

from island_v2.island_environment_summary import (
    deterministic_polygon_points,
    sample_island_environment,
)


def _write_raster(path: Path, values: np.ndarray) -> None:
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=values.shape[0],
        width=values.shape[1],
        count=1,
        dtype="float32",
        crs="EPSG:4326",
        transform=from_origin(0, 2, 1, 1),
        nodata=-9999,
    ) as dataset:
        dataset.write(values.astype("float32"), 1)


def test_polygon_points_are_deterministic_and_inside() -> None:
    polygon = Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    first = deterministic_polygon_points(polygon, max_points=8, grid_size=4)
    second = deterministic_polygon_points(polygon, max_points=8, grid_size=4)
    assert [(point.x, point.y) for point in first] == [(point.x, point.y) for point in second]
    assert len(first) == 8
    assert all(polygon.covers(point) for point in first)


def test_environment_summary_keeps_median_and_range(tmp_path: Path) -> None:
    raster = tmp_path / "bio1.tif"
    _write_raster(raster, np.array([[1, 3], [5, 7]], dtype=float))
    islands = gpd.GeoDataFrame(
        {
            "island_id": ["island_a"],
            "geometry": [Polygon([(0, 0), (2, 0), (2, 2), (0, 2)])],
        },
        crs="EPSG:4326",
    )

    summary, samples = sample_island_environment(
        islands,
        {"bio1": raster},
        max_points_per_island=16,
        grid_size=4,
    )

    assert len(samples) == 16
    assert summary.loc[0, "n_environment_samples"] == 16
    assert summary.loc[0, "bio1"] == summary.loc[0, "bio1_median"]
    assert summary.loc[0, "bio1_p10"] <= summary.loc[0, "bio1_median"]
    assert summary.loc[0, "bio1_median"] <= summary.loc[0, "bio1_p90"]
    assert pd.notna(summary.loc[0, "bio1_median"])
