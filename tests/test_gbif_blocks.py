import geopandas as gpd
from shapely import from_wkt
from shapely.geometry import Polygon

from island_v2.gbif_blocks import block_query_geometries, build_query_catchments


def polygon_parts(geometry):
    return [geometry] if geometry.geom_type == "Polygon" else list(geometry.geoms)


def test_query_catchments_preserve_island_representative_points():
    islands = gpd.GeoDataFrame(
        {
            "island_id": ["a", "b"],
            "area_km2": [10.0, 12.0],
            "geometry_sha256": ["ha", "hb"],
        },
        geometry=[
            Polygon([(0, 0), (0, 0.1), (0.1, 0.1), (0.1, 0), (0, 0)]),
            Polygon([(1, 1), (1, 1.1), (1.1, 1.1), (1.1, 1), (1, 1)]),
        ],
        crs=4326,
    )

    catchments = build_query_catchments(islands, buffer_m=2_000, simplify_m=500)

    for original, catchment in zip(islands.geometry, catchments.geometry, strict=True):
        assert catchment.covers(original.representative_point())
        for polygon in polygon_parts(catchment):
            assert polygon.exterior.is_ccw
            assert all(not ring.is_ccw for ring in polygon.interiors)


def test_regional_blocks_cover_each_island_once_store_valid_ccw_wkt():
    islands = gpd.GeoDataFrame(
        {
            "island_id": ["a", "b", "c"],
            "area_km2": [10.0, 12.0, 14.0],
            "geometry_sha256": ["ha", "hb", "hc"],
        },
        geometry=[
            Polygon([(0, 0), (0, 0.1), (0.1, 0.1), (0.1, 0), (0, 0)]),
            Polygon([(0.5, 0), (0.5, 0.1), (0.6, 0.1), (0.6, 0), (0.5, 0)]),
            Polygon([(35, 0), (35, 0.1), (35.1, 0.1), (35.1, 0), (35, 0)]),
        ],
        crs=4326,
    )
    catchments = build_query_catchments(islands, buffer_m=2_000, simplify_m=500)

    blocks, members = block_query_geometries(
        catchments,
        grid_degrees=30,
        max_wkt_chars=20_000,
        max_islands_per_block=125,
    )

    assert set(members["island_id"]) == {"a", "b", "c"}
    assert not members["island_id"].duplicated().any()
    assert (blocks["query_wkt_chars"] <= 20_000).all()
    for value in blocks["query_geometry_wkt"]:
        geometry = from_wkt(value)
        assert geometry.is_valid
        for polygon in polygon_parts(geometry):
            assert polygon.exterior.is_ccw
            assert all(not ring.is_ccw for ring in polygon.interiors)
