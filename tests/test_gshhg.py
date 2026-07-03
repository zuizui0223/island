import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon

from island_v2.gshhg_source import make_island_units


def test_multipolygon_landmass_is_not_split():
    first = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    second = Polygon([(2, 0), (2, 1), (3, 1), (3, 0), (2, 0)])
    source = gpd.GeoDataFrame(
        {"gshhg_polygon_id": ["x"]},
        geometry=[MultiPolygon([first, second])],
        crs=4326,
    )

    result = make_island_units(
        source,
        source_label="test",
        min_area_km2=20,
        mainland_area_threshold_km2=7_000_000,
    )

    assert len(result) == 1
    assert result.iloc[0].geometry.geom_type == "MultiPolygon"
