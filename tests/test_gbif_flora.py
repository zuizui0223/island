import json

import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon

from island_v2.gbif_flora import gbif_payload, request_rows, stable_island_table


def test_prepare_splits_multipart_and_filters_small_polygons():
    first = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    second = Polygon([(2, 0), (2, 1), (3, 1), (3, 0), (2, 0)])
    tiny = Polygon([(4, 0), (4, 0.001), (4.001, 0.001), (4.001, 0), (4, 0)])
    source = gpd.GeoDataFrame(
        {"name": ["A", "tiny"], "parent": ["p1", "p2"]},
        geometry=[MultiPolygon([first, second]), tiny],
        crs=4326,
    )

    output = stable_island_table(source, "demo", 20, "name", "parent")

    assert len(output) == 2
    assert set(output["parent_feature_id"]) == {"p1"}
    assert output["island_id"].is_unique


def test_payload_uses_exact_geometry_and_quality_filters():
    payload = gbif_payload(
        "POLYGON((0 0, 1 0, 1 1, 0 0))",
        "123",
        "catalogue",
        "SPECIES_LIST",
        None,
        2000,
    )
    predicates = payload["predicate"]["predicates"]

    assert payload["format"] == "SPECIES_LIST"
    assert predicates[0]["type"] == "within"
    assert predicates[1]["key"] == "TAXON_KEY"
    assert any(item.get("key") == "HAS_GEOSPATIAL_ISSUE" for item in predicates)
    assert any(item.get("key") == "YEAR" for item in predicates)


def test_request_rows_never_replaces_polygon_with_bbox():
    triangle = Polygon([(0, 0), (2, 0), (0, 2), (0, 0)])
    islands = gpd.GeoDataFrame(
        {"island_id": ["i1"], "area_km2": [100.0]},
        geometry=[triangle],
        crs=4326,
    )

    rows = request_rows(islands, {"usage_key": "123"}, None, "SPECIES_LIST", None, None)
    payload = json.loads(rows[0]["payload_json"])

    assert "POLYGON ((0 0, 2 0, 0 2, 0 0))" in rows[0]["geometry_wkt"]
    assert payload["predicate"]["predicates"][0]["geometry"] == rows[0]["geometry_wkt"]
