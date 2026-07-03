import geopandas as gpd
from shapely import from_wkt
from shapely.geometry import Polygon

from island_v2.gbif_blocks import (
    block_query_geometries,
    build_query_catchments,
    split_at_antimeridian,
)


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


def test_split_at_antimeridian_leaves_ordinary_geometry_untouched():
    normal = Polygon([(10, 10), (10, 11), (11, 11), (11, 10), (10, 10)])
    parts = split_at_antimeridian(normal)
    assert parts == [normal]


def test_split_at_antimeridian_recovers_true_small_extent():
    # Vertices at 179 and -179 are two degrees apart across the seam, but a
    # naive (non-wrapping) WKT ring between them spans 358 degrees the long
    # way around -- exactly what GBIF rejects as an antimeridian-crossing
    # polygon.
    crossing = Polygon([(179, 60), (179, 61), (-179, 61), (-179, 60), (179, 60)])
    minx, _, maxx, _ = crossing.bounds
    assert maxx - minx > 180.0

    parts = split_at_antimeridian(crossing)

    assert len(parts) == 2
    for part in parts:
        assert part.is_valid
        pminx, _, pmaxx, _ = part.bounds
        assert pmaxx - pminx <= 180.0
    assert abs(sum(part.area for part in parts) - 2.0) < 1e-9


def test_split_at_antimeridian_repairs_invalid_shifted_geometry():
    # Recentring a self-intersecting seam-crossing ring on the seam leaves an
    # invalid polygon; intersecting that directly makes GEOS raise
    # "TopologyException: side location conflict" (this is the real crash that
    # took down the frozen GSHHG submission run on a dateline landmass near
    # 65N). split_at_antimeridian must repair before cutting.
    bowtie = Polygon([(179, 64), (-179, 66), (179, 66), (-179, 64), (179, 64)])
    assert not bowtie.is_valid

    parts = split_at_antimeridian(bowtie)

    assert len(parts) == 2
    for part in parts:
        assert part.is_valid
        pminx, _, pmaxx, _ = part.bounds
        assert pmaxx - pminx <= 180.0


def test_regional_blocks_split_seam_crossing_island_into_two_valid_blocks():
    islands = gpd.GeoDataFrame(
        {
            "island_id": ["dateline_island", "normal_island"],
            "area_km2": [10.0, 12.0],
            "geometry_sha256": ["hd", "hn"],
        },
        geometry=[
            Polygon(
                [(179.9, 60), (179.9, 60.1), (-179.9, 60.1), (-179.9, 60), (179.9, 60)]
            ),
            Polygon([(10, 10), (10, 10.1), (10.1, 10.1), (10.1, 10), (10, 10)]),
        ],
        crs=4326,
    )

    catchments = build_query_catchments(islands, buffer_m=2_000, simplify_m=500)

    assert set(catchments["analysis_island_id"]) == {"dateline_island", "normal_island"}
    assert (catchments["analysis_island_id"] == "dateline_island").sum() == 2
    assert catchments["island_id"].is_unique

    blocks, members = block_query_geometries(
        catchments,
        grid_degrees=30,
        max_wkt_chars=120_000,
        max_islands_per_block=125,
    )

    assert len(members) == 3
    assert set(members["analysis_island_id"]) == {"dateline_island", "normal_island"}
    for value in blocks["query_geometry_wkt"]:
        geometry = from_wkt(value)
        assert geometry.is_valid
        minx, _, maxx, _ = geometry.bounds
        assert maxx - minx <= 180.0


def test_block_id_is_stable_when_an_unrelated_cell_gains_extra_islands():
    # A "frozen" campaign persists ledger entries keyed by block_id and refuses
    # to resubmit anything already accepted by GBIF. If block IDs embedded a
    # manifest-wide sequence number, fixing/changing one regional cell (adding
    # rows there, e.g. from antimeridian splitting) would renumber every other
    # cell's blocks and silently invalidate already-recorded provenance
    # hashes for blocks that were never touched. block_id must depend only on
    # (regional_cell, local_sequence).
    def make_islands(n_extra_in_first_cell: int) -> gpd.GeoDataFrame:
        ids = ["far_island"]
        geoms = [Polygon([(80, 10), (80, 10.1), (80.1, 10.1), (80.1, 10), (80, 10)])]
        for i in range(n_extra_in_first_cell):
            offset = i * 0.2
            ids.append(f"near_origin_{i}")
            geoms.append(
                Polygon(
                    [
                        (offset, 0),
                        (offset, 0.1),
                        (offset + 0.1, 0.1),
                        (offset + 0.1, 0),
                        (offset, 0),
                    ]
                )
            )
        return gpd.GeoDataFrame(
            {
                "island_id": ids,
                "area_km2": [10.0] * len(ids),
                "geometry_sha256": [f"h{i}" for i in range(len(ids))],
            },
            geometry=geoms,
            crs=4326,
        )

    def block_ids_for(n_extra_in_first_cell: int) -> set[str]:
        catchments = build_query_catchments(
            make_islands(n_extra_in_first_cell), buffer_m=2_000, simplify_m=500
        )
        blocks, _ = block_query_geometries(
            catchments, grid_degrees=30, max_wkt_chars=120_000, max_islands_per_block=1
        )
        return set(blocks["block_id"])

    baseline = block_ids_for(1)
    with_extra_island_in_a_different_cell = block_ids_for(2)

    far_island_block_ids = {b for b in baseline if "cell_w+060.0" in b}
    assert far_island_block_ids
    assert far_island_block_ids.issubset(with_extra_island_in_a_different_cell)
