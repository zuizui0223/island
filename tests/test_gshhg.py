import geopandas as gpd
import pytest
from shapely.geometry import MultiPolygon, Polygon
from typer.main import get_command
from typer.testing import CliRunner

from island_v2.gshhg_source import app, dissolve_ids, make_island_units, natural_earth_parts


def test_gshhg_dateline_sibling_parts_are_dissolved_before_area_filtering():
    east = Polygon([(179, 0), (179, 1), (180, 1), (180, 0), (179, 0)])
    west = Polygon([(-180, 0), (-180, 1), (-179, 1), (-179, 0), (-180, 0)])
    source = gpd.GeoDataFrame(
        {"id": ["0-E", "0-W"], "sibling_id": [0, 0]},
        geometry=[east, west],
        crs=4326,
    )

    dissolved = dissolve_ids(source)

    assert len(dissolved) == 1
    assert dissolved.iloc[0]["source_feature_id"] == "0"
    assert dissolved.iloc[0].geometry.geom_type == "MultiPolygon"

    result = make_island_units(
        dissolved,
        source_label="test",
        min_area_km2=1,
        mainland_area_threshold_km2=20_000,
    )

    assert result.empty


def test_gshhg_distinct_sibling_ids_remain_distinct_landmasses():
    first = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    second = Polygon([(2, 0), (2, 1), (3, 1), (3, 0), (2, 0)])
    source = gpd.GeoDataFrame(
        {"id": ["10", "11"], "sibling_id": [110, 111]},
        geometry=[first, second],
        crs=4326,
    )

    result = dissolve_ids(source)

    assert len(result) == 2
    assert set(result["source_feature_id"]) == {"110", "111"}


@pytest.mark.parametrize("invalid_sibling_id", [None, "", -1, 1.5, "not-an-id"])
def test_gshhg_invalid_sibling_ids_fail_closed(invalid_sibling_id):
    polygon = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    source = gpd.GeoDataFrame(
        {"id": ["0-E"], "sibling_id": [invalid_sibling_id]},
        geometry=[polygon],
        crs=4326,
    )

    with pytest.raises(RuntimeError, match="sibling_id contains missing or invalid"):
        dissolve_ids(source)


def test_gshhg_dissolve_falls_back_to_literal_id_without_sibling_linkage():
    first = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    second = Polygon([(2, 0), (2, 1), (3, 1), (3, 0), (2, 0)])
    source = gpd.GeoDataFrame(
        {"ID": ["legacy", "legacy"]},
        geometry=[first, second],
        crs=4326,
    )

    result = dissolve_ids(source)

    assert len(result) == 1
    assert result.iloc[0]["source_feature_id"] == "legacy"
    assert result.iloc[0].geometry.geom_type == "MultiPolygon"


def test_gshhg_multipolygon_landmass_is_not_split():
    first = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    second = Polygon([(2, 0), (2, 1), (3, 1), (3, 0), (2, 0)])
    source = gpd.GeoDataFrame(
        {"source_feature_id": ["x"]},
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


def test_natural_earth_fallback_splits_multipolygon_components():
    first = Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])
    second = Polygon([(2, 0), (2, 1), (3, 1), (3, 0), (2, 0)])
    source = gpd.GeoDataFrame(
        geometry=[MultiPolygon([first, second])],
        crs=4326,
    )

    result = natural_earth_parts(source)

    assert len(result) == 2
    assert result["source_feature_id"].is_unique


def test_build_is_a_cli_subcommand_with_fallback_control():
    # Rich's --help rendering wraps/highlights option cells differently
    # depending on the runner's terminal width and colour support, which
    # made asserting against the rendered text flaky across environments
    # (it passed locally but failed in CI at both the default and a forced
    # COLUMNS=200). Inspect the registered Click command instead: it is
    # rendering-independent and is what actually defines the CLI's options.
    command = get_command(app)
    build_command = command.commands["build"]
    option_names = {opt for param in build_command.params for opt in param.opts}

    assert "--allow-natural-earth-fallback" in option_names

    result = CliRunner().invoke(app, ["build", "--help"])
    assert result.exit_code == 0
