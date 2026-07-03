import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon
from typer.testing import CliRunner

from island_v2.gshhg_source import app, make_island_units, natural_earth_parts


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
    # Force a wide terminal: typer/rich wrap --help onto multiple lines (and
    # can truncate long option names) once COLUMNS drops below its default,
    # which made this assertion flaky under narrow CI/sandbox terminals.
    result = CliRunner().invoke(app, ["build", "--help"], env={"COLUMNS": "200"})

    assert result.exit_code == 0
    assert "--allow-natural-earth-fallback" in result.output
