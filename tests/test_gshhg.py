import geopandas as gpd
import pytest
from shapely.geometry import MultiPolygon, Polygon
from typer.main import get_command
from typer.testing import CliRunner

from island_v2.gshhg_source import app, dissolve_ids, make_island_units, natural_earth_parts


def test_sibling_fragments_are_recombined_before_mainland_filter():
    from shapely.geometry import box

    source = gpd.GeoDataFrame(
        {"id": ["0-E", "0-W", "1"], "sibling_id": [0, 0, 1]},
        geometry=[box(0, 0, 2, 2), box(4, 0, 5, 1), box(8, 0, 8.1, 0.1)],
        crs=4326,
    )
    combined = dissolve_ids(source)
    assert len(combined) == 2
    result = make_island_units(combined, "test", 1, 50_000)
    assert result.parent_feature_id.tolist() == ["1"]


@pytest.mark.parametrize("link", [None, -1, 0.5, "invalid"])
def test_invalid_sibling_linkage_is_rejected(link):
    from shapely.geometry import box

    source = gpd.GeoDataFrame(
        {"id": ["0-W"], "sibling_id": [link]}, geometry=[box(0, 0, 1, 1)], crs=4326
    )
    with pytest.raises(ValueError, match="sibling"):
        dissolve_ids(source)


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
