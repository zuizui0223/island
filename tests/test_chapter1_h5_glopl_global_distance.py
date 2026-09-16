from __future__ import annotations

from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import yaml
from shapely.geometry import Point, Polygon

from island_v2.chapter1_h5_glopl_global_distance import (
    aggregate_measurement_cells,
    assign_context,
    build_context_design,
    build_study_key,
    classify_global_result,
    compute_distance_to_continent_union,
    evaluate_preflight_gate,
    read_preflight_metadata,
    select_seeded_continent_union,
)


def _config() -> dict:
    with open("config/chapter1_h5_glopl_global_distance_v1.yml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def test_preflight_reader_never_materializes_effect_columns(tmp_path: Path) -> None:
    path = tmp_path / "glopl.csv"
    pd.DataFrame(
        {
            "Latitude": [1.0],
            "Longitude": [2.0],
            "DOI": ["10.1/x"],
            "Author": ["A"],
            "Year": ["2001"],
            "PL_Effect_Size": [9.9],
            "PL_Effect_Size_Type1": ["FS"],
        }
    ).to_csv(path, index=False)
    got = read_preflight_metadata(path, _config())
    assert list(got.columns) == _config()["preflight"]["allowed_columns"]
    assert "PL_Effect_Size" not in got.columns


def test_study_key_prefers_doi_and_falls_back_to_author_year() -> None:
    frame = pd.DataFrame(
        {
            "DOI": [" 10.1/ABC ", "", ""],
            "Author": ["X", "Smith", ""],
            "Year": ["1999", "2005", ""],
        }
    )
    assert build_study_key(frame).tolist() == [
        "doi:10.1/abc",
        "author_year:Smith|2005",
        "",
    ]


def test_seeded_continent_distance_is_zero_on_mainland_and_positive_offshore() -> None:
    land = gpd.GeoDataFrame(
        [
            {"geometry": Polygon([(-2, -2), (2, -2), (2, 2), (-2, 2)])},
            {"geometry": Polygon([(8, -1), (10, -1), (10, 1), (8, 1)])},
        ],
        crs="EPSG:4326",
    )
    union = select_seeded_continent_union(
        land,
        {"a": (0.0, 0.0)},
        projected_crs="EPSG:3857",
    )
    points = gpd.GeoSeries([Point(0, 0), Point(5, 0)], crs="EPSG:4326")
    distance = compute_distance_to_continent_union(points, union, projected_crs="EPSG:3857")
    assert distance[0] == pytest.approx(0.0)
    assert distance[1] > 300.0


def test_context_assignment_matches_chapter1_latitude_regimes() -> None:
    got = assign_context(pd.Series([0.0, 30.0, 70.0, -30.0]), _config())
    assert got.tolist() == [
        "tropical",
        "northern_midlatitude",
        "northern_high_latitude",
        "southern_extratropical",
    ]


def test_preflight_gate_requires_global_and_primary_pair_support() -> None:
    cfg = _config()
    summary = {
        "n_unique_sites": 500,
        "n_unique_studies": 200,
        "by_context": {
            "northern_midlatitude": {"n_sites": 100, "n_studies": 50},
            "tropical": {"n_sites": 100, "n_studies": 50},
            "northern_high_latitude": {"n_sites": 5, "n_studies": 3},
            "southern_extratropical": {"n_sites": 50, "n_studies": 20},
        },
    }
    got = evaluate_preflight_gate(summary, cfg)
    assert got["global_gradient_admitted"]
    assert got["north_tropical_pair_admitted"]
    assert not got["four_context_heterogeneity_admitted"]


def test_measurement_cell_aggregation_preserves_publication_total_weight_one() -> None:
    rows = pd.DataFrame(
        {
            "study_key": ["doi:a", "doi:a", "doi:a", "doi:b"],
            "site_key": ["1|1", "1|1", "2|2", "3|3"],
            "analysis_regime": ["tropical", "tropical", "northern_midlatitude", "tropical"],
            "log1p_distance_to_major_continent_km": [0.0, 0.0, 2.0, 1.0],
            "PL_Effect_Size": [1.0, 3.0, 5.0, 7.0],
            "PL_Effect_Size_Type1": ["FS", "FS", "SO", "FS"],
            "PL_Effect_Size_Type2": ["Sup", "Sup", "Sup", "Sup"],
            "Constant_added": ["FALSE"] * 4,
            "Level_of_Supplementation": ["whole_plant"] * 4,
        }
    )
    got = aggregate_measurement_cells(rows, _config())
    a = got.loc[got["study_key"].eq("doi:a")]
    assert len(a) == 2
    assert a.loc[a["site_key"].eq("1|1"), "PL_Effect_Size"].iloc[0] == pytest.approx(2.0)
    assert a["analysis_weight"].sum() == pytest.approx(1.0)


def test_context_design_contains_distance_interactions_and_measurement_controls() -> None:
    frame = pd.DataFrame(
        {
            "analysis_regime": ["northern_midlatitude", "tropical", "southern_extratropical"],
            "z_log1p_distance_to_major_continent_km": [-1.0, 1.0, 0.5],
            "PL_Effect_Size_Type1": ["FS", "SO", "FS"],
            "PL_Effect_Size_Type2": ["Sup", "Sup", "Bagout"],
            "Constant_added": ["FALSE", "FALSE", "TRUE"],
            "Level_of_Supplementation": ["whole_plant", "whole_plant", "partial_plant"],
        }
    )
    design, names = build_context_design(frame, _config())
    assert design.shape[0] == 3
    assert "z_distance" in names
    assert "z_distance:context_tropical" in names
    assert "z_distance:context_southern_extratropical" in names
    assert any(name.startswith("measure:PL_Effect_Size_Type1=") for name in names)
    assert np.isfinite(design).all()


def test_classification_does_not_promote_direction_only_without_frozen_p_values() -> None:
    cfg = _config()
    global_result = {"evaluable": True, "distance_slope": 0.3, "one_sided_positive_p": 0.10}
    pair = {
        "evaluable": True,
        "northern_distance_slope": 0.4,
        "northern_one_sided_positive_p": 0.03,
        "interaction_tropical_minus_northern": -0.2,
        "interaction_one_sided_negative_p": 0.20,
    }
    heterogeneity = {"evaluable": True, "joint_p": 0.04}
    sensitivities = {
        "supplemental_only": {"global_distance_slope": 0.2, "north_slope": 0.2, "interaction": -0.1},
        "no_zero_constant": {"global_distance_slope": 0.4, "north_slope": 0.3, "interaction": -0.2},
    }
    got = classify_global_result(global_result, pair, heterogeneity, sensitivities, cfg)
    assert not got["global_gradient_supported"]
    assert not got["north_tropical_specificity_supported"]
    assert got["four_context_heterogeneity_supported"]
    assert not got["causal_pollination_mechanism_identified"]
