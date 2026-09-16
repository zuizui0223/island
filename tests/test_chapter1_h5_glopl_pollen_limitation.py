from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_h5_glopl_pollen_limitation import (
    aggregate_publication_island,
    build_primary_design,
    classify_h5_result,
    prepare_analysis_rows,
    validate_source_sha256,
)


def _config() -> dict:
    with open("config/chapter1_h5_glopl_pollen_limitation_v1.yml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def test_source_hash_must_match_before_analysis(tmp_path: Path) -> None:
    path = tmp_path / "x.csv"
    path.write_text("a,b\n1,2\n", encoding="utf-8")
    with pytest.raises(ValueError, match="source SHA-256 mismatch"):
        validate_source_sha256(path, _config())


def test_prepare_rows_reuses_exact_preflight_row_id_and_primary_contexts() -> None:
    effects = pd.DataFrame(
        {
            "row_id": [0, 1, 2, 3],
            "PL_Effect_Size": [0.2, 0.4, 0.6, float("nan")],
            "PL_Effect_Size_Type1": ["FS"] * 4,
            "PL_Effect_Size_Type2": ["Sup"] * 4,
            "Constant_added": ["FALSE"] * 4,
            "Level_of_Supplementation": ["whole_plant"] * 4,
        }
    )
    matched = pd.DataFrame(
        {
            "row_id": [0, 1, 2, 3],
            "island_id": ["n1", "t1", pd.NA, "n1"],
            "study_key": ["doi:a", "doi:b", "doi:c", "doi:a"],
            "boundary_only_match": [False, False, False, False],
            "multi_polygon_match": [False, False, False, False],
            "analysis_regime": ["northern_midlatitude", "tropical", "tropical", "northern_midlatitude"],
        }
    )
    cov = pd.DataFrame(
        {
            "island_id": ["n1", "t1"],
            "log_distance_to_continent_km": [1.0, 2.0],
            "log_island_area_km2": [3.0, 4.0],
            "climate_pc1": [0.1, 0.2],
            "climate_pc2": [0.0, 0.1],
            "climate_pc3": [-0.1, 0.0],
            "climate_pc4": [0.2, -0.2],
            "spatial_block": ["b1", "b2"],
        }
    )
    got = prepare_analysis_rows(effects, matched, cov, _config())
    assert got["row_id"].tolist() == [0, 1]
    assert set(got["analysis_regime"]) == {"northern_midlatitude", "tropical"}


def test_publication_island_aggregation_gives_each_publication_total_weight_one() -> None:
    rows = pd.DataFrame(
        {
            "study_key": ["doi:a", "doi:a", "doi:a", "doi:b"],
            "island_id": ["i1", "i1", "i2", "i2"],
            "analysis_regime": ["northern_midlatitude"] * 4,
            "spatial_block": ["b1", "b1", "b2", "b2"],
            "PL_Effect_Size": [1.0, 3.0, 5.0, 7.0],
            "z_log_distance_to_continent_km": [0.0, 0.0, 1.0, 1.0],
            "z_log_island_area_km2": [0.0, 0.0, 1.0, 1.0],
            "z_climate_pc1": [0.0, 0.0, 0.0, 0.0],
            "z_climate_pc2": [0.0, 0.0, 0.0, 0.0],
            "z_climate_pc3": [0.0, 0.0, 0.0, 0.0],
            "z_climate_pc4": [0.0, 0.0, 0.0, 0.0],
        }
    )
    got = aggregate_publication_island(rows)
    a = got.loc[got["study_key"].eq("doi:a")]
    assert len(a) == 2
    assert a.loc[a["island_id"].eq("i1"), "PL_Effect_Size"].iloc[0] == 2.0
    assert pytest.approx(a["analysis_weight"].sum()) == 1.0
    b = got.loc[got["study_key"].eq("doi:b")]
    assert pytest.approx(b["analysis_weight"].sum()) == 1.0


def test_primary_design_has_frozen_distance_context_interaction() -> None:
    frame = pd.DataFrame(
        {
            "analysis_regime": ["northern_midlatitude", "tropical"],
            "z_log_distance_to_continent_km": [-1.0, 1.0],
            "z_log_island_area_km2": [0.0, 0.0],
            "z_climate_pc1": [0.0, 0.0],
            "z_climate_pc2": [0.0, 0.0],
            "z_climate_pc3": [0.0, 0.0],
            "z_climate_pc4": [0.0, 0.0],
        }
    )
    design, names = build_primary_design(frame)
    assert names == [
        "intercept",
        "context_tropical",
        "z_log_distance_to_continent_km",
        "z_log_distance_to_continent_km:context_tropical",
        "z_log_island_area_km2",
        "z_climate_pc1",
        "z_climate_pc2",
        "z_climate_pc3",
        "z_climate_pc4",
    ]
    assert design.shape == (2, 9)
    assert design[0, 3] == 0.0
    assert design[1, 3] == 1.0


def test_classification_separates_north_gradient_from_context_specificity() -> None:
    cfg = _config()
    primary = {
        "evaluable": True,
        "northern_distance_slope": 0.3,
        "northern_one_sided_positive_p": 0.01,
        "interaction_tropical_minus_northern": -0.4,
        "interaction_one_sided_negative_p": 0.20,
    }
    sensitivities = {
        "supplemental_only": {"evaluable": True, "northern_distance_slope": 0.2, "interaction_tropical_minus_northern": -0.1},
        "no_zero_constant": {"evaluable": True, "northern_distance_slope": 0.4, "interaction_tropical_minus_northern": -0.2},
    }
    got = classify_h5_result(primary, sensitivities, cfg)
    assert got["classification"] == "north_pollen_limitation_gradient_supported_context_specificity_not_established"
    assert got["north_service_gradient_supported"]
    assert not got["context_specific_service_gradient_supported"]

    primary["interaction_one_sided_negative_p"] = 0.02
    got = classify_h5_result(primary, sensitivities, cfg)
    assert got["classification"] == "context_specific_pollen_limitation_gradient_supported"
    assert got["full_mechanism_promotion"]
