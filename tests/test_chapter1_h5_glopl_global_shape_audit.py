from __future__ import annotations

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_h5_glopl_global_shape_audit import (
    build_shape_design,
    classify_shape_audit,
    prepare_shape_cells,
)


def _config() -> dict:
    with open("config/chapter1_h5_glopl_global_shape_audit_v1.yml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _cells() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "study_key": ["doi:a", "doi:a", "doi:b", "doi:c"],
            "site_key": ["m1", "i1", "i2", "m2"],
            "analysis_regime": ["northern_midlatitude", "northern_midlatitude", "tropical", "tropical"],
            "log1p_distance_to_major_continent_km": [0.0, 2.0, 4.0, 0.0],
            "PL_Effect_Size": [0.1, 0.4, 0.8, 0.2],
            "analysis_weight": [0.5, 0.5, 1.0, 1.0],
            "PL_Effect_Size_Type1": ["FS", "FS", "SO", "FS"],
            "PL_Effect_Size_Type2": ["Sup", "Sup", "Sup", "Bagout"],
            "Constant_added": ["FALSE", "FALSE", "FALSE", "TRUE"],
            "Level_of_Supplementation": ["whole", "whole", "whole", "partial"],
        }
    )


def test_prepare_shape_cells_reuses_parent_weights_without_renormalizing() -> None:
    got, audit = prepare_shape_cells(_cells(), _config())
    assert got["analysis_weight"].tolist() == [0.5, 0.5, 1.0, 1.0]
    assert audit["n_mainland_cells"] == 2
    assert audit["n_offshore_cells"] == 2
    assert got.loc[got["offshore"].eq(0.0), "offshore_distance_z"].eq(0.0).all()
    assert got.loc[got["offshore"].eq(1.0), "offshore_distance_z"].mean() == pytest.approx(0.0)


def test_offshore_distance_standardization_uses_only_positive_distance_cells() -> None:
    got, audit = prepare_shape_cells(_cells(), _config())
    offshore = got.loc[got["offshore"].eq(1.0), "offshore_distance_z"]
    assert offshore.tolist() == pytest.approx([-1.0, 1.0])
    assert audit["offshore_log_distance_mean"] == pytest.approx(3.0)
    assert audit["offshore_log_distance_sd"] == pytest.approx(1.0)


def test_shape_design_contains_step_and_within_offshore_slope() -> None:
    got, _ = prepare_shape_cells(_cells(), _config())
    design, names = build_shape_design(got)
    assert design.shape[0] == 4
    assert "offshore_indicator" in names
    assert "offshore_distance_z" in names
    assert "context_tropical" in names
    assert any(name.startswith("measure:PL_Effect_Size_Type1=") for name in names)


def test_classification_step_dominated_is_descriptive_only() -> None:
    two_part = {
        "evaluable": True,
        "offshore_step": 0.2,
        "offshore_step_two_sided_p": 0.01,
        "within_offshore_slope": 0.03,
        "within_offshore_slope_two_sided_p": 0.4,
    }
    offshore_only = {
        "evaluable": True,
        "distance_slope": 0.02,
        "distance_slope_two_sided_p": 0.5,
    }
    got = classify_shape_audit(two_part, offshore_only, _config())
    assert got["shape_classification"] == "step_dominated"
    assert not got["can_change_parent_promotion"]


def test_classification_requires_both_within_offshore_tests_for_gradient_label() -> None:
    two_part = {
        "evaluable": True,
        "offshore_step": 0.1,
        "offshore_step_two_sided_p": 0.2,
        "within_offshore_slope": 0.2,
        "within_offshore_slope_two_sided_p": 0.01,
    }
    offshore_only = {
        "evaluable": True,
        "distance_slope": 0.15,
        "distance_slope_two_sided_p": 0.02,
    }
    got = classify_shape_audit(two_part, offshore_only, _config())
    assert got["shape_classification"] == "within_offshore_gradient_present"
    assert not got["can_change_parent_promotion"]
