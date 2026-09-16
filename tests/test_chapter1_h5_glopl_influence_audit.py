from __future__ import annotations

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_h5_glopl_influence_audit import (
    classify_influence,
    recompute_publication_weights,
    summarize_influence,
)


def _config() -> dict:
    with open("config/chapter1_h5_glopl_influence_audit_v1.yml", encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def test_recompute_publication_weights_keeps_total_one_after_deletion() -> None:
    cells = pd.DataFrame(
        {
            "study_key": ["a", "a", "a", "b"],
            "island_id": ["i1", "i2", "i3", "i1"],
        }
    )
    kept = cells.loc[~cells["island_id"].eq("i1")].copy()
    got = recompute_publication_weights(kept)
    assert pytest.approx(got.loc[got["study_key"].eq("a"), "analysis_weight"].sum()) == 1.0
    assert pytest.approx(got.loc[got["study_key"].eq("b"), "analysis_weight"].sum()) == 0.0


def test_summarize_influence_reports_ranges_and_max_delta() -> None:
    refits = pd.DataFrame(
        {
            "removed_unit": ["i1", "i2", "i3"],
            "evaluable": [True, True, True],
            "northern_distance_slope": [0.2, 0.4, 0.5],
            "interaction_tropical_minus_northern": [-0.1, 0.1, -0.2],
        }
    )
    out = summarize_influence(
        refits,
        full_north=0.3,
        full_interaction=-0.05,
    )
    assert out["north_slope_min"] == 0.2
    assert out["north_slope_max"] == 0.5
    assert out["north_slope_sign_positive_fraction"] == 1.0
    assert pytest.approx(out["interaction_sign_negative_fraction"]) == 2 / 3
    assert out["most_influential_removed_unit_for_north_slope"] == "i3"


def test_classification_never_promotes_parent_and_separates_sign_stability() -> None:
    island = {
        "evaluable_fraction": 1.0,
        "north_slope_sign_positive_fraction": 1.0,
        "interaction_sign_negative_fraction": 0.8,
    }
    publication = {
        "evaluable_fraction": 1.0,
        "north_slope_sign_positive_fraction": 1.0,
        "interaction_sign_negative_fraction": 1.0,
    }
    out = classify_influence(island, publication, _config())
    assert out["northern_direction"] == "northern_direction_diffuse"
    assert out["context_interaction"] == "context_interaction_sign_fragile"
    assert out["parent_classification"] == "pollen_limitation_gradient_not_supported"
    assert not out["can_promote_parent_result"]
