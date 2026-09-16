from __future__ import annotations

import pandas as pd

from island_v2.chapter1_h3_observed_taxonomic_depth import (
    build_atomic_taxonomic_residuals,
    build_island_stage_scores,
    classify_depth,
)


def _probability_config():
    return {
        "broad_outcomes": {
            "generalized_form": {
                "trait_name": "floral_form",
                "positive_states": ["open_radial"],
                "negative_states": ["tubular"],
            }
        }
    }


def _depth_config():
    return {
        "atomic_outcomes": ["generalized_form"],
        "species_taxonomic_expectation": {
            "minimum_scored_species_per_group_including_focal": 2
        },
    }


def test_loo_taxonomic_means_exclude_focal_species():
    audit = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "A c", "B a"],
            "trait_name": ["floral_form"] * 4,
            "resolved_for_primary": [True] * 4,
            "canonical_signature": ["open_radial", "tubular", "open_radial", "tubular"],
        }
    )
    taxonomy = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b", "A c", "B a"],
            "family": ["F", "F", "F", "G"],
            "genus": ["A", "A", "A", "B"],
        }
    )
    result = build_atomic_taxonomic_residuals(
        audit, taxonomy, _probability_config(), _depth_config()
    ).set_index("accepted_species")
    # A a has state 1; its LOO genus mean is mean(A b=0, A c=1)=0.5.
    assert result.loc["A a", "genus_loo_mean"] == 0.5
    # A b has state 0; its LOO genus mean is mean(A a=1, A c=1)=1.
    assert result.loc["A b", "genus_loo_mean"] == 1.0
    # Singleton genus B fails common family/genus support and is excluded.
    assert "B a" not in result.index


def test_island_scores_use_identical_species_support_across_stages():
    residuals = pd.DataFrame(
        {
            "accepted_species": ["A a", "A b"],
            "outcome": ["generalized_form", "generalized_form"],
            "family": ["F", "F"],
            "genus": ["A", "A"],
            "observed_score": [1.0, 0.0],
            "family_loo_mean": [0.0, 1.0],
            "genus_loo_mean": [0.0, 1.0],
            "after_family_residual": [1.0, -1.0],
            "after_genus_residual": [1.0, -1.0],
            "family_n": [2, 2],
            "genus_n": [2, 2],
        }
    )
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1"],
            "accepted_species": ["A a", "A b", "X x"],
        }
    )
    score = build_island_stage_scores(flora, residuals).iloc[0]
    assert score["n_species"] == 2
    assert score["observed_score"] == 0.5
    assert score["after_family_residual"] == 0.0
    assert score["after_genus_residual"] == 0.0


def test_classification_requires_observed_common_support_signal():
    omnibus = pd.DataFrame(
        {
            "stage": ["observed_score", "after_family_residual", "after_genus_residual"],
            "p_value": [0.01, 0.02, 0.2],
        }
    )
    assert classify_depth(omnibus, 0.05) == "compatible_with_genus_structuring"
    omnibus.loc[omnibus["stage"].eq("observed_score"), "p_value"] = 0.2
    assert classify_depth(omnibus, 0.05) == "observed_common_support_not_reproduced"
