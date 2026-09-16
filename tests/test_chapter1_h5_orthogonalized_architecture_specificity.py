from __future__ import annotations

import pandas as pd
import yaml

from island_v2.chapter1_h5_orthogonalized_architecture_specificity import (
    focused_configs,
    integrate_scope_decisions,
    summarize_scope,
)


def _config() -> dict:
    return yaml.safe_load(
        open(
            "config/chapter1_h5_orthogonalized_architecture_specificity_v1.yml",
            encoding="utf-8",
        )
    )


def test_focused_configs_use_only_all_observed_north_tropical_and_two_decomposition_families() -> None:
    pattern = {
        "contract": "base",
        "contexts": [
            "northern_midlatitude",
            "northern_high_latitude",
            "tropical",
            "southern_extratropical",
        ],
        "strata": ["all_observed", "all_native"],
        "support_tiers": {"confirmatory": 50, "pilot": 30},
        "minimum_outcomes_per_vector": 2,
    }
    branching = {
        "contract": "base_branching",
        "alpha": 0.05,
        "branch_axes": {"old": {"components": {"old": 1.0}}},
        "axis_sets": {"old": {"axes": ["old"], "role": "old", "classify": True}},
        "context_layers": {
            "analysis_regime": {
                "column": "analysis_regime",
                "contexts": pattern["contexts"],
                "role": "primary_when_where_layer",
            },
            "biogeographic_realm": {
                "column": "biogeographic_realm",
                "contexts": ["Palearctic", "Neotropical"],
                "role": "sensitivity",
            },
        },
    }

    focused_pattern, focused_branching = focused_configs(pattern, branching, _config())

    assert focused_pattern["contexts"] == ["northern_midlatitude", "tropical"]
    assert focused_pattern["strata"] == ["all_observed"]
    assert focused_pattern["support_tiers"] == {"confirmatory": 50}
    assert focused_pattern["minimum_outcomes_per_vector"] == 1
    assert set(focused_branching["axis_sets"]) == {
        "shared_architecture",
        "identity_specific_residuals",
    }
    assert focused_branching["axis_sets"]["shared_architecture"]["axes"] == [
        "shared_architecture_factor"
    ]
    assert focused_branching["axis_sets"]["identity_specific_residuals"]["axes"] == [
        "large_bee_like_residual",
        "butterfly_like_residual",
        "bird_like_residual",
    ]
    assert set(focused_branching["context_layers"]) == {"analysis_regime"}
    assert pattern["strata"] == ["all_observed", "all_native"]
    assert set(branching["axis_sets"]) == {"old"}


def test_summarize_scope_requires_all_three_residual_components() -> None:
    between = pd.DataFrame(
        [
            {
                "context_layer": "analysis_regime",
                "axis_set": "identity_specific_residuals",
                "stratum": "all_observed",
                "support_tier": "confirmatory",
                "context_a": "northern_midlatitude",
                "context_b": "tropical",
                "status": "fit",
                "n_retained_syndromes": 2,
                "retained_syndromes": "large_bee_like_residual|butterfly_like_residual",
                "joint_context_difference_df": 2,
                "p_value": 0.001,
                "n_unique_islands": 100,
                "n_clusters": 20,
            },
            {
                "context_layer": "analysis_regime",
                "axis_set": "shared_architecture",
                "stratum": "all_observed",
                "support_tier": "confirmatory",
                "context_a": "northern_midlatitude",
                "context_b": "tropical",
                "status": "fit",
                "n_retained_syndromes": 1,
                "retained_syndromes": "shared_architecture_factor",
                "joint_context_difference_df": 1,
                "p_value": 0.02,
                "n_unique_islands": 100,
                "n_clusters": 20,
            },
        ]
    )
    slopes = pd.DataFrame()

    result = summarize_scope(
        between,
        slopes,
        evidence_scope="all_analysis_eligible",
        factor_audit={
            "variance_fraction": 0.86,
            "n_complete_scored_species": 1200,
            "n_complete_gift_source_species": 900,
            "source_factor_residual_covariances": {
                "large_bee_like_residual": 0.0,
                "butterfly_like_residual": 0.0,
                "bird_like_residual": 0.0,
            },
        },
        config=_config(),
    )

    assert not result["identity_specific_evaluable"]
    assert result["identity_specific_joint_p"] is None
    assert result["shared_architecture_joint_p"] == 0.02


def test_integrated_promotion_requires_residual_support_in_both_evidence_scopes() -> None:
    primary = {
        "evidence_scope": "all_analysis_eligible",
        "identity_specific_evaluable": True,
        "identity_specific_joint_p": 0.01,
        "identity_specific_joint_df": 3,
    }
    direct = {
        "evidence_scope": "direct_only",
        "identity_specific_evaluable": True,
        "identity_specific_joint_p": 0.20,
        "identity_specific_joint_df": 3,
    }
    not_promoted = integrate_scope_decisions(primary, direct, _config())
    assert not not_promoted["robust_identity_specific_structure"]
    assert not_promoted["classification"] == "all_analysis_only_identity_specific_signal_not_promoted"

    direct["identity_specific_joint_p"] = 0.03
    promoted = integrate_scope_decisions(primary, direct, _config())
    assert promoted["robust_identity_specific_structure"]
    assert promoted["classification"] == "robust_identity_specific_architecture_structure"
    assert not promoted["pollinator_identity_promoted"]
