from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_global_branching import (
    build_branch_scores,
    run_global_branching,
)


def _branching_config() -> dict:
    return {
        "alpha": 0.05,
        "branch_axes": {
            "accessibility_generalization": {
                "components": {"generalized_accessible": 1.0}
            },
            "reproductive_assurance": {"components": {"selfing_core": 1.0}},
            "sampled_large_bee_concordance": {"components": {"large_bee_like": 1.0}},
            "sampled_butterfly_concordance": {"components": {"butterfly_like": 1.0}},
            "sampled_bird_concordance": {"components": {"bird_like": 1.0}},
        },
        "axis_sets": {
            "universal_plant_response": {
                "axes": ["accessibility_generalization", "reproductive_assurance"],
                "role": "primary_when_where_test",
                "classify": True,
            },
            "sampled_guild_concordance": {
                "axes": [
                    "sampled_large_bee_concordance",
                    "sampled_butterfly_concordance",
                    "sampled_bird_concordance",
                ],
                "role": "secondary_discussion_concordance",
                "classify": False,
            },
        },
        "context_layers": {
            "analysis_regime": {
                "column": "analysis_regime",
                "contexts": ["northern_midlatitude", "tropical"],
            },
            "biogeographic_realm": {
                "column": "biogeographic_realm",
                "contexts": ["Palearctic", "Neotropical"],
            },
        },
    }


def _pattern_config() -> dict:
    return {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
    }


def _synthetic(n_per_context: int = 120):
    rng = np.random.default_rng(1382)
    scores = []
    covariates = []
    assignments = []
    for context, realm in [
        ("northern_midlatitude", "Palearctic"),
        ("tropical", "Neotropical"),
    ]:
        cluster_noise = rng.normal(0, 0.025, size=n_per_context // 5)
        for i in range(n_per_context):
            island = f"{context}_{i}"
            distance = -1.0 + 2.0 * i / (n_per_context - 1)
            noise = cluster_noise[i // 5] + rng.normal(0, 0.01)
            if context == "northern_midlatitude":
                access = 0.35 * distance + noise
                selfing = 0.28 * distance + noise
                butterfly = -0.18 * distance + noise
            else:
                access = 0.01 * distance + noise
                selfing = 0.01 * distance + noise
                butterfly = 0.32 * distance + noise
            values = {
                "large_bee_like": -access,
                "generalized_accessible": access,
                "selfing_core": selfing,
                "butterfly_like": butterfly,
                "bird_like": 0.01 * distance + noise,
            }
            for syndrome, value in values.items():
                scores.append(
                    {
                        "island_id": island,
                        "stratum": "all_native",
                        "syndrome": syndrome,
                        "syndrome_score": value,
                        "n_species": 20,
                    }
                )
            covariates.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_b{i // 5}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": np.sin(i / 7),
                    "climate_pc1": np.cos(i / 9),
                }
            )
            assignments.append({"island_id": island, "biogeographic_realm": realm})
    return pd.DataFrame(scores), pd.DataFrame(covariates), pd.DataFrame(assignments)


def test_branch_score_requires_every_predeclared_component():
    scores = pd.DataFrame(
        [
            {
                "island_id": "complete",
                "stratum": "all_native",
                "syndrome": "large_bee_like",
                "syndrome_score": -0.4,
                "n_species": 10,
            },
            {
                "island_id": "complete",
                "stratum": "all_native",
                "syndrome": "generalized_accessible",
                "syndrome_score": 0.6,
                "n_species": 8,
            },
            {
                "island_id": "missing",
                "stratum": "all_native",
                "syndrome": "large_bee_like",
                "syndrome_score": -0.2,
                "n_species": 12,
            },
        ]
    )
    config = {
        "branch_axes": {
            "accessibility_generalization": {
                "components": {"generalized_accessible": 1.0}
            }
        }
    }
    out = build_branch_scores(scores, config)
    assert out["island_id"].tolist() == ["complete"]
    assert out.iloc[0]["syndrome_score"] == 0.6
    assert out.iloc[0]["n_species"] == 8


def test_global_branching_separates_access_selfing_from_alternative_specialization():
    scores, covariates, assignments = _synthetic()
    _, slopes, _, between, states = run_global_branching(
        scores,
        covariates,
        assignments,
        _pattern_config(),
        _branching_config(),
    )
    primary = states.loc[
        states["context_layer"].eq("analysis_regime")
        & states["support_tier"].eq("confirmatory")
    ].set_index("context")
    north = primary.loc["northern_midlatitude"]
    tropical = primary.loc["tropical"]
    assert bool(north["context_vector_supported"])
    assert bool(north["accessibility_generalization_positive"])
    assert bool(north["reproductive_assurance_positive"])
    assert "accessibility_generalization" in north["plant_side_branch_state"]
    assert "reproductive_assurance" in north["plant_side_branch_state"]
    assert not bool(tropical["accessibility_generalization_positive"])

    secondary = slopes.loc[
        slopes["axis_set"].eq("sampled_guild_concordance")
        & slopes["context"].eq("tropical")
        & slopes["syndrome"].eq("sampled_butterfly_concordance")
    ].iloc[0]
    assert bool(secondary["axis_supported"])

    direct = between.loc[
        between["context_layer"].eq("analysis_regime")
        & between["axis_set"].eq("universal_plant_response")
        & between["context_a"].eq("northern_midlatitude")
        & between["context_b"].eq("tropical")
        & between["support_tier"].eq("confirmatory")
    ].iloc[0]
    assert bool(direct["context_vector_difference_supported"])
    fitted_axes = set(
        slopes.loc[
            slopes["context_layer"].eq("analysis_regime"), "syndrome"
        ].astype(str)
    )
    assert "accessibility_generalization" in fitted_axes
    assert "reproductive_assurance" in fitted_axes


def test_under_supported_context_is_not_labeled_as_a_branch():
    scores, covariates, assignments = _synthetic(n_per_context=25)
    _, _, _, _, states = run_global_branching(
        scores,
        covariates,
        assignments,
        _pattern_config(),
        _branching_config(),
    )
    confirmatory = states.loc[states["support_tier"].eq("confirmatory")]
    assert set(confirmatory["plant_side_branch_state"]) == {"not_testable"}
