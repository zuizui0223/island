from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_global_branch_selection import run_selection_sensitivity


def _configs():
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
    }
    branching = {
        "alpha": 0.05,
        "branch_axes": {
            "accessibility_generalization": {
                "components": {"generalized_accessible": 1.0}
            },
            "reproductive_assurance": {"components": {"selfing_core": 1.0}},
        },
        "axis_sets": {
            "universal_plant_response": {
                "axes": ["accessibility_generalization", "reproductive_assurance"],
                "role": "primary_when_where_test",
                "classify": True,
            }
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
    selection = {
        "primary_axis_set": "universal_plant_response",
        "selection_model": {"minimum_observed": 30, "minimum_unobserved": 30},
        "weight_modes": {
            "unweighted": {},
            "stabilized_ipw_clip_0_05": {
                "minimum_propensity": 0.05,
                "maximum_propensity": 0.95,
                "maximum_weight": 20.0,
            },
        },
        "stabilization": {"normalize_within_axis_stratum_context": True},
        "headline_contrasts": {
            "analysis_regime": [["northern_midlatitude", "tropical"]],
            "biogeographic_realm": [["Palearctic", "Neotropical"]],
        },
    }
    return pattern, branching, selection


def _synthetic():
    rng = np.random.default_rng(1383)
    scores = []
    covariates = []
    realms = []
    n_per_context = 180
    for context, realm in [
        ("northern_midlatitude", "Palearctic"),
        ("tropical", "Neotropical"),
    ]:
        for i in range(n_per_context):
            island = f"{context}_{i}"
            distance = -1 + 2 * i / (n_per_context - 1)
            covariates.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_{i // 6}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": np.sin(i / 9),
                    "climate_pc1": np.cos(i / 11),
                }
            )
            realms.append({"island_id": island, "biogeographic_realm": realm})
            selection_probability = 0.35 + 0.25 / (1 + np.exp(-2 * distance))
            if rng.random() >= selection_probability:
                continue
            if context == "northern_midlatitude":
                access = 0.30 * distance + rng.normal(0, 0.03)
                selfing = 0.03 * distance + rng.normal(0, 0.03)
            else:
                access = -0.10 * distance + rng.normal(0, 0.03)
                selfing = 0.30 * distance + rng.normal(0, 0.03)
            for syndrome, value in {
                "generalized_accessible": access,
                "selfing_core": selfing,
            }.items():
                scores.append(
                    {
                        "island_id": island,
                        "stratum": "all_native",
                        "syndrome": syndrome,
                        "syndrome_score": value,
                        "n_species": 10,
                    }
                )
    return pd.DataFrame(scores), pd.DataFrame(covariates), pd.DataFrame(realms)


def test_selection_sensitivity_is_bounded_and_outcome_blind():
    pattern, branching, selection = _configs()
    scores, covariates, realms = _synthetic()
    weights, fits, diagnostics, _, _, between, states, summary = (
        run_selection_sensitivity(
            scores,
            covariates,
            realms,
            pattern,
            branching,
            selection,
        )
    )
    assert set(weights["selection_mode"]) == {
        "unweighted",
        "stabilized_ipw_clip_0_05",
    }
    assert weights["analysis_weight"].gt(0).all()
    assert weights["analysis_weight"].le(20).all()
    assert fits["n_full_islands"].eq(360).all()
    assert diagnostics["effective_sample_fraction"].between(0, 1).all()
    assert set(states["selection_mode"]) == {
        "unweighted",
        "stabilized_ipw_clip_0_05",
    }
    assert not between.empty
    assert len(summary) == 4


def test_missing_score_is_not_inserted_as_a_zero_response():
    pattern, branching, selection = _configs()
    scores, covariates, realms = _synthetic()
    observed_ids = set(scores["island_id"])
    outputs = run_selection_sensitivity(
        scores,
        covariates,
        realms,
        pattern,
        branching,
        selection,
    )
    weights = outputs[0]
    assert set(weights["island_id"]).issubset(observed_ids)
    assert len(observed_ids) < len(covariates)
