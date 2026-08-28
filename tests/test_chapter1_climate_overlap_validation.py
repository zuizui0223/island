from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_climate_overlap_validation import (
    run_climate_overlap_validation,
)


def _configs() -> tuple[dict, dict, dict]:
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
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
    explanation = {
        "validations": {
            "V1_H2_environment_vs_biogeography": {
                "analyses": [
                    {
                        "name": "outcome_blind_context_overlap_weighting",
                        "balance_predictors": [
                            "log_island_area_km2",
                            "climate_pc1",
                            "climate_pc2",
                            "climate_pc3",
                            "climate_pc4",
                        ],
                    }
                ],
                "frozen_implementation": {
                    "primary_axis_set": "universal_plant_response",
                    "headline_contrasts": {
                        "analysis_regime": [
                            ["northern_midlatitude", "tropical"]
                        ],
                        "biogeographic_realm": [["Palearctic", "Neotropical"]],
                    },
                    "overlap_model": {
                        "link": "binomial_logit",
                        "weighting": "overlap_weights",
                        "propensity_clip": 1.0e-8,
                        "normalize_mean_weight_within_context": True,
                        "minimum_context_islands": 50,
                        "poor_overlap_if_either": {
                            "minimum_effective_sample_fraction": 0.25,
                            "maximum_absolute_weighted_standardized_mean_difference": 0.10,
                        },
                    },
                    "distance_x_climate_model": {
                        "climate_interaction_predictors": [
                            "climate_pc1",
                            "climate_pc2",
                            "climate_pc3",
                            "climate_pc4",
                        ],
                        "retain_baseline_main_effects": True,
                        "target": "joint_distance_x_context_vector_after_distance_x_climate",
                    },
                    "multiplicity": {
                        "family": "validation_mode_x_axis_set_x_stratum_x_support_tier",
                        "alpha": 0.05,
                    },
                },
            }
        }
    }
    return pattern, branching, explanation


def _synthetic(
    *,
    n_per_context: int = 240,
    residual_context_slope: float = 0.0,
    seed: int = 1421,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    scores: list[dict] = []
    covariates: list[dict] = []
    realms: list[dict] = []
    for context, realm, climate_mean in [
        ("northern_midlatitude", "Palearctic", -0.75),
        ("tropical", "Neotropical", 0.75),
    ]:
        context_extra = residual_context_slope if context == "tropical" else 0.0
        cluster_noise = rng.normal(
            0, 0.015, size=max(1, (n_per_context + 5) // 6)
        )
        for i in range(n_per_context):
            island_id = f"{context}_{i}"
            distance = rng.uniform(-1.2, 1.2)
            climate_pc1 = rng.normal(climate_mean, 0.75)
            climate_pc2 = rng.normal(0.3 * climate_mean, 1.0)
            climate_pc3 = rng.normal(-0.2 * climate_mean, 1.0)
            climate_pc4 = rng.normal(0.1 * climate_mean, 1.0)
            area = rng.normal(0.4 * climate_mean, 0.9)
            noise = cluster_noise[i // 6] + rng.normal(0, 0.025)
            accessibility = (
                0.75 * distance * climate_pc1 + context_extra * distance + noise
            )
            assurance = (
                -0.60 * distance * climate_pc1
                + 0.8 * context_extra * distance
                + noise
            )
            for syndrome, value in {
                "generalized_accessible": accessibility,
                "selfing_core": assurance,
            }.items():
                scores.append(
                    {
                        "island_id": island_id,
                        "stratum": "all_native",
                        "syndrome": syndrome,
                        "syndrome_score": value,
                        "n_species": 20,
                    }
                )
            covariates.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_{i // 6}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": area,
                    "climate_pc1": climate_pc1,
                    "climate_pc2": climate_pc2,
                    "climate_pc3": climate_pc3,
                    "climate_pc4": climate_pc4,
                }
            )
            realms.append(
                {"island_id": island_id, "biogeographic_realm": realm}
            )
    return pd.DataFrame(scores), pd.DataFrame(covariates), pd.DataFrame(realms)


def test_overlap_weights_are_outcome_blind_and_improve_balance() -> None:
    pattern, branching, explanation = _configs()
    scores, covariates, realms = _synthetic()
    outputs = run_climate_overlap_validation(
        scores, covariates, realms, pattern, branching, explanation
    )
    weights, _, balance, summary, _, _ = outputs
    assert weights["overlap_weight"].gt(0).all()
    assert summary["positivity_pass"].all()
    assert balance["overlap_weighted_standardized_mean_difference"].abs().max() < 0.1
    assert (
        balance["overlap_weighted_standardized_mean_difference"].abs()
        < balance["unweighted_standardized_mean_difference"].abs()
    ).all()

    perturbed = scores.copy()
    perturbed["syndrome_score"] = -17.0 * perturbed["syndrome_score"] + 3.0
    repeated = run_climate_overlap_validation(
        perturbed, covariates, realms, pattern, branching, explanation
    )[0]
    keys = ["context_layer", "context_a", "context_b", "island_id"]
    observed = weights.sort_values(keys).reset_index(drop=True)
    expected = repeated.sort_values(keys).reset_index(drop=True)
    np.testing.assert_allclose(
        observed["overlap_weight"], expected["overlap_weight"], atol=0, rtol=0
    )


def test_distance_x_climate_absorbs_a_pure_climate_contingency() -> None:
    pattern, branching, explanation = _configs()
    scores, covariates, realms = _synthetic(residual_context_slope=0.0)
    vector_tests = run_climate_overlap_validation(
        scores, covariates, realms, pattern, branching, explanation
    )[4]
    confirmatory = vector_tests.loc[
        vector_tests["support_tier"].eq("confirmatory")
        & vector_tests["stratum"].eq("all_native")
    ]
    reference = confirmatory.loc[
        confirmatory["validation_mode"].eq("reference_climate_main_effects")
    ]
    adjusted = confirmatory.loc[
        confirmatory["validation_mode"].eq("distance_x_climate_adjusted")
    ]
    assert reference["v1_vector_difference_supported"].all()
    assert not adjusted["v1_vector_difference_supported"].any()


def test_residual_context_branch_survives_distance_x_climate_adjustment() -> None:
    pattern, branching, explanation = _configs()
    scores, covariates, realms = _synthetic(residual_context_slope=0.55, seed=1422)
    vector_tests = run_climate_overlap_validation(
        scores, covariates, realms, pattern, branching, explanation
    )[4]
    adjusted = vector_tests.loc[
        vector_tests["validation_mode"].eq("distance_x_climate_adjusted")
        & vector_tests["support_tier"].eq("confirmatory")
        & vector_tests["stratum"].eq("all_native")
    ]
    assert adjusted["v1_vector_difference_supported"].all()


def test_v1_does_not_lower_the_frozen_response_support_gate() -> None:
    pattern, branching, explanation = _configs()
    scores, covariates, realms = _synthetic(n_per_context=40)
    vector_tests = run_climate_overlap_validation(
        scores, covariates, realms, pattern, branching, explanation
    )[4]
    assert not vector_tests["status"].eq("fit").any()
    assert vector_tests["threshold"].eq(50).all()
    assert not vector_tests["v1_vector_difference_supported"].any()
