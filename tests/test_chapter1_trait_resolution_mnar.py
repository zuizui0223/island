from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from island_v2.chapter1_trait_resolution_mnar import (
    adjust_island_scores,
    build_tipping_summary,
    build_total_species_counts,
    completed_prevalence_from_resolution_or,
    run_mnar_sensitivity,
)


def test_selection_odds_model_reproduces_baseline_and_is_monotone():
    baseline = completed_prevalence_from_resolution_or(
        total_species=100,
        resolved_species=50,
        resolved_focal_count=25,
        resolution_odds_ratio=1,
    )
    focal_overresolved = completed_prevalence_from_resolution_or(
        total_species=100,
        resolved_species=50,
        resolved_focal_count=25,
        resolution_odds_ratio=3,
    )
    focal_underresolved = completed_prevalence_from_resolution_or(
        total_species=100,
        resolved_species=50,
        resolved_focal_count=25,
        resolution_odds_ratio=1 / 3,
    )
    assert baseline == pytest.approx(0.5)
    assert focal_overresolved < baseline < focal_underresolved
    assert focal_overresolved == pytest.approx(0.375)
    assert focal_underresolved == pytest.approx(0.625)


def _status_flora(islands: list[str], n_species: int = 20) -> pd.DataFrame:
    rows = []
    for island in islands:
        for i in range(n_species):
            rows.append(
                {
                    "island_id": island,
                    "accepted_species": f"{island}_species_{i}",
                    "origin_status": "native",
                    "endemic_status": "nonendemic",
                    "floristic_status": "native_nonendemic",
                }
            )
    return pd.DataFrame(rows)


def test_adjustment_preserves_observed_information_and_unaffected_axis():
    scores = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "stratum": "all_native",
                "syndrome": "selfing_core",
                "syndrome_score": 0.0,
                "n_species": 10,
            },
            {
                "island_id": "i1",
                "stratum": "all_native",
                "syndrome": "generalized_accessible",
                "syndrome_score": 0.4,
                "n_species": 12,
            },
        ]
    )
    flora = _status_flora(["i1"])
    totals = build_total_species_counts(flora, ["all_native"])
    assignments = pd.DataFrame(
        {
            "island_id": ["i1"],
            "analysis_regime": ["northern_midlatitude"],
            "biogeographic_realm": ["Palearctic"],
        }
    )
    scenario = {
        "scenario_id": "shared_or_3",
        "scenario_family": "shared",
        "scenario_type": "selection_grid",
        "context_layer": "all_islands",
        "context_a": "",
        "context_b": "",
        "resolution_odds_ratio": 3.0,
        "log_resolution_odds_ratio": math.log(3),
        "bound_mode": "",
    }
    adjusted, diagnostics = adjust_island_scores(
        scores, totals, assignments, scenario, {"selfing_core"}
    )
    indexed = adjusted.set_index("syndrome")
    assert indexed.loc["selfing_core", "syndrome_score"] < 0
    assert indexed.loc["selfing_core", "n_species"] == 10
    assert indexed.loc["generalized_accessible", "syndrome_score"] == pytest.approx(0.4)
    assert diagnostics.iloc[0]["n_unresolved_species"] == 10


def test_context_contrast_bound_leaves_other_contexts_at_observed_score():
    scores = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "stratum": "all_native",
                "syndrome": "selfing_core",
                "syndrome_score": 0.4,
                "n_species": 10,
            }
        ]
    )
    totals = build_total_species_counts(_status_flora(["i1"]), ["all_native"])
    assignments = pd.DataFrame(
        {
            "island_id": ["i1"],
            "analysis_regime": ["southern_extratropical"],
            "biogeographic_realm": ["Australasian"],
        }
    )
    scenario = {
        "scenario_id": "contrast_bound",
        "scenario_family": "north_tropical",
        "scenario_type": "partial_identification_bound",
        "context_layer": "analysis_regime",
        "context_a": "northern_midlatitude",
        "context_b": "tropical",
        "resolution_odds_ratio": float("nan"),
        "log_resolution_odds_ratio": float("nan"),
        "bound_mode": "first_context_nonfocal_second_context_focal",
    }
    adjusted, _ = adjust_island_scores(
        scores, totals, assignments, scenario, {"selfing_core"}
    )
    assert adjusted.iloc[0]["syndrome_score"] == pytest.approx(0.4)


def _pattern_config() -> dict:
    return {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 10},
    }


def _branching_config() -> dict:
    return {
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


def _explanation_config() -> dict:
    return {
        "contract": "chapter1_explanation_gap_validation_v1",
        "validations": {
            "V5_trait_resolution_MNAR_tipping_point": {
                "freeze_before_execution": {
                    "resolution_odds_ratio_grid": [0.5, 1.0, 2.0],
                    "affected_syndrome_axes": ["selfing_core"],
                    "scenario_families": {
                        "shared": {
                            "context_layer": "all_islands",
                            "contrast_contexts": [],
                        }
                    },
                }
            }
        },
    }


def _synthetic():
    rng = np.random.default_rng(51)
    scores = []
    covariates = []
    assignments = []
    islands = []
    for context, realm, multiplier in (
        ("northern_midlatitude", "Palearctic", 1.0),
        ("tropical", "Neotropical", -0.5),
    ):
        for i in range(20):
            island = f"{context}_{i}"
            islands.append(island)
            distance = -1 + 2 * i / 19
            noise = rng.normal(0, 0.005)
            scores.extend(
                [
                    {
                        "island_id": island,
                        "stratum": "all_native",
                        "syndrome": "generalized_accessible",
                        "syndrome_score": 0.25 * distance + noise,
                        "n_species": 12,
                    },
                    {
                        "island_id": island,
                        "stratum": "all_native",
                        "syndrome": "selfing_core",
                        "syndrome_score": multiplier * 0.3 * distance + noise,
                        "n_species": 8 + i % 5,
                    },
                ]
            )
            covariates.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_{i // 2}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": math.sin(i / 4),
                    "climate_pc1": math.cos(i / 5),
                }
            )
            assignments.append(
                {"island_id": island, "biogeographic_realm": realm}
            )
    return (
        pd.DataFrame(scores),
        _status_flora(islands),
        pd.DataFrame(covariates),
        pd.DataFrame(assignments),
    )


def test_end_to_end_grid_reproduces_or1_and_materializes_tipping_outputs():
    scores, flora, covariates, assignments = _synthetic()
    result = run_mnar_sensitivity(
        island_scores=scores,
        status_flora=flora,
        covariates=covariates,
        realm_assignment=assignments,
        pattern_config=_pattern_config(),
        branching_config=_branching_config(),
        explanation_config=_explanation_config(),
        evidence_scope="synthetic",
    )
    manifest = result["manifest"]
    assert manifest["maximum_OR1_score_difference"] <= 1e-12
    assert manifest["n_scenarios"] == 5
    tipping = result["tipping"]
    assert not tipping.empty
    assert set(tipping["odds_ratio_direction"]) == {"below_one", "above_one"}
    assert set(result["score_diagnostics"]["scenario_type"]) == {
        "selection_grid",
        "partial_identification_bound",
    }


def test_tipping_summary_records_support_loss_without_selecting_gain():
    common = {
        "evidence_scope": "all",
        "scenario_family": "shared",
        "scenario_type": "selection_grid",
        "context_layer": "analysis_regime",
        "axis_set": "universal_plant_response",
        "stratum": "all_native",
        "support_tier": "confirmatory",
        "context": "tropical",
        "syndrome": "reproductive_assurance",
    }
    slopes = pd.DataFrame(
        [
            {
                **common,
                "resolution_odds_ratio": odds_ratio,
                "log_resolution_odds_ratio": math.log(odds_ratio),
                "distance_slope": estimate,
                "q_axis_family": q_value,
                "axis_supported": supported,
            }
            for odds_ratio, estimate, q_value, supported in (
                (0.5, -0.1, 0.2, False),
                (1.0, 0.2, 0.01, True),
                (2.0, 0.1, 0.2, False),
            )
        ]
    )
    tipping = build_tipping_summary(slopes, pd.DataFrame(), pd.DataFrame())
    assert len(tipping) == 2
    assert set(tipping["verdict"]) == {"break_detected"}
    below = tipping.loc[tipping["odds_ratio_direction"].eq("below_one")].iloc[0]
    assert below["tipping_event"] == "support_lost|sign_changed"
