from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_joint_observation_bias import (
    _shared_v5_scenario,
    _summarize_envelope,
    _summarize_primary,
    load_config,
)
from island_v2.chapter1_species_detection_tipping import adjust_species_detection
from island_v2.chapter1_trait_resolution_mnar import adjust_island_scores


def test_frozen_joint_surface_sizes() -> None:
    config = load_config(Path("config/chapter1_joint_observation_bias.yml"))
    v5 = config["parent_contracts"]["V5_trait_resolution"]
    v6 = config["parent_contracts"]["V6_species_detection"]
    n = (
        len(v5["resolution_odds_ratio_grid"])
        * len(v6["median_distance_completeness_grid"])
        * len(v6["distance_completeness_odds_ratio_grid"])
        * len(v6["state_recording_odds_ratio_grid"])
    )
    assert n == config["joint_primary_surface"]["n_surfaces_per_evidence_scope"] == 1575

    env = config["partial_identification_envelope"]
    n_bounds = sum(len(item["modes"]) for item in env["trait_bound_scenarios"])
    corners = env["v6_corner_grid"]
    n_corners = (
        len(corners["median_distance_completeness"])
        * len(corners["distance_completeness_odds_ratio"])
        * len(corners["state_recording_odds_ratio"])
    )
    assert n_bounds * n_corners == env["n_corner_bound_surfaces_per_evidence_scope"] == 48


def _synthetic_scores() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "stratum": ["all_native", "all_native"],
            "syndrome": ["selfing_core", "generalized_accessible"],
            "syndrome_score": [0.0, 0.2],
            "n_species": [5.0, 5.0],
        }
    )


def test_joint_no_differential_bias_is_identity_and_keeps_information_weight() -> None:
    scores = _synthetic_scores()
    totals = pd.DataFrame(
        {
            "island_id": ["i1"],
            "stratum": ["all_native"],
            "n_total_stratum_species": [10.0],
        }
    )
    assignments = pd.DataFrame(
        {
            "island_id": ["i1"],
            "analysis_regime": ["northern_midlatitude"],
            "biogeographic_realm": ["Palearctic"],
        }
    )
    after_v5, _ = adjust_island_scores(
        scores,
        totals,
        assignments,
        _shared_v5_scenario(1.0),
        {"selfing_core"},
    )
    recorded = totals.rename(columns={"n_total_stratum_species": "n_recorded_stratum_species"})
    cov = pd.DataFrame({"island_id": ["i1"], "log_distance_to_continent_km": [1.0]})
    adjusted, _ = adjust_species_detection(
        after_v5,
        recorded,
        cov,
        geography="log_distance_to_continent_km",
        distance_mean=0.0,
        distance_sd=1.0,
        median_completeness=0.5,
        distance_completeness_odds_ratio=0.25,
        state_recording_odds_ratio=1.0,
        affected_syndrome="generalized_accessible",
    )
    np.testing.assert_allclose(adjusted["syndrome_score"], scores["syndrome_score"], atol=1e-12)
    np.testing.assert_allclose(adjusted["n_species"], scores["n_species"], atol=0.0)


def test_joint_processes_act_on_distinct_components() -> None:
    scores = _synthetic_scores()
    totals = pd.DataFrame(
        {
            "island_id": ["i1"],
            "stratum": ["all_native"],
            "n_total_stratum_species": [10.0],
        }
    )
    assignments = pd.DataFrame(
        {
            "island_id": ["i1"],
            "analysis_regime": ["northern_midlatitude"],
            "biogeographic_realm": ["Palearctic"],
        }
    )
    after_v5, _ = adjust_island_scores(
        scores,
        totals,
        assignments,
        _shared_v5_scenario(2.0),
        {"selfing_core"},
    )
    recorded = totals.rename(columns={"n_total_stratum_species": "n_recorded_stratum_species"})
    cov = pd.DataFrame({"island_id": ["i1"], "log_distance_to_continent_km": [1.0]})
    adjusted, _ = adjust_species_detection(
        after_v5,
        recorded,
        cov,
        geography="log_distance_to_continent_km",
        distance_mean=0.0,
        distance_sd=1.0,
        median_completeness=0.5,
        distance_completeness_odds_ratio=0.5,
        state_recording_odds_ratio=0.5,
        affected_syndrome="generalized_accessible",
    )
    original = scores.set_index("syndrome")["syndrome_score"]
    current = adjusted.set_index("syndrome")["syndrome_score"]
    assert not np.isclose(current["selfing_core"], original["selfing_core"])
    assert not np.isclose(current["generalized_accessible"], original["generalized_accessible"])
    np.testing.assert_allclose(adjusted["n_species"], scores["n_species"], atol=0.0)


def test_primary_and_envelope_classification_are_conservative() -> None:
    surface = pd.DataFrame(
        {
            "evidence_scope": ["direct_only"] * 2,
            "target": ["Palearctic_accessibility"] * 2,
            "stratum": ["all_native"] * 2,
            "target_type": ["scalar"] * 2,
            "expected_direction": ["positive"] * 2,
            "status": ["fit", "fit"],
            "robust_cell": [True, False],
            "estimate": [0.2, -0.1],
            "supported": [True, False],
        }
    )
    summary = _summarize_primary(surface).iloc[0]
    assert summary["n_robust_cells"] == 1
    assert not bool(summary["all_fit_cells_robust"])
    assert bool(summary["any_fit_cell_fragile"])

    envelope = _summarize_envelope(surface).iloc[0]
    assert envelope["estimate_lower"] == -0.1
    assert envelope["estimate_upper"] == 0.2
    assert not bool(envelope["expected_sign_identified"])
    assert not bool(envelope["support_identified_across_envelope"])
    assert not bool(envelope["envelope_robust"])


def test_parent_grids_remain_exactly_inherited() -> None:
    joint = load_config(Path("config/chapter1_joint_observation_bias.yml"))
    explanation = yaml.safe_load(Path("config/chapter1_explanation_gap_validation.yml").read_text())
    v6 = yaml.safe_load(Path("config/chapter1_species_detection_tipping.yml").read_text())
    assert joint["parent_contracts"]["V5_trait_resolution"]["resolution_odds_ratio_grid"] == explanation[
        "validations"
    ]["V5_trait_resolution_MNAR_tipping_point"]["freeze_before_execution"][
        "resolution_odds_ratio_grid"
    ]
    assert joint["parent_contracts"]["V6_species_detection"]["median_distance_completeness_grid"] == v6[
        "primary_bias_direction"
    ]["median_distance_completeness_grid"]
    assert joint["parent_contracts"]["V6_species_detection"][
        "distance_completeness_odds_ratio_grid"
    ] == v6["primary_bias_direction"]["distance_completeness_odds_ratio_grid"]
    assert joint["parent_contracts"]["V6_species_detection"]["state_recording_odds_ratio_grid"] == v6[
        "primary_bias_direction"
    ]["positive_state_recording_odds_ratio_grid"]
