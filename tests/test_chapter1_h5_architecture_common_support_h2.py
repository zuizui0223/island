from __future__ import annotations

import pandas as pd
import yaml

from island_v2.chapter1_h5_architecture_common_support_h2 import (
    complete_template_species,
    focused_atomic_config,
    integrate_decisions,
    restrict_common_support,
)


def _config() -> dict:
    with open(
        "config/chapter1_h5_architecture_common_support_h2_v1.yml",
        encoding="utf-8",
    ) as handle:
        return yaml.safe_load(handle)


def test_complete_template_species_requires_all_three_named_templates() -> None:
    scores = pd.DataFrame(
        [
            {"accepted_species": "A a", "syndrome": "large_bee_like", "syndrome_concordance": 0.1},
            {"accepted_species": "A a", "syndrome": "butterfly_like", "syndrome_concordance": 0.2},
            {"accepted_species": "A a", "syndrome": "bird_like", "syndrome_concordance": 0.3},
            {"accepted_species": "B b", "syndrome": "large_bee_like", "syndrome_concordance": 0.1},
            {"accepted_species": "B b", "syndrome": "butterfly_like", "syndrome_concordance": 0.2},
        ]
    )
    assert complete_template_species(scores, _config()) == {"A a"}


def test_restrict_common_support_uses_both_species_and_island_intersection() -> None:
    flora = pd.DataFrame(
        [
            {"island_id": "i1", "accepted_species": "A a", "origin_status": "native", "endemic_status": "nonendemic", "floristic_status": "native"},
            {"island_id": "i1", "accepted_species": "B b", "origin_status": "native", "endemic_status": "nonendemic", "floristic_status": "native"},
            {"island_id": "i2", "accepted_species": "A a", "origin_status": "native", "endemic_status": "nonendemic", "floristic_status": "native"},
        ]
    )
    audit = pd.DataFrame(
        [
            {"accepted_species": "A a", "trait_name": "floral_form", "resolved_for_primary": True, "canonical_signature": "open_radial"},
            {"accepted_species": "B b", "trait_name": "floral_form", "resolved_for_primary": True, "canonical_signature": "open_radial"},
        ]
    )
    covariates = pd.DataFrame(
        [
            {"island_id": "i1", "analysis_regime": "northern_midlatitude"},
            {"island_id": "i2", "analysis_regime": "tropical"},
        ]
    )
    strict_support = pd.DataFrame(
        [
            {"evidence_scope": "all_analysis_eligible", "island_id": "i1", "syndrome": "shared_architecture_factor", "n_species": 50, "stratum": "all_observed"},
        ]
    )
    scores = pd.DataFrame(
        [
            {"accepted_species": "A a", "syndrome": name, "syndrome_concordance": 0.1}
            for name in ("large_bee_like", "butterfly_like", "bird_like")
        ]
    )

    out_flora, out_audit, out_cov = restrict_common_support(
        flora,
        audit,
        covariates,
        strict_support,
        scores,
        config=_config(),
        evidence_scope="all_analysis_eligible",
    )
    assert set(out_flora["island_id"]) == {"i1"}
    assert set(out_flora["accepted_species"]) == {"A a"}
    assert set(out_audit["accepted_species"]) == {"A a"}
    assert set(out_cov["island_id"]) == {"i1"}


def test_focused_atomic_config_keeps_exact_six_outcomes_and_two_contexts() -> None:
    base = {
        "contexts": ["northern_midlatitude", "northern_high_latitude", "tropical", "southern_extratropical"],
        "strata": ["all_observed", "all_native"],
        "support_tiers": {"confirmatory": 50, "pilot": 30},
        "minimum_outcomes_per_vector": 2,
        "broad_outcomes": {name: {"trait_name": name, "positive_states": ["x"], "negative_states": ["y"], "northern_classic_direction": 1} for name in [
            "plain_colour", "generalized_form", "actinomorphic_symmetry", "shallow_open_tube", "small_flower", "self_compatibility", "selfing_mating_system", "autonomous_selfing"
        ]},
    }
    focused = focused_atomic_config(base, _config())
    assert focused["contexts"] == ["northern_midlatitude", "tropical"]
    assert focused["strata"] == ["all_observed"]
    assert focused["support_tiers"] == {"confirmatory": 50}
    assert focused["minimum_outcomes_per_vector"] == 6
    assert list(focused["broad_outcomes"]) == _config()["atomic_outcomes"]


def test_integration_requires_all_six_and_supported_h2_in_both_scopes() -> None:
    all_result = {
        "evidence_scope": "all_analysis_eligible",
        "evaluable": True,
        "n_retained_outcomes": 6,
        "joint_p": 0.01,
    }
    direct = {
        "evidence_scope": "direct_only",
        "evaluable": True,
        "n_retained_outcomes": 6,
        "joint_p": 0.20,
    }
    decision = integrate_decisions(all_result, direct, _config())
    assert decision["classification"] == "common_support_H2_not_robust_across_evidence_scopes"
    assert not decision["architecture_null_informative_about_H2_representation"]

    direct["joint_p"] = 0.03
    decision = integrate_decisions(all_result, direct, _config())
    assert decision["classification"] == "common_support_H2_retained_architecture_null_informative"
    assert decision["architecture_null_informative_about_H2_representation"]
