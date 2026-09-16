from __future__ import annotations

import pandas as pd

from island_v2.chapter1_h5_glopl_floral_architecture_moderation import (
    aggregate_species_measurement_cells,
    build_trait_assignments,
    classify_family_result,
    classify_trait_result,
    normalize_species_name,
    summarize_trait_support,
)


def _config() -> dict:
    return {
        "trait_recodes": {
            "generalized_form": {
                "source_trait": "floral_form",
                "protected_states": ["open_radial", "brush_puff", "composite_head"],
                "unprotected_states": ["tubular", "bell_campanulate", "funnel_trumpet"],
                "excluded_states": "all_other_or_mixed_states",
            },
            "actinomorphic_symmetry": {
                "source_trait": "floral_symmetry",
                "protected_states": ["actinomorphic"],
                "unprotected_states": ["zygomorphic"],
                "excluded_states": ["asymmetric", "mixed_or_variable"],
            },
            "shallow_open_tube": {
                "source_trait": "tube_depth_class",
                "protected_states": ["absent_or_open", "shallow"],
                "unprotected_states": ["deep"],
                "excluded_states": ["intermediate", "mixed_or_variable"],
            },
        },
        "preflight": {
            "support_gates": {
                "global_min_matched_species_per_trait": 2,
                "global_min_publications_per_class": 2,
                "global_min_sites_per_class": 2,
                "global_min_offshore_sites_per_class": 1,
                "context_min_publications_per_class": 1,
                "context_min_sites_per_class": 1,
                "context_min_offshore_sites_per_class": 1,
            }
        },
        "analysis": {"publication_total_weight": 1.0},
    }


def test_normalize_species_name_is_exact_but_format_stable() -> None:
    assert normalize_species_name("  Abelia__chinensis ") == "abelia chinensis"
    assert normalize_species_name("ABELIA chinensis") == "abelia chinensis"
    assert normalize_species_name("") == ""


def test_build_trait_assignments_uses_only_predeclared_atomic_states() -> None:
    ledger = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b", "C c", "D d", "E e", "F f", "G g"],
            "axis": ["floral_structural_complexity"] * 7,
            "trait_name": [
                "floral_form",
                "floral_form",
                "floral_form",
                "floral_symmetry",
                "floral_symmetry",
                "tube_depth_class",
                "tube_depth_class",
            ],
            "normalized_value": [
                "open_radial",
                "tubular",
                "salverform",
                "actinomorphic",
                "asymmetric",
                "shallow",
                "intermediate",
            ],
            "resolution_status": ["resolved"] * 7,
            "quality": ["high", "medium", "high", "high", "medium", "high", "high"],
        }
    )
    out = build_trait_assignments(ledger, _config())
    got = {(r.species_key, r.trait, int(r.trait_state)) for r in out.itertuples()}
    assert got == {
        ("a a", "generalized_form", 1),
        ("b b", "generalized_form", 0),
        ("d d", "actinomorphic_symmetry", 1),
        ("f f", "shallow_open_tube", 1),
    }


def test_summarize_trait_support_requires_both_architecture_classes_and_offshore_support() -> None:
    rows = pd.DataFrame(
        {
            "trait": ["generalized_form"] * 8,
            "trait_state": [0, 0, 0, 0, 1, 1, 1, 1],
            "species_key": ["a", "b", "a", "b", "c", "d", "c", "d"],
            "site_key": ["s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"],
            "study_key": ["p1", "p2", "p1", "p2", "p3", "p4", "p3", "p4"],
            "analysis_regime": ["northern_midlatitude"] * 4 + ["tropical"] * 4,
            "distance_to_major_continent_km": [0, 10, 0, 20, 0, 30, 0, 40],
        }
    )
    summary = summarize_trait_support(rows, _config())
    trait = summary["generalized_form"]
    assert trait["global_admitted"] is True
    assert trait["classes"]["0"]["n_offshore_sites"] == 2
    assert trait["classes"]["1"]["n_offshore_sites"] == 2


def test_species_measurement_cells_keep_total_weight_one_per_publication() -> None:
    rows = pd.DataFrame(
        {
            "study_key": ["p1", "p1", "p1", "p1", "p2"],
            "site_key": ["s1", "s1", "s2", "s2", "s3"],
            "species_key": ["a", "a", "b", "c", "d"],
            "analysis_regime": ["northern_midlatitude"] * 5,
            "z_distance": [0.0, 0.0, 1.0, 1.0, 2.0],
            "trait": ["generalized_form"] * 5,
            "trait_state": [0, 0, 1, 1, 0],
            "PL_Effect_Size": [1.0, 3.0, 2.0, 4.0, 5.0],
            "PL_Effect_Size_Type1": ["A"] * 5,
            "PL_Effect_Size_Type2": ["Sup"] * 5,
            "Constant_added": ["false"] * 5,
            "Level_of_Supplementation": ["full"] * 5,
        }
    )
    out = aggregate_species_measurement_cells(rows, _config())
    sums = out.groupby("study_key")["analysis_weight"].sum().to_dict()
    assert sums == {"p1": 1.0, "p2": 1.0}
    p1_a = out.loc[(out.study_key == "p1") & (out.species_key == "a")].iloc[0]
    assert p1_a.PL_Effect_Size == 2.0


def test_trait_classification_requires_negative_buffering_interaction_and_positive_restricted_slope() -> None:
    primary = {
        "evaluable": True,
        "restricted_distance_slope": 0.30,
        "distance_by_architecture_interaction": -0.20,
        "interaction_one_sided_negative_p": 0.02,
    }
    sensitivities = {
        "supplemental_only": {"distance_by_architecture_interaction": -0.10},
        "no_zero_constant": {"distance_by_architecture_interaction": -0.05},
    }
    assert classify_trait_result(primary, sensitivities)["supported"] is True
    primary["restricted_distance_slope"] = -0.01
    assert classify_trait_result(primary, sensitivities)["supported"] is False


def test_route_b_family_requires_two_supported_atomic_architecture_traits() -> None:
    assert classify_family_result([True, True, False])["route_B_family_supported"] is True
    assert classify_family_result([True, False, False])["route_B_family_supported"] is False
