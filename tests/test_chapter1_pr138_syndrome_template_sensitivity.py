from __future__ import annotations

from island_v2.chapter1_pr138_syndrome_template_sensitivity import build_template_variants


def _config():
    return {
        "score_definition": {"minimum_informative_traits": 2},
        "syndromes": {
            "large_bee_like": {
                "traits": {
                    "flower_primary_color": {"weight": 0.5},
                    "floral_form": {"weight": 1.0},
                    "floral_symmetry": {"weight": 0.75},
                    "tube_depth_class": {"weight": 1.0},
                    "flower_size_class": {"weight": 0.25},
                }
            },
            "generalized_accessible": {
                "traits": {
                    "floral_form": {"weight": 1.0},
                    "floral_symmetry": {"weight": 0.75},
                    "tube_depth_class": {"weight": 1.0},
                }
            },
            "selfing_syndrome": {
                "traits": {
                    "self_incompatibility": {"weight": 1.0},
                    "mating_system": {"weight": 1.0},
                    "autonomous_selfing_capacity": {"weight": 1.0},
                    "flower_size_class": {"weight": 0.5},
                }
            },
            "selfing_core": {
                "traits": {
                    "self_incompatibility": {"weight": 1.0},
                    "mating_system": {"weight": 1.0},
                    "autonomous_selfing_capacity": {"weight": 1.0},
                }
            },
        },
    }


def test_template_variants_are_outcome_blind_and_complete():
    config = _config()
    variants = build_template_variants(config)
    expected_traits = {
        "flower_primary_color",
        "floral_form",
        "floral_symmetry",
        "tube_depth_class",
        "flower_size_class",
        "self_incompatibility",
        "mating_system",
        "autonomous_selfing_capacity",
    }
    assert {"canonical", "equal_weights", "no_colour", "pollination_morphology_only", "minimum_three_traits"}.issubset(variants)
    assert {f"drop_{x}" for x in expected_traits}.issubset(variants)

    assert "flower_primary_color" not in variants["no_colour"]["syndromes"]["large_bee_like"]["traits"]
    assert variants["minimum_three_traits"]["score_definition"]["minimum_informative_traits"] == 3
    assert all(
        spec["weight"] == 1.0
        for syndrome in variants["equal_weights"]["syndromes"].values()
        for spec in syndrome["traits"].values()
    )
    assert set(variants["pollination_morphology_only"]["syndromes"]["large_bee_like"]["traits"]) == {
        "floral_form",
        "floral_symmetry",
        "tube_depth_class",
    }
    assert "flower_size_class" not in variants["drop_flower_size_class"]["syndromes"]["selfing_syndrome"]["traits"]

    # The input configuration itself must not be mutated.
    assert "flower_primary_color" in config["syndromes"]["large_bee_like"]["traits"]
    assert config["syndromes"]["large_bee_like"]["traits"]["flower_primary_color"]["weight"] == 0.5
