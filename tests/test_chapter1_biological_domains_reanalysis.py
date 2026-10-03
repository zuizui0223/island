import math

import pandas as pd
import yaml

from island_v2.chapter1_biological_domains_reanalysis import (
    _bh,
    build_genus_residual_species_scores,
    build_island_scores,
    build_species_response_scores,
    response_specs,
)


def _config():
    return yaml.safe_load(
        """
minimum_species_per_island_score: 1
common_support_responses:
  - reproductive_assurance
  - accessibility_specialization
  - flower_size_reduction
primary_posthoc_H1_responses:
  - reproductive_assurance
  - accessibility_specialization
domains:
  reproductive_assurance_core:
    domain: reproductive_assurance
    role: frozen_PR138_benchmark
    minimum_components_per_species: 2
    components:
      self_incompatibility:
        weight: 1.0
        values: {SI: 0.0, SC: 1.0}
      mating_system:
        weight: 1.0
        values: {predominantly_outcrossing: 0.0, predominantly_selfing: 1.0}
  reproductive_assurance:
    domain: reproductive_assurance
    role: posthoc_extended_domain
    minimum_components_per_species: 2
    components:
      self_incompatibility:
        weight: 1.0
        values: {SI: 0.0, SC: 1.0}
      mating_system:
        weight: 1.0
        values: {predominantly_outcrossing: 0.0, predominantly_selfing: 1.0}
      cleistogamy:
        weight: 1.0
        values: {absent: 0.0, facultative: 1.0, obligate: 1.0}
  accessibility_specialization:
    domain: accessibility_specialization
    role: frozen_PR138_benchmark_and_posthoc_domain
    minimum_components_per_species: 2
    components:
      floral_form:
        weight: 1.0
        values: {open_radial: 1.0, tubular: 0.0}
      floral_symmetry:
        weight: 0.75
        values: {actinomorphic: 1.0, zygomorphic: 0.0}
  signal_display:
    role: exploratory_domain
    directional_component:
      flower_size_reduction:
        source_trait: flower_size_class
        role: frozen_selfing_syndrome_size_recode
        minimum_components_per_species: 1
        weight: 1.0
        values: {very_small: 1.0, small: 1.0, large: 0.0, very_large: 0.0}
    raw_state_traits:
      - flower_primary_color
      - flower_size_class
      - inflorescence_display
"""
    )


def _raw_axis_config():
    return {
        "axes": {
            "flower_colour": {
                "traits": ["flower_primary_color"],
            },
            "floral_structural_complexity": {
                "traits": [
                    "floral_form",
                    "floral_symmetry",
                    "flower_size_class",
                    "inflorescence_display",
                ],
            },
            "reproductive_assurance": {
                "traits": [
                    "self_incompatibility",
                    "mating_system",
                    "cleistogamy",
                ],
            },
        },
        "evidence_scopes": {
            "direct_only": ["high", "medium"],
        },
    }


def _ontology():
    return {
        "traits": {
            "flower_primary_color": {
                "allowed_values": ["red_pink", "white", "unresolved"],
            },
            "floral_form": {
                "allowed_values": ["open_radial", "tubular", "unresolved"],
            },
            "floral_symmetry": {
                "allowed_values": ["actinomorphic", "zygomorphic", "unresolved"],
            },
            "flower_size_class": {
                "allowed_values": ["small", "large", "unresolved"],
            },
            "inflorescence_display": {
                "allowed_values": ["solitary", "umbel_corymb", "unresolved"],
            },
            "self_incompatibility": {
                "allowed_values": ["SC", "SI", "unresolved"],
            },
            "mating_system": {
                "allowed_values": [
                    "predominantly_selfing",
                    "predominantly_outcrossing",
                    "unresolved",
                ],
            },
            "cleistogamy": {
                "allowed_values": ["absent", "facultative", "unresolved"],
            },
        }
    }


def test_signal_display_is_not_one_confirmatory_response():
    specs = response_specs(_config())
    assert "signal_display" not in specs
    assert specs["flower_size_reduction"]["role"] == (
        "frozen_selfing_syndrome_size_recode"
    )
    assert "plain_colour_proxy" not in specs
    assert "inflorescence_display_reduction_proxy" not in specs


def test_accessibility_preserves_frozen_trait_weights():
    specs = response_specs(_config())
    assert specs["accessibility_specialization"]["components"][
        "floral_form"
    ]["weight"] == 1.0
    assert specs["accessibility_specialization"]["components"][
        "floral_symmetry"
    ]["weight"] == 0.75


def test_ambiguous_mixed_states_are_missing_not_midpoint():
    species_axis = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "axis": "floral_structural_complexity",
                "trait_composition": (
                    'floral_form=["open_radial","tubular"]|'
                    'floral_symmetry=["actinomorphic"]'
                ),
                "quality": "high",
            },
            {
                "accepted_species": "Beta two",
                "axis": "floral_structural_complexity",
                "trait_composition": (
                    'floral_form=["open_radial"]|'
                    'floral_symmetry=["actinomorphic"]'
                ),
                "quality": "high",
            },
        ]
    )
    scores, _ = build_species_response_scores(
        species_axis,
        _ontology(),
        _raw_axis_config(),
        _config(),
        evidence_scope="direct_only",
    )
    access = scores.loc[
        scores["response"].eq("accessibility_specialization")
    ].set_index("accepted_species")
    assert "Alpha one" not in access.index
    assert math.isclose(
        float(access.loc["Beta two", "score"]),
        1.0,
    )


def test_common_support_uses_only_species_with_all_three_strict_scores():
    cfg = _config()
    flora = pd.DataFrame(
        {
            "island_id": ["1", "1", "2", "2"],
            "accepted_species": ["A one", "B two", "A one", "C three"],
            "floristic_status": [""] * 4,
            "origin_status": [""] * 4,
        }
    )
    rows = []
    for response, value in (
        ("reproductive_assurance", 0.8),
        ("accessibility_specialization", 0.7),
        ("flower_size_reduction", 0.6),
    ):
        rows.append(
            {
                "accepted_species": "A one",
                "response": response,
                "domain": "x",
                "role": "x",
                "score": value,
            }
        )
    rows.extend(
        [
            {
                "accepted_species": "B two",
                "response": "reproductive_assurance",
                "domain": "x",
                "role": "x",
                "score": 0.1,
            },
            {
                "accepted_species": "C three",
                "response": "flower_size_reduction",
                "domain": "x",
                "role": "x",
                "score": 0.2,
            },
        ]
    )
    scores = pd.DataFrame(rows)
    out = build_island_scores(
        flora,
        scores,
        cfg,
        flora_scope="all_observed",
        support_mode="common_species",
    )
    assert set(out["response"]) == {
        "reproductive_assurance",
        "accessibility_specialization",
        "flower_size_reduction",
    }
    assert set(out["n_scored_species"]) == {1}
    assert out["island_id"].nunique() == 2


def test_genus_residual_is_leave_one_species_out():
    scores = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "response": "r",
                "domain": "d",
                "role": "x",
                "score": 0.2,
            },
            {
                "accepted_species": "Alpha two",
                "response": "r",
                "domain": "d",
                "role": "x",
                "score": 0.8,
            },
            {
                "accepted_species": "Beta one",
                "response": "r",
                "domain": "d",
                "role": "x",
                "score": 0.5,
            },
        ]
    )
    out = build_genus_residual_species_scores(scores)
    alpha = out.set_index("accepted_species")["score"]
    assert math.isclose(alpha["Alpha one"], -0.6)
    assert math.isclose(alpha["Alpha two"], 0.6)
    assert "Beta one" not in set(out["accepted_species"])


def test_bh_is_monotone_and_bounded():
    q = _bh(pd.Series([0.001, 0.02, 0.2, 0.8]))
    assert q.between(0, 1).all()
    assert list(q.sort_values()) == sorted(q.tolist())


def test_repository_config_reuses_frozen_pr138_directional_states():
    config = yaml.safe_load(
        open("config/chapter1_biological_domains_reanalysis.yml", encoding="utf-8")
    )
    frozen = yaml.safe_load(
        open("config/chapter1_pr138_pollination_syndromes.yml", encoding="utf-8")
    )
    selfing = frozen["syndromes"]["selfing_core"]["traits"]
    assurance = config["domains"]["reproductive_assurance"]["components"]
    for trait, spec in selfing.items():
        mapping = assurance[trait]["values"]
        assert set(mapping) == set(spec["preferred"]) | set(spec["opposed"])
        assert all(mapping[state] == 1.0 for state in spec["preferred"])
        assert all(mapping[state] == 0.0 for state in spec["opposed"])

    access = frozen["syndromes"]["generalized_accessible"]["traits"]
    observed = config["domains"]["accessibility_specialization"]["components"]
    for trait, spec in access.items():
        mapping = observed[trait]["values"]
        assert set(mapping) == set(spec["preferred"]) | set(spec["opposed"])
        assert all(mapping[state] == 1.0 for state in spec["preferred"])
        assert all(mapping[state] == 0.0 for state in spec["opposed"])

    frozen_size = frozen["syndromes"]["selfing_syndrome"]["traits"][
        "flower_size_class"
    ]
    observed_size = config["domains"]["signal_display"]["directional_component"][
        "flower_size_reduction"
    ]["values"]
    assert set(observed_size) == set(frozen_size["preferred"]) | set(
        frozen_size["opposed"]
    )
    assert all(
        observed_size[state] == 1.0 for state in frozen_size["preferred"]
    )
    assert all(
        observed_size[state] == 0.0 for state in frozen_size["opposed"]
    )


def test_repository_config_does_not_order_colour_or_inflorescence():
    config = yaml.safe_load(
        open("config/chapter1_biological_domains_reanalysis.yml", encoding="utf-8")
    )
    components = config["domains"]["signal_display"]["directional_component"]
    assert set(components) == {"flower_size_reduction"}
    descriptive = set(
        config["domains"]["signal_display"]["raw_state_traits"]
    )
    assert "flower_primary_color" in descriptive
    assert "inflorescence_display" in descriptive


def test_pairwise_common_support_does_not_require_third_domain():
    cfg = _config()
    cfg["common_support_sets"] = {
        "common_assurance_accessibility": [
            "reproductive_assurance",
            "accessibility_specialization",
        ]
    }
    flora = pd.DataFrame(
        {
            "island_id": ["1", "1"],
            "accepted_species": ["A one", "B two"],
            "floristic_status": ["", ""],
            "origin_status": ["", ""],
        }
    )
    scores = pd.DataFrame(
        [
            {
                "accepted_species": "A one",
                "response": "reproductive_assurance",
                "domain": "x",
                "role": "x",
                "score": 0.8,
            },
            {
                "accepted_species": "A one",
                "response": "accessibility_specialization",
                "domain": "x",
                "role": "x",
                "score": 0.7,
            },
            {
                "accepted_species": "B two",
                "response": "reproductive_assurance",
                "domain": "x",
                "role": "x",
                "score": 0.2,
            },
            {
                "accepted_species": "B two",
                "response": "accessibility_specialization",
                "domain": "x",
                "role": "x",
                "score": 0.3,
            },
        ]
    )
    out = build_island_scores(
        flora,
        scores,
        cfg,
        flora_scope="all_observed",
        support_mode="common_assurance_accessibility",
    )
    assert set(out["response"]) == {
        "reproductive_assurance",
        "accessibility_specialization",
    }
    assert set(out["n_scored_species"]) == {2}
