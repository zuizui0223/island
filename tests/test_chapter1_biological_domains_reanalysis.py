import math

import pandas as pd
import yaml

from island_v2.chapter1_biological_domains_reanalysis import (
    _bh,
    build_genus_residual_species_scores,
    build_island_scores,
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
primary_responses:
  - reproductive_assurance
  - accessibility_specialization
domains:
  reproductive_assurance:
    role: primary_directional
    minimum_components_per_species: 2
    components:
      self_incompatibility: {SI: 0.0, SC: 1.0}
      mating_system: {predominantly_outcrossing: 0.0, predominantly_selfing: 1.0}
  accessibility_specialization:
    role: primary_directional
    minimum_components_per_species: 2
    components:
      floral_form: {open_radial: 1.0, tubular: 0.0}
      floral_symmetry: {actinomorphic: 1.0, zygomorphic: 0.0}
  signal_display:
    role: exploratory_domain
    components:
      flower_size_reduction:
        source_trait: flower_size_class
        role: strict_directional_proxy
        values: {very_small: 1.0, small: 1.0, large: 0.0, very_large: 0.0}
    descriptive_state_traits:
      - flower_primary_color
      - flower_size_class
      - inflorescence_display
"""
    )


def test_signal_display_is_not_one_confirmatory_response():
    specs = response_specs(_config())
    assert "signal_display" not in specs
    assert specs["flower_size_reduction"]["role"] == "strict_directional_proxy"
    assert "plain_colour_proxy" not in specs
    assert "inflorescence_display_reduction_proxy" not in specs


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
