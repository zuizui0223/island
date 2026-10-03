import math
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_biological_domains_reanalysis import (
    _bh,
    build_genus_residual_species_scores,
    build_island_scores,
    h3_summary,
    response_specs,
)


def _config():
    return yaml.safe_load(
        Path("config/chapter1_biological_domains_reanalysis.yml").read_text(
            encoding="utf-8"
        )
    )


def test_domain_contract_has_core_extended_accessibility_and_size_only():
    specs = response_specs(_config())
    assert set(specs) == {
        "reproductive_assurance_core",
        "reproductive_assurance",
        "accessibility_specialization",
        "flower_size_reduction",
    }
    assert specs["reproductive_assurance_core"]["domain"] == (
        "reproductive_assurance"
    )
    assert "cleistogamy" not in specs[
        "reproductive_assurance_core"
    ]["components"]
    assert "cleistogamy" in specs["reproductive_assurance"]["components"]
    assert set(specs["flower_size_reduction"]["components"]) == {
        "flower_size_class"
    }


def test_accessibility_reuses_frozen_pr138_weights_and_states():
    specs = response_specs(_config())
    access = specs["accessibility_specialization"]["components"]
    assert access["floral_form"]["weight"] == 1.0
    assert access["floral_symmetry"]["weight"] == 0.75
    assert access["tube_depth_class"]["weight"] == 1.0
    assert "intermediate" not in access["tube_depth_class"]["values"]
    assert "bell_campanulate" not in access["floral_form"]["values"]
    assert "funnel_trumpet" not in access["floral_form"]["values"]


def test_signal_display_has_no_single_ordinal_domain_score():
    cfg = _config()
    signal = cfg["domains"]["signal_display"]
    assert signal["directional_identification"] == "partial_only"
    assert set(signal["raw_state_traits"]) == {
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    }
    specs = response_specs(cfg)
    assert "signal_display" not in specs
    assert "plain_colour_proxy" not in specs
    assert "inflorescence_display_reduction_proxy" not in specs


def test_common_support_uses_same_species_for_all_three_directional_scores():
    cfg = _config()
    cfg["minimum_species_per_island_score"] = 1
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


def test_contract_marks_reclassification_posthoc_and_blocks_hard_attraction_order():
    cfg = _config()
    assert (
        cfg["inferential_role"]
        == "posthoc_biological_reclassification_stress_test"
    )
    guards = set(cfg["reviewer_guards"])
    assert "colour_hue_is_not_equated_with_visual_attraction_intensity" in guards
    assert (
        "inflorescence_architecture_is_not_equated_with_total_display_area"
        in guards
    )
    assert (
        "posthoc_reclassification_never_replaces_final_confirmatory_H1"
        in guards
    )
