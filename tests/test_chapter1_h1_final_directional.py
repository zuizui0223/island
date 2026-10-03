import math

import numpy as np
import pandas as pd

from island_v2.chapter1_h1_final_directional import (
    _wild_signflip_p,
    intersection_union_summary,
    meta_directional_summary,
    validate_score_weights,
)


def _config():
    return {
        "model_outcomes": [
            "self_compatibility",
            "selfing_mating_system",
            "autonomous_selfing",
            "plain_colour",
            "generalized_form",
            "actinomorphic_symmetry",
            "shallow_open_tube",
        ],
        "response_families": {
            "reproductive_assurance": [
                "self_compatibility",
                "selfing_mating_system",
                "autonomous_selfing",
            ],
            "colour_dulling": ["plain_colour"],
            "accessibility_generalization": [
                "generalized_form",
                "actinomorphic_symmetry",
                "shallow_open_tube",
            ],
        },
        "directional_score": {
            "outcome_weights": {
                "self_compatibility": 1 / 9,
                "selfing_mating_system": 1 / 9,
                "autonomous_selfing": 1 / 9,
                "plain_colour": 1 / 3,
                "generalized_form": 1 / 9,
                "actinomorphic_symmetry": 1 / 9,
                "shallow_open_tube": 1 / 9,
            }
        },
    }


def test_directional_weights_equalize_three_domains():
    cfg = _config()
    weights = validate_score_weights(cfg)
    assert math.isclose(sum(weights.values()), 1.0)
    for outcomes in cfg["response_families"].values():
        assert math.isclose(
            sum(weights[x] for x in outcomes),
            1 / 3,
        )


def test_intersection_union_uses_weakest_required_region():
    contexts = ["a", "b", "c", "d"]
    frame = pd.DataFrame(
        {
            "status": ["fit"] * 4,
            "context": contexts,
            "positive_direction": [True] * 4,
            "p_one_sided_t": [0.001, 0.02, 0.04, 0.03],
            "p_one_sided_wild": [0.002, 0.01, 0.049, 0.02],
        }
    )
    result = intersection_union_summary(
        frame,
        contexts=contexts,
        alpha=0.05,
    )
    assert math.isclose(result["iut_p_one_sided_t"], 0.04)
    assert math.isclose(result["iut_p_one_sided_wild"], 0.049)
    assert result["recurrent_robust_supported"] is True


def test_intersection_union_fails_when_one_region_is_negative():
    contexts = ["a", "b", "c", "d"]
    frame = pd.DataFrame(
        {
            "status": ["fit"] * 4,
            "context": contexts,
            "positive_direction": [True, True, False, True],
            "p_one_sided_t": [0.001, 0.002, 0.9, 0.003],
            "p_one_sided_wild": [0.001, 0.002, 0.8, 0.003],
        }
    )
    result = intersection_union_summary(
        frame,
        contexts=contexts,
        alpha=0.05,
    )
    assert result["all_context_estimates_positive"] is False
    assert result["recurrent_robust_supported"] is False


def test_wild_signflip_is_deterministic_and_directional():
    influences = np.array([0.08, -0.02, 0.05, -0.01, 0.03, 0.02])
    p1 = _wild_signflip_p(
        0.30,
        0.12,
        influences,
        replications=1999,
        seed=17,
    )
    p2 = _wild_signflip_p(
        0.30,
        0.12,
        influences,
        replications=1999,
        seed=17,
    )
    p_negative = _wild_signflip_p(
        -0.30,
        0.12,
        influences,
        replications=1999,
        seed=17,
    )
    assert p1 == p2
    assert p1 < 0.05
    assert p_negative > 0.5


def test_meta_summary_can_support_global_average_when_strict_recurrence_fails():
    contexts = ["north_mid", "north_high", "tropical", "south"]
    frame = pd.DataFrame(
        {
            "status": ["fit"] * 4,
            "context": contexts,
            "estimate": [0.016, 0.114, 0.098, 0.070],
            "cluster_robust_se": [0.0144, 0.0383, 0.0165, 0.0268],
            "positive_direction": [True] * 4,
            "p_one_sided_t": [0.13, 0.002, 1e-7, 0.006],
            "p_one_sided_wild": [0.14, 0.001, 0.001, 0.002],
        }
    )
    strict = intersection_union_summary(
        frame,
        contexts=contexts,
        alpha=0.05,
    )
    pooled = meta_directional_summary(
        frame,
        contexts=contexts,
        alpha=0.05,
    )
    assert strict["recurrent_robust_supported"] is False
    assert pooled["positive_global_average_supported"] is True
    assert pooled["i2"] > 0.5
