import numpy as np
import pandas as pd

from island_v2.chapter1_v13_colour_architecture_audit import (
    build_colour_architecture_counts,
    fit_colour_architecture_models,
    parse_axis_traits,
)


def test_parse_axis_traits_preserves_raw_multistate_values():
    parsed = parse_axis_traits(
        'floral_form=["tubular","funnel_trumpet"]|tube_depth_class=["deep"]'
    )
    assert parsed["floral_form"] == {"tubular", "funnel_trumpet"}
    assert parsed["tube_depth_class"] == {"deep"}


def test_raw_colour_architecture_counts_require_both_colour_and_form():
    axis = pd.DataFrame(
        {
            "accepted_species": ["A a", "A a", "B b", "B b", "C c", "C c"],
            "axis": [
                "flower_colour",
                "floral_structural_complexity",
                "flower_colour",
                "floral_structural_complexity",
                "flower_colour",
                "floral_structural_complexity",
            ],
            "trait_composition": [
                'flower_primary_color=["red_pink"]',
                'floral_form=["tubular"]',
                'flower_primary_color=["yellow_orange"]',
                'floral_form=["bell_campanulate"]',
                'flower_primary_color=["blue_purple"]',
                'floral_symmetry=["actinomorphic"]',
            ],
            "quality": ["high"] * 6,
        }
    )
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i1"],
            "accepted_species": ["A a", "B b", "C c"],
            "origin_status": ["native"] * 3,
            "endemic_status": ["nonendemic"] * 3,
            "floristic_status": ["native_nonendemic"] * 3,
        }
    )
    counts = build_colour_architecture_counts(
        axis,
        flora,
        evidence_scope="all",
        strata=["all_observed"],
    )
    got = counts.set_index("combination")
    # C c lacks floral-form evidence and therefore is outside this pairwise denominator.
    assert got.loc["red_pink__butterfly_form", "trials"] == 2
    assert got.loc["red_pink__butterfly_form", "successes"] == 1
    assert got.loc["yellow_orange__bird_form", "successes"] == 1
    assert got.loc["blue_purple__large_bee_form", "successes"] == 0


def test_colour_architecture_distance_effect_can_persist_after_selfing_adjustment():
    rng = np.random.default_rng(1311)
    rows = []
    scores = []
    covariates = []
    n = 160
    for i in range(n):
        island = f"i{i}"
        distance = -1.0 + 2.0 * i / (n - 1)
        selfing = 0.35 * distance + rng.normal(0, 0.35)
        trials = 70
        p = 1.0 / (1.0 + np.exp(-(-1.0 - 0.8 * distance + 0.05 * selfing)))
        rows.append(
            {
                "island_id": island,
                "stratum": "all_observed",
                "combination": "red_pink__butterfly_form",
                "successes": int(round(trials * p)),
                "trials": trials,
            }
        )
        scores.append(
            {
                "island_id": island,
                "stratum": "all_observed",
                "syndrome": "selfing_core",
                "syndrome_score": selfing,
            }
        )
        covariates.append(
            {
                "island_id": island,
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{i // 4}",
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": np.sin(i / 11),
                "climate_pc1": np.cos(i / 13),
            }
        )
    config = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude"],
        "strata": ["all_observed"],
        "support_tiers": {"confirmatory": 50},
    }
    result = fit_colour_architecture_models(
        pd.DataFrame(rows),
        pd.DataFrame(scores),
        pd.DataFrame(covariates),
        config,
    )
    row = result.loc[
        result["model"].eq("conditional_selfing")
        & result["combination"].eq("red_pink__butterfly_form")
        & result["status"].eq("fit")
    ].iloc[0]
    assert row["distance_estimate"] < 0
    assert row["distance_p"] < 0.05
