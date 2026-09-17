import numpy as np
import pandas as pd

from island_v2.chapter1_v13_raw_colour_audit import (
    build_raw_colour_counts,
    fit_raw_colour_models,
    parse_colour_states,
)


def test_parse_colour_states_preserves_multistate_categories():
    states = parse_colour_states(
        'flower_primary_color=["red_pink","white","yellow_orange"]'
    )
    assert states == {"red_pink", "white", "yellow_orange"}


def test_raw_colour_counts_keep_each_colour_as_an_independent_response():
    axis = pd.DataFrame(
        {
            "accepted_species": ["A a", "B b", "C c"],
            "axis": ["flower_colour"] * 3,
            "trait_composition": [
                'flower_primary_color=["red_pink","white"]',
                'flower_primary_color=["yellow_orange"]',
                'flower_primary_color=["other_described"]',
            ],
            "quality": ["high", "medium", "high"],
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
    counts = build_raw_colour_counts(
        axis,
        flora,
        evidence_scope="all",
        strata=["all_observed"],
        colours=["white", "red_pink", "yellow_orange", "blue_purple"],
    )
    got = counts.set_index("colour")
    assert got.loc["white", "trials"] == 2
    assert got.loc["red_pink", "trials"] == 2
    assert got.loc["yellow_orange", "trials"] == 2
    assert got.loc["white", "successes"] == 1
    assert got.loc["red_pink", "successes"] == 1
    assert got.loc["yellow_orange", "successes"] == 1
    assert got.loc["blue_purple", "successes"] == 0


def test_colour_distance_effect_survives_selfing_adjustment_in_synthetic_data():
    rng = np.random.default_rng(1307)
    rows = []
    scores = []
    covariates = []
    n = 160
    for i in range(n):
        island = f"i{i}"
        distance = -1.0 + 2.0 * i / (n - 1)
        selfing = 0.45 * distance + rng.normal(0, 0.35)
        trials = 80
        red_p = 1.0 / (1.0 + np.exp(-(-0.5 - 0.9 * distance + 0.05 * selfing)))
        white_p = 1.0 / (1.0 + np.exp(-(-0.2 + 0.05 * distance + 0.1 * selfing)))
        rows.extend(
            [
                {
                    "island_id": island,
                    "stratum": "all_observed",
                    "colour": "red_pink",
                    "successes": int(round(trials * red_p)),
                    "trials": trials,
                },
                {
                    "island_id": island,
                    "stratum": "all_observed",
                    "colour": "white",
                    "successes": int(round(trials * white_p)),
                    "trials": trials,
                },
            ]
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
        "colours": ["white", "red_pink"],
    }
    fitted = fit_raw_colour_models(
        pd.DataFrame(rows),
        pd.DataFrame(scores),
        pd.DataFrame(covariates),
        config,
    )
    conditional = fitted.loc[
        fitted["model"].eq("conditional_selfing")
        & fitted["colour"].eq("red_pink")
        & fitted["status"].eq("fit")
    ].iloc[0]
    assert conditional["distance_estimate"] < 0
    assert conditional["distance_p"] < 0.05
