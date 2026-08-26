import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_biogeographic_pattern import (
    build_observed_broad_counts,
    run_observed_pattern,
)


def _config():
    return {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50, "pilot": 30},
        "minimum_outcomes_per_vector": 2,
        "alpha": 0.05,
        "broad_outcomes": {
            "plain_colour": {
                "trait_name": "flower_primary_color",
                "positive_states": ["white", "green_brown_inconspicuous"],
                "negative_states": ["yellow_orange", "red_pink", "blue_purple"],
                "northern_classic_direction": 1,
            },
            "generalized_form": {
                "trait_name": "floral_form",
                "positive_states": ["open_radial", "brush_puff", "composite_head"],
                "negative_states": ["tubular", "spurred"],
                "northern_classic_direction": 1,
            },
            "self_compatibility": {
                "trait_name": "self_incompatibility",
                "positive_states": ["SC"],
                "negative_states": ["SI"],
                "northern_classic_direction": 1,
            },
        },
    }


def _synthetic(*, north_slope=0.8, tropical_slope=0.0, n_per_context=80):
    rows = []
    cov = []
    rng = np.random.default_rng(11)
    for context, slope in [
        ("northern_midlatitude", north_slope),
        ("tropical", tropical_slope),
    ]:
        for i in range(n_per_context):
            island = f"{context}_{i}"
            distance = -1.5 + 3.0 * i / max(1, n_per_context - 1)
            area = rng.normal()
            climate = rng.normal()
            cov.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_b{i // 5}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": area,
                    "climate_pc1": climate,
                }
            )
            for j, outcome in enumerate(
                ["plain_colour", "generalized_form", "self_compatibility"]
            ):
                eta = -0.25 + slope * distance + 0.08 * area - 0.05 * climate + 0.05 * j
                p = 1.0 / (1.0 + np.exp(-eta))
                trials = 80
                successes = int(round(trials * p))
                rows.append(
                    {
                        "island_id": island,
                        "outcome": outcome,
                        "stratum": "all_native",
                        "successes": successes,
                        "trials": trials,
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(cov)


def test_observed_primary_detects_northern_classic_and_between_region_difference():
    counts, covariates = _synthetic(north_slope=0.9, tropical_slope=-0.05)
    _, _, within, between = run_observed_pattern(counts, covariates, _config())
    north = within.loc[
        (within["support_tier"] == "confirmatory")
        & (within["context"] == "northern_midlatitude")
    ].iloc[0]
    contrast = between.loc[
        (between["support_tier"] == "confirmatory")
        & (between["context_a"] == "northern_midlatitude")
        & (between["context_b"] == "tropical")
    ].iloc[0]
    assert bool(north["observed_vector_supported"])
    assert bool(north["classic_direction_supported"])
    assert north["classic_projection_estimate"] > 0
    assert bool(contrast["biogeographic_difference_supported"])


def test_same_slopes_do_not_manufacture_biogeographic_difference():
    counts, covariates = _synthetic(north_slope=0.55, tropical_slope=0.55)
    _, _, _, between = run_observed_pattern(counts, covariates, _config())
    contrast = between.loc[
        (between["support_tier"] == "confirmatory")
        & (between["context_a"] == "northern_midlatitude")
        & (between["context_b"] == "tropical")
    ].iloc[0]
    assert not bool(contrast["biogeographic_difference_supported"])


def test_under_supported_context_is_not_testable():
    counts, covariates = _synthetic(n_per_context=25)
    _, _, within, between = run_observed_pattern(counts, covariates, _config())
    confirmatory = within.loc[within["support_tier"] == "confirmatory"]
    assert set(confirmatory["status"]) == {"not_testable"}
    confirmatory_between = between.loc[between["support_tier"] == "confirmatory"]
    assert set(confirmatory_between["status"]) == {"not_testable"}


def test_broad_builder_fails_closed_for_crossing_multistate_signature():
    status = pd.DataFrame(
        {
            "island_id": ["i1"] * 4,
            "accepted_species": ["plain", "showy", "mixed", "sc"],
            "origin_status": ["native"] * 4,
            "endemic_status": ["nonendemic"] * 4,
            "floristic_status": ["native_nonendemic"] * 4,
        }
    )
    audit = pd.DataFrame(
        {
            "accepted_species": ["plain", "showy", "mixed", "sc"],
            "trait_name": [
                "flower_primary_color",
                "flower_primary_color",
                "flower_primary_color",
                "self_incompatibility",
            ],
            "resolved_for_primary": [True] * 4,
            "canonical_signature": ["white", "red_pink", "white|red_pink", "SC"],
        }
    )
    counts = build_observed_broad_counts(status, audit, _config())
    plain = counts.loc[
        (counts["outcome"] == "plain_colour") & (counts["stratum"] == "all_native")
    ].iloc[0]
    assert plain["trials"] == 2
    assert plain["successes"] == 1
