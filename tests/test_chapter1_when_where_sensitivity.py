import math

import numpy as np
import pandas as pd

from island_v2.chapter1_when_where_sensitivity import run_when_where_sensitivity


def _synthetic():
    rng = np.random.default_rng(29)
    composition = []
    covariates = []
    responses = [
        ("floral_form", "bell_campanulate", 1.0),
        ("floral_form", "papilionaceous", -0.9),
        ("flower_primary_color", "red_pink", -0.8),
    ]
    idx = 0
    for context, latitude in [("northern_midlatitude", 42.0), ("tropical", 8.0)]:
        for j in range(65):
            island_id = f"i{idx}"
            distance = 60.0 + 12.0 * j
            log_distance = math.log1p(distance)
            covariates.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "island_latitude": latitude,
                    "spatial_block": f"b{idx // 3}",
                    "distance_to_continent_km": distance,
                    "log_distance_to_continent_km": log_distance,
                    "log_island_area_km2": 2.5 + rng.normal(0, 0.2),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            z = (j - 32) / 18
            for stratum in ["all_native", "native_nonendemic"]:
                for trait, state, north_slope in responses:
                    slope = north_slope if context == "northern_midlatitude" else -0.2 * north_slope
                    eta = -0.3 + slope * z
                    p = 1.0 / (1.0 + math.exp(-eta))
                    trials = 120
                    composition.append(
                        {
                            "island_id": island_id,
                            "stratum": stratum,
                            "trait_name": trait,
                            "trait_state": state,
                            "successes": int(rng.binomial(trials, p)),
                            "trials": trials,
                        }
                    )
            idx += 1
    return pd.DataFrame(composition), pd.DataFrame(covariates)


def test_when_where_sensitivity_preserves_headline_across_predeclared_scenarios():
    composition, covariates = _synthetic()
    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "cluster_column": "spatial_block",
        "canonical_context_column": "analysis_regime",
        "canonical_isolation_column": "log_distance_to_continent_km",
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native", "native_nonendemic"],
        "support_tiers": {"confirmatory": 50},
        "core_latitude_rule": {
            "tropical_abs_max": 20.0,
            "north_mid_min": 30.0,
            "north_mid_max": 60.0,
            "north_high_min": 70.0,
            "south_mid_min_abs": 30.0,
            "south_mid_max_abs": 60.0,
        },
        "scenarios": [
            {"name": "canonical", "context_rule": "canonical", "isolation_form": "log1p", "distance_filter": "all"},
            {"name": "core", "context_rule": "core_latitudes", "isolation_form": "log1p", "distance_filter": "all"},
            {"name": "sqrt", "context_rule": "canonical", "isolation_form": "sqrt", "distance_filter": "all"},
            {"name": "positive", "context_rule": "canonical", "isolation_form": "log1p", "distance_filter": "positive_only"},
            {"name": "oceanic50", "context_rule": "canonical", "isolation_form": "log1p", "distance_filter": "oceanic_50km"},
        ],
    }
    _, _, _, summary = run_when_where_sensitivity(composition, covariates, config)
    assert len(summary) == 10
    assert summary["headline_replication"].all()
    assert summary["north_n_responses"].ge(3).all()
    assert summary["tropical_n_responses"].ge(3).all()
    assert summary["common_n_responses"].ge(3).all()
