import math

import numpy as np
import pandas as pd

from island_v2.chapter1_trait_probability import run_trait_probability


def test_trait_probability_recovers_positive_and_negative_context_slopes():
    rng = np.random.default_rng(47)
    composition = []
    covariates = []
    idx = 0
    for context, slope in [("northern_midlatitude", 1.0), ("tropical", -0.7)]:
        for j in range(65):
            island_id = f"p{idx}"
            distance = 0.2 + 2.8 * j / 64
            covariates.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{idx // 4}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": 2.0 + rng.normal(0, 0.2),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            p = 1.0 / (1.0 + math.exp(-(-0.4 + slope * (distance - 1.6))))
            trials = 100
            composition.append(
                {
                    "island_id": island_id,
                    "stratum": "all_native",
                    "trait_name": "self_incompatibility",
                    "trait_state": "SC",
                    "successes": int(rng.binomial(trials, p)),
                    "trials": trials,
                }
            )
            idx += 1

    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "distance_gradient": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "min_islands_per_curve": 30,
        "confirmatory_islands_per_curve": 50,
        "distance_grid_quantiles": [0.05, 0.5, 0.95],
    }
    coefficients, curves, fits = run_trait_probability(
        pd.DataFrame(composition), pd.DataFrame(covariates), config
    )

    assert len(fits.loc[fits["support_tier"].eq("confirmatory")]) == 2
    slopes = coefficients.loc[
        coefficients["predictor"].eq("z_log_distance_to_continent_km")
    ].set_index("context")["estimate_log_odds"]
    assert slopes["northern_midlatitude"] > 0
    assert slopes["tropical"] < 0

    north = curves.loc[curves["context"].eq("northern_midlatitude")].sort_values(
        "distance_quantile"
    )
    tropical = curves.loc[curves["context"].eq("tropical")].sort_values(
        "distance_quantile"
    )
    assert north["predicted_probability"].iloc[-1] > north["predicted_probability"].iloc[0]
    assert tropical["predicted_probability"].iloc[-1] < tropical["predicted_probability"].iloc[0]


def test_trait_probability_does_not_create_zero_rows_for_missing_islands():
    composition = pd.DataFrame(
        {
            "island_id": [f"i{x}" for x in range(30)],
            "stratum": ["all_native"] * 30,
            "trait_name": ["self_incompatibility"] * 30,
            "trait_state": ["SC"] * 30,
            "successes": [5] * 30,
            "trials": [10] * 30,
        }
    )
    covariates = pd.DataFrame(
        {
            "island_id": [f"i{x}" for x in range(40)],
            "analysis_regime": ["northern_midlatitude"] * 40,
            "spatial_block": [f"b{x // 3}" for x in range(40)],
            "log_distance_to_continent_km": np.linspace(0.1, 3.0, 40),
            "log_island_area_km2": np.linspace(1.0, 2.0, 40),
            "climate_pc1": np.linspace(-1.0, 1.0, 40),
            "climate_pc2": np.linspace(-0.8, 0.8, 40),
            "climate_pc3": np.linspace(-0.6, 0.6, 40),
            "climate_pc4": np.linspace(-0.4, 0.4, 40),
        }
    )
    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "distance_gradient": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "contexts": ["northern_midlatitude"],
        "strata": ["all_native"],
        "min_islands_per_curve": 30,
        "confirmatory_islands_per_curve": 50,
        "distance_grid_quantiles": [0.05, 0.95],
    }
    _, _, fits = run_trait_probability(composition, covariates, config)
    row = fits.iloc[0]
    assert row["n_unique_islands"] == 30
