import math

import numpy as np
import pandas as pd

from island_v2.chapter1_when_where_omnibus import run_when_where_omnibus


def _synthetic_when_where():
    rng = np.random.default_rng(17)
    composition = []
    covariates = []
    responses = [
        ("floral_form", "bell_campanulate", 1.1),
        ("floral_form", "papilionaceous", -1.0),
        ("flower_primary_color", "red_pink", -0.9),
    ]
    idx = 0
    for context in ["northern_midlatitude", "tropical"]:
        for j in range(60):
            island_id = f"i{idx}"
            isolation = -1.5 + 3.0 * j / 59
            covariates.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{idx // 3}",
                    "log_distance_to_continent_km": isolation,
                    "log_island_area_km2": 2.0 + rng.normal(0, 0.25),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            for stratum in ["all_native", "native_nonendemic"]:
                for trait, state, north_slope in responses:
                    slope = north_slope if context == "northern_midlatitude" else 0.0
                    eta = -0.4 + slope * isolation
                    p = 1.0 / (1.0 + math.exp(-eta))
                    trials = 100
                    successes = int(rng.binomial(trials, p))
                    composition.append(
                        {
                            "island_id": island_id,
                            "stratum": stratum,
                            "trait_name": trait,
                            "trait_state": state,
                            "successes": successes,
                            "trials": trials,
                        }
                    )
            idx += 1
    return pd.DataFrame(composition), pd.DataFrame(covariates)


def test_when_where_omnibus_recovers_northern_vector_and_between_context_difference():
    composition, covariates = _synthetic_when_where()
    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "isolation_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native", "native_nonendemic"],
        "support_tiers": {"confirmatory": 50},
    }
    within, between, classification = run_when_where_omnibus(
        composition, covariates, config
    )
    north = within.loc[
        within["stratum"].eq("all_native")
        & within["context"].eq("northern_midlatitude")
    ].iloc[0]
    assert bool(north["where_supported"])
    contrast = between.loc[between["stratum"].eq("all_native")].iloc[0]
    assert bool(contrast["between_where_supported"])
    north_class = classification.loc[
        classification["context"].eq("northern_midlatitude")
    ].iloc[0]
    assert north_class["when_class"] == "persists_in_native_nonendemic"
