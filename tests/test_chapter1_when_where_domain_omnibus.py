import math

import numpy as np
import pandas as pd

from island_v2.chapter1_when_where_domain_omnibus import run_domain_omnibus


def test_domain_omnibus_recovers_lower_dimensional_when_where_signals():
    rng = np.random.default_rng(41)
    comp = []
    cov = []
    idx = 0
    for context in ["northern_midlatitude", "tropical"]:
        for j in range(65):
            island_id = f"d{idx}"
            isolation = -1.5 + 3.0 * j / 64
            cov.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{idx // 3}",
                    "log_distance_to_continent_km": isolation,
                    "log_island_area_km2": 2.0 + rng.normal(0, 0.2),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            states = {
                "floral_form": [("bell", 1.0), ("salver", -0.9)],
                "flower_primary_color": [("red", -0.8), ("white", 0.7)],
            }
            for trait, values in states.items():
                for state, north_slope in values:
                    slope = north_slope if context == "northern_midlatitude" else -0.35 * north_slope
                    p = 1.0 / (1.0 + math.exp(-(-0.2 + slope * isolation)))
                    comp.append(
                        {
                            "island_id": island_id,
                            "stratum": "all_native",
                            "trait_name": trait,
                            "trait_state": state,
                            "successes": int(rng.binomial(100, p)),
                            "trials": 100,
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
        "isolation_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
    }
    within, between = run_domain_omnibus(pd.DataFrame(comp), pd.DataFrame(cov), config)
    assert set(within["trait_domain"]) == {"floral_form", "flower_primary_color"}
    assert within["domain_where_supported"].all()
    assert set(between["trait_domain"]) == {"floral_form", "flower_primary_color"}
    assert between["domain_between_supported"].all()
