import math

import numpy as np
import pandas as pd

from island_v2.chapter1_when_where_domain_sensitivity import run_domain_sensitivity


def test_domain_sensitivity_keeps_domains_separate():
    rng = np.random.default_rng(31)
    comp = []
    cov = []
    idx = 0
    for context in ["northern_midlatitude", "tropical"]:
        for j in range(60):
            island_id = f"d{idx}"
            isolation = -1.5 + 3.0 * j / 59
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
            for trait, states in {
                "floral_form": [("bell", 1.0), ("salver", -0.9)],
                "flower_primary_color": [("red", -0.8), ("white", 0.7)],
            }.items():
                for state, north_slope in states:
                    slope = north_slope if context == "northern_midlatitude" else 0.0
                    p = 1.0 / (1.0 + math.exp(-(-0.3 + slope * isolation)))
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
    within, between, _ = run_domain_sensitivity(pd.DataFrame(comp), pd.DataFrame(cov), config)
    assert set(within["trait_domain"]) == {"floral_form", "flower_primary_color"}
    north = within.loc[within["context"].eq("northern_midlatitude")]
    assert north["where_supported"].all()
    assert set(between["trait_domain"]) == {"floral_form", "flower_primary_color"}
