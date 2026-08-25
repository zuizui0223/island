import math

import numpy as np
import pandas as pd

from island_v2.chapter1_when_where_block_deletion import run_block_deletion


def test_leave_one_block_headline_is_not_single_block_driven():
    rng = np.random.default_rng(53)
    comp = []
    cov = []
    idx = 0
    for context in ["northern_midlatitude", "tropical"]:
        for block_i in range(10):
            for within_block in range(6):
                island_id = f"i{idx}"
                x = -1.5 + 3.0 * (block_i * 6 + within_block) / 59
                cov.append(
                    {
                        "island_id": island_id,
                        "analysis_regime": context,
                        "spatial_block": f"{context[:2]}{block_i}",
                        "log_distance_to_continent_km": x,
                        "log_island_area_km2": 2.0 + rng.normal(0, 0.1),
                        "climate_pc1": rng.normal(),
                        "climate_pc2": rng.normal(),
                        "climate_pc3": rng.normal(),
                        "climate_pc4": rng.normal(),
                    }
                )
                for state, north_slope in [("bell", 1.5), ("salver", -1.4), ("red", -1.2)]:
                    slope = north_slope if context == "northern_midlatitude" else -0.45 * north_slope
                    p = 1.0 / (1.0 + math.exp(-(-0.2 + slope * x)))
                    trait = "flower_primary_color" if state == "red" else "floral_form"
                    comp.append(
                        {
                            "island_id": island_id,
                            "stratum": "all_native",
                            "trait_name": trait,
                            "trait_state": state,
                            "successes": int(rng.binomial(140, p)),
                            "trials": 140,
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
        "north_context": "northern_midlatitude",
        "tropical_context": "tropical",
        "confirmatory_min_islands": 50,
        "strata": ["all_native"],
    }
    detail, summary = run_block_deletion(pd.DataFrame(comp), pd.DataFrame(cov), config)
    assert detail["deleted_block"].nunique() == 20
    row = summary.iloc[0]
    assert row["n_headline_testable"] == 20
    assert row["headline_replication_fraction"] == 1.0
