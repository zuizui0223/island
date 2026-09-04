import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_equal_island_support import run_equal_island_support


def test_equal_island_support_thresholds_fail_closed_when_support_drops():
    rows = []
    cov = []
    for context in ["northern_midlatitude", "tropical"]:
        for i in range(60):
            island = f"{context}_{i}"
            x = -1 + 2 * i / 59
            cov.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_{i // 5}",
                    "log_distance_to_continent_km": x,
                    "log_island_area_km2": np.sin(i),
                    "climate_pc1": np.cos(i),
                }
            )
            for outcome in ["a", "b", "c"]:
                p = 1 / (1 + np.exp(-0.7 * x))
                trials = 6 if i < 50 else 2
                rows.append(
                    {
                        "island_id": island,
                        "outcome": outcome,
                        "stratum": "all_native",
                        "successes": p * trials,
                        "trials": trials,
                    }
                )
    config = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
        "minimum_outcomes_per_vector": 2,
        "alpha": 0.05,
        "broad_outcomes": {x: {"northern_classic_direction": 1} for x in ["a", "b", "c"]},
    }
    _, _, within, _, summary = run_equal_island_support(
        pd.DataFrame(rows), pd.DataFrame(cov), config, [1, 5, 10]
    )
    assert set(summary["minimum_direct_species_per_island_axis"]) == {1, 5, 10}
    assert summary.loc[summary["minimum_direct_species_per_island_axis"].eq(1), "north_status"].eq("fit").all()
    assert summary.loc[summary["minimum_direct_species_per_island_axis"].eq(5), "north_status"].eq("fit").all()
    assert summary.loc[summary["minimum_direct_species_per_island_axis"].eq(10), "north_status"].eq("not_testable").all()
    assert not within.loc[within["minimum_direct_species_per_island_axis"].eq(10)].empty
