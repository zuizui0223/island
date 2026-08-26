import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_information_weight import run_pr138_information_weight


def test_information_weight_modes_preserve_testability_and_do_not_invent_difference():
    rows = []
    cov = []
    outcomes = ["a", "b", "c"]
    for context in ["northern_midlatitude", "tropical"]:
        for i in range(70):
            island = f"{context}_{i}"
            x = -1.0 + 2.0 * i / 69
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
            for j, outcome in enumerate(outcomes):
                p = 1 / (1 + np.exp(-(-0.3 + 0.4 * x + 0.02 * j)))
                trials = 10 + (i % 7) * 20
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
        "broad_outcomes": {
            x: {"northern_classic_direction": 1} for x in outcomes
        },
    }
    _, between, summary = run_pr138_information_weight(
        pd.DataFrame(rows), pd.DataFrame(cov), config, ["canonical", "equal_island"]
    )
    assert set(summary["information_weight_mode"]) == {"canonical", "equal_island"}
    assert summary["north_testable"].all()
    assert not summary["north_tropical_difference_supported"].any()
    assert set(between["information_weight_mode"]) == {"canonical", "equal_island"}
