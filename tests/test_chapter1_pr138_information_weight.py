import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_information_weight import run_pr138_information_weight


def test_information_weight_modes_preserve_testability_and_do_not_invent_difference():
    rows = []
    cov = []
    outcomes = ["a", "b", "c"]
    rng = np.random.default_rng(138)
    contexts = ["northern_midlatitude", "tropical"]
    n_per_context = 120
    for context in contexts:
        # Independent cluster noise in each region keeps the same true distance slope
        # while avoiding the degenerate zero-residual covariance of a deterministic null.
        cluster_noise = rng.normal(0.0, 0.08, size=n_per_context // 5)
        for i in range(n_per_context):
            island = f"{context}_{i}"
            x = -1.0 + 2.0 * i / (n_per_context - 1)
            block = i // 5
            cov.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"{context}_{block}",
                    "log_distance_to_continent_km": x,
                    "log_island_area_km2": np.sin(i / 7),
                    "climate_pc1": np.cos(i / 11),
                }
            )
            for j, outcome in enumerate(outcomes):
                linear = -0.3 + 0.4 * x + 0.02 * j + cluster_noise[block]
                p = 1 / (1 + np.exp(-linear))
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
        "contexts": contexts,
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
