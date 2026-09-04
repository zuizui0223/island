import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_pathway_decomposition import run_pathway_decomposition


def test_attraction_distance_effect_can_persist_after_selfing_core_adjustment():
    rng = np.random.default_rng(1381)
    score_rows = []
    cov_rows = []
    n = 120
    cluster_self = rng.normal(0, 0.05, size=n // 5)
    cluster_attr = rng.normal(0, 0.05, size=n // 5)
    for i in range(n):
        island = f"i{i}"
        x = -1 + 2 * i / (n - 1)
        block = i // 5
        selfing_core = 0.20 * x + cluster_self[block]
        attraction_shift = 0.28 * x + 0.12 * selfing_core + cluster_attr[block]
        large_bee = -attraction_shift
        generalized = attraction_shift
        cov_rows.append(
            {
                "island_id": island,
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{block}",
                "log_distance_to_continent_km": x,
                "log_island_area_km2": np.sin(i / 9),
                "climate_pc1": np.cos(i / 11),
            }
        )
        for syndrome, score in [
            ("large_bee_like", large_bee),
            ("generalized_accessible", generalized),
            ("selfing_core", selfing_core),
        ]:
            score_rows.append(
                {
                    "island_id": island,
                    "stratum": "all_native",
                    "syndrome": syndrome,
                    "syndrome_score": score,
                }
            )
    config = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
    }
    result = run_pathway_decomposition(
        pd.DataFrame(score_rows), pd.DataFrame(cov_rows), config
    )
    fit = result.loc[result["status"].eq("fit")].set_index("model")
    assert fit.loc["selfing_core_distance", "distance_estimate"] > 0
    assert fit.loc["attraction_unconditional", "distance_estimate"] > 0
    assert fit.loc[
        "attraction_conditional_on_selfing_core", "distance_estimate"
    ] > 0
    assert fit.loc[
        "attraction_conditional_on_selfing_core", "distance_p"
    ] < 0.05
