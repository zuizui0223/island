from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_syndrome_distance_sensitivity import run_distance_sensitivity


def test_raw_sqrt_log1p_keep_strong_synthetic_pattern():
    scores = []
    cov = []
    for ci, context in enumerate(["northern_midlatitude", "tropical"]):
        for i in range(72):
            island = f"{context}_{i}"
            distance = float(i * i + 1)
            x = np.log1p(distance)
            cov.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": f"b{ci}_{i // 6}",
                    "distance_to_continent_km": distance,
                    "log_island_area_km2": np.sin(i / 9),
                    "climate_pc1": np.cos(i / 8),
                }
            )
            for syndrome in ["large_bee_like", "generalized_accessible", "selfing_syndrome"]:
                if context == "northern_midlatitude":
                    slope = {
                        "large_bee_like": -0.18,
                        "generalized_accessible": 0.18,
                        "selfing_syndrome": 0.12,
                    }[syndrome]
                else:
                    slope = {
                        "large_bee_like": 0.13,
                        "generalized_accessible": -0.10,
                        "selfing_syndrome": 0.0,
                    }[syndrome]
                scores.append(
                    {
                        "island_id": island,
                        "syndrome": syndrome,
                        "syndrome_score": slope * x + 0.02 * np.sin(i * 1.7),
                        "n_species": 10,
                        "stratum": "all_native",
                    }
                )
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
    }
    syndrome = {
        "syndromes": {
            "large_bee_like": {"traits": {}},
            "generalized_accessible": {"traits": {}},
            "selfing_syndrome": {"traits": {}},
        }
    }
    _, _, _, headline = run_distance_sensitivity(
        pd.DataFrame(scores), pd.DataFrame(cov), pattern, syndrome
    )
    assert set(headline["distance_form"]) == {"raw", "sqrt", "log1p"}
    assert headline["north_attraction_direction_supported"].all()
    assert headline["north_tropical_vector_difference_supported"].all()
