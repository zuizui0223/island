from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_syndrome_block_deletion import run_syndrome_block_deletion


def _synthetic():
    rows = []
    cov = []
    syndromes = ["large_bee_like", "generalized_accessible", "selfing_syndrome"]
    for context_index, context in enumerate(["northern_midlatitude", "tropical"]):
        for i in range(72):
            island = f"{context}_{i}"
            x = -1.0 + 2.0 * i / 71
            block = f"{context_index}_{i // 6}"
            block_noise = 0.015 * np.sin((i // 6) * 1.3)
            cov.append(
                {
                    "island_id": island,
                    "analysis_regime": context,
                    "spatial_block": block,
                    "log_distance_to_continent_km": x,
                    "log_island_area_km2": np.sin(i / 9),
                    "climate_pc1": np.cos(i / 8),
                }
            )
            for syndrome in syndromes:
                if context == "northern_midlatitude":
                    slope = {
                        "large_bee_like": -0.22,
                        "generalized_accessible": 0.22,
                        "selfing_syndrome": 0.16,
                    }[syndrome]
                else:
                    slope = {
                        "large_bee_like": 0.18,
                        "generalized_accessible": -0.12,
                        "selfing_syndrome": 0.0,
                    }[syndrome]
                rows.append(
                    {
                        "island_id": island,
                        "syndrome": syndrome,
                        "syndrome_score": slope * x + block_noise + 0.01 * np.sin(i * 1.7),
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
        "syndromes": {x: {"traits": {}} for x in syndromes},
        "score_definition": {"minimum_informative_traits": 2},
    }
    return pd.DataFrame(rows), pd.DataFrame(cov), pattern, syndrome


def test_leave_one_block_retains_strong_synthetic_headline():
    scores, cov, pattern, syndrome = _synthetic()
    detail, summary = run_syndrome_block_deletion(scores, cov, pattern, syndrome)
    assert not detail.empty
    assert summary.iloc[0]["n_deleted_blocks"] == 24
    assert summary.iloc[0]["n_headline_testable"] == 24
    assert summary.iloc[0]["headline_replication_fraction"] == 1.0
    assert summary.iloc[0]["north_attraction_direction_fraction"] == 1.0
