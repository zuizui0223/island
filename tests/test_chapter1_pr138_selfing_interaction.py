import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_selfing_interaction import run_interaction_analysis


def _synthetic(interaction_strength: float) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rng = np.random.default_rng(13826)
    score_rows = []
    cov_rows = []
    n = 180
    for i in range(n):
        island = f"i{i}"
        x = -1.5 + 3.0 * i / (n - 1)
        selfing = 0.45 * np.sin(i / 7.0) + 0.10 * x
        cluster = i // 6
        attraction = (
            0.30 * x
            + 0.12 * selfing
            + interaction_strength * x * selfing
            + rng.normal(0, 0.035)
        )
        cov_rows.append(
            {
                "island_id": island,
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{cluster}",
                "log_distance_to_continent_km": x,
                "log_island_area_km2": np.sin(i / 13.0),
                "climate_pc1": np.cos(i / 17.0),
            }
        )
        for syndrome, score in [
            ("large_bee_like", -attraction),
            ("generalized_accessible", attraction),
            ("selfing_core", selfing),
        ]:
            score_rows.append(
                {
                    "source_mode": "geo_k10",
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
    return pd.DataFrame(score_rows), pd.DataFrame(cov_rows), config


def test_positive_distance_by_selfing_interaction_is_recovered():
    scores, covariates, config = _synthetic(0.55)
    result = run_interaction_analysis(
        scores,
        covariates,
        config,
        context_column="analysis_regime",
        contexts=["northern_midlatitude"],
        context_layer="analysis_regime",
    )
    fit = result.loc[result["status"].eq("fit")].iloc[0]
    assert fit["distance_estimate"] > 0
    assert fit["distance_x_selfing_estimate"] > 0
    assert fit["distance_x_selfing_p"] < 0.05


def test_no_selfing_threshold_is_needed_for_null_interaction():
    scores, covariates, config = _synthetic(0.0)
    result = run_interaction_analysis(
        scores,
        covariates,
        config,
        context_column="analysis_regime",
        contexts=["northern_midlatitude"],
        context_layer="analysis_regime",
    )
    fit = result.loc[result["status"].eq("fit")].iloc[0]
    assert fit["distance_estimate"] > 0
    assert abs(fit["distance_x_selfing_estimate"]) < 0.08
