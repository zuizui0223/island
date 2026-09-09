import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_syndrome_analysis import (
    _prepare,
    _trait_concordance,
    _within_context,
    score_species_syndromes,
)


def test_mixed_preferred_and_opposed_state_is_unresolved():
    score = _trait_concordance(
        "red_pink|green_brown_inconspicuous",
        {"red_pink"},
        {"green_brown_inconspicuous"},
    )
    assert np.isnan(score)


def test_species_syndrome_requires_multivariate_support():
    ledger = pd.DataFrame(
        [
            {"accepted_species": "A one", "trait_name": "flower_primary_color", "normalized_value": "red_pink"},
            {"accepted_species": "A one", "trait_name": "tube_depth_class", "normalized_value": "deep"},
            {"accepted_species": "B two", "trait_name": "flower_primary_color", "normalized_value": "red_pink"},
        ]
    )
    config = {
        "score_definition": {"minimum_informative_traits": 2},
        "syndromes": {
            "bird_like": {
                "traits": {
                    "flower_primary_color": {"weight": 1, "preferred": ["red_pink"], "opposed": ["green_brown_inconspicuous"]},
                    "tube_depth_class": {"weight": 1, "preferred": ["deep"], "opposed": ["absent_or_open"]},
                }
            }
        },
    }
    scored = score_species_syndromes(ledger, config)
    assert scored["accepted_species"].tolist() == ["A one"]
    assert scored.iloc[0]["syndrome_concordance"] == 1.0


def test_northern_classic_projection_detects_large_bee_loss_plus_generalist_selfing_gain():
    rows = []
    cov = []
    syndromes = ["large_bee_like", "generalized_accessible", "selfing_syndrome"]
    for i in range(90):
        island = f"n{i}"
        distance = -1.0 + 2.0 * i / 89
        block = i // 5
        block_noise = 0.015 * np.sin(block * 1.3)
        cov.append(
            {
                "island_id": island,
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{block}",
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": np.sin(i / 7),
                "climate_pc1": np.cos(i / 9),
            }
        )
        for syndrome in syndromes:
            slope = -0.22 if syndrome == "large_bee_like" else 0.22
            rows.append(
                {
                    "island_id": island,
                    "syndrome": syndrome,
                    "syndrome_score": slope * distance + block_noise,
                    "n_species": 10,
                    "stratum": "all_native",
                }
            )
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
    }
    data = _prepare(pd.DataFrame(rows), pd.DataFrame(cov), pattern)
    slopes, result = _within_context(
        data,
        stratum="all_native",
        context_value="northern_midlatitude",
        support_tier="confirmatory",
        threshold=50,
        pattern_config=pattern,
        syndrome_config={},
    )
    assert result["status"] == "fit"
    assert result["northern_classic_projection"] > 0
    assert result["northern_classic_projection_one_sided_p"] < 0.01
    direction = slopes.set_index("syndrome")["distance_slope"]
    assert direction["large_bee_like"] < 0
    assert direction["generalized_accessible"] > 0
    assert direction["selfing_syndrome"] > 0
