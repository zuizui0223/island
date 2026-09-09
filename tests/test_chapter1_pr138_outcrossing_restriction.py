import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_outcrossing_restriction import (
    analyse_restricted_attraction,
    reproductive_restriction_sets,
)


def test_reproductive_restrictions_require_exact_unambiguous_states():
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "SI exact",
                "trait_name": "self_incompatibility",
                "normalized_value": "SI",
            },
            {
                "accepted_species": "SI mixed",
                "trait_name": "self_incompatibility",
                "normalized_value": "SC|SI",
            },
            {
                "accepted_species": "Out exact",
                "trait_name": "mating_system",
                "normalized_value": "predominantly_outcrossing",
            },
            {
                "accepted_species": "Out mixed",
                "trait_name": "mating_system",
                "normalized_value": "mixed_mating|predominantly_outcrossing",
            },
            {
                "accepted_species": "Both exact",
                "trait_name": "self_incompatibility",
                "normalized_value": "SI",
            },
            {
                "accepted_species": "Both exact",
                "trait_name": "mating_system",
                "normalized_value": "predominantly_outcrossing",
            },
        ]
    )
    sets = reproductive_restriction_sets(ledger)
    assert sets["si_only"] == {"SI exact", "Both exact"}
    assert sets["predominantly_outcrossing"] == {"Out exact", "Both exact"}
    assert sets["si_and_predominantly_outcrossing"] == {"Both exact"}


def test_distance_effect_is_recovered_in_restricted_attraction_score():
    rng = np.random.default_rng(13827)
    score_rows = []
    cov_rows = []
    n = 140
    for i in range(n):
        island = f"i{i}"
        distance = -1.2 + 2.4 * i / (n - 1)
        attraction = 0.35 * distance + rng.normal(0, 0.04)
        cov_rows.append(
            {
                "island_id": island,
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"b{i // 5}",
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": np.sin(i / 11),
                "climate_pc1": np.cos(i / 13),
            }
        )
        for syndrome, score in [
            ("large_bee_like", -attraction),
            ("generalized_accessible", attraction),
        ]:
            score_rows.append(
                {
                    "restriction": "si_only",
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
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
    }
    result = analyse_restricted_attraction(
        pd.DataFrame(score_rows),
        pd.DataFrame(cov_rows),
        config,
        context_column="analysis_regime",
        contexts=["northern_midlatitude"],
        context_layer="analysis_regime",
        restriction="si_only",
    )
    fit = result.loc[result["status"].eq("fit")].iloc[0]
    assert fit["distance_estimate"] > 0
    assert fit["distance_p"] < 0.05
