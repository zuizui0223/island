import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_palearctic_restricted_block_deletion import (
    FULL_SENTINEL,
    run_palearctic_restricted_block_deletion,
)


def test_restricted_block_deletion_is_outcome_blind_and_preserves_strong_signal():
    rng = np.random.default_rng(13821)
    n = 180
    covariates = []
    realms = []
    scores = []
    blocks = ["lat12_lon20", *[f"b{i}" for i in range(1, 15)]]

    for i in range(n):
        island_id = f"i{i}"
        distance = -1.5 + 3.0 * i / (n - 1)
        block = blocks[i % len(blocks)]
        covariates.append(
            {
                "island_id": island_id,
                "spatial_block": block,
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": rng.normal(),
                "climate_pc1": rng.normal(),
            }
        )
        realms.append(
            {"island_id": island_id, "biogeographic_realm": "Palearctic"}
        )
        attraction = 0.42 * distance + rng.normal(0, 0.045)
        for source_mode, offset in [("geo_k5", 0.0), ("geo_k10", 0.015)]:
            for syndrome, value in [
                ("large_bee_like", -attraction + offset),
                ("generalized_accessible", attraction + offset),
            ]:
                scores.append(
                    {
                        "restriction": "si_only",
                        "source_mode": source_mode,
                        "island_id": island_id,
                        "stratum": "all_native",
                        "syndrome": syndrome,
                        "syndrome_score": value,
                    }
                )

    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "support_tiers": {"confirmatory": 30},
        "alpha": 0.05,
    }
    outputs = run_palearctic_restricted_block_deletion(
        pd.DataFrame(scores),
        pd.DataFrame(covariates),
        pd.DataFrame(realms),
        pattern,
    )

    manifest = outputs["manifest"]
    assert manifest["outcome_values_used_to_select_blocks"] is False
    assert manifest["effect_estimates_used_to_select_blocks"] is False
    assert manifest["p_values_used_to_select_blocks"] is False
    assert manifest["n_palearctic_blocks"] == len(blocks)

    result = outputs["results"]
    assert FULL_SENTINEL in set(result["deleted_block"])
    assert set(blocks).issubset(set(result["deleted_block"]))
    fitted_deletions = result.loc[
        ~result["deleted_block"].eq(FULL_SENTINEL) & result["status"].eq("fit")
    ]
    assert not fitted_deletions.empty
    assert (fitted_deletions["distance_estimate"] > 0).all()

    influence = outputs["influence"]
    aegean = influence.loc[influence["deleted_block"].eq("lat12_lon20")]
    assert not aegean.empty
    assert aegean["is_preidentified_aegean"].all()
    assert np.isfinite(aegean["absolute_deletion_delta"]).all()

    ranking = outputs["block_influence"]
    assert len(ranking) == len(blocks)
    assert set(ranking["overall_influence_rank"]) == set(range(1, len(blocks) + 1))
    assert ranking.loc[
        ranking["deleted_block"].eq("lat12_lon20"), "is_preidentified_aegean"
    ].all()

    summary = outputs["summary"]
    assert len(summary) == 1
    assert summary.iloc[0]["restriction"] == "si_only"
    assert summary.iloc[0]["stratum"] == "all_native"
    assert int(summary.iloc[0]["n_positive"]) == int(
        summary.iloc[0]["n_fitted_block_source_scenarios"]
    )
