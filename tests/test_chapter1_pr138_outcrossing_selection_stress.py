import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_outcrossing_selection_stress import (
    build_restricted_attraction_scores,
    run_stress,
)


def test_build_restricted_attraction_scores_uses_predeclared_contrast():
    rows = []
    for syndrome, score in [
        ("large_bee_like", -0.4),
        ("generalized_accessible", 0.2),
    ]:
        rows.append(
            {
                "restriction": "si_only",
                "source_mode": "geo_k5",
                "island_id": "i1",
                "stratum": "all_native",
                "syndrome": syndrome,
                "syndrome_score": score,
            }
        )
    out = build_restricted_attraction_scores(pd.DataFrame(rows))
    assert len(out) == 1
    assert out.iloc[0]["syndrome"] == "restricted_attraction"
    assert np.isclose(out.iloc[0]["syndrome_score"], 0.3)


def test_joint_stress_refits_after_preidentified_block_deletion():
    rng = np.random.default_rng(1388)
    cov_rows = []
    realm_rows = []
    score_rows = []
    n = 180
    for i in range(n):
        island = f"i{i}"
        x = -1.0 + 2.0 * i / (n - 1)
        palearctic = i < n // 2
        context = "northern_midlatitude" if palearctic else "tropical"
        realm = "Palearctic" if palearctic else "Neotropical"
        block = "lat12_lon20" if i % 17 == 0 else f"b{i // 5}"
        cov_rows.append(
            {
                "island_id": island,
                "analysis_regime": context,
                "spatial_block": block,
                "log_distance_to_continent_km": x,
                "log_island_area_km2": np.sin(i / 13),
                "climate_pc1": np.cos(i / 17),
            }
        )
        realm_rows.append(
            {"island_id": island, "biogeographic_realm": realm}
        )
        # Outcome-blind missingness: independent random thinning, enough observed and unobserved.
        if rng.random() < 0.68:
            attraction = (0.35 if palearctic else -0.30) * x + rng.normal(0, 0.08)
            for syndrome, value in [
                ("large_bee_like", -attraction),
                ("generalized_accessible", attraction),
            ]:
                score_rows.append(
                    {
                        "restriction": "si_only",
                        "source_mode": "geo_k5",
                        "island_id": island,
                        "stratum": "all_native",
                        "syndrome": syndrome,
                        "syndrome_score": value,
                    }
                )

    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "contexts": ["northern_midlatitude", "tropical"],
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 20},
    }
    selection = {
        "primary_axis_set": "universal_plant_response",
        "selection_model": {"minimum_observed": 20, "minimum_unobserved": 20},
        "weight_modes": {
            "unweighted": {"role": "primary"},
            "stabilized_ipw_clip_0_05": {
                "minimum_propensity": 0.05,
                "maximum_propensity": 0.95,
                "maximum_weight": 20.0,
            },
        },
        "stabilization": {"normalize_within_axis_stratum_context": True},
    }
    outputs = run_stress(
        pd.DataFrame(score_rows),
        pd.DataFrame(cov_rows),
        pd.DataFrame(realm_rows),
        pattern,
        selection,
    )
    result = outputs["results"]
    focus = result.loc[
        result["context_layer"].eq("biogeographic_realm")
        & result["context"].eq("Palearctic")
        & result["selection_mode"].eq("stabilized_ipw_clip_0_05")
        & result["status"].eq("fit")
    ]
    assert set(focus["universe_mode"]) == {"full", "drop_aegean_lat12_lon20"}
    assert (focus["distance_estimate"] > 0).all()
    assert (focus["distance_p"] < 0.05).all()
