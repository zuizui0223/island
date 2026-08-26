import numpy as np
import pandas as pd

from island_v2.chapter1_pr138_palearctic_restricted_ipw import (
    build_within_palearctic_weights,
    run_palearctic_restricted_ipw,
)


def test_within_palearctic_weights_use_availability_not_score_values():
    rng = np.random.default_rng(42)
    universe = pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(100)],
            "log_distance_to_continent_km": np.linspace(-1, 1, 100),
            "log_island_area_km2": rng.normal(size=100),
            "climate_pc1": rng.normal(size=100),
            "spatial_block": [f"b{i // 4}" for i in range(100)],
        }
    )
    # Availability must not be deterministically separable by a modeled predictor.
    # Draw a fixed random subset so this test checks outcome blindness rather than
    # creating an artificial distance threshold that forces logistic separation.
    observed = set(rng.choice(universe["island_id"].to_numpy(), size=60, replace=False))
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
    }
    selection = {
        "selection_model": {"minimum_observed": 20, "minimum_unobserved": 20},
        "weight_modes": {
            "unweighted": {},
            "stabilized_ipw_clip_0_05": {
                "minimum_propensity": 0.05,
                "maximum_propensity": 0.95,
                "maximum_weight": 20.0,
            },
        },
        "stabilization": {"normalize_within_axis_stratum_context": True},
    }
    weights, fit, _ = build_within_palearctic_weights(
        observed, universe, pattern, selection
    )
    assert fit["status"] == "fit"
    assert set(weights["island_id"]) == observed
    assert "syndrome_score" not in weights.columns
    weighted = weights.loc[weights["selection_mode"].eq("stabilized_ipw_clip_0_05")]
    assert np.isfinite(weighted["analysis_weight"]).all()
    assert weighted["analysis_weight"].max() <= 20.0


def test_palearctic_restricted_ipw_keeps_positive_signal_with_and_without_aegean():
    rng = np.random.default_rng(1389)
    n = 180
    cov_rows = []
    realm_rows = []
    score_rows = []
    for i in range(n):
        island = f"i{i}"
        x = -1 + 2 * i / (n - 1)
        block = "lat12_lon20" if i % 19 == 0 else f"b{i // 5}"
        cov_rows.append(
            {
                "island_id": island,
                "spatial_block": block,
                "log_distance_to_continent_km": x,
                "log_island_area_km2": np.sin(i / 11),
                "climate_pc1": np.cos(i / 13),
            }
        )
        realm_rows.append(
            {"island_id": island, "biogeographic_realm": "Palearctic"}
        )
        if rng.random() < 0.68:
            attraction = 0.32 * x + rng.normal(0, 0.08)
            for syndrome, score in [
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
                        "syndrome_score": score,
                    }
                )
    pattern = {
        "geography_column": "log_distance_to_continent_km",
        "cluster_column": "spatial_block",
        "baseline_covariates": ["log_island_area_km2", "climate_pc1"],
        "support_tiers": {"confirmatory": 30},
    }
    selection = {
        "selection_model": {"minimum_observed": 30, "minimum_unobserved": 30},
        "weight_modes": {
            "unweighted": {},
            "stabilized_ipw_clip_0_05": {
                "minimum_propensity": 0.05,
                "maximum_propensity": 0.95,
                "maximum_weight": 20.0,
            },
        },
        "stabilization": {"normalize_within_axis_stratum_context": True},
    }
    outputs = run_palearctic_restricted_ipw(
        pd.DataFrame(score_rows),
        pd.DataFrame(cov_rows),
        pd.DataFrame(realm_rows),
        pattern,
        selection,
    )
    result = outputs["results"]
    focus = result.loc[
        result["restriction"].eq("si_only")
        & result["selection_mode"].eq("stabilized_ipw_clip_0_05")
        & result["status"].eq("fit")
    ]
    assert set(focus["universe_mode"]) == {
        "full_palearctic",
        "palearctic_drop_aegean_lat12_lon20",
    }
    assert (focus["distance_estimate"] > 0).all()
    assert (focus["distance_p"] < 0.05).all()
