import math

import numpy as np
import pandas as pd

from island_v2.chapter1_information_weight_sensitivity import (
    reweight_composition,
    run_information_weight_sensitivity,
)


def _synthetic():
    rng = np.random.default_rng(53)
    comp = []
    cov = []
    responses = [("form", "a", 1.0), ("form", "b", -0.8), ("color", "c", 0.7)]
    idx = 0
    for context, scale in [("northern_midlatitude", 1.0), ("tropical", -0.35)]:
        for j in range(65):
            island_id = f"w{idx}"
            distance = -1.5 + 3.0 * j / 64
            cov.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{idx // 3}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": 2.0 + rng.normal(0, 0.2),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            for stratum in ["all_native", "native_nonendemic"]:
                for trait, state, slope in responses:
                    p = 1 / (1 + math.exp(-(-0.2 + scale * slope * distance)))
                    trials = 20 if j % 3 else 500
                    comp.append(
                        {
                            "island_id": island_id,
                            "stratum": stratum,
                            "trait_name": trait,
                            "trait_state": state,
                            "successes": int(rng.binomial(trials, p)),
                            "trials": trials,
                        }
                    )
            idx += 1
    return pd.DataFrame(comp), pd.DataFrame(cov)


def _config():
    return {
        "weighting_modes": ["canonical", "cap_100", "cap_50", "cap_20", "equal_island"],
        "headline_strata": ["all_native", "native_nonendemic"],
        "omnibus": {
            "baseline_covariates": [
                "log_island_area_km2",
                "climate_pc1",
                "climate_pc2",
                "climate_pc3",
                "climate_pc4",
            ],
            "isolation_column": "log_distance_to_continent_km",
            "context_column": "analysis_regime",
            "cluster_column": "spatial_block",
            "contexts": ["northern_midlatitude", "tropical"],
            "strata": ["all_native", "native_nonendemic"],
            "support_tiers": {"confirmatory": 50},
        },
    }


def test_reweighting_preserves_observed_share():
    comp, _ = _synthetic()
    original = comp["successes"] / comp["trials"]
    for mode in ["cap_100", "cap_50", "cap_20", "equal_island"]:
        weighted = reweight_composition(comp, mode)
        assert np.allclose(weighted["successes"] / weighted["trials"], original)
    equal = reweight_composition(comp, "equal_island")
    assert equal["trials"].eq(1.0).all()


def test_headline_is_recovered_under_equal_and_capped_information():
    comp, cov = _synthetic()
    _, _, summary = run_information_weight_sensitivity(comp, cov, _config())
    assert len(summary) == 10
    assert summary["headline_testable"].all()
    assert summary["headline_replication"].all()
