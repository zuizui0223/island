import numpy as np
import pandas as pd

from island_v2.chapter1_v14_h2_decomposition import (
    _clustered_ols,
    build_plant_scores,
    run_decomposition,
)


def test_build_plant_scores_creates_attraction_shift():
    rows = []
    values = {
        "selfing_core": 0.4,
        "generalized_accessible": 0.6,
        "large_bee_like": -0.2,
        "butterfly_like": 0.1,
        "bird_like": -0.1,
    }
    for syndrome, value in values.items():
        rows.append(
            {
                "island_id": "i1",
                "stratum": "all_observed",
                "syndrome": syndrome,
                "syndrome_score": value,
            }
        )
    result = build_plant_scores(pd.DataFrame(rows), "all_observed")
    assert np.isclose(result.loc[0, "attraction_shift"], 0.4)


def test_clustered_ols_recovers_positive_distance_after_selfing():
    n = 120
    distance = np.linspace(-2.0, 2.0, n)
    selfing = np.sin(np.linspace(0.0, 6.0, n))
    response = 0.8 * distance + 0.5 * selfing
    frame = pd.DataFrame(
        {
            "response": response,
            "distance": distance,
            "selfing_core": selfing,
            "block": [f"b{i // 6}" for i in range(n)],
        }
    )
    result = _clustered_ols(
        frame,
        response="response",
        predictors=["distance", "selfing_core"],
        cluster_column="block",
    )
    assert result["status"] == "fit"
    assert result["coefficients"]["z_distance"]["estimate"] > 0



def test_run_decomposition_fits_continuous_and_plain_colour_routes():
    n = 72
    distance = np.linspace(-1.5, 1.5, n)
    island_ids = [f"i{i}" for i in range(n)]
    covariates = pd.DataFrame(
        {
            "island_id": island_ids,
            "log_distance_to_continent_km": distance,
            "analysis_regime": ["test_context"] * n,
            "spatial_block": [f"b{i // 6}" for i in range(n)],
            "log_island_area_km2": 0.4 * np.sin(np.linspace(0.0, 9.0, n)),
            "climate_pc1": np.cos(np.linspace(0.0, 5.0, n)),
            "climate_pc2": np.sin(np.linspace(0.0, 8.0, n)),
            "climate_pc3": np.cos(np.linspace(0.0, 11.0, n)),
            "climate_pc4": np.sin(np.linspace(0.0, 13.0, n)),
        }
    )
    selfing = 0.35 * distance + 0.1 * np.sin(np.linspace(0.0, 5.0, n))
    accessible = 0.45 * distance + 0.2 * selfing
    large_bee = -0.25 * distance
    butterfly = 0.05 * distance
    bird = 0.02 * distance
    score_rows = []
    for island_id, s, a, b, bf, bd in zip(
        island_ids, selfing, accessible, large_bee, butterfly, bird, strict=True
    ):
        for syndrome, value in (
            ("selfing_core", s),
            ("generalized_accessible", a),
            ("large_bee_like", b),
            ("butterfly_like", bf),
            ("bird_like", bd),
        ):
            score_rows.append(
                {
                    "island_id": island_id,
                    "stratum": "all_observed",
                    "syndrome": syndrome,
                    "syndrome_score": value,
                }
            )
    syndrome_scores = pd.DataFrame(score_rows)
    probability = 1.0 / (1.0 + np.exp(-(0.55 * distance + 0.15 * selfing)))
    trials = np.full(n, 30)
    successes = np.rint(probability * trials).astype(int)
    counts = pd.DataFrame(
        {
            "island_id": island_ids,
            "successes": successes,
            "trials": trials,
            "outcome": ["plain_colour"] * n,
            "stratum": ["all_observed"] * n,
        }
    )
    config = {
        "contract": "chapter1_v14_h2_decomposition_v1",
        "geography_column": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "contexts": ["test_context"],
        "primary_stratum": "all_observed",
        "minimum_islands": 50,
    }
    results, summary = run_decomposition(
        counts,
        syndrome_scores,
        covariates,
        config,
        evidence_scope="test",
    )
    assert set(results["response"]) == {
        "selfing_core",
        "generalized_accessible",
        "attraction_shift",
        "large_bee_like",
        "butterfly_like",
        "bird_like",
        "shared_named_architecture",
        "plain_colour",
    }
    assert results["status"].eq("fit").all()
    for response in ("selfing_core", "generalized_accessible", "plain_colour"):
        estimate = results.loc[results["response"].eq(response), "distance_estimate"].iloc[0]
        assert estimate > 0
    assert summary["n_contexts_accessibility_positive_after_selfing"] == 1
