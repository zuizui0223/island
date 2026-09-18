import numpy as np
import pandas as pd

from island_v2.chapter1_v14_h2_decomposition import (
    _clustered_ols,
    build_plant_scores,
)


def test_build_plant_scores_creates_attraction_shift():
    rows = []
    values = {
        "selfing_core": 0.4,
        "generalized_accessible": 0.6,
        "large_bee_like": -0.2,
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
