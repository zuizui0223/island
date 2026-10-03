import math

import numpy as np
import pandas as pd

from island_v2.chapter1_final_reviewer_robust import (
    _holm,
    _jackknife_summary,
    _run_h2_model,
)


def test_holm_is_monotone_in_rank_and_bounded():
    p = pd.Series([0.01, 0.04, 0.03, 0.20])
    q = _holm(p)
    assert np.all((q.dropna() >= 0) & (q.dropna() <= 1))
    order = np.argsort(p.to_numpy())
    ranked_q = q.to_numpy()[order]
    assert np.all(np.diff(ranked_q) >= -1e-12)
    assert math.isclose(float(q.iloc[0]), 0.04, rel_tol=0, abs_tol=1e-12)


def test_jackknife_summary_uses_cluster_df_and_one_sided_direction():
    loo = np.array([0.18, 0.20, 0.21, 0.19, 0.22])
    out = _jackknife_summary(0.20, loo)
    assert out["jackknife_n_clusters"] == 5
    assert out["df"] == 4
    assert out["jackknife_se"] > 0
    assert 0 <= out["p_one_sided_positive"] <= out["p_two_sided"] <= 1
    assert out["loo_nonpositive_count"] == 0


def test_h2_cluster_jackknife_recovers_positive_distance_signal():
    rows = []
    for cluster in range(8):
        for j in range(5):
            distance = cluster + j / 10
            rows.append(
                {
                    "island_id": f"I{cluster}_{j}",
                    "spatial_block": f"C{cluster}",
                    "distance": distance,
                    "area": (j + 1) * 0.2,
                    "response": 0.4 * distance + 0.05 * j,
                }
            )
    frame = pd.DataFrame(rows)
    result = _run_h2_model(
        frame,
        response="response",
        predictors=["distance", "area"],
        geography="distance",
        cluster="spatial_block",
        minimum_clusters=5,
    )
    assert result["status"] == "fit"
    assert result["n_clusters"] == 8
    assert result["estimate"] > 0
    assert result["p_one_sided_positive"] < 0.05


def test_h2_jackknife_refuses_too_few_clusters():
    frame = pd.DataFrame(
        {
            "island_id": ["a", "b", "c", "d"],
            "spatial_block": ["A", "A", "B", "B"],
            "distance": [0.0, 1.0, 2.0, 3.0],
            "response": [0.0, 1.0, 2.0, 3.0],
        }
    )
    result = _run_h2_model(
        frame,
        response="response",
        predictors=["distance"],
        geography="distance",
        cluster="spatial_block",
        minimum_clusters=3,
    )
    assert result["status"] == "not_testable"
