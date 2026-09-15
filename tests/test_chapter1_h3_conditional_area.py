from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_h3_conditional_area import run_conditional_area


def test_conditional_area_contrast_changes_across_fixed_area_points():
    rng = np.random.default_rng(11)
    covariates = []
    counts = []
    for context in ("north", "tropical"):
        for i in range(90):
            island = f"{context}_{i}"
            distance = rng.normal()
            area = rng.normal()
            climate = rng.normal()
            covariates.append(
                {
                    "island_id": island,
                    "distance": distance,
                    "area": area,
                    "climate1": climate,
                    "context": context,
                    "block": f"{context}_{i // 9}",
                }
            )
            for outcome in ("a", "b"):
                distance_effect = 0.5 if context == "north" else 0.1
                moderation = -0.5 if context == "north" else 0.4
                eta = -0.1 + distance_effect * distance + 0.15 * area + moderation * distance * area
                prob = 1.0 / (1.0 + np.exp(-eta))
                trials = 50
                counts.append(
                    {
                        "island_id": island,
                        "successes": rng.binomial(trials, prob),
                        "trials": trials,
                        "outcome": outcome,
                        "stratum": "all_observed",
                    }
                )

    probability = {
        "model_outcomes": ["a", "b"],
        "strata": ["all_observed"],
        "minimum_islands_per_outcome": 20,
        "minimum_outcomes_per_vector": 2,
        "max_iter": 500,
    }
    ladder = {
        "contract": "test",
        "contexts": ["north", "tropical"],
        "primary_between_contexts": [["north", "tropical"]],
        "minimum_islands_per_outcome": 20,
        "minimum_outcomes_per_vector": 2,
        "alpha": 0.05,
        "geography_column": "distance",
        "area_column": "area",
        "context_column": "context",
        "cluster_column": "block",
        "control_columns": ["climate1"],
        "hypotheses": {
            "H3_island_capacity_moderation": {"conditional_area_z_values": [-1.0, 0.0, 1.0]}
        },
    }
    slopes, omnibus = run_conditional_area(
        pd.DataFrame(counts), pd.DataFrame(covariates), probability, ladder
    )
    assert set(omnibus["area_z"]) == {-1.0, 0.0, 1.0}
    assert omnibus["status"].eq("fit").all()
    means = slopes.groupby("area_z")["conditional_distance_difference_b_minus_a"].mean()
    assert means.loc[-1.0] < means.loc[1.0]
