from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_all_data_ladder import run_ladder


def _configs():
    probability = {
        "model_outcomes": ["trait_a", "trait_b"],
        "strata": ["all_observed"],
        "contexts": ["north", "tropical"],
        "between_contexts": [["north", "tropical"]],
        "minimum_islands_per_outcome": 20,
        "minimum_outcomes_per_vector": 2,
        "alpha": 0.05,
        "max_iter": 500,
        "geography_column": "distance",
        "context_column": "context",
        "cluster_column": "block",
        "baseline_covariates": ["area", "climate1"],
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
    }
    return probability, ladder


def _synthetic(seed: int = 7):
    rng = np.random.default_rng(seed)
    islands = []
    counts = []
    island_number = 0
    for context in ("north", "tropical"):
        for i in range(80):
            island_id = f"{context}_{i}"
            distance = rng.normal()
            area = rng.normal()
            climate = rng.normal()
            block = f"{context}_b{i // 8}"
            islands.append(
                {
                    "island_id": island_id,
                    "distance": distance,
                    "area": area,
                    "climate1": climate,
                    "context": context,
                    "block": block,
                }
            )
            for outcome_index, outcome in enumerate(("trait_a", "trait_b")):
                base = -0.2 + 0.15 * outcome_index
                distance_effect = 0.7 if context == "north" else -0.25
                moderation = -0.65 if context == "north" else 0.30
                eta = base + distance_effect * distance + 0.2 * area + moderation * distance * area
                probability = 1.0 / (1.0 + np.exp(-eta))
                trials = 40
                successes = rng.binomial(trials, probability)
                counts.append(
                    {
                        "island_id": island_id,
                        "successes": successes,
                        "trials": trials,
                        "outcome": outcome,
                        "stratum": "all_observed",
                    }
                )
            island_number += 1
    return pd.DataFrame(counts), pd.DataFrame(islands)


def test_ladder_recovers_global_branching_and_area_moderation():
    probability, ladder = _configs()
    counts, covariates = _synthetic()
    result = run_ladder(counts, covariates, probability, ladder)

    h1 = result["h1_global_omnibus"].iloc[0]
    assert h1["status"] == "fit"
    assert bool(h1["universal_response_rejected"])

    h2 = result["h2_between_omnibus"].iloc[0]
    assert h2["status"] == "fit"
    assert h2["p_value"] < 0.05

    h3_within = result["h3_within_omnibus"].set_index("context")
    assert h3_within.loc["north", "status"] == "fit"
    assert h3_within.loc["north", "p_value"] < 0.05
    assert h3_within.loc["north", "directional_classification"] == "small_island_amplification"

    h3_between = result["h3_between_omnibus"].iloc[0]
    assert h3_between["status"] == "fit"
    assert h3_between["p_value"] < 0.05
