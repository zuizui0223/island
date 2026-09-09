from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_pr136_biogeographic_residual import run_biogeographic_residual


def _synthetic(*, opposite: bool, n_per_context: int = 80):
    contexts = ["northern_midlatitude", "tropical"]
    outcomes = ["plain_colour", "generalized_form", "self_compatibility"]
    cov_rows = []
    null_rows = []
    for context_index, context in enumerate(contexts):
        for i in range(n_per_context):
            island_id = f"{context}_{i}"
            distance = (i + 1) / n_per_context
            cluster_index = i // 4
            cov_rows.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{context_index}_{cluster_index}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": np.cos(i / 7),
                    "climate_pc1": np.sin(i / 9),
                    "climate_pc2": np.cos(i / 11),
                    "climate_pc3": np.sin(i / 13),
                    "climate_pc4": np.cos(i / 17),
                }
            )
            cluster_noise = 0.025 * np.sin((cluster_index + 1) * 1.37 + context_index)
            island_noise = 0.008 * np.cos((i + 2) * 1.91)
            for outcome_index, outcome in enumerate(outcomes):
                slope = 0.18 + 0.06 * outcome_index
                if opposite and context == "tropical":
                    slope *= -1
                outcome_noise = 0.004 * np.sin((outcome_index + 1) * (i + 1) * 0.71)
                residual = slope * distance + cluster_noise + island_noise + outcome_noise
                null_rows.append(
                    {
                        "island_id": island_id,
                        "outcome": outcome,
                        "stratum": "all_native",
                        "observed_n_species": 25,
                        "deviation_observed_minus_null": residual,
                    }
                )
    config = {
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
        "contexts": contexts,
        "strata": ["all_native"],
        "support_tiers": {"confirmatory": 50},
        "minimum_outcomes_per_vector": 2,
        "alpha": 0.05,
    }
    return pd.DataFrame(null_rows), pd.DataFrame(cov_rows), config


def test_detects_true_biogeographic_vector_difference():
    genus_null, covariates, config = _synthetic(opposite=True)
    _, _, within, between = run_biogeographic_residual(genus_null, covariates, config)
    assert set(within["context"]) == {"northern_midlatitude", "tropical"}
    assert within["status"].eq("fit").all()
    row = between.iloc[0]
    assert row["status"] == "fit"
    assert bool(row["biogeographic_difference_supported"])
    assert row["q_within_stratum_tier"] < 0.05


def test_same_slopes_do_not_manufacture_biogeographic_difference():
    genus_null, covariates, config = _synthetic(opposite=False)
    _, _, _, between = run_biogeographic_residual(genus_null, covariates, config)
    row = between.iloc[0]
    assert row["status"] == "fit"
    assert not bool(row["biogeographic_difference_supported"])
    assert pd.isna(row["q_within_stratum_tier"]) or row["q_within_stratum_tier"] > 0.05


def test_under_supported_context_is_not_testable_not_null():
    genus_null, covariates, config = _synthetic(opposite=True, n_per_context=25)
    _, _, within, between = run_biogeographic_residual(genus_null, covariates, config)
    assert within["status"].eq("not_testable").all()
    assert between.iloc[0]["status"] == "not_testable"


def test_secondary_input_needs_no_pollinator_columns():
    genus_null, covariates, config = _synthetic(opposite=True)
    assert not any("bombus" in c.lower() or "pollinator" in c.lower() for c in covariates.columns)
    _, _, within, between = run_biogeographic_residual(genus_null, covariates, config)
    assert not within.empty
    assert not between.empty
