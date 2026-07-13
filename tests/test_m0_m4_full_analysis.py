from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.full_analysis_covariates import build_covariates
from island_v2.m0_m4_full_analysis import (
    build_full_analysis_data,
    fit_model_ladder,
    summarize_tier_robustness,
)


def _config() -> dict:
    return {
        "analysis_contract": "test_full_analysis",
        "bombus_channel": {"predictor": "environmental_compatibility_max"},
        "climate": {
            "columns": ["bio1", "bio4", "bio5", "bio6", "bio12", "bio15", "bio18", "bio19", "elevation"],
            "n_pcs": 4,
            "expected_min_cumulative_variance": 0.80,
        },
        "spatial_inference": {"block_degrees": 10.0},
        "baseline_covariates": [
            "log_island_area_km2",
            "abs_latitude",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
            "log_trait_species_richness",
        ],
        "proxies": {
            "floral_phenotype": {
                "trait_name": "pollination_functional_guild",
                "positive_values": ["bumblebees"],
            },
            "autonomous_selfing": {
                "trait_name": "autonomous_selfing_capacity",
                "positive_values": ["autonomous"],
            },
            "self_compatibility": {
                "trait_name": "self_incompatibility",
                "positive_values": ["SC", "likely_SC"],
            },
            "alternative_pollen_vector": {
                "trait_name": "pollen_vector_mode",
                "positive_values": ["wind", "mixed"],
            },
        },
        "reporting": {"minimum_complete_islands": 20, "nominal_z_threshold": 1.96},
    }


def test_build_covariates_fits_deterministic_climate_pcs() -> None:
    config = _config()
    rows = []
    sample_rows = []
    manifest_rows = []
    for index in range(30):
        island_id = f"i{index:03d}"
        manifest_rows.append({"island_id": island_id, "area_km2": 5.0 + index * 2.0})
        climate = {
            "island_id": island_id,
            "bio1": 5 + index * 0.2,
            "bio4": 100 + index * 3,
            "bio5": 20 + index * 0.1,
            "bio6": -5 + index * 0.08,
            "bio12": 300 + index * 8,
            "bio15": 20 + index * 0.5,
            "bio18": 80 + index * 2,
            "bio19": 70 + index * 1.5,
            "elevation": 10 + index * 4,
        }
        rows.append(climate)
        for offset in (-0.05, 0.05):
            sample_rows.append(
                {
                    "island_id": island_id,
                    "decimal_longitude": -170 + index * 10 + offset,
                    "decimal_latitude": -60 + index * 4 + offset,
                }
            )

    covariates, loadings, manifest = build_covariates(
        pd.DataFrame(rows), pd.DataFrame(sample_rows), pd.DataFrame(manifest_rows), config
    )
    assert len(covariates) == 30
    assert not covariates[[f"climate_pc{i}" for i in range(1, 5)]].isna().any().any()
    assert manifest["pca"]["cumulative_explained_variance"] >= 0.80
    assert manifest["n_spatial_blocks"] > 1
    assert set(loadings["variable"]) == set(config["climate"]["columns"])


def _synthetic_tables(n_islands: int = 120) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(42)
    composition_rows = []
    metric_rows = []
    covariate_rows = []
    for index in range(n_islands):
        island_id = f"island_{index:03d}"
        bombus = index / (n_islands - 1)
        richness = 80 + index % 40
        bumblebee_count = int(np.clip(round(10 + 50 * bombus + rng.normal(0, 3)), 1, richness - 1))
        autonomous_count = int(np.clip(round(45 - 20 * bombus + rng.normal(0, 3)), 1, richness - 1))
        sc_count = int(np.clip(round(50 - 15 * bombus + rng.normal(0, 3)), 1, richness - 1))
        wind_count = int(np.clip(round(35 - 10 * bombus + rng.normal(0, 3)), 1, richness - 1))

        for trait_name, positive_value, negative_value, positive_count in (
            ("pollination_functional_guild", "bumblebees", "other_bees", bumblebee_count),
            ("autonomous_selfing_capacity", "autonomous", "dependent", autonomous_count),
            ("self_incompatibility", "SC", "SI", sc_count),
            ("pollen_vector_mode", "wind", "biotic", wind_count),
        ):
            composition_rows.extend(
                [
                    {
                        "island_id": island_id,
                        "trait_name": trait_name,
                        "filled_value": positive_value,
                        "n_species": positive_count,
                    },
                    {
                        "island_id": island_id,
                        "trait_name": trait_name,
                        "filled_value": negative_value,
                        "n_species": richness - positive_count,
                    },
                ]
            )

        metric_rows.append(
            {
                "island_id": island_id,
                "environmental_compatibility_max": bombus,
                "n_trait_species": richness,
            }
        )
        covariate_rows.append(
            {
                "island_id": island_id,
                "log_island_area_km2": np.log1p(10 + index),
                "abs_latitude": 5 + index % 65,
                "climate_pc1": np.sin(index / 10),
                "climate_pc2": np.cos(index / 12),
                "climate_pc3": np.sin(index / 7),
                "climate_pc4": np.cos(index / 9),
                "spatial_block": f"block_{index // 10:02d}",
            }
        )
    return pd.DataFrame(composition_rows), pd.DataFrame(metric_rows), pd.DataFrame(covariate_rows)


def test_full_model_ladder_uses_independent_m4_proxy_and_clustered_binomial() -> None:
    config = _config()
    composition, metrics, covariates = _synthetic_tables()
    data, metadata = build_full_analysis_data(composition, metrics, covariates, config)
    assert metadata["baseline_complete"] is True
    assert not np.allclose(
        data["floral_phenotype_share"].fillna(0),
        1.0 - data["alternative_pollen_vector_share"].fillna(0),
    )

    coefficients, fits, indirect, model_summary = fit_model_ladder(data, config)
    assert set(coefficients["model"]) == {"M0", "M1", "M2", "M3", "M4"}
    assert set(fits["status"]) <= {"fitted", "max_iter_reached"}
    assert (fits["n_clusters"] >= 10).all()
    m1_bombus = coefficients.loc[
        (coefficients["model"] == "M1")
        & (coefficients["outcome"] == "floral_phenotype")
        & (coefficients["predictor"] == "bombus_channel_state")
    ]
    assert len(m1_bombus) == 1
    assert float(m1_bombus.iloc[0]["estimate_log_odds"]) > 0
    m4_predictors = set(coefficients.loc[coefficients["model"] == "M4", "predictor"])
    assert "alternative_pollen_vector" in m4_predictors
    assert "bombus_x_alternative" in m4_predictors
    assert set(indirect["mediator"]) == {
        "autonomous_selfing",
        "self_compatibility",
        "parallel_total",
    }
    assert set(model_summary["model"]) == {"M0", "M1", "M2", "M3", "M4"}


def test_tier_robustness_reports_direction_and_effect_spread() -> None:
    base = pd.DataFrame(
        [
            {
                "model": "M1",
                "equation": 1,
                "outcome": "floral_phenotype",
                "predictor": "bombus_channel_state",
                "estimate_log_odds": 0.4,
                "z_value": 3.0,
            }
        ]
    )
    tables = {
        "primary": base.copy(),
        "broad": base.assign(estimate_log_odds=0.3, z_value=2.5),
        "sensitivity_all": base.assign(estimate_log_odds=0.2, z_value=2.1),
    }
    summary = summarize_tier_robustness(tables)
    assert bool(summary.iloc[0]["direction_consistent_across_tiers"])
    assert bool(summary.iloc[0]["all_tiers_nominally_supported"])
    assert np.isclose(summary.iloc[0]["max_absolute_coefficient_difference"], 0.2)
