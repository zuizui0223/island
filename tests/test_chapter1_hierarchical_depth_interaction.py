from __future__ import annotations

import numpy as np
import pandas as pd

from island_v2.chapter1_hierarchical_depth_interaction import (
    build_paired_stage_delta,
    fit_matched_depth_interaction,
)

COMPONENTS = [
    "shared_architecture_factor",
    "large_bee_like_residual",
    "butterfly_like_residual",
    "bird_like_residual",
]
CONTEXTS = ["northern_midlatitude", "tropical"]
BASELINE = [
    "log_island_area_km2",
    "climate_pc1",
    "climate_pc2",
    "climate_pc3",
    "climate_pc4",
]


def _synthetic_paired(n_per_context: int = 70) -> pd.DataFrame:
    rng = np.random.default_rng(17)
    rows = []
    for context_index, context in enumerate(CONTEXTS):
        for i in range(n_per_context):
            island = f"{context}_{i:03d}"
            distance = (i - n_per_context / 2) / 10
            block = f"b_{context_index}_{i // 4:02d}"
            covariates = {
                "log_distance_to_continent_km": distance,
                "log_island_area_km2": rng.normal(),
                "climate_pc1": rng.normal(),
                "climate_pc2": rng.normal(),
                "climate_pc3": rng.normal(),
                "climate_pc4": rng.normal(),
            }
            for component_index, component in enumerate(COMPONENTS):
                slope = 0.0 if context_index == 0 else 0.35 + 0.05 * component_index
                rows.append(
                    {
                        "island_id": island,
                        "architecture_component": component,
                        "analysis_regime": context,
                        "spatial_block": block,
                        "stage_delta": slope * distance + rng.normal(scale=0.08),
                        **covariates,
                    }
                )
    return pd.DataFrame(rows)


def test_matched_interaction_recovers_context_depth_difference() -> None:
    paired = _synthetic_paired()
    components, summary = fit_matched_depth_interaction(
        paired,
        components=COMPONENTS,
        contexts=CONTEXTS,
        context_column="analysis_regime",
        cluster_column="spatial_block",
        geography_column="log_distance_to_continent_km",
        baseline_covariates=BASELINE,
        minimum_paired_islands=50,
    )
    assert summary["status"] == "fit"
    assert summary["p_value"] < 1e-6
    assert len(components) == 4
    assert (components["context_difference_in_stage_delta_slope"] > 0).all()
    assert (components["n_islands_northern_midlatitude"] == 70).all()
    assert (components["n_islands_tropical"] == 70).all()


def test_matched_interaction_fails_closed_on_low_support() -> None:
    paired = _synthetic_paired(n_per_context=20)
    components, summary = fit_matched_depth_interaction(
        paired,
        components=COMPONENTS,
        contexts=CONTEXTS,
        context_column="analysis_regime",
        cluster_column="spatial_block",
        geography_column="log_distance_to_continent_km",
        baseline_covariates=BASELINE,
        minimum_paired_islands=50,
    )
    assert components.empty
    assert summary["status"] == "not_testable"
    assert "paired support" in summary["failure_reason"]


def test_stage_pairing_uses_only_common_support_and_beyond_genus() -> None:
    covariates = pd.DataFrame(
        {
            "island_id": ["i1", "i2", "i3"],
            "analysis_regime": ["northern_midlatitude", "tropical", "tropical"],
            "spatial_block": ["a", "b", "c"],
            "log_distance_to_continent_km": [1.0, 2.0, 3.0],
            "log_island_area_km2": [2.0, 2.5, 3.0],
            "climate_pc1": [0.0, 0.1, 0.2],
            "climate_pc2": [0.2, 0.1, 0.0],
            "climate_pc3": [1.0, 1.1, 1.2],
            "climate_pc4": [-0.2, -0.1, 0.0],
        }
    )
    pre = pd.DataFrame(
        {
            "island_id": ["i1", "i2", "i3"],
            "stratum": ["all_native"] * 3,
            "source_mode": ["geo_k5"] * 3,
            "syndrome": ["shared_architecture_factor"] * 3,
            "syndrome_score": [0.1, 0.2, 0.3],
        }
    )
    genus = pd.DataFrame(
        {
            "island_id": ["i1", "i2", "i3", "i2"],
            "stratum": ["all_native"] * 4,
            "source_mode": ["geo_k5"] * 4,
            "architecture_component": ["shared_architecture_factor"] * 4,
            "lineage_outcome": [
                "beyond_genus_residual",
                "beyond_genus_residual",
                "entry_enrichment",
                "loading_increment",
            ],
            "lineage_value": [0.4, 0.7, 99.0, 99.0],
            "n_represented_species": [10, 10, 10, 10],
            "n_represented_genera": [7, 7, 7, 7],
        }
    )
    paired = build_paired_stage_delta(
        pre,
        genus,
        covariates,
        source_mode="geo_k5",
        stratum="all_native",
        components=["shared_architecture_factor"],
        contexts=CONTEXTS,
        context_column="analysis_regime",
        cluster_column="spatial_block",
        geography_column="log_distance_to_continent_km",
        baseline_covariates=BASELINE,
    )
    assert paired["island_id"].tolist() == ["i1", "i2"]
    np.testing.assert_allclose(paired["stage_delta"], [0.3, 0.5])
