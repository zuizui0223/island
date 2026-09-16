from __future__ import annotations

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_h5_orthogonalized_architecture_structural_correction import (
    fit_scalar_shared_h2,
    integrate_structural_correction,
    summarize_rank2_residual_h2,
)


def _config() -> dict:
    return yaml.safe_load(
        open(
            "config/chapter1_h5_orthogonalized_architecture_structural_correction_v2.yml",
            encoding="utf-8",
        )
    )


def test_rank2_residual_summary_accepts_three_labels_with_intrinsic_df_two() -> None:
    between = pd.DataFrame(
        [
            {
                "context_layer": "analysis_regime",
                "axis_set": "identity_specific_residuals",
                "stratum": "all_observed",
                "support_tier": "confirmatory",
                "context_a": "northern_midlatitude",
                "context_b": "tropical",
                "status": "fit",
                "n_retained_syndromes": 3,
                "retained_syndromes": "bird_like_residual|butterfly_like_residual|large_bee_like_residual",
                "joint_context_difference_df": 2,
                "p_value": 0.01,
                "n_unique_islands": 154,
                "n_clusters": 57,
            }
        ]
    )
    result = summarize_rank2_residual_h2(between, evidence_scope="all_analysis_eligible", config=_config())
    assert result["evaluable"]
    assert result["joint_df"] == 2
    assert result["joint_p"] == 0.01

    between.loc[0, "joint_context_difference_df"] = 1
    failed = summarize_rank2_residual_h2(between, evidence_scope="all_analysis_eligible", config=_config())
    assert not failed["evaluable"]
    assert failed["joint_p"] is None


def test_scalar_shared_h2_recovers_positive_tropical_distance_interaction() -> None:
    rng = np.random.default_rng(44)
    rows = []
    covariates = []
    n_per_context = 80
    for context_index, context in enumerate(["northern_midlatitude", "tropical"]):
        for i in range(n_per_context):
            island_id = f"{context}_{i}"
            distance = -1.5 + 3.0 * i / (n_per_context - 1)
            area = np.sin(i / 7.0) + rng.normal(0, 0.05)
            climate = [
                np.cos(i / 9.0) + rng.normal(0, 0.05),
                np.sin(i / 11.0) + rng.normal(0, 0.05),
                np.cos(i / 13.0) + rng.normal(0, 0.05),
                np.sin(i / 17.0) + rng.normal(0, 0.05),
            ]
            interaction = 0.65 if context_index == 1 else 0.0
            response = 0.10 * distance + interaction * distance + rng.normal(0, 0.12)
            rows.append(
                {
                    "island_id": island_id,
                    "stratum": "all_observed",
                    "syndrome": "shared_architecture_factor",
                    "syndrome_score": response,
                    "n_species": 70,
                }
            )
            covariates.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"block_{context_index}_{i // 4}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": area,
                    "climate_pc1": climate[0],
                    "climate_pc2": climate[1],
                    "climate_pc3": climate[2],
                    "climate_pc4": climate[3],
                }
            )

    result = fit_scalar_shared_h2(pd.DataFrame(rows), pd.DataFrame(covariates), _config())
    assert result["evaluable"]
    assert result["northern_midlatitude_n_islands"] == 80
    assert result["tropical_n_islands"] == 80
    assert result["distance_by_tropical_estimate"] > 0.4
    assert result["distance_by_tropical_p"] < 0.05


def test_integrated_structural_correction_requires_both_evidence_scopes() -> None:
    all_result = {
        "residual": {"evaluable": True, "joint_p": 0.01},
        "shared": {"evaluable": True, "distance_by_tropical_p": 0.30},
    }
    direct_result = {
        "residual": {"evaluable": True, "joint_p": 0.20},
        "shared": {"evaluable": True, "distance_by_tropical_p": 0.40},
    }
    decision = integrate_structural_correction(all_result, direct_result, _config())
    assert not decision["robust_identity_specific_residual_H2"]
    assert not decision["robust_shared_architecture_H2"]
    assert decision["classification"] == "all_analysis_only_residual_signal_scope_sensitive"
    assert not decision["pollinator_identity_promoted"]
