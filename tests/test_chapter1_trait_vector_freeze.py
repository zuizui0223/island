import pandas as pd

from island_v2.chapter1_trait_vector_freeze import freeze_trait_vectors


def test_trait_vector_freeze_preserves_direction_and_excludes_syndrome_labels():
    slopes = pd.DataFrame(
        {
            "stratum": ["all_native"] * 6,
            "trait_name": ["floral_form"] * 6,
            "trait_state": ["tubular", "tubular", "tubular", "open_radial", "open_radial", "open_radial"],
            "context": [
                "northern_midlatitude", "tropical", "southern_extratropical",
                "northern_midlatitude", "tropical", "southern_extratropical",
            ],
            "reference_context": ["northern_midlatitude"] * 6,
            "isolation_slope_log_odds_per_sd": [0.8, -0.7, 0.05, -0.9, 0.6, 0.02],
            "cluster_robust_se": [0.1] * 6,
            "z_value": [8, -7, 0.5, -9, 6, 0.2],
            "p_value": [1e-8, 1e-7, 0.6, 1e-9, 1e-6, 0.8],
            "n_context_islands": [40] * 6,
        }
    )
    coefficients = pd.DataFrame(
        {
            "stratum": ["all_native", "all_native"],
            "trait_name": ["floral_form", "floral_form"],
            "trait_state": ["tubular", "open_radial"],
            "model": ["M2", "M2"],
            "predictor": [
                "z_isolation:context[tropical]",
                "z_isolation:context[tropical]",
            ],
            "estimate_log_odds": [-1.5, 1.5],
            "cluster_robust_se": [0.2, 0.2],
            "z_value": [-7.5, 7.5],
            "p_value": [1e-12, 1e-12],
        }
    )
    lineage = pd.DataFrame(
        {
            "stratum": ["all_native"],
            "context": ["northern_midlatitude"],
            "outcome": ["generalized_form"],
            "predictor": ["log_distance_to_continent_km"],
            "estimate": [-0.01],
            "cluster_robust_se": [0.002],
            "z_value": [-5.0],
            "p_value": [1e-6],
        }
    )

    vector, interaction, lineage_out, patterns = freeze_trait_vectors(
        slopes, coefficients, lineage
    )
    assert vector["q_within_stratum_trait_context"].notna().all()
    tubular = patterns.loc[patterns["trait_state"].eq("tubular")].iloc[0]
    assert tubular["pattern_class"] == "fdr_directional_reversal"
    assert interaction["interaction_fdr_supported"].all()
    assert lineage_out.iloc[0]["lineage_direction_class"] == "fdr_negative"
    assert not any("bombus" in str(value).lower() for value in vector.columns)
