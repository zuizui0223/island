from __future__ import annotations

import pandas as pd

from island_v2.chapter1_h5_pollination_triangulation import (
    architecture_specificity,
    build_template_consistency,
    integrated_decision,
)


def _config() -> dict:
    return {
        "plant_side_consistency": {
            "primary_evidence_scope": "direct_only",
            "sensitivity_evidence_scope": "all_analysis_eligible",
            "strata": ["all_native", "native_nonendemic"],
            "context_layer": "analysis_regime",
            "predictions": {
                "northern_midlatitude": {
                    "sampled_large_bee_concordance": "negative",
                },
                "tropical": {
                    "sampled_butterfly_concordance": "nonnegative",
                },
            },
            "identity_specificity_guardrail": {
                "architecture_factor_variance_fraction_ceiling": 0.80,
            },
        },
        "independent_existing_checks": {
            "N1_channel_heterogeneity": {"status": "not_promoted"},
            "H5c_biotic_vs_wind": {"status": "not_promoted"},
        },
    }


def _slopes() -> pd.DataFrame:
    rows = []
    for stratum in ["all_native", "native_nonendemic"]:
        rows += [
            {
                "context_layer": "analysis_regime",
                "axis_set": "sampled_guild_concordance",
                "support_tier": "confirmatory",
                "context": "northern_midlatitude",
                "stratum": stratum,
                "syndrome": "sampled_large_bee_concordance",
                "distance_slope": -0.1,
                "q_axis_family": 0.01,
                "n_islands": 100,
            },
            {
                "context_layer": "analysis_regime",
                "axis_set": "sampled_guild_concordance",
                "support_tier": "confirmatory",
                "context": "tropical",
                "stratum": stratum,
                "syndrome": "sampled_butterfly_concordance",
                "distance_slope": 0.1,
                "q_axis_family": 0.01,
                "n_islands": 80,
            },
        ]
    return pd.DataFrame(rows)


def test_triangulation_separates_consistency_from_identity() -> None:
    cfg = _config()
    consistency = build_template_consistency(_slopes(), _slopes(), cfg)
    factors = pd.DataFrame({"variance_fraction": [0.87, 0.87, 0.87]})
    specificity = architecture_specificity(factors, factors, cfg)
    globi = pd.DataFrame({"promoted": [False, False]})
    decision = integrated_decision(consistency, specificity, globi, cfg)
    assert decision["plant_side_pollination_architecture_consistent"] is True
    assert decision["named_pollinator_identity_specificity_passed"] is False
    assert decision["independent_pollination_support_present"] is False
    assert decision["mechanism_status"] == "pollination_plausible_but_not_identified"
