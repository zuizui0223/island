import importlib

import pandas as pd

mod = importlib.import_module("island_v2.flora_status_support_decision")


def _row(ceiling: int, direct: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "regime": ["southern"],
            "stratum": ["all_native"],
            "outcome": ["floral_form"],
            "n_islands_with_stratum": [ceiling],
            "n_direct_eligible_islands": [direct],
        }
    )


def test_status_ceiling_below_50_blocks_trait_acquisition():
    out = mod.classify_support_decision(_row(34, 31)).iloc[0]
    assert out["support_decision"] == "status_ceiling_below_confirmatory"
    assert not bool(out["trait_acquisition_allowed"])
    assert out["status_ceiling_gap_to_confirmatory"] == 16


def test_confirmatory_support_models_before_more_acquisition():
    out = mod.classify_support_decision(_row(241, 235)).iloc[0]
    assert out["support_decision"] == "confirmatory_count_met"
    assert out["next_action"] == "model_before_acquiring_more_traits"
    assert not bool(out["trait_acquisition_allowed"])


def test_only_trait_limited_stratum_allows_trait_acquisition():
    out = mod.classify_support_decision(_row(80, 41)).iloc[0]
    assert out["support_decision"] == "targeted_trait_acquisition_zone"
    assert bool(out["trait_acquisition_allowed"])


def test_threshold_robustness_reports_trial_count_sensitivity():
    support = pd.DataFrame(
        {
            "island_id": ["a", "b", "c"],
            "analysis_regime": ["north"] * 3,
            "stratum": ["endemic"] * 3,
            "outcome": ["flower_colour"] * 3,
            "n_direct_species": [1, 4, 12],
        }
    )
    out = mod.threshold_robustness(support).set_index("min_direct_species")
    assert out.loc[1, "n_islands"] == 3
    assert out.loc[3, "n_islands"] == 2
    assert out.loc[5, "n_islands"] == 1
    assert out.loc[10, "n_islands"] == 1
