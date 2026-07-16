import pandas as pd

from island_v2 import trait_inference_validation as validation
from island_v2 import trait_inference_validation_fast as validation_fast


def _direct_cells() -> pd.DataFrame:
    rows = [
        ("Aaa target", "Aaa", "Fam1", "white"),
        ("Aaa red", "Aaa", "Fam1", "red_pink"),
        ("Aaa blue", "Aaa", "Fam1", "blue_purple"),
        ("Ccc red", "Ccc", "Fam1", "red_pink"),
        ("Ddd red", "Ddd", "Fam1", "red_pink"),
        ("Eee red", "Eee", "Fam1", "red_pink"),
        ("Xxx lone", "Xxx", "Loneaceae", "white"),
    ]
    return pd.DataFrame(
        [
            {
                "accepted_species": species,
                "genus": genus,
                "family": family,
                "trait_name": "flower_primary_color",
                "filled_value": value,
                "support_n": "1",
            }
            for species, genus, family, value in rows
        ]
    )


def test_loo_tie_uses_qualified_family_and_never_global() -> None:
    predictions = validation.leave_one_species_out_predictions(
        _direct_cells(),
        genus_min_support=1,
        family_min_support=3,
        genus_min_consensus=1.0,
        family_min_consensus=0.8,
        family_min_supporting_genera=2,
    ).set_index("accepted_species")

    target = predictions.loc["Aaa target"]
    assert target["genus_top_tie_n"] == 2
    assert target["cascade_prediction_tier"] == "family_inference"
    assert target["cascade_prediction"] == "red_pink"
    assert target["family_support_species"] == 4
    assert target["family_total_support_species"] == 5

    lone = predictions.loc["Xxx lone"]
    assert lone["cascade_prediction_tier"] == "unresolved_no_evidence"
    assert lone["cascade_prediction"] == ""
    assert "global_fallback" not in set(predictions["cascade_prediction_tier"])

    metrics = validation._classification_metrics(
        predictions.reset_index(), "actual_value", "cascade_prediction"
    )
    assert metrics["n"] == int(predictions["cascade_prediction"].ne("").sum())
    assert metrics["n"] < len(predictions)


def test_cached_and_slow_threshold_grids_match_without_global() -> None:
    direct = _direct_cells()
    predictions = validation.leave_one_species_out_predictions(
        direct,
        genus_min_support=1,
        family_min_support=3,
        genus_min_consensus=1.0,
        family_min_consensus=0.8,
        family_min_supporting_genera=2,
    )
    slow = validation.threshold_grid(
        direct,
        genus_thresholds=[1, 2],
        family_thresholds=[3],
        genus_min_consensus=1.0,
        family_min_consensus=0.8,
        family_min_supporting_genera=2,
    ).sort_values(["trait_name", "genus_min_support", "family_min_support"])
    fast = validation_fast.threshold_grid_from_cached_predictions(
        predictions,
        genus_thresholds=[1, 2],
        family_thresholds=[3],
        genus_min_consensus=1.0,
        family_min_consensus=0.8,
        family_min_supporting_genera=2,
    ).sort_values(["trait_name", "genus_min_support", "family_min_support"])

    pd.testing.assert_frame_equal(
        slow.reset_index(drop=True), fast.reset_index(drop=True), check_dtype=False
    )
    assert "n_unresolved_predictions" in slow.columns
    assert "n_global_predictions" not in slow.columns


def test_operational_loo_applies_trait_specific_family_policy() -> None:
    predictions = validation.leave_one_species_out_predictions(
        _direct_cells(),
        genus_min_support=1,
        family_min_support=3,
        genus_min_consensus=1.0,
        family_min_consensus=0.8,
        family_min_supporting_genera=2,
        inference_policies={
            "flower_primary_color": {
                "min_genus_support": 2,
                "family_allowed": False,
            }
        },
    ).set_index("accepted_species")

    target = predictions.loc["Aaa target"]
    assert target["operational_genus_min_support"] == 2
    assert not target["operational_family_allowed"]
    assert target["cascade_prediction_tier"] == "unresolved_no_evidence"
