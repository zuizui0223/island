import pandas as pd

from analysis.build_secondary_probabilistic_genus_low_20260809 import (
    _incremental_predictions,
    _prediction_metrics,
)


def test_probabilistic_low_remains_trait_specific_and_non_strict() -> None:
    predictions = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Alpha one", "Alpha two"],
            "genus": ["Alpha", "Alpha", "Alpha"],
            "axis": [
                "floral_structural_complexity",
                "floral_structural_complexity",
                "flower_colour",
            ],
            "trait_name": ["floral_form", "tube_depth_class", "flower_primary_color"],
            "dominance": [0.8, 0.75, 0.8],
        }
    )
    master = pd.DataFrame(
        {"accepted_species": ["Alpha one", "Alpha two"], "n_islands": [9, 2]}
    )
    result = _incremental_predictions(
        predictions,
        {("Alpha one", "floral_form")},
        {("Alpha one", "floral_structural_complexity")},
        master,
        role="secondary_multiple_imputation_relaxed_min3",
    )

    assert set(zip(result["accepted_species"], result["trait_name"], strict=True)) == {
        ("Alpha one", "tube_depth_class"),
        ("Alpha two", "flower_primary_color"),
    }
    assert not result.loc[
        result["accepted_species"].eq("Alpha one"), "potential_new_species_axis"
    ].item()
    assert result["join_key"].eq("genus_x_trait_name").all()
    assert result["quality"].eq("probabilistic_low").all()
    assert not result["strict_coverage_eligible"].any()
    assert not result["family_inference_used"].any()
    assert not result["global_fallback_used"].any()
    assert _prediction_metrics(result)["potential_new_species_axis"] == 1
