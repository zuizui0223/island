import pandas as pd

from island_v2.trait_bombus_analysis import build_analysis_input


def test_build_analysis_input_preserves_regime_and_trait_shares() -> None:
    traits = pd.DataFrame(
        {
            "island_id": ["a", "a", "a", "a", "b"],
            "accepted_species": ["sp1", "sp2", "sp1", "sp2", "sp3"],
            "trait_name": [
                "flower_primary_color",
                "flower_primary_color",
                "mating_system",
                "mating_system",
                "flower_primary_color",
            ],
            "filled_value": ["white", "red", "selfing", "outcrossing", "white"],
        }
    )
    bombus = pd.DataFrame(
        {
            "island_id": ["a", "b"],
            "n_candidate_bombus_species": [3, 2],
            "n_scored_bombus_species": [2, 2],
            "n_unscored_bombus_species": [1, 0],
            "n_environmentally_compatible_species": [1, 0],
            "environmentally_compatible_species_share": [0.5, 0.0],
            "environmental_compatibility_max": [0.8, 0.2],
            "environmental_compatibility_mean": [0.5, 0.1],
            "environmental_extrapolation_species_share": [0.0, 0.5],
            "bombus_regime": [
                "applicable_environment_compatible",
                "structural_bombus_absence_comparison",
            ],
            "primary_bombus_channel_analysis": [True, False],
            "global_comparison_only": [False, True],
        }
    )

    island_input, joined = build_analysis_input(traits, bombus)

    assert island_input.set_index("island_id").loc["a", "n_trait_species"] == 2
    assert island_input.set_index("island_id").loc["b", "global_comparison_only"]
    color_a = joined.loc[
        (joined["island_id"] == "a")
        & (joined["trait_name"] == "flower_primary_color")
    ]
    assert set(color_a["species_share_within_trait"]) == {0.5}
    assert set(color_a["bombus_regime"]) == {"applicable_environment_compatible"}
