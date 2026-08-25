import pandas as pd

from island_v2.trait_bombus_analysis import build_analysis_input, normalize_bombus_columns


def _traits() -> pd.DataFrame:
    return pd.DataFrame(
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
            "filled_value": ["white", "red", "selfing", "outcrossing", "unresolved"],
            "fill_tier": [
                "species_direct",
                "family_inference",
                "species_direct",
                "family_inference",
                "unresolved_no_evidence",
            ],
            "analysis_tier": ["broad", "broad", "broad", "broad", "sensitivity_all"],
            "analysis_eligible": ["true", "true", "true", "true", "false"],
        }
    )


def _canonical_bombus() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "island_id": ["a", "b"],
            "n_bombus_species_total": [169, 169],
            "n_bombus_species_scored": [142, 142],
            "n_bombus_species_unscored": [27, 27],
            "n_environmentally_compatible_species": [3, 0],
            "share_environmentally_compatible_species": [3 / 142, 0.0],
            "environmental_compatibility_max": [0.8, 0.2],
            "environmental_compatibility_mean": [0.5, 0.1],
            "environmental_extrapolation_share": [0.4, 0.9],
            "mean_univariate_extrapolation_fraction": [0.2, 0.7],
            "extrapolation_quality": ["low", "high"],
            "bombus_regime_set": [
                "applicable_environment_compatible",
                "structural_bombus_absence_comparison",
            ],
            "primary_bombus_channel_analysis": [True, False],
            "global_comparison_only": [False, True],
        }
    )


def test_build_analysis_input_preserves_regime_traits_and_species_master() -> None:
    island_input, joined, species_master = build_analysis_input(_traits(), _canonical_bombus())

    indexed = island_input.set_index("island_id")
    assert indexed.loc["a", "n_trait_species"] == 2
    assert indexed.loc["b", "global_comparison_only"]
    assert indexed.loc["b", "extrapolation_quality"] == "high"

    color_a = joined.loc[
        (joined["island_id"] == "a")
        & (joined["trait_name"] == "flower_primary_color")
    ]
    assert set(color_a["species_share_within_trait"]) == {0.5}
    assert set(color_a["bombus_regime_set"]) == {"applicable_environment_compatible"}
    assert len(species_master) == 4
    assert "unresolved" not in set(species_master["filled_value"])
    assert species_master["bombus_join_status"].eq("both").all()
    assert "fill_tier" in species_master.columns


def test_legacy_bombus_columns_are_normalized() -> None:
    legacy = _canonical_bombus().rename(
        columns={
            "n_bombus_species_total": "n_candidate_bombus_species",
            "n_bombus_species_scored": "n_scored_bombus_species",
            "n_bombus_species_unscored": "n_unscored_bombus_species",
            "share_environmentally_compatible_species": "environmentally_compatible_species_share",
            "environmental_extrapolation_share": "environmental_extrapolation_species_share",
            "bombus_regime_set": "bombus_regime",
        }
    )
    normalized = normalize_bombus_columns(legacy)
    assert normalized["n_bombus_species_total"].eq(169).all()
    assert "bombus_regime_set" in normalized.columns
