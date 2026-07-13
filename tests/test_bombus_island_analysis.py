from __future__ import annotations

import pandas as pd

from island_v2.bombus_island_analysis import aggregate_islands


def test_aggregate_islands_preserves_regime_coverage_and_extrapolation_quality() -> None:
    scores = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2"],
            "bombus_species": ["Bombus alpha", "Bombus beta", "Bombus alpha"],
            "environmental_compatibility": [0.8, None, 0.2],
            "environmental_extrapolation": [False, True, True],
            "univariate_extrapolation_fraction": [2 / 9, None, 7 / 9],
            "bombus_regime": [
                "applicable_environment_compatible",
                "applicable_environment_unresolved",
                "structural_bombus_absence_comparison",
            ],
            "primary_bombus_channel_analysis": [True, True, False],
            "global_comparison_only": [False, False, True],
        }
    )

    result = aggregate_islands(scores).set_index("island_id")

    assert result.loc["i1", "n_bombus_species_total"] == 2
    assert result.loc["i1", "n_bombus_species_scored"] == 1
    assert result.loc["i1", "n_bombus_species_unscored"] == 1
    assert result.loc["i1", "n_environmentally_compatible_species"] == 1
    assert result.loc["i1", "environmental_compatibility_max"] == 0.8
    assert result.loc["i1", "mean_univariate_extrapolation_fraction"] == 2 / 9
    assert result.loc["i1", "extrapolation_quality"] == "low"
    assert result.loc["i1", "primary_bombus_channel_analysis"]
    assert not result.loc["i1", "global_comparison_only"]
    assert result.loc["i2", "extrapolation_quality"] == "high"
    assert result.loc["i2", "global_comparison_only"]
