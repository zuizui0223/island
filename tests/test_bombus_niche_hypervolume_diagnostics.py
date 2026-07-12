from __future__ import annotations

import pandas as pd
import pytest

from island_v2.bombus_niche_hypervolume import score_niche_hypervolumes


def _training() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "bombus_species": "Bombus alpha",
                "bio1": 8.0 + index * 0.2,
                "bio12": 700.0 + index * 12.0 + (index % 3),
            }
            for index in range(30)
        ]
    )


def test_hypervolume_boundary_is_calibrated_to_half_membership() -> None:
    training = _training()
    probe = pd.DataFrame(
        [
            {
                "island_id": "probe",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 880.0,
            }
        ]
    )
    scored = score_niche_hypervolumes(
        training,
        probe,
        environment_columns=["bio1", "bio12"],
    ).iloc[0]

    reconstructed = 2 ** (-float(scored["hypervolume_distance_ratio"]))
    assert scored["environmental_compatibility"] == pytest.approx(reconstructed)
    assert scored["empirical_training_tail_support"] > 0
    assert scored["covariance_condition_number"] >= 1


def test_univariate_extrapolation_is_explicit() -> None:
    targets = pd.DataFrame(
        [
            {
                "island_id": "inside_axes",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 880.0,
            },
            {
                "island_id": "outside_axes",
                "bombus_species": "Bombus alpha",
                "bio1": 30.0,
                "bio12": 2500.0,
            },
        ]
    )
    result = score_niche_hypervolumes(
        _training(),
        targets,
        environment_columns=["bio1", "bio12"],
    ).set_index("island_id")

    assert not bool(result.loc["inside_axes", "environmental_extrapolation"])
    assert result.loc["inside_axes", "univariate_extrapolation_count"] == 0
    assert bool(result.loc["outside_axes", "environmental_extrapolation"])
    assert result.loc["outside_axes", "univariate_extrapolation_count"] == 2
    assert result.loc["outside_axes", "univariate_extrapolation_fraction"] == 1.0
    assert result.loc["outside_axes", "empirical_training_tail_support"] < result.loc[
        "inside_axes", "empirical_training_tail_support"
    ]


def test_unresolved_rows_keep_diagnostic_schema() -> None:
    target = pd.DataFrame(
        [
            {
                "island_id": "unresolved",
                "bombus_species": "Bombus alpha",
                "bio1": 10.0,
                "bio12": 800.0,
            }
        ]
    )
    result = score_niche_hypervolumes(
        _training().head(4),
        target,
        environment_columns=["bio1", "bio12"],
        min_occurrences=10,
    ).iloc[0]

    assert result["model_status"] == "insufficient_occurrences"
    assert pd.isna(result["hypervolume_distance_ratio"])
    assert pd.isna(result["environmental_extrapolation"])
