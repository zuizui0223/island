from __future__ import annotations

import pandas as pd
import pytest
import typer

from island_v2.bombus_niche_hypervolume import score_niche_hypervolumes


def _training() -> pd.DataFrame:
    rows = []
    for index in range(20):
        rows.append(
            {
                "bombus_species": "Bombus alpha",
                "bio1": 10.0 + index * 0.1,
                "bio12": 900.0 + index * 5.0,
            }
        )
    return pd.DataFrame(rows)


def test_species_hypervolume_scores_near_environment_higher_than_far() -> None:
    islands = pd.DataFrame(
        [
            {
                "island_id": "near",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 950.0,
            },
            {
                "island_id": "far",
                "bombus_species": "Bombus alpha",
                "bio1": 30.0,
                "bio12": 2000.0,
            },
        ]
    )
    result = score_niche_hypervolumes(
        _training(),
        islands,
        environment_columns=["bio1", "bio12"],
        min_occurrences=10,
    ).set_index("island_id")

    assert result.loc["near", "model_status"] == "scored"
    assert result.loc["far", "model_status"] == "scored"
    assert result.loc["near", "environmental_compatibility"] > result.loc[
        "far", "environmental_compatibility"
    ]
    assert bool(result.loc["near", "inside_hypervolume"])
    assert not bool(result.loc["far", "inside_hypervolume"])


def test_insufficient_species_occurrences_remain_unresolved() -> None:
    occurrences = _training().head(5)
    islands = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 950.0,
            }
        ]
    )
    result = score_niche_hypervolumes(
        occurrences,
        islands,
        environment_columns=["bio1", "bio12"],
        min_occurrences=10,
    ).iloc[0]

    assert result["model_status"] == "insufficient_occurrences"
    assert pd.isna(result["environmental_compatibility"])


def test_missing_island_environment_is_not_forced_to_zero() -> None:
    islands = pd.DataFrame(
        [
            {
                "island_id": "missing",
                "bombus_species": "Bombus alpha",
                "bio1": "",
                "bio12": 950.0,
            }
        ]
    )
    result = score_niche_hypervolumes(
        _training(),
        islands,
        environment_columns=["bio1", "bio12"],
        min_occurrences=10,
    ).iloc[0]

    assert result["model_status"] == "missing_island_environment"
    assert pd.isna(result["environmental_compatibility"])


def test_scores_are_invariant_to_environment_variable_units() -> None:
    occurrences = _training()
    islands = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "bombus_species": "Bombus alpha",
                "bio1": 11.4,
                "bio12": 970.0,
            }
        ]
    )
    baseline = score_niche_hypervolumes(
        occurrences,
        islands,
        environment_columns=["bio1", "bio12"],
    ).iloc[0]

    rescaled_occurrences = occurrences.copy()
    rescaled_islands = islands.copy()
    rescaled_occurrences["bio1"] *= 1_000_000
    rescaled_islands["bio1"] *= 1_000_000
    rescaled = score_niche_hypervolumes(
        rescaled_occurrences,
        rescaled_islands,
        environment_columns=["bio1", "bio12"],
    ).iloc[0]

    assert rescaled["environmental_compatibility"] == pytest.approx(
        baseline["environmental_compatibility"]
    )
    assert rescaled["mahalanobis_d2"] == pytest.approx(baseline["mahalanobis_d2"])


def test_constant_environment_dimension_is_reported_and_dropped() -> None:
    occurrences = _training()
    occurrences["constant"] = 1.0
    islands = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 950.0,
                "constant": "",
            }
        ]
    )
    result = score_niche_hypervolumes(
        occurrences,
        islands,
        environment_columns=["bio1", "bio12", "constant"],
    ).iloc[0]

    assert result["model_status"] == "scored"
    assert result["n_environmental_dimensions"] == 2
    assert result["n_environmental_dimensions_requested"] == 3
    assert result["dropped_environmental_dimensions"] == "constant"


def test_duplicate_island_species_targets_are_rejected() -> None:
    islands = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "bombus_species": "Bombus alpha",
                "bio1": 11.0,
                "bio12": 950.0,
            },
            {
                "island_id": "i1",
                "bombus_species": "Bombus alpha",
                "bio1": 11.2,
                "bio12": 960.0,
            },
        ]
    )

    with pytest.raises(typer.BadParameter, match="one row per"):
        score_niche_hypervolumes(
            _training(),
            islands,
            environment_columns=["bio1", "bio12"],
        )
