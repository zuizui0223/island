from __future__ import annotations

import pandas as pd

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
