from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import typer

from island_v2.bombus_niche_hypervolume import score_niche_hypervolumes


def _training(n: int = 100) -> pd.DataFrame:
    rng = np.random.default_rng(123)
    x = rng.normal(size=n)
    return pd.DataFrame(
        {
            "bombus_species": ["Bombus alpha"] * n,
            "bio1": 10 + x,
            "bio5": 20 + 2 * x + rng.normal(scale=0.02, size=n),
            "bio12": 900 + rng.normal(scale=80, size=n),
        }
    )


def test_species_ellipsoidal_niche_scores_near_environment_higher_than_far() -> None:
    islands = pd.DataFrame(
        [
            {"island_id": "near", "bombus_species": "Bombus alpha", "bio1": 10.0, "bio5": 20.0, "bio12": 900.0},
            {"island_id": "far", "bombus_species": "Bombus alpha", "bio1": 30.0, "bio5": 60.0, "bio12": 2000.0},
        ]
    )
    result = score_niche_hypervolumes(
        _training(), islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50
    ).set_index("island_id")
    assert result.loc["near", "model_status"] == "scored"
    assert result.loc["far", "model_status"] == "scored"
    assert result.loc["near", "environmental_compatibility"] > result.loc["far", "environmental_compatibility"]
    assert bool(result.loc["near", "inside_ellipsoidal_envelope"])
    assert not bool(result.loc["far", "inside_ellipsoidal_envelope"])
    assert result.loc["near", "model_type"] == "regularized_ellipsoidal_environmental_niche"


def test_insufficient_species_occurrences_remain_unresolved() -> None:
    occurrences = _training().head(10)
    islands = pd.DataFrame(
        [{"island_id": "i1", "bombus_species": "Bombus alpha", "bio1": 10.0, "bio5": 20.0, "bio12": 900.0}]
    )
    result = score_niche_hypervolumes(
        occurrences, islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50
    ).iloc[0]
    assert result["model_status"] == "insufficient_occurrences_after_quality_filtering"
    assert pd.isna(result["environmental_compatibility"])


def test_missing_island_environment_is_not_forced_to_zero() -> None:
    islands = pd.DataFrame(
        [{"island_id": "missing", "bombus_species": "Bombus alpha", "bio1": "", "bio5": 20.0, "bio12": 900.0}]
    )
    result = score_niche_hypervolumes(
        _training(), islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50
    ).iloc[0]
    assert result["model_status"] == "missing_island_environment"
    assert pd.isna(result["environmental_compatibility"])


def test_scores_are_invariant_to_environment_variable_units() -> None:
    occurrences = _training()
    islands = pd.DataFrame(
        [{"island_id": "i1", "bombus_species": "Bombus alpha", "bio1": 10.2, "bio5": 20.4, "bio12": 920.0}]
    )
    baseline = score_niche_hypervolumes(
        occurrences, islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50
    ).iloc[0]
    rescaled_occurrences = occurrences.copy()
    rescaled_islands = islands.copy()
    rescaled_occurrences["bio1"] *= 1_000_000
    rescaled_islands["bio1"] *= 1_000_000
    rescaled = score_niche_hypervolumes(
        rescaled_occurrences, rescaled_islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50
    ).iloc[0]
    assert rescaled["environmental_compatibility"] == pytest.approx(baseline["environmental_compatibility"])
    assert rescaled["mahalanobis_d2"] == pytest.approx(baseline["mahalanobis_d2"])


def test_constant_environment_dimension_is_reported_and_dropped() -> None:
    occurrences = _training()
    occurrences["constant"] = 1.0
    islands = pd.DataFrame(
        [{"island_id": "i1", "bombus_species": "Bombus alpha", "bio1": 10.0, "bio5": 20.0, "bio12": 900.0, "constant": ""}]
    )
    result = score_niche_hypervolumes(
        occurrences, islands, environment_columns=["bio1", "bio5", "bio12", "constant"], min_occurrences=50
    ).iloc[0]
    assert result["model_status"] == "scored"
    assert result["n_environmental_dimensions_active"] == 3
    assert result["n_environmental_dimensions_requested"] == 4
    assert result["dropped_environmental_dimensions"] == "constant"


def test_correlated_dimensions_are_reduced_by_pca() -> None:
    islands = pd.DataFrame(
        [{"island_id": "i1", "bombus_species": "Bombus alpha", "bio1": 10.0, "bio5": 20.0, "bio12": 900.0}]
    )
    result = score_niche_hypervolumes(
        _training(), islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50
    ).iloc[0]
    assert int(result["n_pca_dimensions"]) < int(result["n_environmental_dimensions_active"])


def test_duplicate_island_species_targets_are_rejected() -> None:
    islands = pd.DataFrame(
        [
            {"island_id": "i1", "bombus_species": "Bombus alpha", "bio1": 10.0, "bio5": 20.0, "bio12": 900.0},
            {"island_id": "i1", "bombus_species": "Bombus alpha", "bio1": 10.1, "bio5": 20.2, "bio12": 910.0},
        ]
    )
    with pytest.raises(typer.BadParameter, match="one row per"):
        score_niche_hypervolumes(_training(), islands, environment_columns=["bio1", "bio5", "bio12"], min_occurrences=50)
