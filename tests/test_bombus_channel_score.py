# ruff: noqa: I001
from __future__ import annotations

import pandas as pd
import pytest
import typer

from island_v2 import bombus_channel_score as channel_score


CONFIG = {
    "primary_thresholds": {
        "min_background_records": 50,
        "min_background_spatial_units": 3,
        "min_background_temporal_units": 2,
        "min_distinct_datasets": 2,
    },
    "observation_recency": {
        "reference_year": 2026,
        "max_years_since_latest_background": 20,
    },
    "continuous_score": {"bombus_record_scale": 3.0},
}


def test_low_effort_zero_shrinks_to_neutral() -> None:
    diagnostics = pd.DataFrame(
        [
            {
                "island_id": "low_effort",
                "bombus_record_count": 0,
                "target_group_record_count": 0,
                "target_group_spatial_units": 0,
                "target_group_temporal_units": 0,
                "distinct_dataset_count": 0,
                "latest_record_age_years": "",
                "bombus_occurrence_evidence": "insufficient_effort",
            }
        ]
    )
    result = channel_score.score_observation_evidence(diagnostics, CONFIG).iloc[0]
    assert result["observation_effort_quality"] == 0
    assert result["observation_availability"] == 0.5


def test_strong_effort_zero_moves_toward_deficit() -> None:
    diagnostics = pd.DataFrame(
        [
            {
                "island_id": "searched_zero",
                "bombus_record_count": 0,
                "target_group_record_count": 100,
                "target_group_spatial_units": 5,
                "target_group_temporal_units": 5,
                "distinct_dataset_count": 3,
                "latest_record_age_years": 0,
                "bombus_occurrence_evidence": "adequate_non_detection",
            }
        ]
    )
    result = channel_score.score_observation_evidence(diagnostics, CONFIG).iloc[0]
    assert result["observation_effort_quality"] == 1
    assert result["observation_availability"] == 0
    assert result["observation_deficit"] == 1


def test_single_detection_is_positive_evidence_regardless_of_background_effort() -> None:
    diagnostics = pd.DataFrame(
        [
            {
                "island_id": "strong_background",
                "bombus_record_count": 1,
                "target_group_record_count": 100,
                "target_group_spatial_units": 5,
                "target_group_temporal_units": 5,
                "distinct_dataset_count": 3,
                "latest_record_age_years": 0,
                "bombus_occurrence_evidence": "detected",
            },
            {
                "island_id": "sparse_background",
                "bombus_record_count": 1,
                "target_group_record_count": 1,
                "target_group_spatial_units": 1,
                "target_group_temporal_units": 1,
                "distinct_dataset_count": 1,
                "latest_record_age_years": 20,
                "bombus_occurrence_evidence": "detected",
            },
        ]
    )
    result = channel_score.score_observation_evidence(diagnostics, CONFIG).set_index(
        "island_id"
    )
    assert result.loc["strong_background", "observation_availability"] > 0.5
    assert result.loc["strong_background", "observation_availability"] == result.loc[
        "sparse_background", "observation_availability"
    ]
    assert result.loc["strong_background", "observation_deficit"] < 0.5


def test_repeated_detections_increase_availability() -> None:
    diagnostics = pd.DataFrame(
        [
            {
                "island_id": "detected",
                "bombus_record_count": 6,
                "target_group_record_count": 100,
                "target_group_spatial_units": 5,
                "target_group_temporal_units": 5,
                "distinct_dataset_count": 3,
                "latest_record_age_years": 0,
                "bombus_occurrence_evidence": "detected",
            }
        ]
    )
    result = channel_score.score_observation_evidence(diagnostics, CONFIG).iloc[0]
    assert result["observation_availability"] > 0.9


def test_source_pool_environmental_scores_keep_multiple_summaries() -> None:
    scores = pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "bombus_species": ["Bombus a", "Bombus b"],
            "environmental_compatibility": [0.2, 0.5],
        }
    )
    result = channel_score.aggregate_environmental_compatibility(scores).iloc[0]
    assert result["environmental_compatibility"] == 0.5
    assert result["environmental_compatibility_max"] == 0.5
    assert result["environmental_compatibility_mean"] == 0.35
    assert round(result["environmental_compatibility_noisy_or"], 3) == 0.6


def test_unscored_source_pool_species_remain_in_environmental_audit() -> None:
    scores = pd.DataFrame(
        {
            "island_id": ["partial", "partial", "unresolved"],
            "bombus_species": ["Bombus a", "Bombus b", "Bombus c"],
            "environmental_compatibility": [0.4, "", ""],
        }
    )
    result = channel_score.aggregate_environmental_compatibility(scores).set_index(
        "island_id"
    )

    assert result.loc["partial", "n_source_pool_bombus_species"] == 2
    assert result.loc["partial", "n_scored_source_pool_bombus_species"] == 1
    assert result.loc["partial", "environmental_compatibility_status"] == "scored_partial"
    assert result.loc["unresolved", "environmental_compatibility_status"] == (
        "no_scored_species"
    )
    assert pd.isna(result.loc["unresolved", "environmental_compatibility"])


def test_duplicate_island_species_environmental_scores_are_rejected() -> None:
    scores = pd.DataFrame(
        {
            "island_id": ["i1", "i1"],
            "bombus_species": ["Bombus a", "Bombus a"],
            "environmental_compatibility": [0.2, 0.5],
        }
    )
    with pytest.raises(typer.BadParameter, match="one row per"):
        channel_score.aggregate_environmental_compatibility(scores)


def test_combined_score_is_secondary_geometric_mean() -> None:
    observation = pd.DataFrame(
        {
            "island_id": ["i1"],
            "observation_availability": [0.25],
            "observation_deficit": [0.75],
        }
    )
    environmental = pd.DataFrame(
        {
            "island_id": ["i1"],
            "environmental_compatibility": [1.0],
            "n_source_pool_bombus_species": [1],
            "environmental_compatibility_max": [1.0],
            "environmental_compatibility_mean": [1.0],
        }
    )
    result = channel_score.combine_channel_components(observation, environmental).iloc[0]
    assert result["bombus_channel_availability"] == 0.5
    assert result["bombus_channel_deficit"] == 0.5
