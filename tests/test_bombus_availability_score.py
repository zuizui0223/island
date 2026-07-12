from __future__ import annotations

import pandas as pd

from island_v2.bombus_availability_score import score_occurrence_evidence


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
}


def test_zero_records_with_weak_effort_remain_uncertain() -> None:
    frame = pd.DataFrame(
        [{
            "island_id": "weak",
            "bombus_record_count": 0,
            "target_group_record_count": 0,
            "target_group_spatial_units": 0,
            "target_group_temporal_units": 0,
            "distinct_dataset_count": 0,
            "latest_record_age_years": "",
        }]
    )
    scored = score_occurrence_evidence(frame, CONFIG).iloc[0]
    assert scored["bombus_occurrence_availability_score"] == 0.5
    assert scored["bombus_occurrence_deficit_score"] == 0.5
    assert scored["bombus_occurrence_score_confidence"] == 0.0


def test_zero_records_with_strong_recent_effort_support_deficit() -> None:
    frame = pd.DataFrame(
        [{
            "island_id": "strong_zero",
            "bombus_record_count": 0,
            "target_group_record_count": 100,
            "target_group_spatial_units": 5,
            "target_group_temporal_units": 4,
            "distinct_dataset_count": 3,
            "latest_record_age_years": 0,
        }]
    )
    scored = score_occurrence_evidence(frame, CONFIG).iloc[0]
    assert scored["bombus_observation_coverage"] == 1.0
    assert scored["bombus_recency_support"] == 1.0
    assert scored["bombus_occurrence_availability_score"] == 0.0
    assert scored["bombus_occurrence_deficit_score"] == 1.0
    assert scored["bombus_occurrence_score_confidence"] == 1.0


def test_detections_move_score_above_half() -> None:
    frame = pd.DataFrame(
        [
            {
                "island_id": "one",
                "bombus_record_count": 1,
                "target_group_record_count": 50,
                "target_group_spatial_units": 3,
                "target_group_temporal_units": 2,
                "distinct_dataset_count": 2,
                "latest_record_age_years": 10,
            },
            {
                "island_id": "many",
                "bombus_record_count": 10,
                "target_group_record_count": 50,
                "target_group_spatial_units": 3,
                "target_group_temporal_units": 2,
                "distinct_dataset_count": 2,
                "latest_record_age_years": 10,
            },
        ]
    )
    scored = score_occurrence_evidence(frame, CONFIG).set_index("island_id")
    assert 0.5 < scored.loc["one", "bombus_occurrence_availability_score"] < 1.0
    assert scored.loc["many", "bombus_occurrence_availability_score"] > scored.loc[
        "one", "bombus_occurrence_availability_score"
    ]
