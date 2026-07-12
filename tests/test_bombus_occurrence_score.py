from __future__ import annotations

import pandas as pd

from island_v2.bombus_occurrence_score import score_bombus_occurrence


CONFIG = {
    "reference_year": 2026,
    "thresholds": {
        "bombus_records_saturation": 3,
        "target_group_records": 50,
        "spatial_units": 3,
        "temporal_units": 2,
        "distinct_datasets": 2,
        "max_years_since_latest_background": 20,
    },
}


def _row(island_id: str, bombus: int, effort: int, year: int | str) -> dict[str, object]:
    return {
        "island_id": island_id,
        "bombus_record_count": bombus,
        "target_group_record_count": effort,
        "target_group_spatial_units": 3 if effort else 0,
        "target_group_temporal_units": 2 if effort else 0,
        "distinct_dataset_count": 2 if effort else 0,
        "latest_record_year": year,
    }


def test_zero_detection_with_strong_recent_effort_has_high_deficit() -> None:
    scores = score_bombus_occurrence(pd.DataFrame([_row("x", 0, 50, 2026)]), CONFIG).iloc[0]
    assert scores["effort_strength"] == 1.0
    assert scores["recency_strength"] == 1.0
    assert scores["bombus_occurrence_availability_score"] == 0.0
    assert scores["bombus_occurrence_deficit_score"] == 1.0
    assert scores["measurement_certainty"] == 1.0


def test_zero_detection_without_effort_stays_neutral_and_uncertain() -> None:
    scores = score_bombus_occurrence(pd.DataFrame([_row("x", 0, 0, "")]), CONFIG).iloc[0]
    assert scores["bombus_occurrence_availability_score"] == 0.5
    assert scores["bombus_occurrence_deficit_score"] == 0.5
    assert scores["measurement_certainty"] == 0.0


def test_more_bombus_records_increase_availability() -> None:
    diagnostics = pd.DataFrame(
        [
            _row("one", 1, 50, 2026),
            _row("many", 9, 50, 2026),
        ]
    )
    scores = score_bombus_occurrence(diagnostics, CONFIG).set_index("island_id")
    assert (
        scores.loc["many", "bombus_occurrence_availability_score"]
        > scores.loc["one", "bombus_occurrence_availability_score"]
    )
    assert scores.loc["many", "bombus_occurrence_deficit_score"] < scores.loc[
        "one", "bombus_occurrence_deficit_score"
    ]
