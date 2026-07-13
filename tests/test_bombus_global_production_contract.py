from __future__ import annotations

import pandas as pd

from island_v2.bombus_channel_score import (
    aggregate_environmental_compatibility,
    score_observation_evidence,
)
from island_v2.bombus_regime import classify_regime


def _config() -> dict[str, object]:
    return {
        "primary_thresholds": {
            "min_background_records": 100,
            "min_background_spatial_units": 10,
            "min_background_temporal_units": 5,
            "min_distinct_datasets": 2,
        },
        "observation_recency": {"max_years_since_latest_background": 10},
        "continuous_score": {"bombus_record_scale": 3.0},
    }


def test_detection_is_positive_even_under_low_background_effort() -> None:
    diagnostics = pd.DataFrame(
        {
            "island_id": ["island_a"],
            "bombus_record_count": [1],
            "target_group_record_count": [0],
            "target_group_spatial_units": [0],
            "target_group_temporal_units": [0],
            "distinct_dataset_count": [0],
            "latest_record_age_years": [10],
            "bombus_occurrence_evidence": ["detected"],
        }
    )
    scored = score_observation_evidence(diagnostics, _config())
    assert float(scored.loc[0, "observation_availability"]) > 0.5


def test_environmental_aggregation_retains_unscored_source_pool_species() -> None:
    scores = pd.DataFrame(
        {
            "island_id": ["island_a", "island_a", "island_a"],
            "bombus_species": ["Bombus alpha", "Bombus beta", "Bombus gamma"],
            "environmental_compatibility": [0.25, 0.80, pd.NA],
        }
    )
    aggregated = aggregate_environmental_compatibility(scores)
    assert float(aggregated.loc[0, "environmental_compatibility"]) == 0.80
    assert int(aggregated.loc[0, "n_source_pool_bombus_species"]) == 3
    assert int(aggregated.loc[0, "n_scored_source_pool_bombus_species"]) == 2
    assert int(aggregated.loc[0, "n_unscored_source_pool_bombus_species"]) == 1
    assert aggregated.loc[0, "environmental_compatibility_status"] == "scored_partial"


def test_applicability_is_not_redefined_by_environmental_score() -> None:
    assert (
        classify_regime("structurally_not_applicable", 0.99)
        == "structural_bombus_absence_comparison"
    )
    assert classify_regime("applicable", 0.99) == "applicable_environment_compatible"
    assert classify_regime("applicable", 0.01) == "applicable_environment_mismatch"
