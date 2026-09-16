from __future__ import annotations

import pandas as pd
import yaml

from island_v2.chapter1_h5_identity_matched_two_route import (
    plant_scores,
    strict_channel_rows,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_h5_identity_matched_two_route_v1.yml", encoding="utf-8"))


def test_strict_channel_rows_requires_symmetric_effort_for_retained() -> None:
    rows = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "channel_id": "bombus",
                "observation_state": "detected",
                "background_record_count": 60,
                "background_spatial_units": 4,
                "background_temporal_units": 3,
                "distinct_dataset_count": 3,
                "latest_background_year": 2020,
            },
            {
                "island_id": "i2",
                "channel_id": "bombus",
                "observation_state": "detected",
                "background_record_count": 10,
                "background_spatial_units": 1,
                "background_temporal_units": 1,
                "distinct_dataset_count": 1,
                "latest_background_year": 2025,
            },
            {
                "island_id": "i3",
                "channel_id": "bombus",
                "observation_state": "adequate_non_detection",
                "background_record_count": 60,
                "background_spatial_units": 4,
                "background_temporal_units": 3,
                "distinct_dataset_count": 3,
                "latest_background_year": 2020,
            },
        ]
    )
    out = strict_channel_rows(rows, _config())
    assert set(out["island_id"]) == {"i1", "i3"}
    assert out.set_index("island_id").loc["i1", "strict_state"] == "retained"
    assert out.set_index("island_id").loc["i3", "strict_state"] == "disrupted"


def test_plant_scores_requires_predeclared_named_axes() -> None:
    rows = []
    for syndrome, value in [
        ("selfing_core", 0.1),
        ("large_bee_like", 0.2),
        ("butterfly_like", 0.3),
        ("bird_like", 0.4),
    ]:
        rows.append(
            {
                "island_id": "i1",
                "stratum": "all_observed",
                "syndrome": syndrome,
                "syndrome_score": value,
            }
        )
    out = plant_scores(pd.DataFrame(rows), _config())
    assert out.loc[0, "selfing_core"] == 0.1
    assert out.loc[0, "large_bee_like"] == 0.2
    assert out.loc[0, "butterfly_like"] == 0.3
    assert out.loc[0, "bird_like"] == 0.4
