from __future__ import annotations

import pandas as pd
import yaml

from island_v2.chapter1_h5_multichannel_two_route import (
    build_composite_channel_exposure,
    build_plant_scores,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_h5_multichannel_two_route_v1.yml", encoding="utf-8"))


def test_symmetric_effort_gate_blocks_low_effort_detection() -> None:
    rows = [
        {
            "island_id": "i1",
            "channel_id": "bombus",
            "observation_state": "detected",
            "background_record_count": 10,
            "background_spatial_units": 1,
            "background_temporal_units": 1,
            "distinct_dataset_count": 1,
            "latest_background_year": 2025,
        },
        {
            "island_id": "i1",
            "channel_id": "lepidoptera",
            "observation_state": "adequate_non_detection",
            "background_record_count": 60,
            "background_spatial_units": 4,
            "background_temporal_units": 3,
            "distinct_dataset_count": 3,
            "latest_background_year": 2025,
        },
        {
            "island_id": "i2",
            "channel_id": "bombus",
            "observation_state": "detected",
            "background_record_count": 80,
            "background_spatial_units": 5,
            "background_temporal_units": 4,
            "distinct_dataset_count": 3,
            "latest_background_year": 2024,
        },
    ]
    composite, long = build_composite_channel_exposure(pd.DataFrame(rows), _config())
    i1 = composite.set_index("island_id").loc["i1"]
    assert i1["n_evaluable_channels"] == 1
    assert i1["n_disrupted_channels"] == 1
    assert i1["any_channel_disrupted"] == 1
    low = long.loc[(long["island_id"] == "i1") & (long["channel_id"] == "bombus")].iloc[0]
    assert not bool(low["retained_symmetric"])
    i2 = composite.set_index("island_id").loc["i2"]
    assert i2["n_retained_channels"] == 1


def test_plant_scores_build_two_route_axes() -> None:
    values = {
        "selfing_core": 0.4,
        "selfing_syndrome": 0.2,
        "generalized_accessible": 0.6,
        "large_bee_like": -0.2,
        "butterfly_like": 0.4,
        "bird_like": 0.2,
    }
    rows = [
        {
            "island_id": "i1",
            "stratum": "all_observed",
            "syndrome": syndrome,
            "syndrome_score": value,
        }
        for syndrome, value in values.items()
    ]
    out = build_plant_scores(pd.DataFrame(rows), "all_observed").set_index("island_id").loc["i1"]
    assert abs(out["attraction_shift"] - 0.4) < 1e-12
    assert abs(out["shared_named_architecture"] - ((-0.2 + 0.4 + 0.2) / 3.0)) < 1e-12
