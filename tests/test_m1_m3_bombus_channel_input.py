from __future__ import annotations

import pandas as pd

from island_v2.m1_m3_bombus_channel_input import build_channel_input


def test_build_channel_input_preserves_components_and_evaluability() -> None:
    observation = pd.DataFrame(
        {
            "island_id": ["a", "b", "c"],
            "bombus_occurrence_evidence": ["detected", "adequate_non_detection", "insufficient_effort"],
            "bombus_channel_evaluable": [True, True, False],
            "observation_availability": [0.9, 0.1, 0.5],
            "observation_deficit": [0.1, 0.9, 0.5],
            "observation_effort_quality": [0.2, 0.9, 0.1],
        }
    )
    environmental = pd.DataFrame(
        {
            "island_id": ["a", "b", "c", "d"],
            "environmental_compatibility_max": [0.81, 0.64, 0.49, 0.36],
        }
    )

    result = build_channel_input(observation, environmental).set_index("island_id")

    assert result.loc["a", "environmental_compatibility"] == 0.81
    assert result.loc["a", "observation_availability"] == 0.9
    assert abs(float(result.loc["a", "bombus_channel_availability"]) - 0.853814968) < 1e-6
    assert bool(result.loc["a", "m1_m3_channel_evaluable"])
    assert bool(result.loc["b", "m1_m3_channel_evaluable"])
    assert not bool(result.loc["c", "m1_m3_channel_evaluable"])
    assert result.loc["d", "channel_component_status"] == "missing_observation"
