# ruff: noqa: I001
from __future__ import annotations

import pandas as pd

from island_v2.bombus_absence_evidence import apply_recency_gate


CONFIG = {
    "observation_recency": {
        "reference_year": 2026,
        "max_years_since_latest_background": 20,
    }
}


def _base_rows() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "island_id": "recent_zero",
                "bombus_occurrence_evidence": "adequate_non_detection",
                "bombus_channel_evaluable": True,
                "latest_record_year": 2024,
                "observation_diagnostic_flags": "",
                "occupancy_sensitivity_eligible": True,
            },
            {
                "island_id": "stale_zero",
                "bombus_occurrence_evidence": "adequate_non_detection",
                "bombus_channel_evaluable": True,
                "latest_record_year": 1990,
                "observation_diagnostic_flags": "",
                "occupancy_sensitivity_eligible": True,
            },
            {
                "island_id": "old_detection",
                "bombus_occurrence_evidence": "detected",
                "bombus_channel_evaluable": True,
                "latest_record_year": 1990,
                "observation_diagnostic_flags": "",
                "occupancy_sensitivity_eligible": False,
            },
        ]
    )


def test_recency_gate_only_downgrades_non_detection() -> None:
    result = apply_recency_gate(_base_rows(), CONFIG).set_index("island_id")

    assert result.loc["recent_zero", "bombus_occurrence_evidence"] == "adequate_non_detection"
    assert bool(result.loc["recent_zero", "background_recent_enough"])
    assert result.loc["recent_zero", "latest_record_age_years"] == 2

    assert result.loc["stale_zero", "bombus_occurrence_evidence"] == "insufficient_effort"
    assert not bool(result.loc["stale_zero", "bombus_channel_evaluable"])
    assert "stale_background_effort" in result.loc[
        "stale_zero", "observation_diagnostic_flags"
    ]

    assert result.loc["old_detection", "bombus_occurrence_evidence"] == "detected"
    assert bool(result.loc["old_detection", "bombus_channel_evaluable"])
    assert "stale_background_effort" in result.loc[
        "old_detection", "observation_diagnostic_flags"
    ]
