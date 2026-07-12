# ruff: noqa: I001
from __future__ import annotations

from pathlib import Path

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


def test_production_workflow_passes_occurrences_to_recency_command() -> None:
    workflow = Path(".github/workflows/run-bombus-occurrence-evidence.yml").read_text(
        encoding="utf-8"
    )
    recency_command = workflow.split("island-v2-bombus-absence-evidence run \\", 1)[1]
    recency_command = recency_command.split("\n\n", 1)[0]

    assert (
        "--occurrences-csv "
        "data/v2/staging/bombus/occurrence_evidence/island_pollinator_occurrences.csv"
        in recency_command
    )
    assert "--diagnostics-csv" not in recency_command


def test_production_workflow_can_fit_environmental_scores_before_channel_score() -> None:
    workflow = Path(".github/workflows/run-bombus-occurrence-evidence.yml").read_text(
        encoding="utf-8"
    )

    assert "occurrence_environment_csv:" in workflow
    assert "island_source_pool_environment_csv:" in workflow
    assert "environment_columns:" in workflow
    assert "island-v2-bombus-niche-hypervolume run" in workflow
    assert (
        "niche_hypervolume/bombus_species_environmental_compatibility.csv" in workflow
    )
