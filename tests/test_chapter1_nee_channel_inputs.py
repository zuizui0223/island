from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.chapter1_nee_channel_inputs import (
    OBSERVATION_COLUMNS,
    SOURCE_COLUMNS,
    adapt_bombus_observation,
    adapt_bombus_source_availability,
    load_config,
    project_channel_states,
    qualification_receipt,
    validate_island_observation,
    validate_source_availability,
)


CONFIG = Path("config/chapter1_nee_channel_inputs.yml")


def _config() -> dict[str, object]:
    return load_config(CONFIG)


def _source(rows: list[dict[str, object]]) -> pd.DataFrame:
    defaults = {
        "source_region_id": "src-1",
        "channel_id": "bombus",
        "source_state": "available",
        "evidence_id": "ev-1",
        "evidence_type": "distribution_review",
        "source_citation": "Source citation",
        "source_url": "https://example.org/source",
        "review_status": "accepted",
    }
    return pd.DataFrame([{**defaults, **row} for row in rows], columns=SOURCE_COLUMNS)


def _observation(rows: list[dict[str, object]]) -> pd.DataFrame:
    defaults = {
        "channel_id": "bombus",
        "observation_state": "detected",
        "channel_record_count": 1,
        "background_record_count": 60,
        "background_spatial_units": 4,
        "background_temporal_units": 3,
        "distinct_dataset_count": 2,
        "latest_background_year": 2025,
        "evidence_source": "mock_occurrence",
        "quality_flags": "",
    }
    return pd.DataFrame([{**defaults, **row} for row in rows], columns=OBSERVATION_COLUMNS)


def test_available_detected_projects_to_retained_even_with_sparse_background() -> None:
    config = _config()
    source = _source([{"island_id": "i1"}])
    observation = _observation(
        [
            {
                "island_id": "i1",
                "observation_state": "detected",
                "channel_record_count": 1,
                "background_record_count": 1,
                "background_spatial_units": 1,
                "background_temporal_units": 1,
                "distinct_dataset_count": 1,
            }
        ]
    )
    projected = project_channel_states(source, observation, config)
    assert projected.loc[0, "channel_state"] == "retained"


def test_available_adequate_non_detection_projects_to_disrupted() -> None:
    config = _config()
    source = _source([{"island_id": "i1"}])
    observation = _observation(
        [
            {
                "island_id": "i1",
                "observation_state": "adequate_non_detection",
                "channel_record_count": 0,
            }
        ]
    )
    projected = project_channel_states(source, observation, config)
    assert projected.loc[0, "channel_state"] == "disrupted"


def test_available_insufficient_effort_projects_to_unresolved_not_disrupted() -> None:
    config = _config()
    source = _source([{"island_id": "i1"}])
    observation = _observation(
        [
            {
                "island_id": "i1",
                "observation_state": "insufficient_effort",
                "channel_record_count": 0,
                "background_record_count": 2,
            }
        ]
    )
    projected = project_channel_states(source, observation, config)
    assert projected.loc[0, "channel_state"] == "unresolved"
    assert "observation_not_evaluable" in projected.loc[0, "projection_flags"]


def test_structural_absence_precedes_detected_island_record_and_is_flagged() -> None:
    config = _config()
    source = _source([{"island_id": "i1", "source_state": "structurally_absent"}])
    observation = _observation([{"island_id": "i1"}])
    projected = project_channel_states(source, observation, config)
    assert projected.loc[0, "channel_state"] == "structurally_absent"
    assert "detected_despite_structural_source_absence" in projected.loc[0, "projection_flags"]


def test_missing_island_observation_never_becomes_structural_absence() -> None:
    config = _config()
    source = _source([{"island_id": "i1"}])
    observation = pd.DataFrame(columns=OBSERVATION_COLUMNS)
    projected = project_channel_states(source, observation, config)
    assert projected.loc[0, "source_state"] == "available"
    assert projected.loc[0, "channel_state"] == "unresolved"


def test_decisive_source_state_requires_accepted_review() -> None:
    config = _config()
    source = _source([{"island_id": "i1", "review_status": "pending"}])
    try:
        validate_source_availability(source, config)
    except ValueError as exc:
        assert "require accepted review" in str(exc)
    else:
        raise AssertionError("pending decisive source state should fail validation")


def test_adequate_non_detection_requires_zero_channel_records() -> None:
    config = _config()
    observation = _observation(
        [
            {
                "island_id": "i1",
                "observation_state": "adequate_non_detection",
                "channel_record_count": 1,
            }
        ]
    )
    try:
        validate_island_observation(observation, config)
    except ValueError as exc:
        assert "channel_record_count == 0" in str(exc)
    else:
        raise AssertionError("nonzero adequate non-detection should fail validation")


def test_support_receipt_requires_retained_disrupted_contrast() -> None:
    config = _config()
    source_rows = [
        {"island_id": f"i{x}", "source_region_id": f"src-{1 + (x % 2)}", "evidence_id": f"ev-{x}"}
        for x in range(60)
    ]
    observation_rows = [{"island_id": f"i{x}"} for x in range(60)]
    projected = project_channel_states(_source(source_rows), _observation(observation_rows), config)
    receipt = qualification_receipt(projected, config)
    bombus = receipt.loc[receipt["channel_id"].eq("bombus")].iloc[0]
    assert bombus["n_retained"] == 60
    assert bombus["n_disrupted"] == 0
    assert bombus["support_tier"] == "not_qualified"
    assert "insufficient_disrupted" in bombus["exclusion_reason"]


def test_support_receipt_reaches_confirmatory_with_predeclared_contrast() -> None:
    config = _config()
    source_rows = [
        {"island_id": f"i{x}", "source_region_id": f"src-{1 + (x % 2)}", "evidence_id": f"ev-{x}"}
        for x in range(60)
    ]
    observation_rows = []
    for x in range(60):
        if x < 10:
            observation_rows.append(
                {
                    "island_id": f"i{x}",
                    "observation_state": "adequate_non_detection",
                    "channel_record_count": 0,
                }
            )
        else:
            observation_rows.append({"island_id": f"i{x}"})
    projected = project_channel_states(_source(source_rows), _observation(observation_rows), config)
    receipt = qualification_receipt(projected, config)
    bombus = receipt.loc[receipt["channel_id"].eq("bombus")].iloc[0]
    assert bombus["n_retained"] == 50
    assert bombus["n_disrupted"] == 10
    assert bombus["support_tier"] == "confirmatory"
    assert bool(bombus["N1_gate_eligible"]) is True


def test_existing_bombus_assets_adapt_without_trait_inputs() -> None:
    config = _config()
    applicability = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "source_region_id": "src-1",
                "applicability": "applicable",
                "source_region_evidence_id": "bev-1",
                "source_region_review_status": "accepted",
                "assignment_review_status": "accepted",
            },
            {
                "island_id": "i2",
                "source_region_id": "src-2",
                "applicability": "structurally_not_applicable",
                "source_region_evidence_id": "bev-2",
                "source_region_review_status": "accepted",
                "assignment_review_status": "accepted",
            },
        ]
    )
    evidence = pd.DataFrame(
        [
            {
                "source_region_evidence_id": "bev-1",
                "evidence_type": "distribution_review",
                "source_citation": "Bombus source 1",
                "source_url": "https://example.org/1",
                "review_status": "accepted",
            },
            {
                "source_region_evidence_id": "bev-2",
                "evidence_type": "distribution_review",
                "source_citation": "Bombus source 2",
                "source_url": "https://example.org/2",
                "review_status": "accepted",
            },
        ]
    )
    diagnostics = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "bombus_occurrence_evidence": "adequate_non_detection",
                "bombus_record_count": 0,
                "target_group_record_count": 70,
                "target_group_spatial_units": 5,
                "target_group_temporal_units": 4,
                "distinct_dataset_count": 3,
                "latest_record_year": 2025,
                "observation_diagnostic_flags": "",
            },
            {
                "island_id": "i2",
                "bombus_occurrence_evidence": "detected",
                "bombus_record_count": 2,
                "target_group_record_count": 2,
                "target_group_spatial_units": 1,
                "target_group_temporal_units": 1,
                "distinct_dataset_count": 1,
                "latest_record_year": 2010,
                "observation_diagnostic_flags": "",
            },
        ]
    )
    source = adapt_bombus_source_availability(applicability, evidence, config)
    observation = adapt_bombus_observation(diagnostics, config)
    projected = project_channel_states(source, observation, config).set_index("island_id")
    assert projected.loc["i1", "channel_state"] == "disrupted"
    assert projected.loc["i2", "channel_state"] == "structurally_absent"
    assert "detected_despite_structural_source_absence" in projected.loc["i2", "projection_flags"]
