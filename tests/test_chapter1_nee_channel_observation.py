from __future__ import annotations

import pandas as pd
import yaml

from island_v2.chapter1_nee_channel_observation import (
    classify_channel_observations,
    observation_receipt,
    target_species,
)


def _policy() -> dict:
    return yaml.safe_load(
        open("config/chapter1_nee_channel_observation_policy.yml", encoding="utf-8")
    )


def _catalog() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"channel_id": "lepidoptera", "pollinator_species": "Danaus plexippus", "catalog_tier": "confirmatory"},
            {"channel_id": "lepidoptera", "pollinator_species": "Papilio xuthus", "catalog_tier": "sensitivity"},
            {"channel_id": "diptera", "pollinator_species": "Eristalis tenax", "catalog_tier": "confirmatory"},
        ]
    )


def _row(
    island: str,
    species: str,
    i: int,
    *,
    dataset: str | None = None,
    year: int | None = None,
    establishment: str = "",
) -> dict:
    return {
        "island_id": island,
        "gbif_id": f"gbif-{island}-{i}",
        "species": species,
        "dataset_key": dataset or f"ds-{i % 3}",
        "year": str(year if year is not None else 2024 - (i % 3)),
        "decimal_latitude": str(20.0 + (i % 5) * 0.2),
        "decimal_longitude": str(130.0 + (i % 5) * 0.2),
        "establishment_means": establishment,
    }


def test_confirmatory_catalog_excludes_sensitivity_taxon() -> None:
    assert target_species(_catalog(), "lepidoptera", True) == {"Danaus plexippus"}
    assert target_species(_catalog(), "lepidoptera", False) == {
        "Danaus plexippus",
        "Papilio xuthus",
    }


def test_detection_overrides_low_background_effort() -> None:
    records = pd.DataFrame([_row("a", "Danaus plexippus", 1)])
    result = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy()
    )
    row = result.iloc[0]
    assert row["observation_state"] == "detected"
    assert row["channel_record_count"] == 1
    assert row["background_record_count"] == 1
    assert "below_background_records" in row["quality_flags"]


def test_adequate_non_detection_requires_all_effort_components() -> None:
    records = pd.DataFrame(
        [_row("a", f"Background species {i}", i) for i in range(60)]
    )
    result = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy()
    )
    row = result.iloc[0]
    assert row["channel_record_count"] == 0
    assert row["background_record_count"] == 60
    assert row["background_spatial_units"] >= 3
    assert row["background_temporal_units"] >= 2
    assert row["distinct_dataset_count"] >= 2
    assert row["observation_state"] == "adequate_non_detection"


def test_sparse_zero_is_insufficient_effort_not_loss() -> None:
    records = pd.DataFrame(
        [_row("a", f"Background species {i}", i) for i in range(8)]
    )
    result = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy()
    )
    assert result.iloc[0]["observation_state"] == "insufficient_effort"


def test_stale_background_does_not_qualify_zero() -> None:
    records = pd.DataFrame(
        [_row("a", f"Background species {i}", i, year=1990 + (i % 3)) for i in range(60)]
    )
    result = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy()
    )
    row = result.iloc[0]
    assert row["observation_state"] == "insufficient_effort"
    assert "stale_or_missing_background_recency" in row["quality_flags"]


def test_known_introduced_target_does_not_create_primary_detection() -> None:
    rows = [_row("a", f"Background species {i}", i) for i in range(60)]
    rows.append(_row("a", "Danaus plexippus", 999, establishment="INTRODUCED"))
    result = classify_channel_observations(
        pd.DataFrame(rows), _catalog(), "lepidoptera", _policy()
    )
    row = result.iloc[0]
    assert row["observation_state"] == "adequate_non_detection"
    assert row["channel_record_count"] == 0
    assert "introduced_target_rows_excluded=1" in row["quality_flags"]


def test_liberal_effort_is_sensitivity_only_but_classifies() -> None:
    records = pd.DataFrame(
        [_row("a", f"Background species {i}", i) for i in range(30)]
    )
    primary = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy(), effort_tier="primary"
    )
    liberal = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy(), effort_tier="liberal"
    )
    assert primary.iloc[0]["observation_state"] == "insufficient_effort"
    assert liberal.iloc[0]["observation_state"] == "adequate_non_detection"


def test_receipt_keeps_claim_ceiling() -> None:
    records = pd.DataFrame([_row("a", "Danaus plexippus", 1)])
    result = classify_channel_observations(
        records, _catalog(), "lepidoptera", _policy()
    )
    receipt = observation_receipt(result, "lepidoptera", "primary")
    assert receipt["uses_focal_plant_traits"] is False
    assert receipt["detection_claim"] == "potential_partner_channel_present_not_realized_service"
