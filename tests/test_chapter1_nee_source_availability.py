from __future__ import annotations

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_nee_source_availability import (
    build_island_source_availability,
    confirmatory_species_by_channel,
    source_availability_receipt,
    validate_structural_absence,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_nee_source_availability.yml", encoding="utf-8"))


def _catalog() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"channel_id": "bombus", "pollinator_species": "Bombus ardens", "catalog_tier": "confirmatory"},
            {"channel_id": "lepidoptera", "pollinator_species": "Danaus plexippus", "catalog_tier": "confirmatory"},
            {"channel_id": "diptera", "pollinator_species": "Eristalis tenax", "catalog_tier": "confirmatory"},
            {"channel_id": "flower_visiting_birds", "pollinator_species": "Birdus example", "catalog_tier": "sensitivity"},
        ]
    )


def _assignments() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"island_id": "i1", "source_region_id": "r1", "review_status": "accepted"},
            {"island_id": "i2", "source_region_id": "r2", "review_status": "accepted"},
            {"island_id": "i3", "source_region_id": "r3", "review_status": "pending"},
        ]
    )


def _positive() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source_region_id": "r1",
                "pollinator_species": "Bombus ardens",
                "evidence_id": "p1",
                "evidence_type": "curated_native_distribution",
                "source_citation": "Bombus source",
                "source_url": "https://example.org/p1",
                "review_status": "accepted",
            },
            {
                "source_region_id": "r1",
                "pollinator_species": "Danaus plexippus",
                "evidence_id": "p2",
                "evidence_type": "quality_filtered_source_region_occurrence",
                "source_citation": "Butterfly occurrence",
                "source_url": "https://example.org/p2",
                "review_status": "accepted",
            },
            # Sensitivity-only catalog taxon must not establish primary availability.
            {
                "source_region_id": "r1",
                "pollinator_species": "Birdus example",
                "evidence_id": "p3",
                "evidence_type": "curated_regional_checklist",
                "source_citation": "Bird list",
                "source_url": "https://example.org/p3",
                "review_status": "accepted",
            },
        ]
    )


def _structural() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source_region_id": "r2",
                "channel_id": "bombus",
                "evidence_id": "a1",
                "evidence_type": "authoritative_native_range_exclusion",
                "source_citation": "Bombus absent from source region",
                "source_url": "https://example.org/a1",
                "review_status": "accepted",
            }
        ]
    )


def test_confirmatory_catalog_only_defines_source_target_taxa() -> None:
    targets = confirmatory_species_by_channel(_catalog(), _config())
    assert targets["bombus"] == {"Bombus ardens"}
    assert targets["lepidoptera"] == {"Danaus plexippus"}
    assert targets["flower_visiting_birds"] == set()


def test_positive_source_evidence_establishes_available_only_for_matching_channel() -> None:
    result, conflicts = build_island_source_availability(
        _assignments(), _catalog(), _positive(), _structural(), _config()
    )
    i1 = result.loc[result["island_id"].eq("i1")].set_index("channel_id")
    assert i1.loc["bombus", "source_state"] == "available"
    assert i1.loc["lepidoptera", "source_state"] == "available"
    assert i1.loc["flower_visiting_birds", "source_state"] == "unresolved"
    assert conflicts.empty


def test_explicit_structural_absence_is_preserved_and_zero_evidence_is_unresolved() -> None:
    result, _ = build_island_source_availability(
        _assignments(), _catalog(), _positive(), _structural(), _config()
    )
    i2 = result.loc[result["island_id"].eq("i2")].set_index("channel_id")
    assert i2.loc["bombus", "source_state"] == "structurally_absent"
    assert i2.loc["lepidoptera", "source_state"] == "unresolved"
    assert i2.loc["diptera", "source_state"] == "unresolved"


def test_pending_source_assignment_forces_all_channels_unresolved() -> None:
    result, _ = build_island_source_availability(
        _assignments(), _catalog(), _positive(), _structural(), _config()
    )
    i3 = result.loc[result["island_id"].eq("i3")]
    assert set(i3["source_state"]) == {"unresolved"}
    assert set(i3["review_status"]) == {"pending"}


def test_positive_and_structural_absence_conflict_has_no_precedence() -> None:
    positive = pd.concat(
        [
            _positive(),
            pd.DataFrame(
                [
                    {
                        "source_region_id": "r2",
                        "pollinator_species": "Bombus ardens",
                        "evidence_id": "p-conflict",
                        "evidence_type": "curated_native_distribution",
                        "source_citation": "Contradictory Bombus presence",
                        "source_url": "https://example.org/conflict",
                        "review_status": "accepted",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    result, conflicts = build_island_source_availability(
        _assignments(), _catalog(), positive, _structural(), _config()
    )
    row = result.loc[
        result["island_id"].eq("i2") & result["channel_id"].eq("bombus")
    ].iloc[0]
    assert row["source_state"] == "unresolved"
    assert row["review_status"] == "pending"
    assert len(conflicts) == 1


def test_zero_gbif_records_is_prohibited_as_structural_absence() -> None:
    structural = pd.DataFrame(
        [
            {
                "source_region_id": "r2",
                "channel_id": "bombus",
                "evidence_id": "bad",
                "evidence_type": "zero_GBIF_records",
                "source_citation": "none",
                "source_url": "",
                "review_status": "accepted",
            }
        ]
    )
    with pytest.raises(ValueError, match="prohibited evidence type"):
        validate_structural_absence(structural, _config())


def test_receipt_preserves_independence_boundary() -> None:
    result, conflicts = build_island_source_availability(
        _assignments(), _catalog(), _positive(), _structural(), _config()
    )
    receipt = source_availability_receipt(result, conflicts)
    assert receipt["uses_focal_plant_traits"] is False
    assert receipt["uses_island_channel_observation"] is False
    assert receipt["zero_occurrence_can_create_structural_absence"] is False
    assert receipt["n_rows"] == 3 * 5
