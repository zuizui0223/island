from __future__ import annotations

import pandas as pd
import pytest
import yaml

import island_v2.chapter1_nee_source_positive_scan as scan
from island_v2.chapter1_nee_source_positive_scan import (
    confirmatory_target_keys,
    find_confirmatory_positive,
    normalized_binomial_key,
    scan_entity_channel,
)


def _config() -> dict:
    return yaml.safe_load(
        open("config/chapter1_nee_source_positive_scan.yml", encoding="utf-8")
    )


def test_normalized_binomial_key_is_conservative() -> None:
    assert normalized_binomial_key(" Bombus terrestris (Linnaeus, 1758) ") == "bombus terrestris"
    assert normalized_binomial_key("Danaus plexippus") == "danaus plexippus"
    assert normalized_binomial_key("Bombus sp.") == ""
    assert normalized_binomial_key("") == ""


def test_confirmatory_target_keys_exclude_sensitivity_taxa() -> None:
    catalog = pd.DataFrame(
        [
            {"channel_id": "lepidoptera", "pollinator_species": "Danaus plexippus", "catalog_tier": "confirmatory"},
            {"channel_id": "lepidoptera", "pollinator_species": "Papilio xuthus", "catalog_tier": "sensitivity"},
            {"channel_id": "diptera", "pollinator_species": "Eristalis tenax", "catalog_tier": "confirmatory"},
        ]
    )
    assert confirmatory_target_keys(catalog, "lepidoptera") == {"danaus plexippus"}


def test_positive_hit_accepts_confirmatory_species_and_flags_unknown_establishment() -> None:
    records = [
        {
            "key": 123,
            "species": "Danaus plexippus",
            "family": "Nymphalidae",
            "occurrenceStatus": "PRESENT",
            "basisOfRecord": "HUMAN_OBSERVATION",
            "establishmentMeans": "",
        }
    ]
    hit, examined, flags = find_confirmatory_positive(records, {"danaus plexippus"})
    assert hit is not None
    assert examined == 1
    assert "matched_target_unknown_establishment" in flags


def test_positive_hit_rejects_introduced_captive_and_explicit_absent() -> None:
    targets = {"bombus terrestris"}
    records = [
        {"species": "Bombus terrestris", "establishmentMeans": "INTRODUCED", "occurrenceStatus": "PRESENT", "basisOfRecord": "HUMAN_OBSERVATION"},
        {"species": "Bombus terrestris", "establishmentMeans": "", "occurrenceStatus": "PRESENT", "basisOfRecord": "LIVING_SPECIMEN"},
        {"species": "Bombus terrestris", "establishmentMeans": "", "occurrenceStatus": "ABSENT", "basisOfRecord": "HUMAN_OBSERVATION"},
    ]
    hit, examined, _ = find_confirmatory_positive(records, targets)
    assert hit is None
    assert examined == 3


def test_non_bombus_family_filter_rejects_non_bee_hymenoptera() -> None:
    records = [
        {"species": "Target species", "family": "Vespidae", "occurrenceStatus": "PRESENT", "basisOfRecord": "HUMAN_OBSERVATION", "establishmentMeans": ""},
        {"species": "Target species", "family": "Halictidae", "occurrenceStatus": "PRESENT", "basisOfRecord": "HUMAN_OBSERVATION", "establishmentMeans": ""},
    ]
    hit, examined, _ = find_confirmatory_positive(
        records,
        {"target species"},
        allowed_families={"Halictidae", "Apidae"},
    )
    assert hit is not None
    assert hit["family"] == "Halictidae"
    assert examined == 2


class _FakeClient:
    pass


def test_capped_no_hit_is_unresolved_never_structural_absence(monkeypatch: pytest.MonkeyPatch) -> None:
    config = _config()

    def fake_get_json(*args, **kwargs):
        return {
            "results": [
                {
                    "key": 999,
                    "species": "Irrelevant species",
                    "family": "Nymphalidae",
                    "occurrenceStatus": "PRESENT",
                    "basisOfRecord": "HUMAN_OBSERVATION",
                    "establishmentMeans": "",
                }
            ],
            "endOfRecords": True,
        }

    monkeypatch.setattr(scan, "_get_json_with_retry", fake_get_json)
    row = scan_entity_channel(
        _FakeClient(),
        entity_id="42",
        geometry_wkt="POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
        channel_id="lepidoptera",
        taxon_key="797",
        target_keys={"danaus plexippus"},
        config=config,
    )
    assert row["source_state"] == "unresolved"
    assert row["review_status"] == "pending"
    assert row["evidence_id"] == ""
    assert "no_confirmatory_hit" in row["quality_flags"]


def test_positive_entity_scan_returns_available(monkeypatch: pytest.MonkeyPatch) -> None:
    config = _config()

    def fake_get_json(*args, **kwargs):
        return {
            "results": [
                {
                    "key": 321,
                    "species": "Danaus plexippus",
                    "family": "Nymphalidae",
                    "occurrenceStatus": "PRESENT",
                    "basisOfRecord": "HUMAN_OBSERVATION",
                    "establishmentMeans": "NATIVE",
                }
            ],
            "endOfRecords": True,
        }

    monkeypatch.setattr(scan, "_get_json_with_retry", fake_get_json)
    row = scan_entity_channel(
        _FakeClient(),
        entity_id="42",
        geometry_wkt="POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
        channel_id="lepidoptera",
        taxon_key="797",
        target_keys={"danaus plexippus"},
        config=config,
    )
    assert row["source_state"] == "available"
    assert row["review_status"] == "accepted"
    assert row["evidence_id"] == "GBIF:321"


def test_no_channel_targets_fails_closed() -> None:
    catalog = pd.DataFrame(
        [{"channel_id": "diptera", "pollinator_species": "Eristalis tenax", "catalog_tier": "confirmatory"}]
    )
    with pytest.raises(ValueError, match="no confirmatory"):
        confirmatory_target_keys(catalog, "lepidoptera")
