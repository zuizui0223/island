from __future__ import annotations

import pandas as pd
import pytest
import yaml
from shapely.geometry import MultiPolygon, Polygon
from shapely import wkt as shapely_wkt

import island_v2.chapter1_nee_source_positive_scan as scan
from island_v2.chapter1_nee_source_positive_scan import (
    _component_record_quotas,
    confirmatory_target_keys,
    find_confirmatory_positive,
    gift_entity_query_wkts,
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


def test_component_record_quotas_share_one_entity_channel_budget() -> None:
    quotas = _component_record_quotas(72, 900)
    assert len(quotas) == 72
    assert sum(quotas) == 900
    assert max(quotas) == 13
    assert min(quotas) == 12
    assert quotas[:36] == [13] * 36
    assert quotas[36:] == [12] * 36


def test_multipart_scan_never_exceeds_shared_900_record_budget(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _config()
    calls: list[dict[str, object]] = []

    def fake_get_json(*args, **kwargs):
        params = dict(kwargs["params"])
        calls.append(params)
        limit = int(params["limit"])
        return {
            "results": [
                {
                    "key": index,
                    "species": "Irrelevant species",
                    "family": "Nymphalidae",
                    "occurrenceStatus": "PRESENT",
                    "basisOfRecord": "HUMAN_OBSERVATION",
                    "establishmentMeans": "",
                }
                for index in range(limit)
            ],
            "endOfRecords": False,
        }

    monkeypatch.setattr(scan, "_get_json_with_retry", fake_get_json)
    wkts = [f"POLYGON (({i} 0, {i + 0.1} 0, {i + 0.1} 0.1, {i} 0.1, {i} 0))" for i in range(72)]
    row = scan_entity_channel(
        _FakeClient(),
        entity_id="345",
        geometry_wkts=wkts,
        channel_id="lepidoptera",
        taxon_key="797",
        target_keys={"danaus plexippus"},
        config=config,
    )
    assert row["source_state"] == "unresolved"
    assert row["n_records_examined"] == 900
    assert sum(int(call["limit"]) for call in calls) == 900
    assert len(calls) == 72
    assert "multipart_exact_component_query" in row["quality_flags"]
    assert "structurally_absent" not in row["quality_flags"]


def test_component_query_error_remains_unresolved_not_absent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _config()
    call_count = 0

    def fake_get_json(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        if call_count == 1:
            raise RuntimeError("synthetic transport failure")
        return {"results": [], "endOfRecords": True}

    monkeypatch.setattr(scan, "_get_json_with_retry", fake_get_json)
    row = scan_entity_channel(
        _FakeClient(),
        entity_id="345",
        geometry_wkts=[
            "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
            "POLYGON ((2 0, 3 0, 3 1, 2 1, 2 0))",
        ],
        channel_id="lepidoptera",
        taxon_key="797",
        target_keys={"danaus plexippus"},
        config=config,
    )
    assert row["source_state"] == "unresolved"
    assert row["evidence_id"] == ""
    assert "source_scan_error=RuntimeError" in row["quality_flags"]
    assert "n_component_query_errors=1" in row["quality_flags"]


def test_exact_component_wkts_are_ordered_by_area_without_simplification(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    large = Polygon([(0, 0), (4, 0), (4, 2), (0, 2), (0, 0)])
    small = Polygon([(10, 0), (11, 0), (11, 1), (10, 1), (10, 0)])
    geometry = MultiPolygon([small, large])
    monkeypatch.setattr(scan, "_gift_entity_geometry", lambda *args, **kwargs: geometry)

    wkts = gift_entity_query_wkts(_FakeClient(), "345", _config())
    recovered = [shapely_wkt.loads(value) for value in wkts]
    assert [part.area for part in recovered] == [8.0, 1.0]
    assert recovered[0].equals(large)
    assert recovered[1].equals(small)


def test_no_channel_targets_fails_closed() -> None:
    catalog = pd.DataFrame(
        [{"channel_id": "diptera", "pollinator_species": "Eristalis tenax", "catalog_tier": "confirmatory"}]
    )
    with pytest.raises(ValueError, match="no confirmatory"):
        confirmatory_target_keys(catalog, "lepidoptera")
