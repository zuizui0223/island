from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml
from shapely.geometry import MultiPolygon, Polygon

import island_v2.chapter1_nee_island_search_observation as search
from island_v2.chapter1_nee_island_search_observation import (
    SearchMeta,
    eligible_islands,
    exact_query_wkts,
    finalize_observation,
    scan_island,
)
from island_v2.chapter1_nee_channel_observation import load_policy

CONFIG = Path("config/chapter1_nee_island_search_acquisition.yml")
POLICY = Path("config/chapter1_nee_channel_observation_policy.yml")


def _config() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def _policy() -> dict:
    return load_policy(POLICY)


def _catalog() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "channel_id": "lepidoptera",
                "pollinator_species": "Danaus plexippus",
                "catalog_tier": "confirmatory",
            }
        ]
    )


def _adequate_records(target: bool = False) -> pd.DataFrame:
    rows = []
    for index in range(60):
        rows.append(
            {
                "island_id": "i1",
                "species": "Danaus plexippus" if target and index == 0 else f"Species {index}",
                "dataset_key": "d1" if index % 2 == 0 else "d2",
                "year": 2025 if index % 2 == 0 else 2024,
                "decimal_latitude": 10.0 + (index % 3) * 0.2,
                "decimal_longitude": 20.0 + (index % 3) * 0.2,
                "establishment_means": "NATIVE" if target and index == 0 else "",
                "gbif_id": str(index + 1),
            }
        )
    return pd.DataFrame(rows)


def test_source_available_scope_is_outcome_blind() -> None:
    table = pd.DataFrame(
        [
            {"island_id": "i1", "channel_id": "lepidoptera", "source_state": "available"},
            {"island_id": "i2", "channel_id": "lepidoptera", "source_state": "unresolved"},
            {"island_id": "i3", "channel_id": "diptera", "source_state": "available"},
        ]
    )
    assert eligible_islands(table, "lepidoptera") == ["i1"]


def test_exact_components_are_not_simplified_and_are_area_sorted() -> None:
    small = Polygon([(10, 0), (11, 0), (11, 1), (10, 1), (10, 0)])
    large = Polygon([(0, 0), (4, 0), (4, 2), (0, 2), (0, 0)])
    wkts = exact_query_wkts(MultiPolygon([small, large]))
    assert len(wkts) == 2
    assert wkts[0].startswith("POLYGON")
    from shapely import wkt

    recovered = [wkt.loads(value) for value in wkts]
    assert recovered[0].equals(large)
    assert recovered[1].equals(small)


def test_complete_search_can_yield_adequate_non_detection() -> None:
    row = finalize_observation(
        _adequate_records(target=False),
        SearchMeta(True, False, "", 60, 1, 1, 1),
        island_id="i1",
        channel_id="lepidoptera",
        catalog=_catalog(),
        policy=_policy(),
    )
    assert row["observation_state"] == "adequate_non_detection"
    assert row["search_complete"] is True


def test_truncated_search_can_never_yield_adequate_non_detection() -> None:
    row = finalize_observation(
        _adequate_records(target=False),
        SearchMeta(False, True, "", 900, 0, 1, 1),
        island_id="i1",
        channel_id="lepidoptera",
        catalog=_catalog(),
        policy=_policy(),
    )
    assert row["observation_state"] == "insufficient_effort"
    assert "search_truncated_fixed_budget" in row["quality_flags"]


def test_detection_remains_detection_even_when_background_is_incomplete() -> None:
    row = finalize_observation(
        _adequate_records(target=True),
        SearchMeta(False, True, "", 900, 0, 1, 1),
        island_id="i1",
        channel_id="lepidoptera",
        catalog=_catalog(),
        policy=_policy(),
    )
    assert row["observation_state"] == "detected"


def test_api_error_without_target_is_unresolved() -> None:
    row = finalize_observation(
        _adequate_records(target=False),
        SearchMeta(False, False, "RuntimeError: synthetic", 30, 0, 1, 1),
        island_id="i1",
        channel_id="lepidoptera",
        catalog=_catalog(),
        policy=_policy(),
    )
    assert row["observation_state"] == "unresolved"
    assert "search_error" in row["quality_flags"]


class _FakeClient:
    pass


def test_shared_900_record_budget_across_components_and_subqueries(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config = _config()
    calls: list[dict[str, object]] = []

    def fake_get_json(*args, **kwargs):
        params = dict(kwargs["params"])
        calls.append(params)
        limit = int(params["limit"])
        offset = int(params["offset"])
        return {
            "results": [
                {
                    "key": offset + index + 1,
                    "species": "Irrelevant species",
                    "datasetKey": "d1",
                    "year": 2025,
                    "decimalLatitude": 10.0,
                    "decimalLongitude": 20.0,
                    "establishmentMeans": "",
                    "basisOfRecord": "HUMAN_OBSERVATION",
                    "occurrenceStatus": "PRESENT",
                }
                for index in range(limit)
            ],
            "endOfRecords": False,
        }

    monkeypatch.setattr(search, "_get_json_with_retry", fake_get_json)
    records, meta = scan_island(
        _FakeClient(),
        island_id="i1",
        geometry_wkts=[
            "POLYGON ((0 0, 1 0, 1 1, 0 1, 0 0))",
            "POLYGON ((2 0, 3 0, 3 1, 2 1, 2 0))",
        ],
        background_queries=[("A", "1"), ("B", "2")],
        targets={"Danaus plexippus"},
        config=config,
    )
    assert len(records) == 900
    assert meta.raw_records_examined == 900
    assert meta.truncated is True
    assert meta.complete is False
    assert sum(int(call["limit"]) for call in calls) == 900


def test_empty_complete_search_is_insufficient_not_absence() -> None:
    row = finalize_observation(
        pd.DataFrame(),
        SearchMeta(True, False, "", 0, 1, 1, 1),
        island_id="i1",
        channel_id="lepidoptera",
        catalog=_catalog(),
        policy=_policy(),
    )
    assert row["observation_state"] == "insufficient_effort"
