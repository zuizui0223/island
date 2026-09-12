from __future__ import annotations

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_nee_source_proxy import (
    aggregate_source_mode,
    receipt,
    source_mode_spec,
    validate_assignments,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_nee_source_proxy.yml", encoding="utf-8"))


def _assignments() -> pd.DataFrame:
    rows = []
    for island, entities in {
        "i1": ["e1", "e2", "e3", "e4", "e5"],
        "i2": ["e6", "e7", "e8", "e9", "e10"],
    }.items():
        for rank, entity in enumerate(entities, start=1):
            rows.append(
                {
                    "island_id": island,
                    "source_mode": "geo_k5",
                    "source_rank": rank,
                    "entity_ID": entity,
                }
            )
    return pd.DataFrame(rows)


def _state(entity: str, channel: str, state: str, evidence: str) -> dict:
    decisive = state in {"available", "structurally_absent"}
    return {
        "entity_ID": entity,
        "channel_id": channel,
        "source_state": state,
        "evidence_id": evidence if decisive else "",
        "evidence_type": (
            "positive"
            if state == "available"
            else "absence"
            if state == "structurally_absent"
            else ""
        ),
        "source_citation": f"citation-{evidence}" if decisive else "",
        "source_url": f"https://example.org/{evidence}" if decisive else "",
        "review_status": "accepted" if decisive else "pending",
    }


def _set_state(
    table: pd.DataFrame,
    entity: str,
    channel: str,
    state: str,
    evidence: str,
) -> None:
    """Replace exactly one synthetic entity x channel row without Series index alignment."""
    mask = table["entity_ID"].eq(entity) & table["channel_id"].eq(channel)
    assert int(mask.sum()) == 1
    replacement = _state(entity, channel, state, evidence)
    for column, value in replacement.items():
        table.loc[mask, column] = value


def _entity_states() -> pd.DataFrame:
    rows = []
    channels = _config()["channels"]
    for entity in [f"e{i}" for i in range(1, 11)]:
        for channel in channels:
            rows.append(_state(entity, channel, "unresolved", f"u-{entity}-{channel}"))
    table = pd.DataFrame(rows)
    # i1: one positive Bombus source is enough for source-set availability.
    _set_state(table, "e3", "bombus", "available", "p1")
    # i2: all five source entities explicitly lack Bombus -> structural absence.
    for entity in ["e6", "e7", "e8", "e9", "e10"]:
        _set_state(table, entity, "bombus", "structurally_absent", f"a-{entity}")
    assert not table[["entity_ID", "channel_id"]].duplicated().any()
    return table


def test_primary_mode_is_frozen_geo_k5() -> None:
    config = _config()
    assert config["primary_source_proxy"]["source_mode"] == "geo_k5"
    assert source_mode_spec(config, "geo_k5") == {"expected_k": 5, "role": "primary"}
    assert source_mode_spec(config, "geo_k10") == {"expected_k": 10, "role": "sensitivity"}


def test_available_if_any_selected_source_entity_is_available() -> None:
    config = _config()
    output, audit = aggregate_source_mode(_assignments(), _entity_states(), config, "geo_k5")
    row = output.loc[output["island_id"].eq("i1") & output["channel_id"].eq("bombus")].iloc[0]
    arow = audit.loc[audit["island_id"].eq("i1") & audit["channel_id"].eq("bombus")].iloc[0]
    assert row["source_state"] == "available"
    assert row["source_region_id"] == "e1"
    assert row["review_status"] == "accepted"
    assert arow["n_selected_entities"] == 5
    assert arow["n_available_entities"] == 1
    assert arow["selected_entity_IDs"] == "e1|e2|e3|e4|e5"


def test_structural_absence_requires_all_selected_source_entities() -> None:
    config = _config()
    output, audit = aggregate_source_mode(_assignments(), _entity_states(), config, "geo_k5")
    row = output.loc[output["island_id"].eq("i2") & output["channel_id"].eq("bombus")].iloc[0]
    arow = audit.loc[audit["island_id"].eq("i2") & audit["channel_id"].eq("bombus")].iloc[0]
    assert row["source_state"] == "structurally_absent"
    assert arow["n_structurally_absent_entities"] == 5
    assert arow["n_unresolved_entities"] == 0


def test_partial_absence_without_positive_stays_unresolved() -> None:
    states = _entity_states()
    _set_state(states, "e1", "diptera", "structurally_absent", "a1")
    output, _ = aggregate_source_mode(_assignments(), states, _config(), "geo_k5")
    row = output.loc[output["island_id"].eq("i1") & output["channel_id"].eq("diptera")].iloc[0]
    assert row["source_state"] == "unresolved"
    assert row["review_status"] == "pending"


def test_incomplete_source_rank_set_fails_closed() -> None:
    assignments = _assignments()
    assignments = assignments.loc[
        ~(assignments["island_id"].eq("i1") & assignments["source_rank"].eq(5))
    ]
    with pytest.raises(ValueError, match="incomplete/unexpected"):
        validate_assignments(assignments, _config(), "geo_k5")


def test_source_context_is_rank1_not_selected_after_channel_state() -> None:
    states = _entity_states()
    output, _ = aggregate_source_mode(_assignments(), states, _config(), "geo_k5")
    i1 = output.loc[output["island_id"].eq("i1")]
    # Even though Bombus evidence is at rank 3, every channel uses the predeclared rank-1 context.
    assert set(i1["source_region_id"]) == {"e1"}


def test_receipt_keeps_primary_proxy_and_independence_boundary() -> None:
    config = _config()
    output, audit = aggregate_source_mode(_assignments(), _entity_states(), config, "geo_k5")
    report = receipt(output, audit, "geo_k5", config)
    assert report["role"] == "primary"
    assert report["expected_k"] == 5
    assert report["uses_pollinator_retention_outcomes"] is False
    assert report["uses_focal_plant_traits"] is False
    assert report["all_selected_source_entities_preserved_in_audit"] is True
