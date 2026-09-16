from __future__ import annotations

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_h5_bombus_upstream_bridge import (
    overlap_gate,
    project_bombus_state,
    state_support,
)


def _config() -> dict:
    return yaml.safe_load(open("config/chapter1_h5_bombus_upstream_bridge_v1.yml", encoding="utf-8"))


def test_state_projection_is_fail_closed() -> None:
    obs = pd.DataFrame(
        {
            "island_id": ["a", "b", "c", "d"],
            "channel_id": ["bombus"] * 4,
            "observation_state": [
                "detected",
                "adequate_non_detection",
                "insufficient_effort",
                "unresolved",
            ],
        }
    )
    got = project_bombus_state(obs).set_index("island_id")
    assert got.loc["a", "bombus_state"] == "retained"
    assert got.loc["b", "bombus_state"] == "disrupted"
    assert pd.isna(got.loc["c", "bombus_state"])
    assert pd.isna(got.loc["d", "bombus_state"])
    assert got.loc["a", "bombus_disrupted"] == 0.0
    assert got.loc["b", "bombus_disrupted"] == 1.0


def test_overlap_gate_requires_both_states_in_both_primary_contexts() -> None:
    cfg = _config()
    rows = []
    for context, retained, disrupted in [
        ("northern_midlatitude", 100, 9),
        ("tropical", 12, 20),
    ]:
        rows.extend(
            {"analysis_regime": context, "bombus_state": "retained"} for _ in range(retained)
        )
        rows.extend(
            {"analysis_regime": context, "bombus_state": "disrupted"} for _ in range(disrupted)
        )
    data = pd.DataFrame(rows)
    support = state_support(data, cfg)
    gate = overlap_gate(support, cfg)
    assert gate["n_primary_contexts_passing"] == 1
    assert gate["overlap_gate_passed"] is False


def test_overlap_gate_passes_when_both_primary_contexts_have_overlap() -> None:
    cfg = _config()
    rows = []
    for context in ["northern_midlatitude", "tropical"]:
        rows.extend(
            {"analysis_regime": context, "bombus_state": "retained"} for _ in range(12)
        )
        rows.extend(
            {"analysis_regime": context, "bombus_state": "disrupted"} for _ in range(11)
        )
    support = state_support(pd.DataFrame(rows), cfg)
    assert overlap_gate(support, cfg)["overlap_gate_passed"] is True


def test_no_record_count_is_needed_for_projection() -> None:
    obs = pd.DataFrame(
        {
            "island_id": ["x", "y"],
            "channel_id": ["bombus", "bombus"],
            "observation_state": ["detected", "adequate_non_detection"],
            "channel_record_count": [900, 0],
        }
    )
    got = project_bombus_state(obs)
    assert np.allclose(got["bombus_disrupted"].to_numpy(), [0.0, 1.0])
