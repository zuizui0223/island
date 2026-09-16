from __future__ import annotations

import pandas as pd
import yaml

from island_v2.chapter1_h5_channel_state_overlap_audit import audit_channel, summarize


def cfg() -> dict:
    return yaml.safe_load(open("config/chapter1_h5_channel_state_overlap_audit_v1.yml", encoding="utf-8"))


def test_audit_counts_strict_states_only() -> None:
    table = pd.DataFrame(
        {
            "island_id": ["a", "b", "c", "d"],
            "channel_id": ["bombus"] * 4,
            "observation_state": ["detected", "adequate_non_detection", "insufficient_effort", "unresolved"],
        }
    )
    cov = pd.DataFrame(
        {"island_id": ["a", "b", "c", "d"], "analysis_regime": ["tropical"] * 4}
    )
    out = audit_channel(table, cov, "bombus", cfg()).iloc[0]
    assert out["retained"] == 1
    assert out["disrupted"] == 1
    assert out["strict_evaluable"] == 2
    assert out["insufficient_effort"] == 1
    assert out["unresolved"] == 1


def test_summary_requires_both_primary_contexts() -> None:
    rows = []
    for channel in cfg()["channels"]:
        rows.extend(
            [
                {"channel": channel, "context": "northern_midlatitude", "passes_reference_overlap": True},
                {"channel": channel, "context": "tropical", "passes_reference_overlap": channel == "bombus"},
            ]
        )
    result = summarize(pd.DataFrame(rows), cfg())
    assert result["n_channels_passing_both_primary_contexts"] == 1
    assert result["channel_primary_overlap_pass"]["bombus"] is True
    assert result["mechanism_promoted"] is False
