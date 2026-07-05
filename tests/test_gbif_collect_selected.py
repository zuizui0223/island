import importlib
import json
from pathlib import Path

import pandas as pd
import pytest
import typer

selected = importlib.import_module("island_v2.gbif_collect_selected")


def _campaign() -> dict[str, object]:
    return {
        "ledger": [
            {"block_id": "north_temperate", "request_status": "succeeded", "download_key": "north"},
            {"block_id": "south", "request_status": "succeeded", "download_key": "south"},
            {"block_id": "running", "request_status": "running", "download_key": "running"},
        ]
    }


def test_select_succeeded_block_accepts_only_named_succeeded_uncollected_entry():
    entry = selected.select_succeeded_block(_campaign(), pd.DataFrame(), "north_temperate")

    assert entry["block_id"] == "north_temperate"


def test_select_succeeded_block_rejects_running_unknown_and_completed_ids():
    with pytest.raises(typer.BadParameter, match="not a succeeded"):
        selected.select_succeeded_block(_campaign(), pd.DataFrame(), "running")

    with pytest.raises(typer.BadParameter, match="not a succeeded"):
        selected.select_succeeded_block(_campaign(), pd.DataFrame(), "missing")

    manifest = pd.DataFrame({"block_id": ["south"], "collection_status": ["collected"]})
    with pytest.raises(typer.BadParameter, match="already collected"):
        selected.select_succeeded_block(_campaign(), manifest, "south")


def test_selected_block_request_is_geographic_operational_wiring_only():
    request = json.loads(
        Path("config/gbif_selected_collection_request_001.json").read_text(encoding="utf-8")
    )

    assert request["block_id"] == "gbif_block_cell_w+120.0_s+30.0_003"
    assert request["rationale_category"] == "geographic_source_region_coverage"
    assert request["expected_workflow"] == "collect_gbif_selected_block"
    assert "model_based_priority_decision" in request["safety_exclusions"]
    assert "analysis_inclusion_decision" in request["safety_exclusions"]
