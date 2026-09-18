from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path("scripts/freeze_chapter1_h4_pollimcrop_preflight_result.py")
SPEC = importlib.util.spec_from_file_location("h4_pollimcrop_freezer", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

PREFLIGHT_CONFIG = Path("config/chapter1_h4_pollimcrop_transportability_preflight_v1.yml")
MAPPING = Path("config/chapter1_h4_pollimcrop_response_mapping_v1.yml")


def _base_preflight() -> dict:
    mapping_sha = MODULE._sha256(MAPPING)
    return {
        "contract": MODULE.PREFLIGHT_CONTRACT,
        "outcomes_read": False,
        "response_mapping": {
            "sha256": mapping_sha,
            "response_column": "PL_effect_size",
            "outcome_values_read": False,
        },
        "validated_expected_scope": {
            "studies": 294,
            "experiments": 1169,
            "crop_species": 108,
            "countries": 62,
        },
        "source_format": {
            "delimiter": ";",
            "decimal_mark": ",",
            "encoding": "utf-8-sig",
        },
        "figshare": {
            "collection_id": 6640595,
            "article_id": 1,
            "file_id": 2,
            "file_name": "PolLimCrop_dataset.zip",
            "file_sha256": "a" * 64,
            "csv_member": "PolLimCrop_dataset.csv",
            "csv_sha256": "b" * 64,
        },
        "metadata_columns_read": ["article_code", "species"],
        "support": {
            "H4a_reproductive_assurance": {"evaluable": True},
            "H4b_accessibility_generalization": {"evaluable": False},
        },
    }


def _freeze(preflight: dict) -> dict:
    return MODULE.build_result_lock(
        preflight,
        preflight_config_path=PREFLIGHT_CONFIG,
        response_mapping_path=MAPPING,
        workflow_run_id=123,
        workflow_head_sha="abc123",
        artifact_id=456,
        artifact_digest="sha256:" + "c" * 64,
    )


def test_freezer_authorizes_only_support_admitted_hypothesis() -> None:
    result = _freeze(_base_preflight())
    assert result["status"] == "outcome_blind_support_gate_passed"
    assert result["decision"]["admitted_hypotheses"] == [
        "H4a_reproductive_assurance"
    ]
    assert result["decision"]["support_failed_hypotheses"] == [
        "H4b_accessibility_generalization"
    ]
    assert result["decision"]["outcome_extraction_authorized"] is True
    assert result["outcomes_read"] is False
    assert result["blinding"]["preflight_outcomes_read"] is False


def test_freezer_records_no_support_without_authorizing_outcomes() -> None:
    preflight = _base_preflight()
    for value in preflight["support"].values():
        value["evaluable"] = False
    result = _freeze(preflight)
    assert result["status"] == "outcome_blind_support_gate_not_met"
    assert result["decision"]["admitted_hypotheses"] == []
    assert result["decision"]["outcome_extraction_authorized"] is False


def test_freezer_rejects_response_mapping_hash_mismatch() -> None:
    preflight = _base_preflight()
    preflight["response_mapping"]["sha256"] = "0" * 64
    with pytest.raises(Exception, match="response-mapping SHA-256 mismatch"):
        _freeze(preflight)


def test_freezer_rejects_scope_mismatch() -> None:
    preflight = _base_preflight()
    preflight["validated_expected_scope"]["crop_species"] = 107
    with pytest.raises(Exception, match="validated scope mismatch"):
        _freeze(preflight)
