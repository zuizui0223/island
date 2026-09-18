from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path("scripts/freeze_chapter1_h4_pollimcrop_result.py")
SPEC = importlib.util.spec_from_file_location("h4_pollimcrop_result_freezer", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

MAPPING = Path("config/chapter1_h4_pollimcrop_response_mapping_v1.yml")


def _write_preflight_lock(tmp_path: Path, admitted: list[str]) -> Path:
    all_names = list(MODULE.HYPOTHESES)
    failed = [name for name in all_names if name not in admitted]
    value = {
        "lock_contract": MODULE.PREFLIGHT_LOCK_CONTRACT,
        "decision": {
            "outcome_extraction_authorized": bool(admitted),
            "admitted_hypotheses": admitted,
            "support_failed_hypotheses": failed,
        },
        "frozen_inputs": {
            "response_mapping_sha256": MODULE._sha256(MAPPING),
        },
    }
    path = tmp_path / "preflight_lock.json"
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def _result(admitted: list[str], supported: list[str]) -> dict:
    primary = {}
    for name in MODULE.HYPOTHESES:
        if name in admitted:
            primary[name] = {
                "evaluable": True,
                "supported": name in supported,
                "outcomes_used": True,
                "estimate": -0.2 if name in supported else 0.1,
                "one_sided_negative_p": 0.01 if name in supported else 0.8,
            }
        else:
            primary[name] = {
                "evaluable": False,
                "reason": "frozen_preflight_support_gate_failed",
                "outcomes_used": False,
            }
    return {
        "contract": MODULE.RESULT_CONTRACT,
        "primary_results": primary,
        "supported_primary_hypotheses": supported,
        "n_supported_primary_hypotheses": len(supported),
        "claim_ceiling": {"wild_plant_temporal_replication": False},
    }


def _freeze(
    tmp_path: Path,
    *,
    admitted: list[str],
    supported: list[str],
) -> dict:
    preflight_path = _write_preflight_lock(tmp_path, admitted)
    result_path = tmp_path / "RESULT.json"
    result_path.write_text(
        json.dumps(_result(admitted, supported)),
        encoding="utf-8",
    )
    return MODULE.build_result_lock(
        _result(admitted, supported),
        json.loads(preflight_path.read_text(encoding="utf-8")),
        result_json_path=result_path,
        preflight_lock_path=preflight_path,
        response_mapping_path=MAPPING,
        workflow_run_id=123,
        workflow_head_sha="abc123",
        artifact_id=456,
        artifact_digest="sha256:" + "d" * 64,
    )


def test_result_freezer_preserves_supported_admitted_hypothesis(tmp_path: Path) -> None:
    lock = _freeze(
        tmp_path,
        admitted=["H4a_reproductive_assurance"],
        supported=["H4a_reproductive_assurance"],
    )
    assert (
        lock["classification"]["H4a_reproductive_assurance"]
        == "transportability_supported"
    )
    assert (
        lock["classification"]["H4b_accessibility_generalization"]
        == "not_evaluable_support_gate_failed"
    )
    assert lock["boundary"]["wild_temporal_replication_repaired"] is False


def test_result_freezer_retains_admitted_negative_result(tmp_path: Path) -> None:
    lock = _freeze(
        tmp_path,
        admitted=["H4a_reproductive_assurance"],
        supported=[],
    )
    assert (
        lock["classification"]["H4a_reproductive_assurance"]
        == "transportability_not_supported"
    )


def test_result_freezer_rejects_outcome_use_for_failed_support(tmp_path: Path) -> None:
    preflight_path = _write_preflight_lock(
        tmp_path,
        ["H4a_reproductive_assurance"],
    )
    result = _result(["H4a_reproductive_assurance"], [])
    result["primary_results"]["H4b_accessibility_generalization"]["outcomes_used"] = True
    result_path = tmp_path / "RESULT.json"
    result_path.write_text(json.dumps(result), encoding="utf-8")
    with pytest.raises(Exception, match="support-failed hypothesis used outcomes"):
        MODULE.build_result_lock(
            result,
            json.loads(preflight_path.read_text(encoding="utf-8")),
            result_json_path=result_path,
            preflight_lock_path=preflight_path,
            response_mapping_path=MAPPING,
            workflow_run_id=123,
            workflow_head_sha="abc123",
            artifact_id=456,
            artifact_digest="sha256:" + "d" * 64,
        )


def test_result_freezer_rejects_mapping_change(tmp_path: Path) -> None:
    preflight_path = _write_preflight_lock(
        tmp_path,
        ["H4a_reproductive_assurance"],
    )
    result = _result(["H4a_reproductive_assurance"], [])
    result_path = tmp_path / "RESULT.json"
    result_path.write_text(json.dumps(result), encoding="utf-8")
    lock = json.loads(preflight_path.read_text(encoding="utf-8"))
    lock["frozen_inputs"]["response_mapping_sha256"] = "0" * 64
    with pytest.raises(Exception, match="response mapping changed"):
        MODULE.build_result_lock(
            result,
            lock,
            result_json_path=result_path,
            preflight_lock_path=preflight_path,
            response_mapping_path=MAPPING,
            workflow_run_id=123,
            workflow_head_sha="abc123",
            artifact_id=456,
            artifact_digest="sha256:" + "d" * 64,
        )
