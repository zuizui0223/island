from __future__ import annotations

import json

import pandas as pd
import yaml

from island_v2.schemas import TraitExtractionResponse
from island_v2.trait_extraction import (
    candidate_columns,
    is_retriable_error,
    strict_json_schema,
    write_run_artifacts,
)


class _StatusError(Exception):
    def __init__(self, status_code: int):
        self.status_code = status_code
        super().__init__(f"status={status_code}")


def _walk(node: object):
    if isinstance(node, dict):
        yield node
        for value in node.values():
            yield from _walk(value)
    elif isinstance(node, list):
        for value in node:
            yield from _walk(value)


def test_strict_schema_removes_sdk_metadata_and_keeps_required_objects() -> None:
    schema = strict_json_schema(TraitExtractionResponse)
    for node in _walk(schema):
        assert "default" not in node
        assert "title" not in node
        assert "examples" not in node
        if isinstance(node.get("properties"), dict):
            assert node["additionalProperties"] is False
            assert set(node["required"]) == set(node["properties"])


def test_retriable_status_codes_are_distinguished_from_invalid_requests() -> None:
    assert is_retriable_error(_StatusError(429))
    assert is_retriable_error(_StatusError(503))
    assert not is_retriable_error(_StatusError(400))


def test_failure_snapshot_has_a_valid_empty_candidate_csv(tmp_path) -> None:
    manifest = {
        "run_id": "run_test",
        "status": "partial_failure",
        "n_completed_batches": 0,
        "n_failed_batches": 1,
        "candidate_row_count": 0,
        "batches": [],
        "failures": [{"exception_type": "BadRequestError", "message": "schema rejected"}],
    }
    candidate_path, manifest_path = write_run_artifacts(
        output_dir=tmp_path,
        run_id="run_test",
        rows=[],
        manifest=manifest,
    )
    table = pd.read_csv(candidate_path, dtype=str)
    assert list(table.columns) == candidate_columns()
    assert table.empty
    assert json.loads(manifest_path.read_text(encoding="utf-8"))["status"] == "partial_failure"


def test_empty_candidate_snapshot_can_enter_logic_audit(tmp_path) -> None:
    manifest = {"run_id": "run_test", "status": "partial_failure"}
    candidate_path, _ = write_run_artifacts(
        output_dir=tmp_path,
        run_id="run_test",
        rows=[],
        manifest=manifest,
    )
    table = pd.read_csv(candidate_path, dtype=str).fillna("")
    ontology = yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))
    from island_v2.trait_candidate_logic import audit_table

    audited = audit_table(table, ontology)
    assert audited.empty
    assert "logic_status" in audited.columns
