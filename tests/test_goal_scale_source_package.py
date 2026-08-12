from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2.merge_reviewed_source_packages import merge_packages
from island_v2.open_web_common import reviewed_source_package_evidence

ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/goal_scale_source_package_20260811"
)
MANIFEST_PATH = ROOT / "goal_scale_source_package_20260811_manifest.json"


def _canonical_hash(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in {".csv", ".json", ".md", ".txt"}:
        payload = payload.replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


def test_goal_scale_package_passes_the_shared_trait_gate() -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    evidence = pd.read_csv(
        ROOT / "goal_scale_source_package_20260811_evidence.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "goal_scale_source_package_20260811_audit.csv.gz", dtype=str
    ).fillna("")

    selected, trait_audit, gate = reviewed_source_package_evidence(evidence, audit)

    assert len(evidence) == 31_262
    assert len(audit) == 2_836
    assert len(selected) == 31_227
    assert selected[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 30_984
    assert evidence["source_record_id"].is_unique
    assert audit["candidate_id"].is_unique
    assert gate["package_gate_passed"] is True
    assert gate["reviewed"] == 2_836
    assert gate["accepted_correct"] == 2_801
    assert gate["precision"] == pytest.approx(2_801 / 2_836)
    assert gate["cultivar_contamination_rate"] == 0.0
    assert gate["explicitly_rejected_candidate_ids"] == 35
    assert gate["approved_trait_count"] == 10
    assert trait_audit["production_approved"].all()
    assert selected["source_group"].value_counts().to_dict() == {
        "wang_2023_floral_symmetry": 28_048,
        "zell_2025_bsdb": 1_451,
        "latest_public_web": 666,
        "two_backbone_bulk_synonym_recovery": 439,
        "gobotany_completion_checkpoint": 333,
        "proteus": 290,
    }
    assert not set(selected["trait_name"]).intersection(
        {"reward_type", "pollen_vector_mode"}
    )
    assert manifest["source_package_gate"] == gate

    for filename, expected_hash in manifest["output_hashes"].items():
        assert _canonical_hash(ROOT / filename) == expected_hash


def test_goal_scale_package_rebuild_is_byte_reproducible(tmp_path: Path) -> None:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    input_packages = manifest["input_packages"]
    rebuilt = merge_packages(
        [Path(item["evidence_path"]) for item in input_packages],
        [Path(item["audit_path"]) for item in input_packages],
        tmp_path,
        "goal_scale_source_package_20260811",
        manifest["generated_at_utc"],
    )

    assert rebuilt["candidate_evidence_rows"] == 31_262
    assert rebuilt["audit_rows"] == 2_836
    assert rebuilt["selected_evidence_rows"] == 31_227
    assert rebuilt["selected_species_trait"] == 30_984
    assert rebuilt["output_hashes"] == manifest["output_hashes"]
    for item in input_packages:
        assert _canonical_hash(Path(item["evidence_path"])) == item["evidence_sha256"]
        assert _canonical_hash(Path(item["audit_path"])) == item["audit_sha256"]
