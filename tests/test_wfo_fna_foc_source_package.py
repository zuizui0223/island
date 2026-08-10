import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from analysis.prepare_wfo_fna_foc_source_package_20260810 import run
from island_v2.open_web_common import reviewed_source_package_evidence

ROOT = Path("data/v2/staging/traits/direct_llm_pilot")
FIRST = ROOT / "20260810_wfo_high_yield_source_acquisition"
PACKAGE = ROOT / "20260810_wfo_fna_foc_africa_source_acquisition"


def test_combined_package_approves_only_audited_high_precision_traits() -> None:
    evidence = pd.read_csv(
        PACKAGE / "wfo_combined_high_yield_evidence.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        PACKAGE / "wfo_combined_high_yield_manual_audit_320.csv", dtype=str
    ).fillna("")
    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert len(evidence) == 10367
    assert len(audit) == 320
    assert summary["accepted_correct"] == 305
    assert summary["precision"] == 0.953125
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"] is True
    assert set(summary["approved_traits"]) == {
        "floral_form",
        "flower_primary_color",
        "inflorescence_display",
    }
    size = scopes.loc[scopes["trait_name"].eq("flower_size_class")].iloc[0]
    assert size["precision"] == 0.8
    assert bool(size["production_approved"]) is False
    assert len(selected) == 9337
    assert selected["accepted_species"].nunique() == 7617
    assert len(selected[["accepted_species", "trait_name"]].drop_duplicates()) == 9049
    assert len(selected[["accepted_species", "axis"]].drop_duplicates()) == 8484


def test_committed_package_is_exactly_reproducible(tmp_path: Path) -> None:
    summary = run(
        argparse.Namespace(
            extracted_evidence=PACKAGE
            / "wfo_fna_foc_africa_extracted_evidence.csv.gz",
            audit_seed=PACKAGE / "wfo_fna_foc_africa_audit_seed_120.csv",
            manual_audit=PACKAGE / "wfo_fna_foc_africa_manual_audit_120.csv",
            first_evidence=FIRST / "wfo_high_yield_incremental_evidence.csv.gz",
            first_audit=FIRST / "wfo_high_yield_manual_audit_200.csv",
            output_dir=tmp_path,
        )
    )
    assert summary["incremental_candidate_rows"] == 3972
    assert summary["incremental_species_trait"] == 3920
    assert summary["combined_candidate_rows"] == 10367
    assert summary["combined_reviewed"] == 320
    regenerated = pd.read_csv(
        tmp_path / "wfo_combined_high_yield_evidence.csv.gz", dtype=str
    ).fillna("")
    committed = pd.read_csv(
        PACKAGE / "wfo_combined_high_yield_evidence.csv.gz", dtype=str
    ).fillna("")
    assert regenerated["source_record_id"].tolist() == committed[
        "source_record_id"
    ].tolist()


def test_source_manifest_pins_official_archives_and_backbone() -> None:
    manifest = json.loads((PACKAGE / "source_run_manifest.json").read_text())
    assert manifest["fixed_universe_species"] == 106295
    assert manifest["matched_species_descriptions"] == 14819
    assert manifest["backbone"]["sha256"] == (
        "ccfc7e3ca85c8b80aa96eb52d4e50e07e1811d28916b544280d8a4f1a043ba72"
    )
    assert {row["provider"] for row in manifest["raw_archives"]} == {
        "fna",
        "flora_of_china",
        "flora_d_afrique",
    }


def test_committed_file_manifest_matches_bytes_and_hashes() -> None:
    manifest = json.loads((PACKAGE / "file_manifest.json").read_text())
    for row in manifest["files"]:
        payload = (PACKAGE / row["file"]).read_bytes()
        canonical = payload if row["file"].endswith(".gz") else payload.replace(
            b"\r\n", b"\n"
        )
        assert len(canonical) == row["bytes"]
        assert hashlib.sha256(canonical).hexdigest() == row["sha256"]
