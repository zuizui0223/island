"""Fail-closed check of generated manuscript/SI proof packaging.

These tests verify source identity and the proof-only status; they do not
authorize a real Ecology Letters submission or redistribute trait records.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from zipfile import ZipFile

from pypdf import PdfReader

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "submission/chapter1_current"
PROOF = SOURCE / "editorial_proofs"


def _checksum(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_editorial_main_and_si_are_source_locked_and_not_submission_ready() -> None:
    audit = json.loads((PROOF / "EDITORIAL_PROOF_MANIFEST.json").read_text(encoding="utf-8"))
    assert audit["contract"] == "chapter1_q1_main_si_editorial_composition_v1"
    assert audit["status"] == "EDITORIAL_PROOF_NOT_SUBMISSION_READY"
    assert audit["scope"] == "layout_proof_only_no_numerical_reanalysis"
    for field in (
        "submission_authorized",
        "author_approval",
        "trait_redistribution_rights_cleared",
        "final_release_doi_minted",
    ):
        assert audit[field] is False
    assert audit["source_audit"]["original_main_body_words"] == 4675
    assert audit["source_audit"]["main_figure_captions"] == 5
    for name, expected_sha in audit["source_sha256"].items():
        assert _checksum(SOURCE / name) == expected_sha


def test_editorial_word_and_pdf_outputs_can_be_reopened() -> None:
    audit = json.loads((PROOF / "EDITORIAL_PROOF_MANIFEST.json").read_text(encoding="utf-8"))
    for family in ("main", "si"):
        files = audit[family]
        word = PROOF / files["docx"]
        pdf = PROOF / files["pdf"]
        assert _checksum(word) == files["docx_sha256"]
        assert _checksum(pdf) == files["pdf_sha256"]
        with ZipFile(word) as contents:
            assert contents.testzip() is None
            assert "word/document.xml" in contents.namelist()
        reader = PdfReader(pdf)
        assert len(reader.pages) == files["pdf_pages"]
        page_text = " ".join((p.extract_text() or "") for p in reader.pages[:3])
        assert "EDITORIAL PROOF - NOT FOR SUBMISSION" in page_text
    for family, count in (("main_plus_figures", 5), ("si_plus_figures", 7)):
        record = audit[family]
        combined = PROOF / record["file"]
        assert _checksum(combined) == record["sha256"]
        pages = len(PdfReader(combined).pages)
        assert pages == record["pages"]
        assert pages == audit["main" if count == 5 else "si"]["pdf_pages"] + count


def test_every_page_has_text_and_no_painted_outside_page_geometry() -> None:
    for role, expected in (("MAIN", 23), ("SI", 14)):
        report_path = PROOF / f"Q1_{role}_PAGES_AUDIT.json"
        report = json.loads(report_path.read_text(encoding="utf-8"))
        assert report["contract"] == "chapter1_q1_editorial_all_pages_audit_v1"
        assert report["every_page_nonempty_and_on_canvas"] is True
        assert report["total_pages"] == expected
        assert len(report["pages"]) == expected
        assert len(report["tiles"]) >= 3
        assert all(page["passes"] for page in report["pages"])
        assert all(page["text_characters"] >= 35 for page in report["pages"])
        assert all((PROOF / t).stat().st_size > 1000 for t in report["tiles"])
