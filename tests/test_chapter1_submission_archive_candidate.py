import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "zenodo" / "chapter1_submission_zenodo_metadata.json"
RELEASE = ROOT / "zenodo" / "CHAPTER1_SUBMISSION_RELEASE.md"


def test_submission_archive_metadata_is_draft_and_rights_safe() -> None:
    meta = json.loads(META.read_text(encoding="utf-8"))
    release = RELEASE.read_text(encoding="utf-8")

    assert meta["upload_type"] == "software"
    assert meta["title"].startswith("island: Chapter 1 corrected submission")
    assert any(
        item["identifier"] == "10.5281/zenodo.22704973"
        for item in meta["related_identifiers"]
    )
    assert "DRAFT METADATA ONLY" in meta["notes"]

    assert "prepared but not published" in release
    assert "full 222,688-cell Chapter 1 trait ledger" in release
    assert "Exclude from this software/results archive" in release
    assert "chapter1-submission-v1" in release
    assert "explicitly approves the external release" in release


def test_submission_archive_candidate_keeps_dataset_and_software_roles_separate() -> None:
    meta = json.loads(META.read_text(encoding="utf-8"))
    release = RELEASE.read_text(encoding="utf-8")

    assert meta["upload_type"] != "dataset"
    assert "rights-filtered trait derivative" in meta["description"]
    assert "10.5281/zenodo.22704973" in meta["description"]
    assert "software + corrected derived results" in release
    assert "rights-filtered trait derivative" in release
    assert "GloPL source data" in release
