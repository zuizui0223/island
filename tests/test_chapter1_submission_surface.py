import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V14_LOCK = ROOT / "config/chapter1_v14_canonical_result_lock.json"
CURRENT = ROOT / "config/chapter1_submission_current.json"
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"
DATABASE = ROOT / "docs/DATABASE_BUILD.md"
V14_MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md"
V14_FREEZE = ROOT / "docs/chapter1_submission_freeze_v14_20260918.md"
V13_MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
V13_FREEZE = ROOT / "docs/chapter1_submission_freeze_v13_20260917.md"


def test_v14_lock_is_preserved_as_verified_historical_provenance():
    lock = json.loads(V14_LOCK.read_text(encoding="utf-8"))
    assert lock["contract"] == "chapter1_v14_canonical_result_lock_v1"
    assert lock["status"] == "canonical_v14_reproduced"
    assert lock["canonical"] is True

    provenance = lock["promotion_provenance"]
    assert provenance["head_sha"] == "7793c6b2a4568a609df2095cc796c367252ea68f"
    assert provenance["ci_run_id"] == 35314955780
    assert provenance["ci_artifact_id"] == 10535020072
    assert provenance["verification"]["verified"] is True


def test_corrected_submission_is_the_only_active_chapter1_surface():
    current = json.loads(CURRENT.read_text(encoding="utf-8"))
    readme = README.read_text(encoding="utf-8")
    pipeline = PIPELINE.read_text(encoding="utf-8")

    assert current["status"] == "primary_submission_baseline"
    assert current["primary"] is True
    assert current["population"]["universe"] == 8264
    assert current["submission_package"] == "submission/chapter1_current"
    assert current["manuscript"] == "submission/chapter1_current/MANUSCRIPT.md"
    assert current["supersedes"] == "config/chapter1_v14_canonical_result_lock.json"

    assert readme.startswith("# Island — Chapter 1 corrected submission baseline")
    assert "submission/chapter1_current/MANUSCRIPT.md" in readme
    assert "8,264 island units" in readme
    assert "Superseded v14 surface (provenance)" in readme
    assert pipeline.startswith("# Current submission pipeline — corrected geography")
    assert "The **only active Chapter 1 submission selector**" in pipeline


def test_historical_surfaces_are_explicitly_marked_and_not_current():
    manuscript = V14_MANUSCRIPT.read_text(encoding="utf-8")
    freeze = V14_FREEZE.read_text(encoding="utf-8")
    database = DATABASE.read_text(encoding="utf-8")
    v13_manuscript = V13_MANUSCRIPT.read_text(encoding="utf-8")
    v13_freeze = V13_FREEZE.read_text(encoding="utf-8")

    assert manuscript.startswith("> **SUPERSEDED FOR SUBMISSION")
    assert "submission/chapter1_current/MANUSCRIPT.md" in manuscript
    assert freeze.startswith("> **SUPERSEDED SUBMISSION FREEZE")
    assert v13_manuscript.startswith("> **SUPERSEDED FOR SUBMISSION")
    assert v13_freeze.startswith("> **SUPERSEDED SUBMISSION FREEZE")
    assert "8,264 = current corrected Chapter 1 analysis universe" in database
    assert "8,265 = historical alpha1 / frozen provenance universe" in database
