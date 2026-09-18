import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LOCK = ROOT / "config/chapter1_v14_canonical_result_lock.json"
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"
MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md"
FREEZE = ROOT / "docs/chapter1_submission_freeze_v14_20260918.md"


def test_v14_canonical_lock_has_verified_ci_provenance():
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    assert lock["contract"] == "chapter1_v14_canonical_result_lock_v1"
    assert lock["status"] == "canonical_v14_reproduced"
    assert lock["canonical"] is True

    provenance = lock["promotion_provenance"]
    assert provenance["head_sha"] == "7793c6b2a4568a609df2095cc796c367252ea68f"
    assert provenance["ci_run_id"] == 35314955780
    assert provenance["ci_artifact_id"] == 10535020072
    assert provenance["ci_artifact_digest"] == (
        "sha256:4665e68341cea35bfb16c33afeb30ccf2705b5a47982f81e409e8507204be811"
    )
    verification = provenance["verification"]
    assert verification["verified"] is True
    assert verification["H1_reproduced"] is True
    assert verification["H2_reproduced"] is True
    assert verification["H4_exact_H2_score_bridge_reproduced"] is True
    assert verification["H4_atomic_reconstruction_reproduced"] is True


def test_v14_is_publication_facing_and_v13_is_parent_provenance():
    readme = README.read_text(encoding="utf-8")
    pipeline = PIPELINE.read_text(encoding="utf-8")
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    freeze = FREEZE.read_text(encoding="utf-8")

    assert readme.startswith("# Island — Chapter 1 v14 paper repository")
    assert "current publication-facing surface is Chapter 1 v14" in readme
    assert "Frozen v13 parent paper surface" in readme
    assert "Chapter 1 v14 candidate reanalysis" not in readme

    assert pipeline.startswith("# Chapter 1 v14 paper pipeline")
    assert "Frozen v13 parent paper pipeline" in pipeline
    assert "Chapter 1 v14 candidate pipeline" not in pipeline

    assert "canonical reordered H1–H4 analysis" in manuscript
    assert "v14 candidate" not in manuscript
    assert "canonical v14 reproduced analysis surface" in freeze


def test_v14_surface_points_to_canonical_lock_and_freeze():
    expected = (
        "config/chapter1_v14_canonical_result_lock.json",
        "docs/chapter1_submission_freeze_v14_20260918.md",
    )
    for path in (README, PIPELINE, MANUSCRIPT):
        text = path.read_text(encoding="utf-8")
        for token in expected:
            assert token in text
