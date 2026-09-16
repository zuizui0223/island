from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"
MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
FIGURE_SYNC = ROOT / "docs/chapter1_v13_submission_figure_sync_20260917.md"
LOCK = ROOT / "config/chapter1_v13_unified_island_syndrome_result_lock.json"

V13_MANUSCRIPT = "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
V13_LOCK = "config/chapter1_v13_unified_island_syndrome_result_lock.json"
V13_FIGURES = "docs/chapter1_v13_submission_figure_sync_20260917.md"
V11_MANUSCRIPT = "docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md"


def test_readme_promotes_v13_and_retains_v11_as_historical() -> None:
    text = README.read_text(encoding="utf-8")
    assert V13_MANUSCRIPT in text
    assert V13_LOCK in text
    assert V13_FIGURES in text
    assert V11_MANUSCRIPT in text
    assert text.index(V13_MANUSCRIPT) < text.index(V11_MANUSCRIPT)
    assert "historical" in text.casefold()


def test_pipeline_promotes_v13_first() -> None:
    text = PIPELINE.read_text(encoding="utf-8")
    assert V13_MANUSCRIPT in text
    assert V13_LOCK in text
    assert V13_FIGURES in text
    assert V11_MANUSCRIPT in text
    assert text.index(V13_MANUSCRIPT) < text.index(V11_MANUSCRIPT)


def test_manuscript_contains_locked_headlines() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    assert "0.0794" in manuscript
    assert "-0.4467" in manuscript
    assert "4/4 -> 4/4 -> 0/4" in manuscript
    assert str(lock["parents"]["functional_bridge"]["artifact_id"]) in manuscript
    assert "post-hoc functional triangulation" in manuscript.casefold()


def test_manuscript_does_not_make_prohibited_positive_causal_claims() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8").casefold()
    prohibited = [
        "pollen limitation caused trait evolution",
        "pollen limitation caused the observed trait evolution",
        "globi proves",
        "global pollinator abundance declines with isolation",
        "genus attenuation proves dispersal",
    ]
    for phrase in prohibited:
        assert phrase not in text


def test_figure_sync_uses_evidence_role_labels() -> None:
    text = FIGURE_SYNC.read_text(encoding="utf-8")
    for label in (
        "frozen primary / confirmatory parent result",
        "frozen sensitivity / robustness",
        "post-hoc functional triangulation",
        "descriptive synthesis",
    ):
        assert label in text
