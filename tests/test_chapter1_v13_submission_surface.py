from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"
MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
FIGURE_SYNC = ROOT / "docs/chapter1_v13_submission_figure_sync_20260917.md"
HYPOTHESIS = ROOT / "docs/chapter1_unified_hypothesis_20260917.md"
FREEZE = ROOT / "docs/chapter1_submission_freeze_v13_20260917.md"
LOCK = ROOT / "config/chapter1_v13_unified_island_syndrome_result_lock.json"

V13_MANUSCRIPT = "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
V13_LOCK = "config/chapter1_v13_unified_island_syndrome_result_lock.json"
V13_FIGURES = "docs/chapter1_v13_submission_figure_sync_20260917.md"
V11_MANUSCRIPT = "docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md"


def _section(text: str, start: str, end: str) -> str:
    a = text.index(start)
    b = text.index(end, a)
    return text[a:b]


def test_readme_promotes_global_only_v13_and_retains_v11_as_historical() -> None:
    text = README.read_text(encoding="utf-8")
    assert V13_MANUSCRIPT in text
    assert V13_LOCK in text
    assert V13_FIGURES in text
    assert V11_MANUSCRIPT in text
    assert text.index(V13_MANUSCRIPT) < text.index(V11_MANUSCRIPT)
    current = _section(text, "## 6. Current submission surface", "## 7. Chapter 1 / Chapter 2 division of labour")
    assert "H1–H4" in current or "H1-H4" in current
    for phrase in ("H5 — taxonomic realization", "Palearctic", "North–Tropical", "regional modification", "Figure 1–4"):
        assert phrase not in current


def test_pipeline_promotes_global_only_v13() -> None:
    text = PIPELINE.read_text(encoding="utf-8")
    assert V13_MANUSCRIPT in text
    assert V13_LOCK in text
    assert V13_FIGURES in text
    assert V11_MANUSCRIPT in text
    current = _section(text, "## 13. Canonical paper surface", "## 14. Chapter 1 / Chapter 2 handoff")
    assert "H1–H4" in current or "H1-H4" in current
    for phrase in ("H5: flora-layer-specific taxonomic realization", "North–Tropical", "Palearctic", "regional modification"):
        assert phrase not in current


def test_manuscript_contains_locked_global_headlines() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    assert "0.0794" in manuscript
    assert "-0.4467" in manuscript
    assert str(lock["parents"]["functional_bridge"]["artifact_id"]) in manuscript
    assert "post-hoc functional triangulation" in manuscript.casefold()
    assert "all four" in manuscript.casefold()


def test_current_submission_surface_has_no_regional_or_taxonomic_branch() -> None:
    files = (MANUSCRIPT, FIGURE_SYNC, HYPOTHESIS, FREEZE)
    prohibited = (
        "North–Tropical",
        "Palearctic",
        "taxonomic realization",
        "regional modification",
        "4/4 -> 4/4 -> 0/4",
        "Figure 4",
        "H5 —",
        "### H5",
    )
    for path in files:
        text = path.read_text(encoding="utf-8")
        for phrase in prohibited:
            assert phrase not in text, f"{phrase!r} remains in {path}"


def test_manuscript_does_not_make_prohibited_positive_causal_claims() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8").casefold()
    prohibited = [
        "pollen limitation caused trait evolution",
        "pollen limitation caused the observed trait evolution",
        "globi proves",
        "global pollinator abundance declines with isolation",
    ]
    for phrase in prohibited:
        assert phrase not in text


def test_figure_sync_has_three_main_figures_and_evidence_roles() -> None:
    text = FIGURE_SYNC.read_text(encoding="utf-8")
    assert "Main Figure 1" in text
    assert "Main Figure 2" in text
    assert "Main Figure 3" in text
    assert "Main Figure 4" not in text
    for label in (
        "frozen primary / confirmatory parent result",
        "frozen sensitivity / robustness",
        "post-hoc functional triangulation",
        "descriptive synthesis",
    ):
        assert label in text
