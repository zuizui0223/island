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
PROSPECTIVE_SUPPORT_LOCK = ROOT / "config/chapter1_h4_prospective_temporal_support_decision_lock.json"

V13_MANUSCRIPT = "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
V13_LOCK = "config/chapter1_v13_unified_island_syndrome_result_lock.json"
V13_FIGURES = "docs/chapter1_v13_submission_figure_sync_20260917.md"
V11_MANUSCRIPT = "docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md"
LEGACY_DIR = "legacy/chapter1-pre-v13/"


def test_readme_is_v13_only_and_points_pre_v13_to_legacy() -> None:
    text = README.read_text(encoding="utf-8")
    assert V13_MANUSCRIPT in text
    assert V13_LOCK in text
    assert V13_FIGURES in text
    assert LEGACY_DIR in text
    assert V11_MANUSCRIPT not in text
    assert "H1–H4" in text or "H1-H4" in text
    for phrase in (
        "frozen H1-H5 hypothesis contract",
        "P1 assembly",
        "P2 component",
        "P3 joint",
        "chapter1_v10_figure2_p2_result_lock.json",
        "North–Tropical",
        "Palearctic",
    ):
        assert phrase not in text


def test_pipeline_is_v13_only_and_points_pre_v13_to_legacy() -> None:
    text = PIPELINE.read_text(encoding="utf-8")
    assert V13_MANUSCRIPT in text
    assert V13_LOCK in text
    assert V13_FIGURES in text
    assert LEGACY_DIR in text
    assert V11_MANUSCRIPT not in text
    assert "H1–H4" in text or "H1-H4" in text
    for phrase in (
        "H1-H5 scientific contract",
        "P1 assembly-inference defense",
        "P2 component-support defense",
        "P3 joint observation-bias defense",
        "chapter1_v10_figure2_p2_result_lock.json",
        "North–Tropical",
        "Palearctic",
    ):
        assert phrase not in text


def test_current_handoff_is_global_only() -> None:
    readme = README.read_text(encoding="utf-8").casefold()
    pipeline = PIPELINE.read_text(encoding="utf-8").casefold()
    for text in (readme, pipeline):
        for phrase in ("genus structuring", "lineage-assembly", "h5d", "at what lineage-assembly level"):
            assert phrase not in text
        assert "pollen limitation" in text
        assert "reproductive assurance" in text
        assert "floral accessibility" in text


def test_manuscript_contains_locked_global_headlines() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    lock = json.loads(LOCK.read_text(encoding="utf-8"))
    assert "0.0794" in manuscript
    assert "-0.4467" in manuscript
    assert str(lock["parents"]["functional_bridge"]["artifact_id"]) in manuscript
    assert "post-hoc functional triangulation" in manuscript.casefold()
    assert "all four" in manuscript.casefold()


def test_manuscript_reports_prospective_h4_support_stop_without_relabelling() -> None:
    text = MANUSCRIPT.read_text(encoding="utf-8").casefold()
    support = json.loads(PROSPECTIVE_SUPPORT_LOCK.read_text(encoding="utf-8"))

    assert "prospective post-2015 temporal validation audit" in text
    assert "9 matched species across 6 publications" in text
    assert "9 matched species across 4 publications" in text
    assert "outcomes remained unopened" in text
    assert "support insufficiency" in text
    assert "post-hoc" in text

    decision = support["decision"]
    ceiling = support["claim_ceiling"]
    assert decision["outcome_extraction_authorized"] is False
    assert decision["threshold_relaxation_authorized"] is False
    assert ceiling["existing_v13_H4_remains_posthoc"] is True
    assert ceiling["support_insufficient"] is True
    assert ceiling["post2015_outcomes_remain_unopened_for_primary_analysis"] is True

    freeze = FREEZE.read_text(encoding="utf-8").casefold()
    assert "post-freeze prospective h4 audit" in freeze
    assert "support insufficiency, not a biological null" in freeze
    assert "post-2015 primary outcomes remain unopened" in freeze


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
