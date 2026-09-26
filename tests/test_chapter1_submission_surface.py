import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
V14_LOCK = ROOT / "config/chapter1_v14_canonical_result_lock.json"
CURRENT = ROOT / "config/chapter1_submission_current.json"
README = ROOT / "README.md"
PIPELINE = ROOT / "docs/PAPER_PIPELINE.md"
DATABASE = ROOT / "docs/DATABASE_BUILD.md"
THESIS = ROOT / "THESIS_CHAPTER_POSITIONING.md"
ANALYSIS_V2 = ROOT / "analysis/v2/README.md"
PROGRESSIVE_CONFIG = ROOT / "config/chapter1_progressive_analysis.yml"
ALL_DATA_PROGRESSIVE_CONFIG = ROOT / "config/chapter1_all_data_progressive_analysis.yml"
DATABASE_RELEASE = ROOT / "docs/CHAPTER1_DATABASE_RELEASE.md"
DATABASE_VERSIONS = ROOT / "config/chapter1_database_versions/README.md"
DATA_V2 = ROOT / "data/v2/README.md"
V14_MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v14_reordered_hypotheses_20260918.md"
V14_FREEZE = ROOT / "docs/chapter1_submission_freeze_v14_20260918.md"
V13_MANUSCRIPT = ROOT / "docs/chapter1_manuscript_full_v13_global_pollination_constraint_20260917.md"
V13_FREEZE = ROOT / "docs/chapter1_submission_freeze_v13_20260917.md"

HISTORICAL_WORKING_DOCS = [
    ROOT / "docs/chapter1_v9_submission_figure_sync_20260915.md",
    ROOT / "docs/chapter1_submission_freeze_20260915_p1_defended.md",
    ROOT / "docs/chapter1_submission_freeze_20260915_p1_p2_defended.md",
    ROOT / "docs/chapter1_manuscript_v8_reframing_20260913.md",
    ROOT / "docs/chapter1_manuscript_full_v7_submission_order_20260909.md",
    ROOT / "docs/chapter1_manuscript_full_v8_hierarchical_syndrome_20260914.md",
    ROOT / "docs/chapter1_manuscript_full_v9_island_first_p1_defended_20260915.md",
    ROOT / "docs/chapter1_manuscript_full_v10_island_first_p1_p2_defended_20260915.md",
    ROOT / "docs/chapter1_manuscript_full_v11_island_first_p1_p2_p3_defended_20260915.md",
    ROOT / "docs/chapter1_latest_trait_reanalysis_20260908.md",
    ROOT / "docs/chapter1_submission_freeze_20260909.md",
    ROOT / "docs/chapter1_submission_freeze_20260915_p1_p2_p3_defended.md",
    ROOT / "docs/chapter1_v8_submission_figure_sync_20260915.md",
    ROOT / "docs/chapter1_v10_submission_figure_sync_20260915.md",
    ROOT / "docs/chapter1_v11_submission_figure_sync_20260915.md",
    ROOT / "docs/chapter1_v13_submission_figure_sync_20260917.md",
]


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


def test_current_scientific_guidance_matches_corrected_h1_h4_surface():
    thesis = THESIS.read_text(encoding="utf-8")
    assert "Current scientific surface: corrected geography baseline" in thesis
    assert "### H1 — recurrent multivariate island response" in thesis
    assert "### H2 — two partially separable plant-response components" in thesis
    assert "### H3 — independent ecological-pressure correlate" in thesis
    assert "### H4 — functional compatibility" in thesis
    assert "What Chapter 1 no longer claims" in thesis

    for path in [ANALYSIS_V2, PROGRESSIVE_CONFIG, ALL_DATA_PROGRESSIVE_CONFIG, DATABASE_RELEASE, DATABASE_VERSIONS]:
        prefix = path.read_text(encoding="utf-8")[:600].lower()
        assert "historical" in prefix, path
        assert "chapter1_submission_current.json" in prefix, path


def test_historical_execution_docs_do_not_present_stale_canonical_routes():
    analysis = ANALYSIS_V2.read_text(encoding="utf-8")
    database_release = DATABASE_RELEASE.read_text(encoding="utf-8")
    database_versions = DATABASE_VERSIONS.read_text(encoding="utf-8")
    data_v2 = DATA_V2.read_text(encoding="utf-8")

    assert "Canonical workflow:" not in analysis
    assert "run-chapter1-context-main.yml" not in analysis
    assert "There is **no active canonical Chapter 1 workflow selected from this directory**" in analysis

    assert "The analysis contract remains:" not in database_release
    assert "run-chapter1-progressive-trait-analysis.yml" not in database_release
    assert "config/chapter1_submission_current.json" in database_release

    assert "active execution pointer" not in database_versions
    assert "legacy database-version pointer" in database_versions
    assert "config/chapter1_submission_current.json" in database_versions

    assert "does **not** make Bombus a predictor or mechanism" in data_v2
    assert "config/chapter1_submission_current.json" in data_v2[:800]


def test_pre_corrected_working_documents_are_visibly_historical():
    for path in HISTORICAL_WORKING_DOCS:
        text = path.read_text(encoding="utf-8")
        assert text.startswith("> **HISTORICAL / SUPERSEDED"), path
        assert "config/chapter1_submission_current.json" in text[:700], path
        assert "submission/chapter1_current/MANUSCRIPT.md" in text[:700], path
