"""GEB submission-specific editorial validation; no scientific refits."""
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GEB = ROOT / "submission/geb_current"


def test_geb_is_first_shot_without_changing_science():
    config = json.loads((ROOT / "config/chapter1_geb_submission.json").read_text())
    science = json.loads((ROOT / config["science_selector"]).read_text())
    assert config["target_journal"] == "Global Ecology and Biogeography"
    assert config["status"] == "GEB_FIRST_SHOT_DRAFT_NOT_SUBMISSION_READY"
    assert config["science_lock_unchanged"] is True
    assert science["status"] == "primary_submission_baseline"
    assert science["population"]["universe"] == 8264
    assert science["population"]["broad_H1_union"] == 4379
    assert science["final_H1_inference"]["all_112_converged"] is True


def test_geb_structured_abstract_keywords_and_blind_text():
    text = (GEB / "MANUSCRIPT_BLINDED.md").read_text()
    assert text.startswith("# Island isolation")
    abstract = text.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0]
    headings = ["Aim", "Location", "Time period", "Major taxa studied",
                "Methods", "Results", "Main conclusions"]
    for title in headings:
        assert abstract.count("**" + title + ":**") == 1
    assert len(abstract.split()) <= 300
    kw = text.split("**Keywords:**", 1)[1].split("\n", 1)[0]
    labels = [x.strip() for x in kw.split(",")]
    assert 6 <= len(labels) <= 10
    assert labels == sorted(labels, key=str.casefold)
    body = "## Introduction" + text.split("## Introduction", 1)[1].split("## References")[0]
    assert len(body.split()) <= 5000
    running = text.split("**Running title:**", 1)[1].split("\n", 1)[0].strip()
    assert len(running) < 40
    assert not re.search(r"zuizui0223|github\.com/zuizui|@gmail", text, re.I)
    assert "Data and Code Availability Statement" in text
    assert "176,414" in text and "46,274" in text
    assert "PENDING" in text
    assert "β = 0.22031" in text and "β = 0.09191" in text
    assert "cannot be directly compared" in text
    assert "not experimentally assigned" in text


def test_geb_anonymity_and_archive_gates():
    metadata = (GEB / "TITLE_PAGE_IDENTIFYING_TEMPLATE.md").read_text()
    assert "ONE CORRESPONDING AUTHOR" not in metadata  # template uses ONE name placeholder
    assert "exactly one" in metadata.lower()
    assert "[AUTHOR NAMES AND ORDER" in metadata
    data = (GEB / "DATA_CODE_AVAILABILITY.md").read_text()
    assert "editorial discretion" in data.lower()
    assert "rights review" in data.lower()
    assert "not" in data.lower()
    assert "10.5281/zenodo.22704973" in data
    config = json.loads((ROOT / "config/chapter1_geb_submission.json").read_text())
    assert all(x == "pending" for x in config["submission_gates"].values())


def test_geb_figure_captions_and_map_source():
    captions = (GEB / "FIGURE_CAPTIONS_GEB.md").read_text()
    assert captions.count("## Figure ") == 6
    assert "Mollweide equal-area" in captions
    assert "(a)" in captions and "(b)" in captions
    src = (ROOT / "scripts/submission/render_chapter1_geb_figures.py").read_text()
    assert 'projection="mollweide"' in src
    assert "EQUATOR_5KKM_DEGREES" in src
    assert "scientific_refit" in src


def test_geb_package_claim_limits_preserved():
    checklist = (GEB / "SUBMISSION_CHECKLIST_GEB.md").read_text()
    assert "No external Wiley submission" in checklist
    manuscript = (GEB / "MANUSCRIPT_BLINDED.md").read_text()
    assert "historical pollinator declines and causal mediation" in manuscript
    assert "post-hoc" in manuscript
    assert "WCVP" in manuscript
