from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "submission" / "chapter1_current"
MANUSCRIPT = PACKAGE / "MANUSCRIPT.md"
CAPTIONS = PACKAGE / "FIGURE_CAPTIONS.md"
COVER = PACKAGE / "COVER_LETTER_DRAFT.md"
TITLE_PAGE = PACKAGE / "TITLE_PAGE_TEMPLATE.md"
TARGET = PACKAGE / "ECOLOGY_LETTERS_TARGET.md"
NOVELTY = PACKAGE / "NOVELTY_STATEMENT.md"
DATA_ACCESS = PACKAGE / "DATA_ACCESSIBILITY_DRAFT.md"
SUPPLEMENT = PACKAGE / "SUPPLEMENT_PLAN.md"
GRAPHICAL = PACKAGE / "GRAPHICAL_ABSTRACT_BRIEF.md"
DATA_GATE = PACKAGE / "ECOLOGY_LETTERS_DATA_GATE.md"
DATA_INQUIRY = PACKAGE / "ECOLOGY_LETTERS_DATA_POLICY_INQUIRY.md"
FIGURE3 = PACKAGE / "figures" / "Figure3_constraint_response_triangle.svg"


def _words(text: str) -> int:
    return len(text.strip().split())


def test_ecology_letters_letter_limits_and_title_sync() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    title = manuscript.splitlines()[0].removeprefix("# ").strip()

    abstract_start = manuscript.index("## Abstract") + len("## Abstract")
    abstract_end = manuscript.index("**Keywords:**")
    intro_start = manuscript.index("## Introduction")
    refs_start = manuscript.index("## References")

    assert _words(manuscript[abstract_start:abstract_end]) <= 150
    assert _words(manuscript[intro_start:refs_start]) <= 5000

    captions = CAPTIONS.read_text(encoding="utf-8")
    assert captions.count("## Figure ") <= 6

    assert title in COVER.read_text(encoding="utf-8")
    assert title in TITLE_PAGE.read_text(encoding="utf-8")


def test_ecology_letters_title_page_metadata() -> None:
    text = TITLE_PAGE.read_text(encoding="utf-8")
    running = next(
        line.split(":", 1)[1].strip()
        for line in text.splitlines()
        if line.startswith("**Running title:**")
    )
    keywords = next(
        line.split(":", 1)[1].strip()
        for line in text.splitlines()
        if line.startswith("**Keywords:**")
    )

    assert len(running) < 45
    assert len([item for item in keywords.split(";") if item.strip()]) <= 10
    assert "**Article type:** Letter" in text
    assert "**Figures:** 6" in text
    assert "**Tables:** 0" in text
    assert "**Text boxes:** 0" in text


def test_ecology_letters_supporting_submission_files_exist() -> None:
    for path in [TARGET, NOVELTY, DATA_ACCESS, SUPPLEMENT, GRAPHICAL, DATA_GATE, DATA_INQUIRY, FIGURE3]:
        assert path.is_file(), path
        assert path.read_text(encoding="utf-8").strip(), path

    target = TARGET.read_text(encoding="utf-8")
    assert "maximum 5,000 words" in target
    assert "maximum 150 words" in target
    assert "Global Ecology and Biogeography" in target

    novelty = NOVELTY.read_text(encoding="utf-8")
    assert "constraint" in novelty.lower()
    assert "historical pollen limitation" in novelty

    data_access = DATA_ACCESS.read_text(encoding="utf-8")
    assert "cannot be redistributed wholesale" in data_access
    assert "46,274" in data_access


def test_ecology_letters_data_policy_gate_is_fail_closed() -> None:
    gate = DATA_GATE.read_text(encoding="utf-8")
    data_access = DATA_ACCESS.read_text(encoding="utf-8")
    target = TARGET.read_text(encoding="utf-8")

    assert "222,688" in gate
    assert "46,274" in gate
    assert "176,414" in gate
    assert "Do not submit to Ecology Letters while this gate is unresolved." in gate
    assert "Not submission-ready for Ecology Letters" in data_access
    assert "conditional on data-policy clearance" in target.lower()


def test_constraint_response_triangle_matches_current_inference() -> None:
    svg = FIGURE3.read_text(encoding="utf-8")

    for label in (
        "Geographic isolation",
        "Plant response (H1–H2)",
        "Pollen limitation (H3)",
        "H4: β &lt; 0",
        "Historical causal edge not identified",
        "not a mediation model",
    ):
        assert label in svg
