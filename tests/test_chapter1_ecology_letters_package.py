import xml.etree.ElementTree as ET
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
GRAPHICAL_SHORT = PACKAGE / "GRAPHICAL_ABSTRACT_SHORT_TEXT.md"
GEB_FALLBACK = PACKAGE / "GEB_FALLBACK.md"
SI_DRAFT = PACKAGE / "SUPPLEMENTARY_INFORMATION_DRAFT.md"
H1_CONVERGENCE_AUDIT = ROOT / "results" / "geography_20260924" / "h1_direct_northern_high_convergence_audit.json"
H1_DOMAIN_TABLE = PACKAGE / "supplement" / "Table_S2e_H1_domain_descriptive.csv"
H1_RAW_AXIS_TABLE = PACKAGE / "supplement" / "Table_S2f_H1_three_axis_primary.csv"
H1_ORIGIN_TABLE = PACKAGE / "supplement" / "Table_S2g_H1_floristic_origin.csv"
H1_AXIS_AUDIT = PACKAGE / "supplement" / "Table_S2h_H1_axis_cell_audit.csv"


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


def test_h1_uses_final_directional_score_and_retains_three_axis_audit() -> None:
    import pandas as pd

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    captions = CAPTIONS.read_text(encoding="utf-8")
    novelty = NOVELTY.read_text(encoding="utf-8")

    assert "222,688 resolved cells" in manuscript
    assert "seven binary traits is fitted separately" in manuscript
    assert "no pooled floral-syndrome score" in manuscript
    assert "Seven traits under final finite-cluster inference" in captions
    assert "reproductive function" in novelty.lower()
    assert "phenotype" in novelty.lower()

    final_h1 = pd.read_csv(
        ROOT / "results/h1_final_traitwise_t_20261004/traitwise_results.csv"
    )
    assert len(final_h1) == 112
    assert final_h1["optimizer_success"].all()

    primary = final_h1.loc[
        final_h1["evidence_scope"].eq("all")
        & final_h1["flora_scope"].eq("broad")
        & final_h1["outcome"].eq("self_compatibility")
    ]
    assert len(primary) == 4
    assert (primary["estimate"] > 0).all()
    assert (primary["p_two_sided"] < 0.05).all()

    wcvp = final_h1.loc[
        final_h1["evidence_scope"].eq("all")
        & final_h1["flora_scope"].eq("wcvp")
        & final_h1["outcome"].eq("self_compatibility")
    ]
    assert len(wcvp) == 4
    assert (wcvp["estimate"] > 0).all()
    assert (wcvp["p_two_sided"] < 0.05).all()

def test_title_page_word_counts_match_manuscript() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    title_page = TITLE_PAGE.read_text(encoding="utf-8")
    abstract_start = manuscript.index("## Abstract") + len("## Abstract")
    abstract_end = manuscript.index("**Keywords:**")
    intro_start = manuscript.index("## Introduction")
    refs_start = manuscript.index("## References")
    abstract_words = _words(manuscript[abstract_start:abstract_end])
    main_words = _words(manuscript[intro_start:refs_start])

    assert f"**Abstract word count:** {abstract_words} " in title_page
    assert f"**Main-text word count:** {main_words:,} " in title_page


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
    for path in [TARGET, NOVELTY, DATA_ACCESS, SUPPLEMENT, GRAPHICAL, DATA_GATE, DATA_INQUIRY, FIGURE3, GRAPHICAL_SHORT, GEB_FALLBACK, SI_DRAFT]:
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
    assert "10.5281/zenodo.22704973" in data_access


def test_ecology_letters_data_policy_gate_is_fail_closed() -> None:
    gate = DATA_GATE.read_text(encoding="utf-8")
    data_access = DATA_ACCESS.read_text(encoding="utf-8")
    target = TARGET.read_text(encoding="utf-8")

    assert "222,688" in gate
    assert "46,274" in gate
    assert "176,414" in gate
    assert "10.5281/zenodo.22704973" in gate
    assert "Do not submit to Ecology Letters while this gate is unresolved." in gate
    assert "Not submission-ready for Ecology Letters" in data_access
    assert "conditional on data-policy clearance" in target.lower()

    inquiry = DATA_INQUIRY.read_text(encoding="utf-8")
    assert "ecolets2@cefe.cnrs.fr" in inquiry
    assert "ecolets@cefe.cnrs.fr" in inquiry
    assert "do not explicitly state whether third-party legal/licensing restrictions qualify" in inquiry

    assert "do not explicitly guarantee an exception for third-party licensing restrictions" in gate.lower()


def test_constraint_response_triangle_matches_current_inference() -> None:
    tree = ET.parse(FIGURE3)
    assert tree.getroot().tag.endswith("svg")
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


def test_graphical_abstract_short_text_and_geb_fallback() -> None:
    short = GRAPHICAL_SHORT.read_text(encoding="utf-8")
    body = short.split("\n\n", 1)[1].strip()
    assert len(body) <= 500

    geb = GEB_FALLBACK.read_text(encoding="utf-8")
    for heading in (
        "**Aim:**",
        "**Location:**",
        "**Time period:**",
        "**Major taxa studied:**",
        "**Methods:**",
        "**Results:**",
        "**Main conclusions:**",
    ):
        assert heading in geb
    assert "double anonymous" in geb.lower()
    assert "legal requirements" in geb


def test_supplementary_information_draft_is_bound_to_corrected_outputs() -> None:
    text = SI_DRAFT.read_text(encoding="utf-8")
    required = [
        "corrected 8,264-unit",
        "results/geography_20260924/all/beta_binomial_within_slopes.csv",
        "results/geography_20260924/direct/h2_decomposition_models.csv",
        "results/geography_20260924/h3_original_corrected_comparison.json",
        "results/geography_20260924/h4_exact_corrected.csv",
        "ECOLOGY_LETTERS_DATA_GATE.md",
    ]
    for marker in required:
        assert marker in text


def test_ecology_letters_reference_list_and_glopl_data_citation() -> None:
    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    title_page = TITLE_PAGE.read_text(encoding="utf-8")
    refs = manuscript.split("## References", 1)[1]

    entries = [
        line
        for line in refs.splitlines()
        if line.strip() and line[0].isalpha() and "(" in line and ")." in line
    ]
    assert len(entries) == 18
    assert "**References:** 18" in title_page
    assert "Bennett et al. 2018a,b" in manuscript
    assert "10.5061/dryad.dt437" in refs
    assert "Bennett, J.M., Steets, J.A., Burns, J.H., Durka, W., Vamosi, J.C., Arceo-Gómez, G. et al. (2018a)." in refs
    assert "Fenster, C.B., Armbruster, W.S., Wilson, P., Dudash, M.R. & Thomson, J.D. (2004)." in refs

    assert "Zell et al. 2025" in manuscript
    assert "Dawson-Glass & Hargreaves 2022" in manuscript
    assert "10.1111/nph.20234" not in manuscript  # journal style omits DOI for article refs


def test_h1_direct_northern_high_convergence_audit_closes_warning() -> None:
    import json
    import pandas as pd

    final_h1 = pd.read_csv(
        ROOT / "results/h1_final_traitwise_t_20261004/traitwise_results.csv"
    )
    assert len(final_h1) == 112
    assert final_h1["optimizer_success"].all()

    audit = json.loads(H1_CONVERGENCE_AUDIT.read_text(encoding="utf-8"))
    assert audit["decision"]["audit_pass"] is True
    assert audit["independent_multistart_confirmation"]["status"] == "pass"

    manuscript = MANUSCRIPT.read_text(encoding="utf-8")
    assert "All 112 fits converged" in manuscript
    assert "one-dimensional score" not in manuscript
    assert "2.93 × 10^-6" not in manuscript
    assert "1.43 × 10^-8" not in manuscript

