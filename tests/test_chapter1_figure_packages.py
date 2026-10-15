import json
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "submission" / "chapter1_current"
MAIN = PACKAGE / "figures"
SUPP = PACKAGE / "supplement" / "figures"


def _assert_svg(path: Path) -> None:
    assert path.is_file(), path
    assert path.stat().st_size > 1000, path
    assert ET.parse(path).getroot().tag.endswith("svg")


def _assert_pdf(path: Path) -> None:
    assert path.is_file(), path
    assert path.stat().st_size > 1000, path
    assert path.read_bytes().startswith(b"%PDF"), path


def test_main_figure_package_is_complete() -> None:
    manifest = json.loads((MAIN / "FIGURE_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["contract"] == "chapter1_main_figures_corrected_v1"
    assert (
        manifest["source_surface"]
        == "corrected_geography_20260924_plus_final_traitwise_H1_20261004"
    )
    assert [item["figure"] for item in manifest["figures"]] == [1, 2, 3, 4, 5, 6]

    for item in manifest["figures"]:
        svg = MAIN / item["svg"]
        pdf = MAIN / item["pdf"]
        _assert_svg(svg)
        _assert_pdf(pdf)
        assert svg.stat().st_size == item["svg_bytes"]
        assert pdf.stat().st_size == item["pdf_bytes"]

    figure3 = (MAIN / "Figure3_constraint_response_triangle.svg").read_text(
        encoding="utf-8"
    )
    assert "Historical causal edge not identified" in figure3
    assert "not a mediation model" in figure3

    assert (
        (MAIN / "Figure4_H1_recurrent_multivariate_response.svg").read_bytes()
        == (ROOT / "results/h1_final_traitwise_t_20261004/traitwise_H1.svg").read_bytes()
    )
    assert (
        (MAIN / "Figure4_H1_recurrent_multivariate_response.pdf").read_bytes()
        == (ROOT / "results/h1_final_traitwise_t_20261004/traitwise_H1.pdf").read_bytes()
    )


def test_supplement_figure_package_is_complete() -> None:
    manifest = json.loads((SUPP / "FIGURE_MANIFEST.json").read_text(encoding="utf-8"))
    assert manifest["contract"] == "chapter1_supplement_figures_v1"
    assert manifest["source_surface"] == "corrected_geography_20260924"
    assert manifest["n_figures"] == 7

    names = [item["file"] for item in manifest["files"]]
    for idx in range(1, 8):
        svg = next(SUPP.glob(f"Figure_S{idx}_*.svg"))
        pdf = next(SUPP.glob(f"Figure_S{idx}_*.pdf"))
        _assert_svg(svg)
        _assert_pdf(pdf)
        assert svg.name in names
        assert pdf.name in names


def test_current_q1_five_main_figures_and_source_six_replay() -> None:
    current = MAIN / "q1_main5"
    current_manifest = json.loads((current / "FIGURE_MANIFEST.json").read_text(encoding="utf-8"))
    assert current_manifest["contract"] == "chapter1_q1_five_main_figures_v1"
    assert [item["figure"] for item in current_manifest["figures"]] == [1, 2, 3, 4, 5]
    for item in current_manifest["figures"]:
        _assert_svg(current / item["svg"])
        _assert_pdf(current / item["pdf"])
        assert (current / item["svg"]).stat().st_size == item["svg_bytes"]
        assert (current / item["pdf"]).stat().st_size == item["pdf_bytes"]
    # The print-facing H1 arrangement is newly faceted from the same frozen
    # 112 source fits; the historical six-source artwork remains byte-exact.
    q1_figure = (current / "Q1_Fig3_seven_atomic_trait_responses.svg").read_text(encoding="utf-8")
    assert "Seven traits in four geographic regions" in q1_figure
    assert "Self-compatibility" in q1_figure
    print_audit = json.loads((current / "PRINT_WIDTH_AUDIT.json").read_text(encoding="utf-8"))
    assert print_audit["contract"] == "chapter1_q1_print_width_text_audit_v1"
    assert print_audit["all_pass"] is True
    assert all(item["min_effective_text_pt"] >= 6.5 for item in print_audit["figures"])
    captions = (PACKAGE / "FIGURE_CAPTIONS.md").read_text(encoding="utf-8")
    assert captions.count("## Figure ") == 5
    assert (PACKAGE / "history/FIGURE_CAPTIONS_SOURCE6_20261015.md").is_file()
    assert len(json.loads((MAIN / "FIGURE_MANIFEST.json").read_text())["figures"]) == 6
    proof_dir = current / "qa_print"
    proof_manifest = json.loads((proof_dir / "PRINT_PROOF_MANIFEST.json").read_text())
    assert proof_manifest["contract"] == "chapter1_q1_a4_print_width_proof_v1"
    assert proof_manifest["figure_width_mm"] == 178
    assert len(proof_manifest["figures"]) == 5
    assert all(item["fits_a4_at_print_width"] for item in proof_manifest["figures"])
    _assert_pdf(proof_dir / proof_manifest["proof_pdf"])
    from pypdf import PdfReader
    assert len(PdfReader(proof_dir / proof_manifest["proof_pdf"]).pages) == 5
