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


def test_main_caption_count_matches_six_figure_manifest() -> None:
    captions = (PACKAGE / "FIGURE_CAPTIONS.md").read_text(encoding="utf-8")
    assert captions.count("## Figure ") == 6
