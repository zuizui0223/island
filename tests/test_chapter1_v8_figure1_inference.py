import json
from pathlib import Path

from PIL import Image
import pytest

from island_v2.chapter1_v8_figure1_inference import (
    ANCHORS,
    FigureInputError,
    render_figure,
    validate_spec,
)


def test_validate_canonical_figure1_spec() -> None:
    text = validate_spec(Path("docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md"))
    assert "family→genus attenuation visually central" in text
    assert "Do not show a pollinator icon as the cause" in text


def test_spec_validation_fails_closed_on_missing_anchor(tmp_path: Path) -> None:
    path = tmp_path / "bad_spec.md"
    path.write_text("4/4 -> 4/4 -> 0/4\n", encoding="utf-8")
    with pytest.raises(FigureInputError):
        validate_spec(path)


def test_render_v8_figure1(tmp_path: Path) -> None:
    output = tmp_path / "output"
    manifest = render_figure(
        spec_path=Path("docs/chapter1_figure1_hierarchical_syndrome_spec_20260914.md"),
        output_dir=output,
    )
    assert manifest["contract"] == "chapter1_v8_figure1_hierarchical_inference_v1"
    assert manifest["new_biological_models_fitted"] is False
    assert manifest["new_p_values_generated"] is False
    assert manifest["numerical_anchors"] == ANCHORS

    png = output / "chapter1_v8_figure1_hierarchical_inference.png"
    svg = output / "chapter1_v8_figure1_hierarchical_inference.svg"
    pdf = output / "chapter1_v8_figure1_hierarchical_inference.pdf"
    manifest_path = output / "chapter1_v8_figure1_manifest.json"
    for path in (png, svg, pdf, manifest_path):
        assert path.is_file()
        assert path.stat().st_size > 0

    width, height = Image.open(png).size
    assert width >= 3000
    assert height >= 2000
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert payload["dashed_arrow_meaning"].startswith("mechanistic or cross-scale")
    assert "pollinator-specific causal mechanism" in payload["claim_boundary"]
