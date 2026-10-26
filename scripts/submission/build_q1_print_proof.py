"""Compose five Q1 vector figures at a true 178-mm width on A4 proof pages.

This is an editorial print-width layout proof, NOT a final journal submission
file and NOT a numerical analysis. It does not modify any source figure PDF.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation

PT_PER_MM = 72 / 25.4
A4_WIDTH_PT = 210 * PT_PER_MM
A4_HEIGHT_PT = 297 * PT_PER_MM
FIGURE_WIDTH_MM = 178
FIGURE_WIDTH_PT = FIGURE_WIDTH_MM * PT_PER_MM
SAFE_TOP_BOTTOM_MM = 12


def build(figure_dir: Path, output_dir: Path) -> dict[str, object]:
    index = json.loads((figure_dir / "FIGURE_MANIFEST.json").read_text())
    if index["contract"] != "chapter1_q1_five_main_figures_v1":
        raise ValueError("Unexpected current Q1 figure manifest")
    rows = sorted(index["figures"], key=lambda item: item["figure"])
    if [x["figure"] for x in rows] != [1, 2, 3, 4, 5]:
        raise ValueError("Require five numbered figures")
    writer = PdfWriter()
    metrics = []
    for item in rows:
        src = PdfReader(figure_dir / item["pdf"])
        if len(src.pages) != 1:
            raise ValueError(f"Figure {item['figure']} is not single-page")
        figure = src.pages[0]
        w, h = float(figure.mediabox.width), float(figure.mediabox.height)
        scale = FIGURE_WIDTH_PT / w
        new_w, new_h = w * scale, h * scale
        if new_h > A4_HEIGHT_PT - 2 * SAFE_TOP_BOTTOM_MM * PT_PER_MM:
            raise ValueError(f"Figure {item['figure']} is too tall at true 178-mm width")
        page = writer.add_blank_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
        x, y = (A4_WIDTH_PT - new_w) / 2, (A4_HEIGHT_PT - new_h) / 2
        page.merge_transformed_page(
            figure,
            Transformation().scale(scale).translate(x, y),
        )
        metrics.append({
            "figure": item["figure"],
            "source_pdf": item["pdf"],
            "native_pdf_width_pt": round(w, 3),
            "native_pdf_height_pt": round(h, 3),
            "print_width_mm": FIGURE_WIDTH_MM,
            "print_height_mm": round(new_h / PT_PER_MM, 2),
            "a4_left_margin_mm": round(x / PT_PER_MM, 2),
            "a4_bottom_margin_mm": round(y / PT_PER_MM, 2),
            "fits_a4_at_print_width": True,
        })
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf = output_dir / "Q1_Main5_178mm_A4_print_proof.pdf"
    with pdf.open("wb") as out:
        writer.write(out)
    verify = PdfReader(pdf)
    if len(verify.pages) != 5:
        raise RuntimeError("Print proof did not contain five pages")
    result = {
        "contract": "chapter1_q1_a4_print_width_proof_v1",
        "status": "visual_layout_proof_not_final_journal_composition",
        "a4_size_mm": [210, 297],
        "figure_width_mm": FIGURE_WIDTH_MM,
        "proof_pdf": pdf.name,
        "figures": metrics,
    }
    (output_dir / "PRINT_PROOF_MANIFEST.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = build(args.figure_dir, args.output_dir)
    for row in result["figures"]:
        print(f"Figure {row['figure']}: 178 mm wide, {row['print_height_mm']} mm tall, A4 fit")


if __name__ == "__main__":
    main()
