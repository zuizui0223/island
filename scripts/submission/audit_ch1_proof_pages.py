"""Inspect every page of the Q1 Main/SI text PDFs and produce legible contact tiles.

This is a *proof-quality* audit: it catches missing/empty/out-of-bounds pages
and provides small inspection PNGs. It does not establish author approval,
third-party redistribution rights, or journal acceptance.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import fitz
from PIL import Image, ImageDraw


def run(docpath: Path, outdir: Path, role: str) -> dict[str, object]:
    outdir.mkdir(parents=True, exist_ok=True)
    pdf = fitz.open(docpath)
    page_data: list[dict[str, object]] = []
    sheet_w, sheet_h = 330, 456
    gap, title_h, cols, per_sheet = 14, 24, 3, 6
    tiles: list[str] = []
    for start in range(0, len(pdf), per_sheet):
        slots = min(per_sheet, len(pdf) - start)
        rows = math.ceil(slots / cols)
        sheet = Image.new(
            "RGB",
            (cols * (sheet_w + gap) + gap, rows * (sheet_h + title_h + gap) + gap),
            "white",
        )
        pen = ImageDraw.Draw(sheet)
        for offset in range(slots):
            idx = start + offset
            page = pdf[idx]
            text = page.get_text()
            chars = len(text.strip())
            blocks = page.get_text("blocks")
            outside = []
            for x0, y0, x1, y1, *_ in blocks:
                if x0 < -0.5 or y0 < -0.5 or x1 > page.rect.width + 0.5 or y1 > page.rect.height + 0.5:
                    outside.append([round(float(x0), 2), round(float(y0), 2),
                                    round(float(x1), 2), round(float(y1), 2)])
            page_data.append({
                "page": idx + 1,
                "text_characters": chars,
                "text_blocks": len(blocks),
                "out_of_page_boxes": outside,
                "page_width_pt": round(float(page.rect.width), 3),
                "page_height_pt": round(float(page.rect.height), 3),
                "passes": chars >= 35 and not outside,
            })
            pix = page.get_pixmap(matrix=fitz.Matrix(0.53, 0.53), alpha=False)
            im = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
            im.thumbnail((sheet_w, sheet_h))
            x = gap + (offset % cols) * (sheet_w + gap)
            y = gap + (offset // cols) * (sheet_h + title_h + gap)
            sheet.paste(im, (x, y + title_h))
            pen.text((x, y + 4), f"{role.upper()} {idx + 1} / {len(pdf)}", fill="#142737")
        name = f"Q1_{role.upper()}_PAGES_{start + 1:02d}_{start + slots:02d}.png"
        output = outdir / name
        sheet.save(output, optimize=True)
        tiles.append(name)
    info = {
        "contract": "chapter1_q1_editorial_all_pages_audit_v1",
        "role": role,
        "source_pdf": docpath.name,
        "total_pages": len(pdf),
        "every_page_nonempty_and_on_canvas": all(p["passes"] for p in page_data),
        "tiles": tiles,
        "pages": page_data,
    }
    path = outdir / f"Q1_{role.upper()}_PAGES_AUDIT.json"
    path.write_text(json.dumps(info, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if not info["every_page_nonempty_and_on_canvas"]:
        bad = [p["page"] for p in page_data if not p["passes"]]
        raise ValueError(f"{role} has empty or out-of-bound pages: {bad}")
    return info


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proof-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.proof_dir
    counts = {}
    for role in ("main", "si"):
        result = run(root / f"Q1_{'Main' if role == 'main' else 'SI'}_EDIT_PROOF.pdf", root, role)
        counts[role] = len(result["pages"])
    print(f"Checked every page (Main {counts['main']}; SI {counts['si']})")


if __name__ == "__main__":
    main()
