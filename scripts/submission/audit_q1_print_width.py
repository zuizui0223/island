"""Check Q1 PDF companion SVG text size after resizing to 178-mm figure width.

Matplotlib embeds text as glyph paths with a translate(...) scale(s -s)
wrapper. Source text size (pt) = 100*s for the fixed font path metrics.
This calculation is only a geometry/text audit; it does not assert
scientific validity or substitute for full-page visual proof review.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

TARGET_WIDTH_PT = 178 * 72 / 25.4
TEXT_RE = re.compile(r'transform="translate\([^)]+\) scale\(([0-9.]+) -\1\)"')
WIDTH_RE = re.compile(r'<svg[^>]+width="([0-9.]+)pt"')


def audit_file(path: Path, threshold: float) -> dict[str, object]:
    body = path.read_text(encoding="utf-8")
    width = WIDTH_RE.search(body)
    if not width:
        raise ValueError(f"Missing SVG width in points: {path}")
    native_width = float(width.group(1))
    sizes_pt = [100 * float(m.group(1)) for m in TEXT_RE.finditer(body)]
    if not sizes_pt:
        raise ValueError(f"No Matplotlib text-glyph transforms found: {path}")
    minimum_native_pt = min(sizes_pt)
    effective_min_pt = minimum_native_pt * TARGET_WIDTH_PT / native_width
    return {
        "filename": path.name,
        "native_width_pt": native_width,
        "print_width_mm": 178,
        "text_label_count": len(sizes_pt),
        "min_native_text_pt": round(minimum_native_pt, 3),
        "min_effective_text_pt": round(effective_min_pt, 3),
        "threshold_pt": threshold,
        "pass": effective_min_pt >= threshold,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--figure-dir", type=Path, required=True)
    parser.add_argument("--min-pt", type=float, default=6.5)
    args = parser.parse_args()
    figures = sorted(args.figure_dir.glob("Q1_Fig*.svg"))
    if len(figures) != 5:
        raise SystemExit(f"Expected five Q1 Main SVGs, found {len(figures)}")
    results = [audit_file(path, args.min_pt) for path in figures]
    report = {
        "contract": "chapter1_q1_print_width_text_audit_v1",
        "target_width_mm": 178,
        "scaling_assumption": "uniform scaling to 178 mm; Matplotlib text path sizes",
        "minimum_effective_text_pt": min(x["min_effective_text_pt"] for x in results),
        "all_pass": all(x["pass"] for x in results),
        "figures": results,
    }
    out = args.figure_dir / "PRINT_WIDTH_AUDIT.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for row in results:
        print(f'{row["filename"]}: min={row["min_effective_text_pt"]:.2f} pt @178 mm; pass={row["pass"]}')
    if not report["all_pass"]:
        raise SystemExit("Some Q1 figures contain text below print-width font threshold")


if __name__ == "__main__":
    main()
