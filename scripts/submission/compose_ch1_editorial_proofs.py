"""Compose non-submittable Chapter 1 Main/SI review proofs from frozen files.

Reads only the current manuscript, SI, captions and pre-existing figure PDFs.
Neither scientific numbers nor the submission baseline are recomputed here.
Rights clearance, author approval and final journal composition are *not* implied.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from pypdf import PdfReader, PdfWriter, Transformation
from weasyprint import CSS, HTML

A4_WIDTH_PT = 210 * 72 / 25.4
A4_HEIGHT_PT = 297 * 72 / 25.4
FIGURE_WIDTH_PT = 178 * 72 / 25.4
SOURCE_ROOT = Path("submission/chapter1_current")
PROOF_BANNER = (
    '<div class="draft-banner"><strong>EDITORIAL PROOF - NOT FOR SUBMISSION</strong>'
    "<span>Authorship, rights clearance, DOI and journal production are pending.</span></div>\n\n"
)
CSS_BODY = """
@page { size: A4; margin: 18mm 19mm 18mm 19mm;
  @bottom-center { content: "Q1 - editorial proof | Page " counter(page);
    font-family: DejaVu Sans, sans-serif; font-size: 7.5pt; color: #6b7280; } }
html, body {font-family: "DejaVu Serif", "Liberation Serif", serif; color: #141a20;}
body {font-size: 9.3pt; line-height: 1.36; overflow-wrap: anywhere;}
p {margin: 0 0 0.53em 0; text-align: left;}
h1, h2, h3, h4 {font-family: "DejaVu Sans", sans-serif;
  break-after: avoid; page-break-after: avoid; color: #13293b;}
h1 {font-size: 15pt; margin: 0.35em 0 0.7em;}
h2 {font-size: 11.7pt; margin: 1.0em 0 0.45em;}
h3 {font-size: 10pt; margin: 0.82em 0 0.35em;}
h4 {font-size: 9.4pt; margin: 0.65em 0 0.28em;}
ol, ul {padding-left: 1.45em; margin: 0.25em 0 0.6em;}
li {margin-bottom: 0.15em;}
table {width: 100%; border-collapse: collapse; table-layout: fixed; margin: 0.65em 0 0.85em; break-inside: avoid; page-break-inside: avoid;}
thead {display: table-header-group;}
td, th {border: 0.55px solid #b8c2ca; padding: 0.26em 0.33em;
  vertical-align: top; overflow-wrap: anywhere; font-size: 8pt;}
th {background: #f0f3f6; font-weight: bold;}
tr {break-inside: avoid;}
blockquote {padding: 0 0 0 0.8em; border-left: 2px solid #9eabb4; margin: 0.3em 0 0.65em;}
a {color: #25527a; text-decoration: none;}
code {font-family: "DejaVu Sans Mono", monospace; font-size: 0.88em; white-space: normal;}
pre {white-space: pre-wrap;}
hr {border: 0; border-top: 1px solid #adb9c1; margin: 1em 0;}
.draft-banner {font: 8pt "DejaVu Sans", sans-serif; border: 1px solid #9d521b;
  background: #fff8ea; padding: 0.5em; margin-bottom: 1em;}
.draft-banner strong {display: block;}
.draft-banner span {display: block;}
.page-break {break-before: page; page-break-before: always;}
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(args: list[str]) -> None:
    subprocess.run(args, check=True, text=True)


def _strip_caps_heading(value: str) -> str:
    """Retain the five caption entries, discard editor-only preamble and history link."""
    lines = value.splitlines()
    try:
        first = next(i for i, line in enumerate(lines) if line.startswith("## Figure 1 "))
    except StopIteration as exc:
        raise ValueError("Five caption entries are missing") from exc
    lines = lines[first:]
    while lines and not lines[-1].strip():
        lines.pop()
    if lines and lines[-1].startswith("Original source-six caption provenance:"):
        lines.pop()
    result = "\n".join(lines).rstrip() + "\n"
    if len(re.findall(r"^## Figure [1-5] \|", result, re.MULTILINE)) != 5:
        raise ValueError("Exactly five current Main captions must be available")
    if re.search(r"^## Figure 6 \|", result, re.MULTILINE):
        raise ValueError("Superseded Figure 6 cannot enter the Main proof")
    return result


def compose_markdown(root: Path) -> tuple[str, str, dict[str, object]]:
    src = root / SOURCE_ROOT
    main = (src / "MANUSCRIPT.md").read_text(encoding="utf-8")
    si = (src / "SUPPLEMENTARY_INFORMATION.md").read_text(encoding="utf-8")
    caps = _strip_caps_heading((src / "FIGURE_CAPTIONS.md").read_text(encoding="utf-8"))
    titlepage = (src / "TITLE_PAGE_TEMPLATE.md").read_text(encoding="utf-8")
    if "# Do pollen-supply constraints underlie floral island syndromes?" not in main:
        raise ValueError("Manuscript title mismatch")
    if "## References" not in main or "## Introduction" not in main:
        raise ValueError("Manuscript section layout changed")
    if "Appendix S7" not in si:
        raise ValueError("SI is missing its final appendix")
    if re.search(r"\bFigure 6\b", main):
        raise ValueError("Main contains superseded Figure 6 reference")
    start = main.index("## Introduction")
    stop = main.index("## References", start)
    word_count = len(main[start:stop].strip().split())
    if f"**Main-text word count:** {word_count:,} " not in titlepage:
        raise ValueError("Title-page word count does not match frozen Main heading-inclusive count")
    if not 4500 < word_count <= 5000:
        raise ValueError(f"Unexpected scientific Main count: {word_count}")
    # Put a physical page break before science and before legends; the source files stay unchanged.
    main_markdown = (
        PROOF_BANNER + titlepage.strip() +
        '\n\n<div class="page-break"></div>\n\n' +
        main.strip() +
        '\n\n<div class="page-break"></div>\n\n' +
        "# Figure legends\n\n" + caps
    )
    si_markdown = PROOF_BANNER + si.strip() + "\n"
    return main_markdown, si_markdown, {
        "original_main_body_words": word_count,
        "word_count_method": "whitespace_split_from_Introduction_heading_to_References_heading",
        "original_main_text_chars": len(main),
        "original_si_chars": len(si),
        "main_figure_captions": 5,
        "historical_causal_mediation_claimed": False,
    }


def _pandoc(md: Path, dst: Path, fmt: str) -> None:
    _run(["pandoc", "--standalone", "--from=gfm", f"--to={fmt}", str(md), "--output", str(dst)])


def _render_one(markdown: str, dest: Path, name: str, style: str) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="ch1-compose-") as tmp:
        tmpdir = Path(tmp)
        input_md = tmpdir / f"{name}.md"
        html_file = tmpdir / f"{name}.html"
        input_md.write_text(markdown, encoding="utf-8")
        docx = dest / f"{name}.docx"
        pdf = dest / f"{name}.pdf"
        _pandoc(input_md, docx, "docx")
        _pandoc(input_md, html_file, "html5")
        HTML(filename=str(html_file), base_url=str(html_file.parent)).write_pdf(
            pdf, stylesheets=[CSS(string=style)]
        )
        page_count = len(PdfReader(pdf).pages)
        if page_count < 2:
            raise ValueError(f"Unexpectedly short PDF {pdf}: {page_count} pages")
        if not docx.read_bytes().startswith(b"PK"):
            raise ValueError("Pandoc did not create a valid DOCX package")
        if not pdf.read_bytes().startswith(b"%PDF"):
            raise ValueError("WeasyPrint did not create PDF")
        return {
            "docx": docx.name, "docx_sha256": sha256(docx),
            "pdf": pdf.name, "pdf_sha256": sha256(pdf),
            "pdf_pages": page_count,
        }


def _append_pdf(writer: PdfWriter, path: Path) -> int:
    reader = PdfReader(path)
    for page in reader.pages:
        writer.add_page(page)
    return len(reader.pages)


def _supplement_figure_pdf(root: Path, dest: Path) -> int:
    source = root / SOURCE_ROOT / "supplement/figures"
    manifest = json.loads((source / "FIGURE_MANIFEST.json").read_text())
    if manifest["contract"] != "chapter1_supplement_figures_v1":
        raise ValueError("Unexpected supplement figure contract")
    files = {f["file"]: f["bytes"] for f in manifest["files"]}
    writer = PdfWriter()
    for i in range(1, 8):
        candidates = [name for name in files if name.startswith(f"Figure_S{i}_") and name.endswith(".pdf")]
        if len(candidates) != 1:
            raise ValueError(f"Supplementary figure S{i} ambiguous")
        source_pdf = source / candidates[0]
        if source_pdf.stat().st_size != files[candidates[0]]:
            raise ValueError(f"Supplementary figure {i} size changed")
        reader = PdfReader(source_pdf)
        if len(reader.pages) != 1:
            raise ValueError(f"Supplementary figure {i} must be one page")
        old = reader.pages[0]
        width, height = float(old.mediabox.width), float(old.mediabox.height)
        scale = min(FIGURE_WIDTH_PT / width, (A4_HEIGHT_PT - 42 * 72 / 25.4) / height)
        display_w, display_h = scale * width, scale * height
        page = writer.add_blank_page(width=A4_WIDTH_PT, height=A4_HEIGHT_PT)
        page.merge_transformed_page(
            old,
            Transformation().scale(scale).translate(
                (A4_WIDTH_PT - display_w) / 2,
                (A4_HEIGHT_PT - display_h) / 2,
            ),
        )
    with dest.open("wb") as target:
        writer.write(target)
    if len(PdfReader(dest).pages) != 7:
        raise ValueError("Supplementary figures proof page count mismatch")
    return 7


def _combined(target: Path, components: list[Path]) -> int:
    writer = PdfWriter()
    for path in components:
        _append_pdf(writer, path)
    with target.open("wb") as out:
        writer.write(out)
    return len(PdfReader(target).pages)


def _text_audit(path: Path, must_contain: list[str]) -> None:
    content = " ".join(
        page.extract_text() or "" for page in PdfReader(path).pages
    )
    squashed = re.sub(r"\s+", " ", content)
    for token in must_contain:
        if token not in squashed:
            raise ValueError(f"Proof PDF {path.name} missing text token: {token}")


def build(root: Path, out: Path) -> dict[str, object]:
    if shutil.which("pandoc") is None:
        raise RuntimeError("pandoc is required to render the source Markdown")
    config = json.loads((root / "config/chapter1_submission_current.json").read_text())
    if config["status"] != "primary_submission_baseline":
        raise ValueError("Submission baseline is not frozen")
    out.mkdir(parents=True, exist_ok=True)
    main_md, si_md, audit = compose_markdown(root)
    main = _render_one(main_md, out, "Q1_Main_EDIT_PROOF", CSS_BODY)
    si_style = CSS_BODY + "body { font-size: 8.8pt; } table td, table th { font-size: 7.65pt; }"
    si = _render_one(si_md, out, "Q1_SI_EDIT_PROOF", si_style)
    _text_audit(out / main["pdf"], [
        "pollen", "4,379", "106,295", "2,969", "References", "Figure legends"
    ])
    _text_audit(out / si["pdf"], [
        "Supplementary Information", "Appendix S1", "Appendix S7", "222,688"
    ])
    main_figures = root / SOURCE_ROOT / "figures/q1_main5/qa_print/Q1_Main5_178mm_A4_print_proof.pdf"
    if len(PdfReader(main_figures).pages) != 5:
        raise ValueError("Main figure proof does not contain five pages")
    main_proof = out / "Q1_Main_plus_5_figures_REVIEW.pdf"
    main_proof_pages = _combined(main_proof, [out / main["pdf"], main_figures])
    supplement_pdf = out / "Q1_Seven_SI_figures_A4.pdf"
    supp_count = _supplement_figure_pdf(root, supplement_pdf)
    si_proof = out / "Q1_SI_plus_7_figures_REVIEW.pdf"
    si_proof_pages = _combined(si_proof, [out / si["pdf"], supplement_pdf])
    result = {
        "contract": "chapter1_q1_main_si_editorial_composition_v1",
        "status": "EDITORIAL_PROOF_NOT_SUBMISSION_READY",
        "scope": "layout_proof_only_no_numerical_reanalysis",
        "main": main,
        "si": si,
        "main_plus_figures": {
            "file": main_proof.name, "pages": main_proof_pages,
            "sha256": sha256(main_proof), "figure_pages": 5
        },
        "si_plus_figures": {
            "file": si_proof.name, "pages": si_proof_pages,
            "sha256": sha256(si_proof), "figure_pages": supp_count
        },
        "supplement_figures_only": supplement_pdf.name,
        "source_sha256": {
            name: sha256(root / SOURCE_ROOT / name)
            for name in [
                "MANUSCRIPT.md", "SUPPLEMENTARY_INFORMATION.md",
                "FIGURE_CAPTIONS.md", "TITLE_PAGE_TEMPLATE.md"
            ]
        },
        "author_approval": False,
        "trait_redistribution_rights_cleared": False,
        "final_release_doi_minted": False,
        "submission_authorized": False,
        "source_audit": audit,
    }
    (out / "EDITORIAL_PROOF_MANIFEST.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path,
                        default=Path("submission/chapter1_current/editorial_proofs"))
    args = parser.parse_args()
    result = build(args.root.resolve(), args.out.resolve())
    print(json.dumps({
        "status": result["status"],
        "Main PDF pages": result["main"]["pdf_pages"],
        "SI PDF pages": result["si"]["pdf_pages"],
        "Main review pages": result["main_plus_figures"]["pages"],
        "SI review pages": result["si_plus_figures"]["pages"],
    }, indent=2))


if __name__ == "__main__":
    main()
