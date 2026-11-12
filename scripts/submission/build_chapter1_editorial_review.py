"""Build a RIGHTS-SAFE, non-submittable Chapter 1 editorial review packet.

No trait-ledger/raw-data inputs are opened. The packet is an author-review aid;
journal submission, licensing clearance and software release are separate gates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import zipfile
from pathlib import Path

import fitz
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor

REVIEW_WARNING = (
    "EDITORIAL REVIEW DRAFT — NOT FOR JOURNAL SUBMISSION. "
    "DATA-RIGHTS AND AUTHOR APPROVAL GATES REMAIN OPEN."
)
DOC_SOURCE = "submission/chapter1_current"
DOC_FILES = (
    "MANUSCRIPT.md",
    "SUPPLEMENTARY_INFORMATION.md",
    "TITLE_PAGE_TEMPLATE.md",
    "COVER_LETTER_DRAFT.md",
    "FIGURE_CAPTIONS.md",
    "ECOLOGY_LETTERS_DATA_GATE.md",
    "DATA_ACCESSIBILITY_DRAFT.md",
    "SUBMISSION_CHECKLIST.md",
)
TABLE_FILES = (
    "Table_S1_data_summary.csv",
    "Table_S2_H1_traitwise.csv",
    "Table_S3a_H2_all_finite_cluster.csv",
    "Table_S3b_H2_direct_finite_cluster.csv",
    "Table_S5_H3_pollen_limitation.csv",
    "Table_S6a_H4_scores.csv",
)
MANUSCRIPT_MAX = 5000
ABSTRACT_MAX = 150

def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def _read_utf8(path: Path) -> str:
    return path.read_text(encoding="utf-8")

def _md_words(text: str) -> int:
    return len(text.strip().split())

def _word_counts(text: str) -> dict[str, int]:
    if "## Abstract" not in text or "**Keywords:**" not in text:
        raise ValueError("Manuscript abstract markers are missing")
    if "## Introduction" not in text or "## References" not in text:
        raise ValueError("Manuscript body markers are missing")
    abstract = text.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0]
    main = text.split("## Introduction", 1)[1].split("## References", 1)[0]
    return {"abstract": _md_words(abstract), "main": _md_words("## Introduction" + main)}

def _set_cell_border(cell) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    borders = OxmlElement("w:tcBorders")
    for edge in ("top", "left", "bottom", "right"):
        el = OxmlElement("w:" + edge)
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), "2")
        el.set(qn("w:color"), "CBD5E1")
        borders.append(el)
    tc_pr.append(borders)

def _add_inline(paragraph, text: str) -> None:
    # Keep content verbatim and style basic Markdown; protect substantive text.
    tick = chr(96)
    token = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*|" + tick +
                       r"[^" + tick + r"]+" + tick +
                       r"|\[[^\]]+\]\([^)]+\))")
    pos = 0
    for m in token.finditer(text):
        if m.start() > pos:
            paragraph.add_run(text[pos:m.start()])
        x = m.group(0)
        if x.startswith("**") and x.endswith("**"):
            run = paragraph.add_run(x[2:-2]); run.bold = True
        elif x.startswith("*") and x.endswith("*"):
            run = paragraph.add_run(x[1:-1]); run.italic = True
        elif x.startswith(chr(96)) and x.endswith(chr(96)):
            run = paragraph.add_run(x[1:-1])
            run.font.name = "Consolas"; run.font.size = Pt(8.3)
        elif x.startswith("[") and "](" in x:
            label, link = x[1:-1].split("](", 1)
            run = paragraph.add_run(label + " (" + link + ")")
            run.font.color.rgb = RGBColor(30, 80, 140)
        else:
            paragraph.add_run(x)
        pos = m.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])

def _add_table(doc: Document, lines: list[str]) -> None:
    rows: list[list[str]] = []
    for line in lines:
        cells = [v.strip() for v in line.strip().strip("|").split("|")]
        if all(re.fullmatch(r":?-{3,}:?", c or "") for c in cells):
            continue
        rows.append(cells)
    if not rows:
        return
    width = max(len(row) for row in rows)
    if width > 10:
        para = doc.add_paragraph(" | ".join(rows[0]))
        para.style = "Normal"
        for row in rows[1:]:
            doc.add_paragraph(" | ".join(row))
        return
    table = doc.add_table(rows=len(rows), cols=width)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    for i, line in enumerate(rows):
        for j, value in enumerate(line):
            cell = table.cell(i, j)
            cell.text = value
            _set_cell_border(cell)
            for p in cell.paragraphs:
                for r in p.runs:
                    r.font.size = Pt(7.5)
                    if i == 0:
                        r.bold = True
    doc.add_paragraph("")

def _add_markdown(doc: Document, text: str) -> None:
    # Combine wrapped Markdown lines into paragraphs while preserving all content.
    pending: list[str] = []
    table_lines: list[str] = []
    code_fence = False

    def flush() -> None:
        if pending:
            p = doc.add_paragraph(style="Normal")
            _add_inline(p, " ".join(pending))
            pending.clear()
    def flush_table() -> None:
        if table_lines:
            _add_table(doc, table_lines)
            table_lines.clear()

    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith(chr(96) * 3):
            flush(); flush_table()
            code_fence = not code_fence
            continue
        if code_fence:
            p = doc.add_paragraph()
            p.style = "Normal"
            r = p.add_run(raw)
            r.font.name = "Consolas"
            r.font.size = Pt(8)
            continue
        if not line:
            flush(); flush_table()
            continue
        if line.startswith("|") and line.endswith("|"):
            flush(); table_lines.append(line)
            continue
        flush_table()
        head = re.match(r"^(#{1,6})\s+(.*)$", line)
        if head:
            flush()
            level = min(len(head.group(1)), 3)
            p = doc.add_heading(level=level)
            _add_inline(p, head.group(2))
            continue
        if line in ("---", "***"):
            flush()
            doc.add_paragraph("")
            continue
        list_match = re.match(r"^([-*]|\d+\.)\s+(.*)$", line)
        if list_match:
            flush()
            p = doc.add_paragraph(style="List Bullet" if list_match.group(1) in ("-", "*") else "List Number")
            _add_inline(p, list_match.group(2))
            continue
        pending.append(line)
    flush(); flush_table()

def _setup_document(kind: str, source: str) -> Document:
    doc = Document()
    sec = doc.sections[0]
    sec.top_margin = Inches(0.83)
    sec.bottom_margin = Inches(0.77)
    sec.left_margin = Inches(0.89)
    sec.right_margin = Inches(0.83)
    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Cambria"
    normal.font.size = Pt(10)
    normal.paragraph_format.space_after = Pt(5)
    normal.paragraph_format.line_spacing = 1.14
    for key, size, colour in (
        ("Title", 16, RGBColor(18, 48, 83)),
        ("Heading 1", 13, RGBColor(18, 48, 83)),
        ("Heading 2", 11, RGBColor(18, 48, 83)),
        ("Heading 3", 10, RGBColor(18, 48, 83)),
    ):
        style = styles[key]; style.font.name = "Cambria"
        style.font.size = Pt(size); style.font.bold = True
        style.font.color.rgb = colour
    header = sec.header.paragraphs[0]
    header.text = "CHAPTER 1 / " + kind.upper() + " / AUTHOR REVIEW ONLY"
    header.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    for run in header.runs:
        run.font.size = Pt(7)
        run.font.color.rgb = RGBColor(130, 48, 48)
    footer = sec.footer.paragraphs[0]
    footer.text = "Source: " + source + "  •  NOT FOR SUBMISSION"
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    for run in footer.runs:
        run.font.size = Pt(7)
    doc.add_heading("EDITORIAL REVIEW COPY", 0)
    doc.add_paragraph(REVIEW_WARNING)
    doc.add_paragraph("The source manuscript is not changed by this conversion. Author "
                      "details, data-sharing authorization and permanent software DOI "
                      "remain pending.")
    doc.add_page_break()
    return doc

def _figure_pdf_to_png(pdf: Path, out_png: Path) -> None:
    with fitz.open(pdf) as source:
        if len(source) < 1:
            raise ValueError("Figure PDF contains zero pages: " + str(pdf))
        pix = source[0].get_pixmap(matrix=fitz.Matrix(1.4, 1.4), alpha=False)
        pix.save(out_png)
    if out_png.stat().st_size < 1500:
        raise ValueError("Rasterized figure too small: " + str(out_png))

def _append_figures(doc: Document, images: list[tuple[str, Path]]) -> None:
    if not images:
        return
    doc.add_page_break()
    doc.add_heading("Figures — review copy", level=1)
    for i, (name, path) in enumerate(images):
        if i:
            doc.add_page_break()
        doc.add_heading(name, level=2)
        doc.add_picture(str(path), width=Inches(6.65))

def _build_doc(markdown: Path, out: Path, label: str,
               images: list[tuple[str, Path]] | None = None,
               captions: Path | None = None) -> None:
    doc = _setup_document(label, markdown.name)
    _add_markdown(doc, _read_utf8(markdown))
    if captions:
        doc.add_page_break()
        _add_markdown(doc, _read_utf8(captions))
    _append_figures(doc, images or [])
    out.parent.mkdir(parents=True, exist_ok=True)
    doc.save(out)

def build(root: Path, out_dir: Path, *, source_commit: str) -> dict:
    root = root.resolve()
    out_dir = out_dir.resolve()
    package = root / DOC_SOURCE
    selector = json.loads(_read_utf8(root / "config/chapter1_submission_current.json"))
    if selector.get("status") != "primary_submission_baseline":
        raise ValueError("Only the current fixed submission baseline may be compiled")
    if selector.get("manuscript") != DOC_SOURCE + "/MANUSCRIPT.md":
        raise ValueError("Unexpected manuscript selector")
    manuscript = _read_utf8(package / "MANUSCRIPT.md")
    counts = _word_counts(manuscript)
    if counts["abstract"] > ABSTRACT_MAX or counts["main"] > MANUSCRIPT_MAX:
        raise ValueError("Letter word limits violated: " + repr(counts))
    title_page = _read_utf8(package / "TITLE_PAGE_TEMPLATE.md")
    if "**Abstract word count:** " + str(counts["abstract"]) not in title_page:
        raise ValueError("Title-page abstract count does not match")
    if "**Main-text word count:** " + format(counts["main"], ",") not in title_page:
        raise ValueError("Title-page body count does not match")

    import_dir = out_dir / "REVIEW_ONLY"
    import_dir.mkdir(parents=True, exist_ok=True)
    sources = import_dir / "source_markdown"
    sources.mkdir(exist_ok=True)
    for name in DOC_FILES:
        path = package / name
        if not path.is_file():
            raise FileNotFoundError(path)
        shutil.copyfile(path, sources / name)

    main_manifest = json.loads(_read_utf8(package / "figures/FIGURE_MANIFEST.json"))
    if [rec["figure"] for rec in main_manifest["figures"]] != list(range(1, 7)):
        raise ValueError("Main Figure 1–6 sequence is incomplete")
    main_pdf_out = import_dir / "main_figures"
    main_pdf_out.mkdir(exist_ok=True)
    png_dir = out_dir / "_generated_png"
    png_dir.mkdir(parents=True, exist_ok=True)
    main_images: list[tuple[str, Path]] = []
    for record in main_manifest["figures"]:
        pdf = package / "figures" / record["pdf"]
        svg = package / "figures" / record["svg"]
        if not pdf.is_file() or not svg.is_file():
            raise FileNotFoundError(str(pdf) + " / " + str(svg))
        if pdf.stat().st_size != record["pdf_bytes"] or svg.stat().st_size != record["svg_bytes"]:
            raise ValueError("Figure manifest byte-size mismatch: " + pdf.name)
        shutil.copyfile(pdf, main_pdf_out / pdf.name)
        png = png_dir / ("Main_Fig_" + str(record["figure"]) + ".png")
        _figure_pdf_to_png(pdf, png)
        main_images.append(("Figure " + str(record["figure"]), png))

    si_fig_dir = package / "supplement/figures"
    supplement_pdfs = sorted(si_fig_dir.glob("Figure_S[1-7]_*.pdf"))
    if len(supplement_pdfs) != 7 or {
        int(re.search(r"Figure_S(\d+)_", f.name).group(1))
        for f in supplement_pdfs
    } != set(range(1, 8)):
        raise ValueError("Expected exactly seven supplementary figure PDFs S1–S7")
    sup_pdf_out = import_dir / "supplementary_figures"
    sup_pdf_out.mkdir(exist_ok=True)
    sup_images: list[tuple[str, Path]] = []
    for figure in supplement_pdfs:
        shutil.copyfile(figure, sup_pdf_out / figure.name)
        no = int(re.search(r"Figure_S(\d+)_", figure.name).group(1))
        png = png_dir / ("Supp_Fig_S" + str(no) + ".png")
        _figure_pdf_to_png(figure, png)
        sup_images.append(("Supplementary Figure S" + str(no), png))

    table_out = import_dir / "supplementary_tables"
    table_out.mkdir(exist_ok=True)
    for name in TABLE_FILES:
        table = package / "supplement" / name
        if not table.is_file():
            raise FileNotFoundError(table)
        shutil.copyfile(table, table_out / name)

    word_dir = import_dir / "word_review"
    _build_doc(package / "MANUSCRIPT.md",
               word_dir / "Manuscript_WITH_FIGURES_REVIEW.docx",
               "manuscript", main_images, package / "FIGURE_CAPTIONS.md")
    _build_doc(package / "SUPPLEMENTARY_INFORMATION.md",
               word_dir / "Supplementary_Information_REVIEW.docx",
               "supplementary information", sup_images)
    _build_doc(package / "TITLE_PAGE_TEMPLATE.md",
               word_dir / "Title_Page_INCOMPLETE_REVIEW.docx",
               "title page")
    _build_doc(package / "COVER_LETTER_DRAFT.md",
               word_dir / "Cover_Letter_UNAPPROVED_REVIEW.docx",
               "cover letter")

    # PDF is a QA/export convenience only; print-ready journal files require author sign-off.
    import subprocess
    executable = shutil.which("libreoffice") or shutil.which("soffice")
    if executable:
        docs = sorted(word_dir.glob("*.docx"))
        cmd = [executable, "-env:UserInstallation=file:///tmp/island-ch1-lo-review",
               "--headless", "--convert-to", "pdf", "--outdir", str(word_dir),
               *[str(p) for p in docs]]
        result = subprocess.run(cmd, check=False, capture_output=True, text=True, timeout=240)
        if result.returncode != 0:
            raise RuntimeError("LibreOffice PDF conversion failed: " + result.stderr)
        for doc in docs:
            pdf = doc.with_suffix(".pdf")
            if not pdf.is_file():
                raise RuntimeError("Missing rendered PDF: " + str(pdf))
            with fitz.open(pdf) as check:
                if len(check) < 1 or not "".join(page.get_text() for page in check).strip():
                    raise RuntimeError("Empty review PDF: " + str(pdf))
    else:
        raise RuntimeError("LibreOffice is required for the review PDF bundle")

    readme = (
        "# Island Chapter 1 — editorial review packet (NOT a submission)\n\n"
        + REVIEW_WARNING + "\n\n"
        "Contains editable review DOCX files, Markdown sources, 6 main-figure PDFs, "
        "7 supplementary-figure PDFs, and six deterministic aggregate tables. "
        "The publication-use figure captions are appended to the review manuscript.\n\n"
        "EXCLUDED: full trait ledger, unlicensed excerpts, raw third-party datasets, "
        "taxon/source row-level tables, release DOI, confirmed author identities.\n\n"
        "OPEN GATES:\n"
        "1. Editor guidance on 176,414 restricted/review-required trait cells.\n"
        "2. All authors' final metadata and approval.\n"
        "3. Rights-safe software/results release tag and DOI (requires explicit approval).\n"
        "4. Final journal-typeset audit. This is a layout review, not a submission export.\n\n"
        "Scientific baseline: config/chapter1_submission_current.json\n"
        "Source commit: " + source_commit + "\n"
        "Word counts: abstract " + str(counts["abstract"])
        + ", main " + str(counts["main"]) + "\n"
    )
    (import_dir / "README_REVIEW_ONLY.txt").write_text(readme, encoding="utf-8")
    manifest = {
        "contract": "chapter1_rights_safe_editorial_review_v1",
        "ready_for_submission": False,
        "contains_restricted_trait_ledger": False,
        "source_commit": source_commit,
        "selector": "config/chapter1_submission_current.json",
        "word_counts": counts,
        "main_figures": 6,
        "supplement_figures": 7,
        "rights_gate": "UNRESOLVED",
        "author_metadata_gate": "UNRESOLVED",
        "software_doi_gate": "UNRESOLVED",
        "files": {},
    }
    for path in sorted(import_dir.rglob("*")):
        if path.is_file():
            rel = str(path.relative_to(import_dir))
            manifest["files"][rel] = {"size": path.stat().st_size,
                                       "sha256": _sha256(path)}
    (import_dir / "PACKAGE_MANIFEST_REVIEW.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    archive = out_dir / "Chapter1_Editorial_REVIEW_ONLY.zip"
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED,
                         compresslevel=7) as handle:
        for path in sorted(import_dir.rglob("*")):
            if path.is_file():
                handle.write(path, path.relative_to(out_dir))
    return {"archive": str(archive), "manifest": manifest,
            "docs": [str(p) for p in sorted(word_dir.glob("*.docx"))]}

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--source-commit", default=os.environ.get("GITHUB_SHA", "UNFROZEN"))
    args = ap.parse_args()
    print(json.dumps(build(args.root, args.output_dir,
                           source_commit=args.source_commit), indent=2))

if __name__ == "__main__":
    main()
