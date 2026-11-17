"""GEB first-shot: double-anonymous author-review packet from frozen Chapter 1.

Presentation only. Does not refit scientific analyses, publish data, create a DOI,
or submit to the journal. Does not redistribute source-restricted trait cells.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import fitz
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_chapter1_editorial_review import _add_markdown, _setup_document
from render_chapter1_geb_figures import render as render_geb_figures

ABSTRACT_TITLES = ("Aim", "Location", "Time period", "Major taxa studied",
                   "Methods", "Results", "Main conclusions")
WARN = "GEB EDITORIAL REVIEW COPY. NOT SUBMISSION-READY. AUTHOR AND DATA-RIGHTS APPROVALS PENDING."


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_markdown(root: Path, path: str) -> str:
    return (root / path).read_text(encoding="utf-8")


def counts(text: str) -> tuple[int, int]:
    abstract = text.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0]
    body = "## Introduction" + text.split("## Introduction", 1)[1].split("## References", 1)[0]
    return len(abstract.split()), len(body.split())


def ensure_content(root: Path) -> tuple[dict, str]:
    contract = json.loads(load_markdown(root, "config/chapter1_geb_submission.json"))
    if contract["status"] != "GEB_FIRST_SHOT_DRAFT_NOT_SUBMISSION_READY":
        raise ValueError("Submission gate must remain closed")
    science = json.loads(load_markdown(root, contract["science_selector"]))
    if science.get("status") != "primary_submission_baseline":
        raise ValueError("Scientific baseline not locked")
    text = load_markdown(root, contract["blinded_manuscript"])
    abstract = text.split("## Abstract", 1)[1].split("**Keywords:**", 1)[0]
    if len(abstract.split()) > 300:
        raise ValueError("GEB abstract exceeds 300 words")
    for title in ABSTRACT_TITLES:
        if abstract.count("**" + title + ":**") != 1:
            raise ValueError("Missing/duplicate GEB abstract section " + title)
    kws = text.split("**Keywords:**", 1)[1].split("\n", 1)[0].strip()
    keywords = [v.strip() for v in kws.split(",")]
    if not 6 <= len(keywords) <= 10:
        raise ValueError("GEB requires 6–10 keywords")
    if keywords != sorted(keywords, key=str.casefold):
        raise ValueError("Keywords must be alphabetized")
    if re.search(r"(?i)github\\.com/zuizui0223|@gmail\\.com|@tohoku|zhang ruiqi|\brei?chel\b", text):
        raise ValueError("Possible author identity in blinded manuscript")
    if "PENDING" not in text or "222,688" not in text or "46,274" not in text:
        raise ValueError("The required access limitations are missing")
    if "**Running title:**" not in text:
        raise ValueError("Running title missing")
    running = text.split("**Running title:**", 1)[1].split("\n", 1)[0].strip()
    if len(running) >= 40:
        raise ValueError("Running title exceeds 39 characters")
    return contract, text


def anonymous_document(kind: str, source: str):
    doc = _setup_document(kind, source)
    props = doc.core_properties
    props.author = "Anonymous"
    props.last_modified_by = "Anonymous"
    props.comments = ""
    props.identifier = ""
    props.keywords = ""
    props.subject = ""
    # GEB requires continuous line numbering in the blinded research manuscript.
    sect = doc.sections[0]._sectPr
    ln = OxmlElement("w:lnNumType")
    ln.set(qn("w:countBy"), "1")
    ln.set(qn("w:start"), "1")
    ln.set(qn("w:restart"), "continuous")
    sect.append(ln)
    return doc


def add_figures(doc, figure_files: list[Path], captions: str | None = None) -> None:
    caption_by_index: dict[int, str] = {}
    if captions:
        for group in re.split(r"(?=^## Figure [1-6] \|)", captions, flags=re.MULTILINE):
            match = re.match(r"^## Figure ([1-6]) \|([^\n]*)", group)
            if match:
                caption_by_index[int(match.group(1))] = group.strip()
    for idx, path in enumerate(figure_files, 1):
        doc.add_page_break()
        if captions:
            caption = caption_by_index.get(idx)
            if caption is None:
                raise ValueError("Figure " + str(idx) + " missing journal legend")
            _add_markdown(doc, caption)
        else:
            doc.add_heading("Supplementary Figure S" + str(idx), 2)
        with fitz.open(path) as pdf:
            if len(pdf) < 1:
                raise ValueError("Empty figure PDF " + str(path))
            pix = pdf[0].get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
            tmp = path.parent / ("_inline_" + str(idx) + ".png")
            pix.save(tmp)
            doc.add_picture(str(tmp), width=Inches(6.45))
            tmp.unlink()


def convert_review_pdfs(folder: Path) -> None:
    exe = shutil.which("libreoffice") or shutil.which("soffice")
    if not exe:
        raise RuntimeError("LibreOffice is required for the GEB review PDF")
    docs = sorted(folder.rglob("*.docx"))
    for d in docs:
        profile = "file:///tmp/island-geb-" + d.stem[:30]
        cmd = [exe, "-env:UserInstallation=" + profile,
               "--headless", "--convert-to", "pdf", "--outdir",
               str(d.parent), str(d)]
        result = subprocess.run(cmd, check=False, capture_output=True,
                                text=True, timeout=180)
        if result.returncode != 0 or not d.with_suffix(".pdf").is_file():
            raise RuntimeError("GEB PDF conversion failed for " + d.name + ": " + result.stderr)
        with fitz.open(d.with_suffix(".pdf")) as pdf:
            if len(pdf) < 1 or not pdf[0].get_text().strip():
                raise ValueError("GEB review PDF is empty: " + d.name)


def build(root: Path, out: Path, sha: str) -> dict:
    root = root.resolve()
    out = out.resolve()
    contract, blind = ensure_content(root)
    outfile = out / "GEB_REVIEW_ONLY"
    blind_dir = outfile / "BLINDED_REVIEW"
    identify_dir = outfile / "IDENTIFYING_ONLY"
    blind_dir.mkdir(parents=True, exist_ok=True)
    identify_dir.mkdir(parents=True, exist_ok=True)
    figure_dir = out / "_working_figures"
    map_audit = render_geb_figures(root, figure_dir)
    fmanifest = json.loads((figure_dir / "FIGURE_MANIFEST.json").read_text())
    figure_pdfs = [figure_dir / entry["pdf"] for entry in fmanifest["figures"]]
    if len(figure_pdfs) != 6 or map_audit["figure1_equal_area"] is not True:
        raise ValueError("GEB Figure 1 equal-area map is not verified")
    doc = anonymous_document("GEB blinded manuscript", "MANUSCRIPT_BLINDED.md")
    _add_markdown(doc, blind)
    captions = load_markdown(root, contract["figure_captions"])
    add_figures(doc, figure_pdfs, captions)
    doc.save(blind_dir / "Manuscript_GEB_BLINDED_REVIEW.docx")

    si = anonymous_document("GEB supporting information", "SUPPLEMENTARY_INFORMATION.md")
    _add_markdown(si, load_markdown(root, "submission/geb_current/SUPPLEMENTARY_INFORMATION.md"))
    si_dir = root / "submission/chapter1_current/supplement/figures"
    supplemental = sorted(si_dir.glob("Figure_S[1-7]_*.pdf"))
    if len(supplemental) != 7:
        raise ValueError("Missing one or more SI figure PDFs")
    add_figures(si, supplemental)
    si.save(blind_dir / "Supporting_Information_GEB_BLINDED_REVIEW.docx")

    # Identifying title page and editorial cover letter remain segregated.
    for src, filename in (
        ("submission/geb_current/TITLE_PAGE_IDENTIFYING_TEMPLATE.md", "Title_Page_IDENTIFYING_INCOMPLETE.docx"),
        ("submission/geb_current/COVER_LETTER_GEB.md", "Cover_Letter_UNAPPROVED.docx"),
    ):
        d = _setup_document("GEB author-only", Path(src).name)
        _add_markdown(d, load_markdown(root, src))
        d.save(identify_dir / filename)
    convert_review_pdfs(outfile)

    fig_out = blind_dir / "FIGURES_GEB"
    fig_out.mkdir(exist_ok=True)
    for p in figure_pdfs:
        shutil.copy2(p, fig_out / p.name)
    shutil.copy2(figure_dir / "GEB_FIGURE_AUDIT.json", fig_out / "GEB_FIGURE_AUDIT.json")

    for name in ["MANUSCRIPT_BLINDED.md", "FIGURE_CAPTIONS_GEB.md",
                 "SUPPLEMENTARY_INFORMATION.md"]:
        shutil.copy2(root / "submission/geb_current" / name, blind_dir / name)
    for name in ["TITLE_PAGE_IDENTIFYING_TEMPLATE.md", "COVER_LETTER_GEB.md",
                 "DATA_CODE_AVAILABILITY.md", "SUBMISSION_CHECKLIST_GEB.md"]:
        shutil.copy2(root / "submission/geb_current" / name, identify_dir / name)
    for name in ["Table_S1_data_summary.csv", "Table_S2_H1_traitwise.csv",
                 "Table_S3a_H2_all_finite_cluster.csv",
                 "Table_S3b_H2_direct_finite_cluster.csv",
                 "Table_S5_H3_pollen_limitation.csv", "Table_S6a_H4_scores.csv"]:
        table = root / "submission/chapter1_current/supplement" / name
        if not table.is_file():
            raise FileNotFoundError(table)
        tables = blind_dir / "AGGREGATE_RESULTS"
        tables.mkdir(exist_ok=True)
        shutil.copy2(table, tables / name)

    readme = (
        "# GEB FIRST-SHOT EDITORIAL REVIEW PACKAGE — NOT FOR JOURNAL SUBMISSION\n"
        + WARN + "\n"
        "BLINDED_REVIEW contains the double-anonymous manuscript/SI in DOCX and PDF "
        "and corrected maps/results; IDENTIFYING_ONLY must remain separate at "
        "Wiley upload. The identifying metadata is currently incomplete.\n"
        "H1–H4 inference is unchanged. GEB structured abstract, keywords and equal-area "
        "world map are journal-specific editorial adaptations.\n"
        "The 176,414 third-party trait cells are source-restricted or rights-review "
        "required; no blanket exception has been received.\n"
        "No software/results DOI or anonymous reviewer link is yet available.\n"
        "SOURCE COMMIT: " + sha + "\n"
    )
    (outfile / "README_DO_NOT_SUBMIT.txt").write_text(readme)
    inventory = {
        "contract": contract["contract"],
        "ready_to_submit": False,
        "source_commit": sha,
        "scientific_results_unchanged": True,
        "contains_restricted_trait_ledger": False,
        "double_anonymous_review_only": True,
        "abstract_main_word_counts": dict(zip(("abstract", "body"), counts(blind), strict=True)),
        "figure1_equal_area": True,
        "figure_count": 6,
        "supplement_figure_count": 7,
        "files": {
            str(p.relative_to(outfile)): {
                "sha256": sha256(p), "bytes": p.stat().st_size,
            } for p in sorted(outfile.rglob("*")) if p.is_file()
        },
    }
    (outfile / "MANIFEST_REVIEW_ONLY.json").write_text(json.dumps(inventory, indent=2)+"\n")
    zip_path = out / "Chapter1_GEB_DoubleAnonymous_REVIEW_ONLY.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=7) as z:
        for p in sorted(outfile.rglob("*")):
            if p.is_file():
                z.write(p, arcname=str(p.relative_to(out)))

    # Render first-page PNGs for visual inspection.
    previews = out / "PREVIEW_QA"
    previews.mkdir(exist_ok=True)
    for inp in (
        blind_dir / "Manuscript_GEB_BLINDED_REVIEW.pdf",
        blind_dir / "Supporting_Information_GEB_BLINDED_REVIEW.pdf",
    ):
        with fitz.open(inp) as pdf:
            pdf[0].get_pixmap(matrix=fitz.Matrix(1.6, 1.6)).save(
                previews / (inp.stem + "_PAGE1.png"))
    with fitz.open(figure_pdfs[0]) as pdf:
        pdf[0].get_pixmap(matrix=fitz.Matrix(1.5, 1.5)).save(
            previews / "GEB_Figure1_Mollweide.png")
    return {"zip": str(zip_path), "manifest": inventory}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--source-commit", default=os.environ.get("GITHUB_SHA", "UNFROZEN"))
    args = ap.parse_args()
    print(json.dumps(build(args.root, args.output_dir, args.source_commit), indent=2))


if __name__ == "__main__":
    main()
