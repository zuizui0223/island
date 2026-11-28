"""Render the integrated GEB candidate without changing active submission files.

Author review only: original H1–H4 remain frozen. Explicit allowlist of
anonymous aggregate tables, figure derivatives and approved candidate prose;
no restricted source ledger, no external publication or submission.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

import cairosvg
import fitz
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from PIL import Image, ImageDraw

ROOT=Path(__file__).resolve().parents[2]
R=ROOT/"research/chapter1_fd_integration_20261125/geb_review"
P=R.parent
sys.path.insert(0,str(ROOT/"scripts/submission"))
from build_chapter1_editorial_review import _add_markdown, _setup_document
from build_chapter1_geb_review import anonymous_document, convert_review_pdfs, counts
from render_chapter1_geb_figures import render as render_geb

BANNED=re.compile(r"(?i)zuizui0223|zhang ruiqi|@gmail\.com|@tohoku|10\.5281/zenodo\.22704973")
ORIGINAL_TABLES=("Table_S1_data_summary.csv","Table_S2_H1_traitwise.csv","Table_S3a_H2_all_finite_cluster.csv","Table_S3b_H2_direct_finite_cluster.csv","Table_S5_H3_pollen_limitation.csv","Table_S6a_H4_scores.csv")
NEW_TABLES=("Table_S8a_flower_pair_isolation.csv","Table_S8b_reproductive_conditioning.csv","Table_S8c_depth_specific_components.csv","Table_S8d_grade_matched_diagnostic.csv","Table_S8e_GloPL_joint_support.csv")
def read(p):return p.read_text(encoding="utf-8")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def clear_para(p):
    for el in list(p._p):
        if el.tag!=qn("w:pPr"):p._p.remove(el)
def configure(doc,kind,anonymous):
    # The original helper adds a review-only preamble. Remove only that
    # presentation text; keep the actual candidate source intact.
    for p in list(doc.paragraphs):
        p._element.getparent().remove(p._element)
    for section in doc.sections:
        section.top_margin=Inches(.82);section.bottom_margin=Inches(.8)
        section.left_margin=Inches(.95);section.right_margin=Inches(.8)
        h=section.header.paragraphs[0];clear_para(h);h.alignment=WD_ALIGN_PARAGRAPH.RIGHT
        h.add_run(kind+" · AUTHOR REVIEW ONLY").font.size=Pt(8)
        f=section.footer.paragraphs[0];clear_para(f);f.alignment=WD_ALIGN_PARAGRAPH.CENTER
        f.add_run("Not for submission  |  ").font.size=Pt(8)
        field=OxmlElement("w:fldSimple");field.set(qn("w:instr"),"PAGE");f._p.append(field)
    if anonymous:
        props=doc.core_properties
        props.author="Anonymous";props.last_modified_by="Anonymous";props.title=kind
        props.subject="";props.keywords="";props.comments="";props.identifier=""
    return doc
def repeat_headers(doc):
    for table in doc.tables:
        if not table.rows:continue
        h=OxmlElement("w:tblHeader");h.set(qn("w:val"),"true")
        table.rows[0]._tr.get_or_add_trPr().append(h)
        for row in table.rows:
            # Prevent a single row splitting across pages. Existing cell
            # widths remain renderer-controlled.
            row._tr.get_or_add_trPr().append(OxmlElement("w:cantSplit"))
def legend_groups(text):
    parts=re.split(r"(?=^## Figure (?:[1-6]|S8) \|)",text,flags=re.MULTILINE)
    groups={}
    for s in parts:
        a=re.match(r"^## Figure ([1-6]|S8) \|",s)
        if a:groups[a.group(1)]=s.strip()
    assert set(groups)==set("123456")|{"S8"}
    return groups
def append_figure(doc,pdf,title,work,caption=None):
    doc.add_page_break()
    if caption:_add_markdown(doc,caption)
    else:doc.add_heading(title,level=2)
    with fitz.open(pdf) as d:
        assert len(d)==1
        p=d[0];img=work/(pdf.stem+".png")
        p.get_pixmap(matrix=fitz.Matrix(2.6,2.6),alpha=False).save(img)
        # Fit both width and height, leaving room for even the longer legends.
        w=min(6.25,5.7*p.rect.width/p.rect.height)
        par=doc.add_paragraph();par.alignment=WD_ALIGN_PARAGRAPH.CENTER
        par.paragraph_format.space_after=Pt(0)
        par.add_run().add_picture(str(img),width=Inches(w))
    img.unlink()
def save_blind_doc(source,path,kind,figs,work,legends=None):
    text=read(source)
    assert not BANNED.search(text),source.name
    d=configure(anonymous_document(kind,source.name),kind,True)
    _add_markdown(d,text)
    repeat_headers(d)
    for ident,p in figs:
        append_figure(d,p,"Supplementary Figure "+ident if ident.startswith("S") else "Figure "+ident,work,legends.get(ident) if legends else None)
    d.save(path)
def safe_pdf_copy(src,dst):
    # Copy only figure pages, not original PDF metadata or attachments.
    with fitz.open(src) as original, fitz.open() as d:
        assert len(original)==1
        d.insert_pdf(original,links=False,annots=False)
        d.set_metadata({"author":"Anonymous","creator":"Figure renderer","title":dst.stem})
        d.save(dst,garbage=4,deflate=True)
def raster_qa(pdf,folder,anonymous):
    folder.mkdir(parents=True,exist_ok=True)
    entries=[]
    with fitz.open(pdf) as d:
        if anonymous:
            assert not BANNED.search(json.dumps(d.metadata,ensure_ascii=False))
        assert len(d)>0
        for i,p in enumerate(d,1):
            t=p.get_text()
            if anonymous:assert not BANNED.search(t),(pdf.name,i)
            # The raster is retained for all-page human visual proof. Bounds
            # checks are a structural preflight, not proof of zero overlap.
            bad=[]
            for block in p.get_text("blocks"):
                x0,y0,x1,y1=block[:4]
                if min(x0,y0)<-1 or x1>p.rect.width+1 or y1>p.rect.height+1:bad.append("text_outside_page")
            for img in p.get_image_info():
                x0,y0,x1,y1=img["bbox"]
                if min(x0,y0)<-1 or x1>p.rect.width+1 or y1>p.rect.height+1:bad.append("image_outside_page")
            # Ignore the review-only header/footer and page number when checking
            # for blank content pages. Flag rather than silently remove.
            main=p.get_text(clip=fitz.Rect(45,48,p.rect.width-32,p.rect.height-45)).strip()
            if not main and not p.get_image_info():bad.append("possible_blank_content_page")
            png=folder/f"page-{i:03d}.png"
            p.get_pixmap(matrix=fitz.Matrix(1.65,1.65),alpha=False).save(png)
            entries.append({"page":i,"characters":len(t),"image_count":len(p.get_image_info()),"flags":sorted(set(bad))})
    # Contact sheets ease review; individual PNGs are available at larger size.
    imgs=sorted(folder.glob("page-*.png"))
    for start in range(0,len(imgs),12):
        sheet=Image.new("RGB",(4*300,3*415),"#e7ebef");draw=ImageDraw.Draw(sheet)
        for k,png in enumerate(imgs[start:start+12]):
            with Image.open(png) as im:
                im.thumbnail((284,387));x=(k%4)*300+(300-im.width)//2;y=(k//4)*415+8
                sheet.paste(im,(x,y));draw.text(((k%4)*300+12,(k//4)*415+396),png.stem,fill="#1c2530")
        sheet.save(folder/f"contact-{start//12+1:02d}.png")
    return entries
def check_docx(path,anonymous):
    d=Document(path)
    if anonymous:
        assert d.core_properties.author==d.core_properties.last_modified_by=="Anonymous"
        with zipfile.ZipFile(path) as z:
            for n in z.namelist():
                if n.endswith((".xml",".rels")):
                    assert not BANNED.search(z.read(n).decode("utf-8",errors="replace")),(path.name,n)
        assert d.sections[0]._sectPr.find(qn("w:lnNumType")) is not None
def build(out,commit):
    subprocess.run([sys.executable,str(P/"verify_review_package.py")],check=True)
    subprocess.run([sys.executable,str(R/"check_geb_candidate.py")],check=True)
    out.mkdir(parents=True,exist_ok=True)
    bundle=out/"GEB_INTEGRATED_AUTHOR_REVIEW_ONLY";blind=bundle/"BLINDED_REVIEW";ident=bundle/"IDENTIFYING_ONLY";work=out/"working"
    for p in (blind,ident,work):p.mkdir(parents=True,exist_ok=True)
    main_fig=blind/"MAIN_FIGURES";supp_fig=blind/"SUPPLEMENTARY_FIGURES";tables=blind/"AGGREGATE_RESULTS"
    for p in (main_fig,supp_fig,tables):p.mkdir()
    g=render_geb(ROOT,work/"main")
    assert g["figure1_equal_area"] is True and g["scientific_refit"] is False
    manifest=json.loads(read(work/"main/FIGURE_MANIFEST.json"))
    assert len(manifest["figures"])==6
    mains=[]
    for i,r in enumerate(manifest["figures"],1):
        dst=main_fig/f"Figure_{i}.pdf";safe_pdf_copy(work/"main"/r["pdf"],dst);mains.append((str(i),dst))
    supp=[]
    for i in range(1,8):
        paths=list((ROOT/"submission/chapter1_current/supplement/figures").glob(f"Figure_S{i}_*.pdf"))
        assert len(paths)==1
        dst=supp_fig/f"Figure_S{i}.pdf";safe_pdf_copy(paths[0],dst);supp.append((f"S{i}",dst))
    svg=read(R/"Figure_S8.svg");assert not BANNED.search(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(),write_to=str(work/"S8_raw.pdf"))
    dst=supp_fig/"Figure_S8.pdf";safe_pdf_copy(work/"S8_raw.pdf",dst);supp.append(("S8",dst))
    legends=legend_groups(read(R/"FIGURE_CAPTIONS_REVIEW.md"))
    main_doc=blind/"Manuscript_GEB_Integrated_REVIEW.docx"
    si_doc=blind/"Supporting_Information_GEB_Integrated_REVIEW.docx"
    save_blind_doc(R/"MANUSCRIPT_BLINDED_REVIEW.md",main_doc,"GEB manuscript",mains,work,legends)
    save_blind_doc(R/"SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW.md",si_doc,"GEB supporting information",supp,work,{"S8":legends["S8"]})
    # Metadata placeholders remain segregated and are never inferred.
    for source,filename in [(R/"COVER_LETTER_REVIEW.md","Cover_Letter_UNAPPROVED.docx"),(ROOT/"submission/geb_current/TITLE_PAGE_IDENTIFYING_TEMPLATE.md","Title_Page_IDENTIFYING_INCOMPLETE.docx")]:
        d=configure(_setup_document("GEB author-only",source.name),"Author-only",False)
        _add_markdown(d,read(source));d.save(ident/filename)
    for name in ("MANUSCRIPT_BLINDED_REVIEW.md","SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW.md","FIGURE_CAPTIONS_REVIEW.md"):
        shutil.copy2(R/name,blind/name)
    for name in ("COVER_LETTER_REVIEW.md","README.md"):
        shutil.copy2(R/name,ident/name)
    for name in ("TITLE_PAGE_IDENTIFYING_TEMPLATE.md","DATA_CODE_AVAILABILITY.md","SUBMISSION_CHECKLIST_GEB.md"):
        shutil.copy2(ROOT/"submission/geb_current"/name,ident/name)
    for name in ORIGINAL_TABLES:shutil.copy2(ROOT/"submission/chapter1_current/supplement"/name,tables/name)
    for name in NEW_TABLES:shutil.copy2(P/"supplement"/name,tables/name)
    assert len(list(tables.glob("*.csv")))==11
    convert_review_pdfs(bundle)
    page_audit={}
    for docx in sorted(bundle.rglob("*.docx")):
        anon=docx.parent==blind;check_docx(docx,anon)
        pdf=docx.with_suffix(".pdf");assert pdf.is_file()
        page_audit[docx.stem]=raster_qa(pdf,out/"PAGE_QA"/docx.stem,anon)
    # Validate the actual rendered text includes the new result and SI.
    with fitz.open(main_doc.with_suffix(".pdf")) as d:
        t=" ".join(p.get_text() for p in d)
        for marker in ("0.02762","0.01967","0.02463","0.02089","0.000221"):assert marker in t,marker
    with fitz.open(si_doc.with_suffix(".pdf")) as d:
        assert "S8" in " ".join(p.get_text() for p in d)
    # Only put anonymous aggregate derivatives in the reviewer-facing tree.
    for p in blind.rglob("*"):
        if p.suffix in {".md",".csv",".txt"} and p.is_file():assert not BANNED.search(read(p)),p.name
    flags=[{"document":k,"page":p["page"],"flags":p["flags"]} for k,a in page_audit.items() for p in a if p["flags"]]
    aw,bw=counts(read(R/"MANUSCRIPT_BLINDED_REVIEW.md"))
    audit={"source_commit":commit,"ready_to_submit":False,"raw_restricted_cells_exported":False,"original_science_selector_changed":False,"new_H4_fit":False,"figure1_equal_area":True,"main_figures":6,"supplementary_figures":8,"aggregate_tables":11,"abstract_words":aw,"body_words":bw,"page_counts":{k:len(v) for k,v in page_audit.items()},"page_flags":flags,"all_pages_rasterized":True,"human_all_page_visual_proof":False,"automated_anonymity_scan_passed":True}
    (out/"PAGE_QA"/"PAGE_AUDIT.json").write_text(json.dumps({"summary":audit,"pages":page_audit},indent=2)+"\n")
    (bundle/"README_DO_NOT_SUBMIT.txt").write_text("GEB integrated candidate — author review only. Do not send or submit. Separate BLINDED_REVIEW and IDENTIFYING_ONLY. Rights, anonymous access, coauthor approvals, archive authorization and visual proof remain open. No raw restricted data included.\n")
    (bundle/"REVIEW_MANIFEST.json").write_text(json.dumps({"audit":audit,"files":{str(p.relative_to(bundle)):{"sha256":sha(p),"bytes":p.stat().st_size} for p in sorted(bundle.rglob("*")) if p.is_file()}},indent=2)+"\n")
    z=out/"Chapter1_GEB_Integrated_AUTHOR_REVIEW_ONLY.zip"
    with zipfile.ZipFile(z,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as a:
        for p in sorted(bundle.rglob("*")):
            if p.is_file():a.write(p,p.relative_to(bundle.parent))
    print("GEB_INTEGRATED_RENDER_RESULT "+json.dumps(audit,allow_nan=False),flush=True)
    return audit
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);p.add_argument("--commit",required=True)
    a=p.parse_args();build(a.out,a.commit)
