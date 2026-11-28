"""Versioned transport author-review render, retaining original figure and document helpers."""
import argparse,json,shutil,subprocess,sys,zipfile
from pathlib import Path
import cairosvg
import fitz
from build_chapter1_geb_integrated_candidate import (
    ROOT,R,P,BANNED,ORIGINAL_TABLES,NEW_TABLES,read,sha,render_geb,
    safe_pdf_copy,legend_groups,save_blind_doc,configure,_setup_document,
    _add_markdown,convert_review_pdfs,check_docx,raster_qa,counts,
)

def build(out,commit):
    subprocess.run([sys.executable,str(P/"verify_review_package.py")],check=True)
    subprocess.run([sys.executable,str(R/"check_geb_candidate.py")],check=True)
    out.mkdir(parents=True,exist_ok=True)
    bundle=out/"GEB_TRANSPORT_AUTHOR_REVIEW_ONLY";blind=bundle/"BLINDED_REVIEW";ident=bundle/"IDENTIFYING_ONLY";work=out/"working"
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
    main_doc=blind/"Manuscript_GEB_Transport_REVIEW.docx"
    si_doc=blind/"Supporting_Information_GEB_Transport_REVIEW.docx"
    save_blind_doc(R/"MANUSCRIPT_BLINDED_REVIEW_TRANSPORT_20261128.md",main_doc,"GEB manuscript",mains,work,legends)
    save_blind_doc(R/"SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW_TRANSPORT_20261128.md",si_doc,"GEB supporting information",supp,work,{"S8":legends["S8"]})
    # Metadata placeholders remain segregated and are never inferred.
    for source,filename in [(R/"COVER_LETTER_REVIEW.md","Cover_Letter_UNAPPROVED.docx"),(ROOT/"submission/geb_current/TITLE_PAGE_IDENTIFYING_TEMPLATE.md","Title_Page_IDENTIFYING_INCOMPLETE.docx")]:
        d=configure(_setup_document("GEB author-only",source.name),"Author-only",False)
        _add_markdown(d,read(source));d.save(ident/filename)
    for name in ("MANUSCRIPT_BLINDED_REVIEW_TRANSPORT_20261128.md","SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW_TRANSPORT_20261128.md","FIGURE_CAPTIONS_REVIEW.md"):
        shutil.copy2(R/name,blind/name)
    for name in ("COVER_LETTER_REVIEW.md","README.md"):
        shutil.copy2(R/name,ident/name)
    for name in ("TITLE_PAGE_IDENTIFYING_TEMPLATE.md","DATA_CODE_AVAILABILITY.md","SUBMISSION_CHECKLIST_GEB.md"):
        shutil.copy2(ROOT/"submission/geb_current"/name,ident/name)
    for name in ORIGINAL_TABLES:shutil.copy2(ROOT/"submission/chapter1_current/supplement"/name,tables/name)
    for name in NEW_TABLES:shutil.copy2(P/"supplement"/name,tables/name)
    subprocess.run([sys.executable,str(R/"check_transport_candidate.py"),"--out",str(tables)],check=True)
    assert len(list(tables.glob("*.csv")))==14
    convert_review_pdfs(bundle)
    page_audit={}
    for docx in sorted(bundle.rglob("*.docx")):
        anon=docx.parent==blind;check_docx(docx,anon)
        pdf=docx.with_suffix(".pdf");assert pdf.is_file()
        page_audit[docx.stem]=raster_qa(pdf,out/"PAGE_QA"/docx.stem,anon)
    # Validate the actual rendered text includes the new result and SI.
    with fitz.open(main_doc.with_suffix(".pdf")) as d:
        t=" ".join(p.get_text() for p in d)
        for marker in ("0.02762","0.01967","0.02463","0.02089","0.000221","0.00397","0.00371"):assert marker in t,marker
    with fitz.open(si_doc.with_suffix(".pdf")) as d:
        text=" ".join(p.get_text() for p in d)
        for marker in ("S8.7","S8.8","S8.9","0.05087","0.01298","0.01066"):assert marker in text,marker
    # Only put anonymous aggregate derivatives in the reviewer-facing tree.
    for p in blind.rglob("*"):
        if p.suffix in {".md",".csv",".txt"} and p.is_file():assert not BANNED.search(read(p)),p.name
    flags=[{"document":k,"page":p["page"],"flags":p["flags"]} for k,a in page_audit.items() for p in a if p["flags"]]
    aw,bw=counts(read(R/"MANUSCRIPT_BLINDED_REVIEW_TRANSPORT_20261128.md"))
    audit={"source_commit":commit,"ready_to_submit":False,"raw_restricted_cells_exported":False,"original_science_selector_changed":False,"new_H4_fit":False,"figure1_equal_area":True,"main_figures":6,"supplementary_figures":8,"aggregate_tables":14,"transport_table_rows":20,"source_receipts_replayed":True,"abstract_words":aw,"body_words":bw,"page_counts":{k:len(v) for k,v in page_audit.items()},"page_flags":flags,"all_pages_rasterized":True,"human_all_page_visual_proof":False,"automated_anonymity_scan_passed":True}
    (out/"PAGE_QA"/"PAGE_AUDIT.json").write_text(json.dumps({"summary":audit,"pages":page_audit},indent=2)+"\n")
    shutil.copy2(P/"MORPHOLOGY_ADOPTION_TRANSPORT_20261128.md",ident/"MORPHOLOGY_ADOPTION_TRANSPORT_20261128.md")
    (bundle/"README_DO_NOT_SUBMIT.txt").write_text("GEB transport candidate — author review only. Do not send or submit. Separate BLINDED_REVIEW and IDENTIFYING_ONLY. Rights, anonymous access, coauthor approvals, archive authorization and visual proof remain open. No raw restricted data included.\n")
    (bundle/"REVIEW_MANIFEST.json").write_text(json.dumps({"audit":audit,"files":{str(p.relative_to(bundle)):{"sha256":sha(p),"bytes":p.stat().st_size} for p in sorted(bundle.rglob("*")) if p.is_file()}},indent=2)+"\n")
    z=out/"Chapter1_GEB_Transport_AUTHOR_REVIEW_ONLY.zip"
    with zipfile.ZipFile(z,"w",compression=zipfile.ZIP_DEFLATED,compresslevel=6) as a:
        for p in sorted(bundle.rglob("*")):
            if p.is_file():a.write(p,p.relative_to(bundle.parent))
    print("GEB_TRANSPORT_RENDER_RESULT "+json.dumps(audit,allow_nan=False),flush=True)
    if flags:raise ValueError("Page structural preflight flags; candidate not verified")
    return audit

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);p.add_argument("--commit",required=True)
    a=p.parse_args();build(a.out,a.commit)
