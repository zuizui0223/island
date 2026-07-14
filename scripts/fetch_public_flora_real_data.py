#!/usr/bin/env python3
"""Fetch substantive public flora data and extract auditable eFloras trait candidates."""
from __future__ import annotations
import argparse, csv, hashlib, html, json, os, re, ssl, time, urllib.parse, urllib.request
from html.parser import HTMLParser
from pathlib import Path

UA={"User-Agent":"island-v2-trait-harvester/1.2 (+https://github.com/zuizui0223/island)"}
MASTER=Path("data/v2/staging/gbif/collected/island_taxa.csv")

class Links(HTMLParser):
    def __init__(self): super().__init__(); self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a":
            href=dict(attrs).get("href")
            if href: self.links.append(href)

class VisibleText(HTMLParser):
    def __init__(self): super().__init__(); self.parts=[]; self.skip=0
    def handle_starttag(self,tag,attrs):
        if tag.lower() in {"script","style","noscript"}: self.skip+=1
        if tag.lower() in {"p","br","div","tr","li","h1","h2","h3","td"}: self.parts.append("\n")
    def handle_endtag(self,tag):
        if tag.lower() in {"script","style","noscript"} and self.skip: self.skip-=1
        if tag.lower() in {"p","div","tr","li","h1","h2","h3"}: self.parts.append("\n")
    def handle_data(self,data):
        if not self.skip: self.parts.append(data)

COLOR_WORDS=("white","cream","yellow","green","orange","red","pink","purple","violet","blue","brown","black")
FORM_WORDS=("campanulate","tubular","funnelform","funnel-shaped","rotate","salverform","urceolate","bilabiate","two-lipped","papilionaceous","spurred")
SYMMETRY={"actinomorphic":"radial","radially symmetric":"radial","radially symmetrical":"radial","zygomorphic":"bilateral","bilaterally symmetric":"bilateral","bilaterally symmetrical":"bilateral"}
SEX_SYSTEM=("bisexual","hermaphroditic","unisexual","monoecious","dioecious","andromonoecious","gynomonoecious")
FLORAL_CONTEXT=re.compile(r"\b(flower|flowers|corolla|petal|petals|perianth|tepals?|calyx)\b",re.I)
SENTENCE_SPLIT=re.compile(r"(?<=[.!?;])\s+|\n+")
BINOMIAL=re.compile(r"\b([A-Z][a-zA-Z.-]+\s+[a-z][a-zA-Z.-]+)\b")


def get(url: str, attempts: int=4, timeout: int=180) -> tuple[bytes,str,str]:
    last=None
    for attempt in range(attempts):
        try:
            req=urllib.request.Request(url,headers=UA)
            with urllib.request.urlopen(req,timeout=timeout,context=ssl.create_default_context()) as r:
                return r.read(), r.geturl(), r.headers.get("content-type","")
        except Exception as exc:
            last=exc
            if attempt+1<attempts: time.sleep(2**attempt)
    raise RuntimeError(f"GET failed after {attempts} attempts: {url}: {last}")


def decode(b: bytes) -> str:
    for enc in ("utf-8-sig","utf-8","cp1252","latin-1"):
        try: return b.decode(enc)
        except UnicodeDecodeError: pass
    return b.decode("utf-8","replace")


def visible_text(raw_html: str) -> str:
    parser=VisibleText(); parser.feed(raw_html)
    return re.sub(r"[ \t]+"," ",html.unescape("".join(parser.parts))).strip()


def master_names() -> set[str]:
    with MASTER.open(encoding="utf-8-sig",newline="") as f:
        return {" ".join((r.get("accepted_species") or "").split()) for r in csv.DictReader(f) if r.get("accepted_species")}


def stable_bucket(value: str, count: int) -> int:
    return int.from_bytes(hashlib.sha1(value.encode("utf-8")).digest()[:8],"big") % count


def focal_species_for_page(text: str, page_matches: set[str]) -> str:
    """Choose the treatment taxon from the page header, not incidental names in keys/references."""
    if not page_matches: return ""
    head=" ".join(text[:4000].split())
    ranked=[]
    for name in page_matches:
        pos=head.find(name)
        if pos>=0: ranked.append((pos,-len(name),name))
    if ranked: return min(ranked)[2]
    return next(iter(page_matches)) if len(page_matches)==1 else ""


def treatment_text_present(text: str, focal: str) -> bool:
    if not focal: return False
    low=text.lower()
    return focal in text[:4000] and any(term in low for term in ("flowers", "corolla", "petals", "perianth", "flowering")) and len(text)>800


def extract_trait_candidates(text: str, focal_species: str, source_url: str, source_file: str) -> list[dict[str,str]]:
    rows=[]; seen=set()
    if not focal_species: return rows
    for sentence in SENTENCE_SPLIT.split(text):
        sentence=" ".join(sentence.split())
        if len(sentence)<12 or len(sentence)>1000: continue
        low=sentence.lower(); hits=[]
        if FLORAL_CONTEXT.search(sentence):
            colors=[word for word in COLOR_WORDS if re.search(rf"\b{re.escape(word)}(?:ish)?\b",low)]
            if colors: hits.append(("flower_primary_color","|".join(colors),"floral colour phrase"))
            forms=[word for word in FORM_WORDS if word in low]
            if forms: hits.append(("floral_form","|".join(forms),"floral form phrase"))
            sym=[value for phrase,value in SYMMETRY.items() if phrase in low]
            if sym: hits.append(("floral_symmetry","|".join(sorted(set(sym))),"floral symmetry phrase"))
        sex=[word for word in SEX_SYSTEM if re.search(rf"\b{re.escape(word)}\b",low)]
        if sex: hits.append(("sex_system","|".join(sex),"sexual system phrase"))
        if "self-incompat" in low or "self incompatible" in low:
            hits.append(("self_incompatibility","self_incompatible","self-incompatibility phrase"))
        elif "self-compat" in low or "self compatible" in low:
            hits.append(("self_incompatibility","self_compatible","self-compatibility phrase"))
        if re.search(r"\b(autogam|autonomous self|spontaneous self)\w*",low):
            hits.append(("autonomous_selfing_capacity","reported_present","autonomous selfing phrase"))
        pollinators=[]
        for phrase,value in (("bee","bee"),("bumblebee","bumblebee"),("butterfl","butterfly"),("moth","moth"),("bird","bird"),("bat","bat"),("fly","fly"),("wind-pollinat","wind")):
            if phrase in low and ("pollinat" in low or "visited by" in low): pollinators.append(value)
        if pollinators: hits.append(("pollination_functional_guild","|".join(sorted(set(pollinators))),"pollinator phrase"))
        for trait,value,rule in hits:
            key=(focal_species,trait,value,sentence)
            if key in seen: continue
            seen.add(key)
            rows.append({"accepted_species":focal_species,"trait_name":trait,"candidate_value":value,"source_name":"eFloras","source_trait":rule,"source_url":source_url,"source_record_id":source_file,"source_excerpt":sentence[:1000],"evidence_scope":"species_direct","evidence_status":"source_backed_rule_extracted_unreviewed"})
    return rows


def fetch_wfo_plant_list(out: Path) -> dict:
    landing="https://wfoplantlist.org/"; body,final,_=get(landing); p=Links(); p.feed(decode(body))
    links=[]
    for href in p.links:
        u=urllib.parse.urljoin(final,href); low=u.lower()
        if any(x in low for x in ("snapshot","download","archive")) or low.endswith((".zip",".csv",".tsv",".txt",".gz")): links.append(u)
    files=[]
    for i,u in enumerate(dict.fromkeys(links[:30]),1):
        try: b,f,c=get(u)
        except Exception: continue
        if "html" in c.lower() and b.lstrip().startswith(b"<"): continue
        name=Path(urllib.parse.urlparse(f).path).name or f"wfo_data_{i}.bin"; out.mkdir(parents=True,exist_ok=True); (out/name).write_bytes(b); files.append({"url":f,"file":name,"bytes":len(b)})
    if not files or sum(x["bytes"] for x in files)<10000: raise RuntimeError("WFO Plant List: no substantive downloadable data")
    report={"source":"WFO Plant List","real_data_acquired":True,"files_saved":len(files),"bytes_saved":sum(x["bytes"] for x in files),"files":files}; (out/"coverage_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8"); return report


def fetch_efloras(out: Path) -> dict:
    landing="http://www.efloras.org/"; shard_index=int(os.environ.get("EFLORAS_SHARD_INDEX","0")); shard_count=max(1,int(os.environ.get("EFLORAS_SHARD_COUNT","1"))); max_pages=max(1,int(os.environ.get("EFLORAS_MAX_PAGES","1000")))
    body,final,_=get(landing); p=Links(); p.feed(decode(body)); content=[]
    for href in p.links:
        u=urllib.parse.urljoin(final,href)
        if urllib.parse.urlparse(u).netloc.endswith("efloras.org") and any(k in u.lower() for k in ("florataxon.aspx","browse.aspx","flora_page.aspx")): content.append(u)
    content=sorted(dict.fromkeys(content))
    if not content: raise RuntimeError("eFloras: no flora/taxon content links discovered")
    seeds=[u for u in content if stable_bucket(u,shard_count)==shard_index] or content[shard_index:shard_index+1]
    out.mkdir(parents=True,exist_ok=True); master=master_names(); matched=set(); focal_matched=set(); saved=[]; candidates=[]; queue=seeds[:100]; seen=set(); treatment_pages=0
    while queue and len(saved)<max_pages:
        u=queue.pop(0)
        if u in seen: continue
        seen.add(u)
        try: b,f,_=get(u)
        except Exception: continue
        if len(b)<1000: continue
        raw=decode(b); text=visible_text(raw); taxa=set(BINOMIAL.findall(text)); page_matches=taxa & master; matched |= page_matches
        focal=focal_species_for_page(text,page_matches)
        is_treatment="florataxon.aspx" in f.lower() and treatment_text_present(text,focal)
        if is_treatment:
            treatment_pages+=1; focal_matched.add(focal)
        name=f"page_{len(saved)+1:05d}.html"; (out/name).write_bytes(b)
        page_candidates=extract_trait_candidates(text,focal if is_treatment else "",f,name); candidates.extend(page_candidates)
        saved.append({"url":f,"file":name,"bytes":len(b),"candidate_taxa":len(taxa),"exact_master_matches":len(page_matches),"focal_species":focal,"is_species_treatment":is_treatment,"trait_candidates":len(page_candidates)})
        q=Links(); q.feed(raw)
        for href in q.links:
            v=urllib.parse.urljoin(f,href)
            if urllib.parse.urlparse(v).netloc.endswith("efloras.org") and "florataxon.aspx" in v.lower() and v not in seen and stable_bucket(v,shard_count)==shard_index and len(queue)<20000: queue.append(v)
    if not saved: raise RuntimeError("eFloras: no substantive flora content pages acquired")
    fields=["accepted_species","trait_name","candidate_value","source_name","source_trait","source_url","source_record_id","source_excerpt","evidence_scope","evidence_status"]
    with (out/"efloras_species_direct_candidates.csv").open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields); writer.writeheader(); writer.writerows(candidates)
    report={"source":"eFloras","real_data_acquired":True,"shard_index":shard_index,"shard_count":shard_count,"pages_saved":len(saved),"species_treatment_pages":treatment_pages,"bytes_saved":sum(x["bytes"] for x in saved),"unique_exact_master_matches_anywhere":len(matched),"unique_focal_species_treatments":len(focal_matched),"species_direct_trait_candidate_rows":len(candidates),"species_with_trait_candidates":len({r['accepted_species'] for r in candidates}),"trait_candidate_cells":len({(r['accepted_species'],r['trait_name']) for r in candidates}),"trait_cells_by_trait":{t:len({(r['accepted_species'],r['trait_name']) for r in candidates if r['trait_name']==t}) for t in sorted({r['trait_name'] for r in candidates})},"pages":saved}
    (out/"coverage_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8"); return report


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True,choices=["wfo-plant-list","efloras"]); args=ap.parse_args()
    if args.source=="efloras": out=Path("data/v2/staging/traits/bulk/efloras")/f"shard-{os.environ.get('EFLORAS_SHARD_INDEX','0')}"; report=fetch_efloras(out)
    else: out=Path("data/v2/staging/traits/bulk")/args.source; report=fetch_wfo_plant_list(out)
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
