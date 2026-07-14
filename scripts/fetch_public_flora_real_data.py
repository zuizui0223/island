#!/usr/bin/env python3
"""Fetch real public data from WFO Plant List or eFloras.

A run is successful only when non-trivial data files/pages are actually written.
Landing-page reachability alone is not success.
"""
from __future__ import annotations
import argparse, csv, json, re, ssl, time, urllib.parse, urllib.request
from html.parser import HTMLParser
from pathlib import Path

UA={"User-Agent":"island-v2-trait-harvester/1.0 (+https://github.com/zuizui0223/island)"}
MASTER=Path("data/v2/staging/gbif/collected/island_taxa.csv")

class Links(HTMLParser):
    def __init__(self): super().__init__(); self.links=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a":
            d=dict(attrs); href=d.get("href")
            if href: self.links.append(href)

def get(url: str, attempts: int=4, timeout: int=180) -> tuple[bytes,str,str]:
    last=None
    for attempt in range(attempts):
        try:
            req=urllib.request.Request(url,headers=UA)
            ctx=ssl.create_default_context()
            with urllib.request.urlopen(req,timeout=timeout,context=ctx) as r:
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

def master_names() -> set[str]:
    with MASTER.open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f)
        return {" ".join((r.get("accepted_species") or "").split()) for r in rd if r.get("accepted_species")}

def fetch_wfo_plant_list(out: Path) -> dict:
    landing="https://wfoplantlist.org/"
    body,final,ctype=get(landing)
    p=Links(); p.feed(decode(body))
    links=[]
    for href in p.links:
        u=urllib.parse.urljoin(final,href)
        low=u.lower()
        if any(x in low for x in ("snapshot","download","archive")) or low.endswith((".zip",".csv",".tsv",".txt",".gz")):
            if u not in links: links.append(u)
    # Follow likely snapshot/download pages one level to discover actual files.
    file_links=[]
    for u in links[:30]:
        if u.lower().endswith((".zip",".csv",".tsv",".txt",".gz")):
            file_links.append(u); continue
        try:
            b,f,c=get(u)
            if "html" not in c.lower() and not b.lstrip().startswith(b"<"):
                file_links.append(u); continue
            q=Links(); q.feed(decode(b))
            for href in q.links:
                v=urllib.parse.urljoin(f,href)
                if v.lower().split("?",1)[0].endswith((".zip",".csv",".tsv",".txt",".gz")) and v not in file_links:
                    file_links.append(v)
        except Exception:
            continue
    if not file_links:
        raise RuntimeError(f"WFO Plant List: no downloadable snapshot/data file discovered from {len(links)} candidate links")
    out.mkdir(parents=True,exist_ok=True)
    saved=[]
    # Download a bounded number of discovered official files; prefer archive/data files.
    for i,u in enumerate(file_links[:8],1):
        b,f,c=get(u)
        name=Path(urllib.parse.urlparse(f).path).name or f"wfo_data_{i}.bin"
        path=out/name; path.write_bytes(b)
        saved.append({"url":f,"file":name,"bytes":len(b),"content_type":c})
    if not saved or sum(x["bytes"] for x in saved)<10000:
        raise RuntimeError("WFO Plant List: discovered files were trivial/empty")
    report={"source":"WFO Plant List","real_data_acquired":True,"files_saved":len(saved),"bytes_saved":sum(x["bytes"] for x in saved),"files":saved}
    (out/"coverage_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    return report

def fetch_efloras(out: Path) -> dict:
    landing="http://www.efloras.org/"
    body,final,ctype=get(landing)
    p=Links(); p.feed(decode(body))
    # Follow actual flora/taxon content links on the official site.
    content=[]
    for href in p.links:
        u=urllib.parse.urljoin(final,href)
        if urllib.parse.urlparse(u).netloc.endswith("efloras.org") and any(k in u.lower() for k in ("florataxon.aspx","browse.aspx","flora_page.aspx")):
            if u not in content: content.append(u)
    if not content:
        raise RuntimeError("eFloras: no flora/taxon content links discovered")
    out.mkdir(parents=True,exist_ok=True)
    master=master_names(); matched=set(); saved=[]
    binomial=re.compile(r"\b([A-Z][a-zA-Z.-]+\s+[a-z][a-zA-Z.-]+)\b")
    # Acquire real content pages. Keep bounded for CI pilot but require actual HTML bodies.
    queue=content[:20]
    seen=set()
    while queue and len(saved)<100:
        u=queue.pop(0)
        if u in seen: continue
        seen.add(u)
        try: b,f,c=get(u)
        except Exception: continue
        text=decode(b)
        if len(b)<1000: continue
        name=f"page_{len(saved)+1:04d}.html"
        (out/name).write_bytes(b)
        taxa=set(binomial.findall(re.sub(r"<[^>]+>"," ",text)))
        matched |= taxa & master
        saved.append({"url":f,"file":name,"bytes":len(b),"candidate_taxa":len(taxa),"exact_master_matches":len(taxa & master)})
        q=Links(); q.feed(text)
        for href in q.links:
            v=urllib.parse.urljoin(f,href)
            if urllib.parse.urlparse(v).netloc.endswith("efloras.org") and "florataxon.aspx" in v.lower() and v not in seen and len(queue)<500:
                queue.append(v)
    if not saved:
        raise RuntimeError("eFloras: no substantive flora content pages acquired")
    report={"source":"eFloras","real_data_acquired":True,"pages_saved":len(saved),"bytes_saved":sum(x["bytes"] for x in saved),"unique_exact_master_matches":len(matched),"pages":saved}
    (out/"coverage_report.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    return report

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--source",required=True,choices=["wfo-plant-list","efloras"]); args=ap.parse_args()
    out=Path("data/v2/staging/traits/bulk")/args.source
    report=fetch_wfo_plant_list(out) if args.source=="wfo-plant-list" else fetch_efloras(out)
    print(json.dumps(report,indent=2))
if __name__=="__main__": main()
