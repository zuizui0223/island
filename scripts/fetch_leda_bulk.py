#!/usr/bin/env python3
"""Download real LEDA Traitbase TXT files from the official data-files page and audit master coverage."""
from __future__ import annotations
import csv, json, re, time, urllib.parse, urllib.request
from html.parser import HTMLParser
from pathlib import Path

LANDING="https://uol.de/en/landeco/research/leda/data-files"
MASTER=Path("data/v2/staging/gbif/collected/island_taxa.csv")
OUT=Path("data/v2/staging/traits/bulk/leda")
UA={"User-Agent":"island-v2-trait-harvester/1.0 (+https://github.com/zuizui0223/island)"}
BINOMIAL=re.compile(r"^[A-Z][A-Za-z.-]+\s+[a-z][A-Za-z.-]+(?:\s+(?:subsp\.|ssp\.|var\.)\s+[a-z][A-Za-z.-]+)?$")

class Anchors(HTMLParser):
    def __init__(self): super().__init__(); self.items=[]; self.href=None; self.text=[]
    def handle_starttag(self,tag,attrs):
        if tag.lower()=="a": self.href=dict(attrs).get("href"); self.text=[]
    def handle_data(self,data):
        if self.href is not None: self.text.append(data)
    def handle_endtag(self,tag):
        if tag.lower()=="a" and self.href is not None:
            self.items.append((" ".join("".join(self.text).split()),self.href)); self.href=None; self.text=[]

def get(url,limit=None,attempts=5):
    last=None
    for attempt in range(attempts):
        try:
            req=urllib.request.Request(url,headers=UA)
            with urllib.request.urlopen(req,timeout=180) as r:
                return r.read() if limit is None else r.read(limit)
        except Exception as exc:
            last=exc
            if attempt+1<attempts: time.sleep(2**attempt)
    raise RuntimeError(f"GET failed after {attempts} attempts: {url}: {last}")

def norm_name(s): return " ".join((s or "").replace("×","x").split()).strip()

def master_names():
    with MASTER.open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f); return {norm_name(r.get("accepted_species","")) for r in rd if r.get("accepted_species")}

def decode(blob):
    for enc in ("utf-8-sig","cp1252","latin-1"):
        try: return blob.decode(enc)
        except UnicodeDecodeError: pass
    return blob.decode("utf-8","replace")

def candidate_names(text):
    found=set()
    for line in text.splitlines():
        if not line.strip(): continue
        cells=re.split(r"\t|;",line)
        for cell in cells[:8]:
            x=norm_name(cell.strip(' "\''))
            if BINOMIAL.match(x): found.add(x)
    return found

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    html=decode(get(LANDING))
    p=Anchors(); p.feed(html)
    links=[]; seen=set()
    for label,href in p.items:
        url=urllib.parse.urljoin(LANDING,href)
        low=url.lower().split("?",1)[0]
        if label.upper()=="TXT" or low.endswith((".txt",".csv",".tsv")):
            if url not in seen: seen.add(url); links.append(url)
    if not links: raise RuntimeError("LEDA official page exposed no TXT/CSV/TSV data links")
    master=master_names(); reports=[]; union=set(); failures=[]
    for i,url in enumerate(links,1):
        name=Path(urllib.parse.urlparse(url).path).name or f"leda_{i}.txt"
        try:
            blob=get(url); (OUT/name).write_bytes(blob)
            text=decode(blob); taxa=candidate_names(text); matches=taxa & master; union |= matches
            reports.append({"file":name,"url":url,"bytes":len(blob),"candidate_taxa":len(taxa),"exact_master_matches":len(matches),"error":""})
            print(f"LEDA {i}/{len(links)} {name}: bytes={len(blob)} taxa={len(taxa)} matches={len(matches)}")
        except Exception as exc:
            failures.append({"file":name,"url":url,"error":str(exc)})
            reports.append({"file":name,"url":url,"bytes":0,"candidate_taxa":0,"exact_master_matches":0,"error":str(exc)})
            print(f"LEDA {i}/{len(links)} {name}: FAILED: {exc}")
    successful=sum(1 for x in reports if x["bytes"]>0)
    if successful==0: raise RuntimeError(f"LEDA: all {len(links)} downloads failed")
    report={"source":"LEDA Traitbase","official_landing":LANDING,"real_data_acquired":True,"files_discovered":len(links),"files_downloaded":successful,"files_failed":len(failures),"downloaded_bytes":sum(x["bytes"] for x in reports),"unique_exact_master_matches":len(union),"failures":failures,"files":reports}
    (OUT/"coverage_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    with (OUT/"file_coverage.csv").open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=reports[0].keys()); w.writeheader(); w.writerows(reports)
    print(json.dumps(report,indent=2,ensure_ascii=False))
if __name__=="__main__": main()
