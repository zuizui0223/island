#!/usr/bin/env python3
"""Audit verified official plant trait/flora sources without guessing package names.

The registry contains only official landing pages verified manually. Each source can be run
independently so GitHub Actions can probe sources in parallel and preserve source-specific
artifacts even when another source fails.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.error
import urllib.parse
import urllib.request
from html.parser import HTMLParser
from pathlib import Path

BASE_OUT = Path("data/v2/staging/traits/bulk/verified_source_audit")

SOURCES = [
    {"key":"leda","source":"LEDA Traitbase","kind":"trait_db","url":"https://uol.de/en/landeco/research/leda","access":"public_download"},
    {"key":"gift","source":"GIFT","kind":"flora_trait_db","url":"https://gift.uni-goettingen.de/home","access":"public_web"},
    {"key":"wfo-plant-list","source":"World Flora Online Plant List","kind":"global_flora_taxonomy","url":"https://wfoplantlist.org/","access":"public_snapshots"},
    {"key":"wfo","source":"World Flora Online","kind":"global_flora_descriptions","url":"https://www.worldfloraonline.org/","access":"public_web"},
    {"key":"powo","source":"Plants of the World Online","kind":"global_flora_descriptions","url":"https://powo.science.kew.org/","access":"public_web"},
    {"key":"efloras","source":"eFloras","kind":"regional_floras_descriptions","url":"http://www.efloras.org/","access":"public_web"},
    {"key":"bien","source":"BIEN","kind":"plant_distribution_traits","url":"https://bien.nceas.ucsb.edu/bien/","access":"public_web"},
    {"key":"try","source":"TRY Plant Trait Database","kind":"trait_db","url":"https://www.try-db.org/","access":"request_or_terms"},
    {"key":"sid","source":"Kew Seed Information Database","kind":"seed_traits","url":"https://ser-sid.org/","access":"public_web"},
]

BULK_RE = re.compile(r"(?:download|data|dataset|archive|snapshot|api|\.csv(?:$|\?)|\.tsv(?:$|\?)|\.zip(?:$|\?)|\.xlsx?(?:$|\?)|\.json(?:$|\?)|\.xml(?:$|\?))", re.I)

class Links(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a": return
        d=dict(attrs); href=d.get("href")
        if href: self.links.append((href, d.get("download")))

def fetch(url: str):
    req=urllib.request.Request(url, headers={"User-Agent":"island-v2-source-audit/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        body=r.read(5_000_000)
        return r.status, r.geturl(), r.headers.get("content-type", ""), body

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source", action="append", dest="sources", help="Registry key to audit; repeatable. Default: all")
    args=parser.parse_args()
    selected=set(args.sources or [])
    known={s["key"] for s in SOURCES}
    unknown=selected-known
    if unknown:
        raise SystemExit("Unknown source key(s): "+", ".join(sorted(unknown)))
    targets=[s for s in SOURCES if not selected or s["key"] in selected]
    if not targets:
        raise SystemExit("No sources selected")

    rows=[]; discovered=[]
    for src in targets:
        out_dir=BASE_OUT/src["key"]
        out_dir.mkdir(parents=True, exist_ok=True)
        row={**src,"ok":False,"http_status":"","final_url":"","content_type":"","error":"","discovered_links":0}
        try:
            status, final_url, ctype, body=fetch(src["url"])
            row.update(ok=True,http_status=status,final_url=final_url,content_type=ctype)
            if "html" in ctype.lower() or body.lstrip().startswith(b"<"):
                p=Links(); p.feed(body.decode("utf-8","replace"))
                seen=set()
                for href, download_attr in p.links:
                    absolute=urllib.parse.urljoin(final_url, href)
                    if absolute in seen: continue
                    if download_attr is not None or BULK_RE.search(absolute):
                        seen.add(absolute)
                        discovered.append({"key":src["key"],"source":src["source"],"landing_url":final_url,"discovered_url":absolute})
                row["discovered_links"]=len(seen)
        except urllib.error.HTTPError as e:
            row.update(http_status=e.code,error=f"HTTPError: {e}")
        except Exception as e:
            row["error"]=f"{type(e).__name__}: {e}"
        rows.append(row)
        (out_dir/"source_status.json").write_text(json.dumps(row,indent=2,ensure_ascii=False),encoding="utf-8")
        with (out_dir/"discovered_machine_links.csv").open("w",newline="",encoding="utf-8") as f:
            fields=["key","source","landing_url","discovered_url"]
            w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows([x for x in discovered if x["key"]==src["key"]])

    BASE_OUT.mkdir(parents=True, exist_ok=True)
    with (BASE_OUT/"source_status.csv").open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
    with (BASE_OUT/"discovered_machine_links.csv").open("w",newline="",encoding="utf-8") as f:
        fields=["key","source","landing_url","discovered_url"]
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(discovered)
    (BASE_OUT/"source_status.json").write_text(json.dumps(rows,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(rows,indent=2,ensure_ascii=False))
    print(f"discovered_machine_links={len(discovered)}")

if __name__ == "__main__": main()
