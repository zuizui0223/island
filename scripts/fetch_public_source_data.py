#!/usr/bin/env python3
"""Fetch real public data files from verified official plant/flora sources.

A source is successful only when at least one non-HTML machine-readable payload is actually
saved. Landing-page reachability alone is not success.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import mimetypes
import os
import re
import subprocess
import sys
import urllib.parse
from html.parser import HTMLParser
from pathlib import Path

SOURCES = {
    "gift": {
        "name": "GIFT",
        "start": "https://gift.uni-goettingen.de/home",
        "hosts": {"gift.uni-goettingen.de"},
        "max_pages": 80,
    },
    "wfo": {
        "name": "World Flora Online (taxonomic backbone download)",
        # The WFO homepage is a JS app with no bulk file. Its Darwin Core
        # backbone archives are published on the download page and served from
        # the files.* host, so start there and allow that host.
        "start": "https://www.worldfloraonline.org/downloadData",
        "hosts": {"www.worldfloraonline.org", "worldfloraonline.org", "files.worldfloraonline.org"},
        "max_pages": 60,
    },
    "powo": {
        "name": "Plants of the World Online (WCVP bulk backbone)",
        # POWO's website exposes no bulk payload to a crawler. Its authoritative
        # machine-readable checklist is published as the World Checklist of
        # Vascular Plants (WCVP) in Kew's data repository, so pull the real bulk
        # archive from there instead of scraping the POWO web app.
        "start": "https://sftp.kew.org/pub/data-repositories/WCVP/",
        "hosts": {"sftp.kew.org"},
        "max_pages": 40,
    },
    "try": {
        "name": "TRY Plant Trait Database (public datasets)",
        # The try-db.org root exposes no file links at depth 1, so the earlier
        # crawl found nothing. Start on the public Data page, which links the
        # openly released datasets (e.g. the Categorical Traits Dataset), and
        # allow the whole try-db.org host so the crawler can follow to the files.
        "start": "https://www.try-db.org/TryWeb/Data.php",
        "hosts": {"www.try-db.org", "try-db.org"},
        "max_pages": 120,
    },
    "sid": {
        "name": "Kew Seed Information Database",
        "start": "https://ser-sid.org/",
        "hosts": {"ser-sid.org", "www.ser-sid.org"},
        "max_pages": 80,
    },
}

DATA_EXT = re.compile(r"\.(?:csv|tsv|txt|json|xml|zip|gz|tgz|xlsx?|parquet|dwca)(?:$|[?#])", re.I)
DATA_HINT = re.compile(r"(?:download|export|dataset|archive|snapshot|bulk|api|data)", re.I)
HTML_TYPES = ("text/html", "application/xhtml+xml")
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/131 Safari/537.36 island-v2/1.0"

class LinkParser(HTMLParser):
    def __init__(self):
        super().__init__(); self.links=[]
    def handle_starttag(self, tag, attrs):
        d=dict(attrs)
        if tag.lower() in {"a","link","script"}:
            href=d.get("href") or d.get("src")
            if href: self.links.append(href)


def curl(url: str, head: bool=False) -> tuple[bytes, str, str]:
    cmd=["curl","--fail","--location","--retry","4","--retry-delay","2","--connect-timeout","30","--max-time","180","-A",UA,"-sS"]
    if head:
        cmd += ["-I"]
    cmd += ["-D","-",url]
    p=subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    raw=p.stdout
    # Split final header block from body. curl -D - may emit multiple redirect blocks.
    marker=b"\r\n\r\n"
    parts=raw.split(marker)
    body=parts[-1] if not head else b""
    headers=parts[-2].decode("latin-1","replace") if len(parts)>=2 else ""
    ctype=""
    final=url
    for line in headers.splitlines():
        if line.lower().startswith("content-type:"):
            ctype=line.split(":",1)[1].strip().lower()
    return body, ctype, final


def is_html(ctype: str, body: bytes) -> bool:
    if any(x in ctype for x in HTML_TYPES): return True
    s=body[:500].lstrip().lower()
    return s.startswith(b"<!doctype html") or s.startswith(b"<html")


def safe_name(url: str, ctype: str) -> str:
    path=urllib.parse.urlparse(url).path
    name=Path(path).name
    if not name or "." not in name:
        ext=mimetypes.guess_extension(ctype.split(";",1)[0].strip()) or ".bin"
        name=hashlib.sha256(url.encode()).hexdigest()[:16]+ext
    return re.sub(r"[^A-Za-z0-9._-]+","_",name)[:180]


def main() -> int:
    ap=argparse.ArgumentParser(); ap.add_argument("--source", required=True, choices=sorted(SOURCES)); args=ap.parse_args()
    cfg=SOURCES[args.source]
    out=Path("data/v2/staging/traits/bulk/public_sources")/args.source
    raw_dir=out/"raw"; raw_dir.mkdir(parents=True, exist_ok=True)
    queue=[cfg["start"]]; seen=set(); candidates=[]; pages=[]; saved=[]; errors=[]

    while queue and len(pages) < cfg["max_pages"]:
        url=queue.pop(0)
        if url in seen: continue
        seen.add(url)
        try:
            body, ctype, _=curl(url)
        except Exception as exc:
            errors.append({"url":url,"error":str(exc)})
            continue
        if not is_html(ctype, body):
            candidates.append(url)
            continue
        pages.append(url)
        p=LinkParser(); p.feed(body.decode("utf-8","replace"))
        for href in p.links:
            absolute=urllib.parse.urljoin(url,href)
            parsed=urllib.parse.urlparse(absolute)
            if parsed.scheme not in {"http","https"}: continue
            if DATA_EXT.search(absolute) or DATA_HINT.search(absolute):
                candidates.append(absolute)
            if parsed.hostname in cfg["hosts"] and absolute not in seen and len(queue) < cfg["max_pages"]*4:
                queue.append(absolute)

    for url in dict.fromkeys(candidates):
        try:
            body, ctype, _=curl(url)
            if is_html(ctype, body):
                continue
            if len(body) < 20:
                continue
            name=safe_name(url,ctype)
            target=raw_dir/name
            if target.exists():
                target=raw_dir/(hashlib.sha256(url.encode()).hexdigest()[:10]+"_"+name)
            target.write_bytes(body)
            saved.append({"url":url,"file":str(target),"bytes":len(body),"content_type":ctype})
        except Exception as exc:
            errors.append({"url":url,"error":str(exc)})

    report={
        "source":cfg["name"],
        "source_key":args.source,
        "official_start":cfg["start"],
        "pages_crawled":len(pages),
        "candidate_machine_links":len(set(candidates)),
        "files_downloaded":len(saved),
        "downloaded_bytes":sum(x["bytes"] for x in saved),
        "real_data_acquired":bool(saved),
        "files":saved,
        "errors":errors[:100],
    }
    (out/"coverage_report.json").write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding="utf-8")
    print(json.dumps(report,indent=2,ensure_ascii=False))
    if not saved:
        print(f"No public machine-readable payload acquired for {args.source}", file=sys.stderr)
        return 3
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
