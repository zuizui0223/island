#!/usr/bin/env python3
"""Retrieve publicly distributed island geology reference tables for provenance audit.

Files are *not* committed; keep them only as ephemeral GitHub Actions artifacts.
No geological classification is inferred from them until their fields are checked.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import urllib.request

FILES = {
    "pnas_island_attributes_05.xlsx": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9478640/bin/pnas.2208629119.sd05.xlsx",
    "pnas_island_attributes_06.xlsx": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9478640/bin/pnas.2208629119.sd06.xlsx",
    "pnas_supplement_01.xlsx": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9478640/bin/pnas.2208629119.sd01.xlsx",
    "weigelt_2013_islanddata.csv": "https://zenodo.org/records/4971536/files/Weigelt_etal_2013_PNAS_islanddata.csv?download=1",
}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    status={}
    for name,url in FILES.items():
        path=args.output/name
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"fd-origin-research/1.0"})
            with urllib.request.urlopen(req,timeout=75) as src,path.open("wb") as dst:
                size=0
                while part:=src.read(1024*1024):
                    size+=len(part)
                    if size>20_000_000:raise ValueError("Input >20MB")
                    dst.write(part)
            magic=path.read_bytes()[:4]
            if name.endswith(".xlsx") and magic!=b"PK\x03\x04":
                raise ValueError(f"Not an XLSX: {magic!r}")
            if name.endswith(".csv") and b"," not in path.open("rb").readline():
                raise ValueError("Not CSV")
            status[name]={"status":"ok","size":size,"source_url":url}
        except Exception as e:
            path.unlink(missing_ok=True)
            status[name]={"status":"unavailable","detail":str(e),"source_url":url}
    (args.output/"download_audit.json").write_text(json.dumps(status,indent=2)+"\n")
    print(json.dumps(status))
if __name__=="__main__":main()
