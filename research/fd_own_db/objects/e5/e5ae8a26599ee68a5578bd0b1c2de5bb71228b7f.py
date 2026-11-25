#!/usr/bin/env python3
"""Inspect public geography and island-geological-origin sources.

This is a source audit, not a classifier. It does not publish any complete
third-party input archive or assert geological origin from coastline distance.
"""
from __future__ import annotations
import argparse
import csv
import io
import json
import re
import tempfile
import urllib.request
import zipfile
from pathlib import Path

ISLAND_COMMIT = "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c"
GEO_URL = (f"https://raw.githubusercontent.com/zuizui0223/island/{ISLAND_COMMIT}/"
           "results/geography_20260924/corrected_geography_covariates.csv")
ZENODO_URL = "https://zenodo.org/records/6325640/files/input.zip?download=1"
WEIGELT_2013_URL_UNUSED = "https://zenodo.org/records/4971536/files/Weigelt_etal_2013_PNAS_islanddata.csv?download=1"


def download(url: str, path: Path, max_bytes: int) -> int:
    req = urllib.request.Request(url, headers={"User-Agent": "fd-origin-metadata-audit/1.0"})
    size = 0
    with urllib.request.urlopen(req, timeout=120) as src, path.open("wb") as out:
        while chunk := src.read(1024 * 1024):
            size += len(chunk)
            if size > max_bytes:
                raise ValueError(f"Public source exceeds limit {max_bytes}: {url}")
            out.write(chunk)
    return size


def csv_sample(content: str) -> dict:
    rows = csv.DictReader(io.StringIO(content))
    entries = []
    for idx, row in enumerate(rows):
        if idx >= 3: break
        entries.append(row)
    return {"columns": rows.fieldnames, "sample_rows": entries}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = {"source_commit": ISLAND_COMMIT, "geography": {}, "zenodo": {}, "warnings": []}
    with tempfile.TemporaryDirectory() as d:
        directory = Path(d)
        geo = directory / "geo.csv"
        report["geography"]["bytes"] = download(GEO_URL, geo, 60_000_000)
        with geo.open(encoding="utf-8") as handle:
            report["geography"].update(csv_sample("".join(handle.readline() for _ in range(5))))
        with geo.open(encoding="utf-8") as handle:
            report["geography"]["n_rows"] = max(0, sum(1 for _ in handle) - 1)
        origin = directory / "input.zip"
        try:
            report["zenodo"]["bytes"] = download(ZENODO_URL, origin, 110_000_000)
            with zipfile.ZipFile(origin) as z:
                names = z.namelist()
                report["zenodo"]["archive_entries"] = len(names)
                report["zenodo"]["matching_paths"] = [
                    {"filename": n, "bytes": z.getinfo(n).file_size}
                    for n in names if re.search("island|geo|model|predict|environment|origin", n, re.I)
                ][:160]
                report["zenodo"]["non_raster_sources"] = [
                    {"name": n, "size": z.getinfo(n).file_size}
                    for n in names if not n.endswith("/")
                    and not n.startswith("input/PHYLACINE")
                    and n.lower().endswith((".csv", ".rds", ".rda", ".rdata", ".dbf", ".shp", ".gpkg", ".xlsx"))
                ][:160]
                dbf="input/geoentities_simple.dbf"
                if dbf in names:
                    raw=z.read(dbf)
                    header_len=int.from_bytes(raw[8:10],"little")
                    row_len=int.from_bytes(raw[10:12],"little")
                    count=int.from_bytes(raw[4:8],"little")
                    fields=[]
                    for pos in range(32, header_len-1,32):
                        block=raw[pos:pos+32]
                        if block[0]==13:break
                        fields.append({"name":block[:11].split(bytes([0]),1)[0].decode("ascii",errors="replace"),"width":block[16]})
                    report["zenodo"]["geoentities_dbf"]={"records":count,"record_len":row_len,"fields":fields}
                candidates=[]
                for n in names:
                    if n.endswith("/") or not re.search(r"\.(csv|tsv|txt)$", n, re.I):
                        continue
                    info=z.getinfo(n)
                    if info.file_size>8_000_000:continue
                    with z.open(n) as source:
                        excerpt=source.read(min(2800,info.file_size)).decode("utf-8-sig",errors="replace")
                    lines=excerpt.splitlines()[:3]
                    candidates.append({"name":n,"size":info.file_size,"first_lines":lines})
                report["zenodo"]["csv_candidates"]=[x for x in candidates if not x["name"].startswith("input/PHYLACINE")][:100]
        except Exception as exc:
            report["warnings"].append(f"zenodo: {type(exc).__name__}: {exc}")
    (args.output/"origin_source_audit.json").write_text(json.dumps(report,indent=2,ensure_ascii=False))
    print("GEOGRAPHY",json.dumps(report["geography"],ensure_ascii=False)[:4500])
    print("NON_RASTER",json.dumps(report["zenodo"].get("non_raster_sources",[]),ensure_ascii=False)[:15000])
    print("DBF",json.dumps(report["zenodo"].get("geoentities_dbf",{}),ensure_ascii=False)[:3500])
    print("CSV",json.dumps(report["zenodo"].get("csv_candidates",[]),ensure_ascii=False)[:10000])
    print("WARNINGS",report["warnings"])


if __name__ == "__main__":
    main()
