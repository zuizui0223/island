#!/usr/bin/env python3
"""Probe publicly released Roeble et al. (2024) Supplementary Data 3 geological origin.

Only save header/schema and first few example records for verification.
Do not publish a copy of the source XLSX in this repository or artifact.
"""
import argparse
import csv
import io
import json
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

BASENAME="41467_2024_51556_MOESM6_ESM.xlsx"
URLS=[
 "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-024-51556-7/MediaObjects/"+BASENAME,
 "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-51556-7/MediaObjects/"+BASENAME,
 "https://pmc.ncbi.nlm.nih.gov/articles/PMC11343744/bin/"+BASENAME,
 "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11343744/supplementaryFiles",
]
NS={"x":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

def get_sheets(blob):
    with zipfile.ZipFile(io.BytesIO(blob)) as archive:
        names=archive.namelist()
        if "xl/workbook.xml" not in names:
            return {"is_zip":True,"files":names[:30],"sheets":[]}
        strings=[]
        if "xl/sharedStrings.xml" in names:
            root=ET.fromstring(archive.read("xl/sharedStrings.xml"))
            strings=["".join(t.text or "" for t in si.findall(".//x:t",NS))
                     for si in root.findall("x:si",NS)]
        sheets=[]
        for name in names:
            if not name.startswith("xl/worksheets/sheet") or not name.endswith(".xml"):
                continue
            root=ET.fromstring(archive.read(name))
            records=[]
            for row in root.findall(".//x:sheetData/x:row",NS)[:8]:
                cells={}
                for cell in row.findall("x:c",NS):
                    col="".join(ch for ch in cell.attrib.get("r","") if ch.isalpha())
                    index=0
                    for char in col:
                        index=index*26+(ord(char)-64)
                    value=cell.find("x:v",NS)
                    inline=cell.find("x:is",NS)
                    v=value.text if value is not None else None
                    if cell.attrib.get("t")=="s" and v is not None:
                        v=strings[int(v)]
                    if inline is not None:
                        v="".join(t.text or "" for t in inline.findall(".//x:t",NS))
                    cells[index-1]=v
                records.append([cells.get(i) for i in range(min(36,max(cells,default=-1)+1))])
            sheets.append({"path":name,"rows":records})
        return {"is_zip":True,"files":len(names),"shared_strings":len(strings),"sheets":sheets}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    report={"source":"Roeble et al 2024 Supplementary Data 3; 510 islands; 272 oceanic 238 continental",
            "downloads":[]}
    for url in URLS:
        try:
            req=urllib.request.Request(url,headers={"User-Agent":"Mozilla/5.0 (compatible; fd-research/1.0)","Accept":"*/*"})
            with urllib.request.urlopen(req,timeout=50) as resp:
                blob=resp.read(24_000_000)
                if len(blob)>=24_000_000:
                    raise ValueError("Exceeded download size cap")
                if not blob.startswith(b"PK"):
                    raise ValueError("Not a ZIP/XLSX byte stream")
                if url.endswith("/supplementaryFiles"):
                    with zipfile.ZipFile(io.BytesIO(blob)) as z:
                        matches=[n for n in z.namelist() if BASENAME.lower() in n.lower()]
                        if not matches:
                            matches=[n for n in z.namelist() if n.lower().endswith(".xlsx")]
                        if not matches:raise ValueError("Supplementary ZIP has no xlsx")
                        report["epmc_entries"]=z.namelist()[:40]
                        blob=z.read(matches[0])
                schema=get_sheets(blob)
                report["downloads"].append({"url":url,"status":"ok","bytes":len(blob),"schema":schema})
                break
        except Exception as exc:
            report["downloads"].append({"url":url,"status":"error","error":f"{type(exc).__name__}: {exc}"[:300]})
    (args.output/"roeble_supplement_schema.json").write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(report,ensure_ascii=False)[:13000])

if __name__=="__main__":main()
