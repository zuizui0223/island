#!/usr/bin/env python3
"""Verify the transport review sources and export only anonymous aggregate tables."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import math
import re

R = Path(__file__).resolve().parent
P = R.parent
ROOT = R.parents[2]
MANUSCRIPT = "MANUSCRIPT_BLINDED_REVIEW_TRANSPORT_20261128.md"
SI = "SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW_TRANSPORT_20261128.md"
TABLES = ("Table_S8f_local_genus_margins.csv", "Table_S8g_genus_transport.csv",
          "Table_S8h_genus_transport_coverage.csv")
PINS = {
    "FLORAL_AXIS_STRATUM_RESULT_20261128.json": ("b6dac06abf334a9edd6e1ed82e0684e1ed018ed9",46105670014,4),
    "FLORAL_GENUS_TRANSPORT_RESULT_20261128.json": ("ecbb76cdcd2f505e693f554b863d1c21c3aed333",46106237400,8),
    "FLORAL_GENUS_TRANSPORT_COVERAGE_RESULT_20261128.json": ("8ae333af1a7a60c65163f79db3814f4e9f360997",46113724276,8),
}
BANNED = re.compile(r"(?i)zuizui0223|zhang ruiqi|@gmail\.com|@tohoku|10\.5281/zenodo\.22704973")

def read(p):
    return p.read_text(encoding="utf-8")

def blob(p):
    b=p.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()

def load(name):
    pin,run,n=PINS[name]
    p=P/name
    assert blob(p)==pin,name
    x=json.loads(read(p))
    assert x["run"]==run and x["owner"]=="island" and len(x["panels"])==n,name
    return x["panels"]

def equal(a,b):
    if isinstance(b,dict):
        for k,v in b.items():equal(a[k],v)
    elif isinstance(b,list):
        assert len(a)==len(b)
        for x,y in zip(a,b):equal(x,y)
    elif isinstance(b,(int,float)) and not isinstance(b,bool):
        assert math.isclose(a,b,abs_tol=1e-10,rel_tol=0),(a,b)
    else:
        assert a==b,(a,b)

def flatten(obj,prefix=""):
    out={}
    for k,v in obj.items():
        if k in ("warning","interval_note"):continue
        key=prefix+str(k)
        if isinstance(v,dict):out.update(flatten(v,key+"__"))
        elif isinstance(v,list):out.update(flatten(dict(enumerate(v)),key+"__"))
        else:out[key]=v
    return out

def write_table(path,rows):
    flat=[flatten(r) for r in rows]
    fields=list(dict.fromkeys(k for r in flat for k in r))
    assert not any(k in ("island_id","accepted_species","genus_name","latitude","longitude") for k in fields)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(flat)
    assert not BANNED.search(read(path))
    with path.open(newline="",encoding="utf-8") as f:
        assert len(list(csv.DictReader(f)))==len(rows)

def check(out):
    m=read(R/MANUSCRIPT);si=read(R/SI)
    base_si=read(R/"SUPPLEMENTARY_INFORMATION_BLINDED_REVIEW.md").rstrip()
    assert si.startswith(base_si+"\n"),"original SI changed"
    old=read(ROOT/"submission/geb_current/MANUSCRIPT_BLINDED.md")
    oldres=old.split("## Results\n",1)[1].split("## Discussion\n",1)[0]
    newres=m.split("## Results\n",1)[1].split("## Discussion\n",1)[0]
    assert all(p in newres for p in oldres.strip().split("\n\n") if p.strip()),"original Results changed"
    a=m.split("## Abstract",1)[1].split("**Keywords:**",1)[0]
    b=m.split("## Introduction",1)[1].split("## References",1)[0]
    assert len(a.split())<=300 and len(b.split())<=5400
    for h in ("Aim","Location","Time period","Major taxa studied","Methods","Results","Main conclusions"):
        assert a.count("**"+h+":**")==1
    k=[x.strip() for x in m.split("**Keywords:**",1)[1].splitlines()[0].split(",")]
    assert 6<=len(k)<=10 and k==sorted(k,key=str.casefold)
    assert "provisional Appendix" not in m and "PENDING" in m and "176,414" in m
    assert not BANNED.search(m) and not BANNED.search(si)
    for marker in ("0.02762","0.01967","0.02463","0.02089","0.000221","0.00397","0.00371","32 global","eight tropical"):
        assert marker in m,marker
    for n in TABLES:assert n in si,n
    for n in ("S8.7","S8.8","S8.9"):assert n in si
    assert "S8.7–S8.9" in m
    g,prior,latest=[load(n) for n in PINS]
    assert sum(x["status"]=="fit" for x in g)==3
    assert sum(x["status"]=="fit" for x in latest)==5
    for x in g:
        if x["status"]=="fit":
            assert x["metrics"]["observed"]["beta"]==x["metrics"]["genus_expected"]["beta"]
            assert x["species_coherent_null"]["informative_strata"]==0
            assert x["species_coherent_null"]["p"] is None
    cover=[]
    for x,y in zip(latest,prior,strict=True):
        equal(x,y)
        row={k:x[k] for k in ("scope","nmin","mode","islands","blocks","status")}
        for k in ("scored_islands","failed_folds"):
            if k in x:row[k]=x[k]
        if x["status"]=="fit":
            c=x["coverage_audit"];s=c["conditional_residual_slopes"]
            assert math.isclose(s["supported_genus"]+s["unsupported_genus"],s["total_residual"],abs_tol=1e-10)
            assert math.isclose(s["total_residual"],x["slopes"]["residual"],abs_tol=1e-10)
            assert math.isclose(sum(c["donor_support_mass_mean"].values()),1,abs_tol=1e-10)
            assert math.isclose(sum(c["category_squared_error_signed_projections"].values()),c["category_squared_error_total"],abs_tol=1e-10)
            row["coverage_audit"]=c
        cover.append(row)
    out.mkdir(parents=True,exist_ok=True)
    for n,rows in zip(TABLES,(g,prior,cover),strict=True):write_table(out/n,rows)
    print(json.dumps({"status":"GEB_TRANSPORT_TEXT_TABLES_PASS","abstract_words":len(a.split()),"body_words":len(b.split()),"new_table_rows":[4,8,8],"original_results_preserved":True,"original_SI_preserved":True,"new_biological_fit":False,"ready_to_submit":False}))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    check(p.parse_args().out)
