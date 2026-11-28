#!/usr/bin/env python3
"""Europe PMC article-level license metadata triage for frozen Chapter 1.

Read only pinned source archives, query the official REST search API using
exact historical publication identifiers, keep no article abstracts or full
texts. A permissive metadata label is a review candidate, NOT a permission
grant. No scientific data or source-policy changes. Public CI emits aggregates
only; an optional local metadata receipt is for an authorized private custodian.
"""
from __future__ import annotations
import argparse
import collections as C
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
import time
from urllib.parse import urlsplit
import zipfile

import pandas as pd
import requests

ROOT="https://www.ebi.ac.uk/europepmc/webservices/rest/search"
FAMILY="domain:europepmc.org"
PINS={
 "ledger":("10042377593-source-scale-batch-integration-34191508045.zip","5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd"),
 "provenance":("10154261753-chapter1-database-v1-rights-audit-34482428564.zip","5d3e5ea80e1cf9757cc7b3bbb64db2ca67f0b8ecce7f7f52668472511f07ad36"),
 "rights":("10186293306-chapter1-database-v1-cell-rights-34566271397.zip","6c5f2577e85e95e841b8659427e27c72ce631f59a3da5a335f54f505de5ea32a"),
}
SRC_COUNTS={"MED":1302,"AGR":537,"PPR":128,"CBA":25,"ETH":21,"PAT":10,"PMC":1}
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def read(z,n):
    return pd.read_csv(z.open(n),compression="gzip" if n.endswith(".gz") else None,dtype=str,keep_default_na=False)
def exact_id(token):
    u=urlsplit(token.removeprefix("url:"))
    m=re.fullmatch(r"/article/([A-Za-z]+)/([A-Za-z0-9_.:-]+)/?",u.path)
    if u.scheme!="https" or u.hostname not in ("europepmc.org","www.europepmc.org") or not m or u.query or u.fragment:
        raise ValueError("unrecognized historical publication URL")
    s,k=m.groups()
    return s.upper(),k.upper()
def normalize_license(v):
    # Strict allowlist, with optional explicit version. Never use substring
    # matching that would mistake CC BY-NC or BY-ND for CC BY.
    s=str(v or "").strip().lower()
    m=re.fullmatch(r"https?://(?:www\.)?creativecommons\.org/(licenses/(by(?:-sa|-nc|-nd|-nc-sa|-nc-nd)?)/([1-4]\.0)(?:/[^?]*)?|publicdomain/zero/1\.0(?:/[^?]*)?)",s)
    if m:
        if "publicdomain/zero/" in s:return "cc0","1.0"
        return "cc-"+m.group(2),m.group(3)
    s=re.sub(r"[ _-]+"," ",s)
    if s in ("cc0","cc 0","cc0 1.0","cc 0 1.0"):return "cc0","1.0" if "1.0" in s else ""
    m=re.fullmatch(r"cc by( nc nd| nc sa| nc| nd| sa)?(?: ([1-4]\.0))?",s)
    if m:return "cc-by"+(m.group(1) or "").replace(" ","-"),m.group(2) or ""
    return "missing" if not s else "other",""
def classify(rec,key):
    if rec is None:return "not_returned"
    src,k=key
    lic,version=normalize_license(rec.get("license"))
    if lic in {"cc-by-nc","cc-by-nc-sa","cc-by-nd","cc-by-nc-nd"}:return "nc_or_nd_metadata"
    if lic=="missing":return "oa_license_missing" if rec.get("isOpenAccess")=="Y" else "no_explicit_license"
    if lic not in {"cc0","cc-by","cc-by-sa"}:return "other_license_review"
    if src not in {"MED","PMC"}:return "version_or_source_review"
    if rec.get("isOpenAccess")!="Y" or not rec.get("pmcid") or not rec.get("doi"):
        return "permissive_label_scope_unverified"
    return "permissive_metadata_candidate"
def exact_match(rec,key):
    s,k=key
    if s=="PMC":return str(rec.get("pmcid","")).upper()==k or (str(rec.get("source","")).upper()==s and str(rec.get("id","")).upper()==k)
    return str(rec.get("source","")).upper()==s and str(rec.get("id","")).upper()==k
def get_metadata(keys):
    records={};failures=0;requests_ok=0;duplicate_conflicts=0;unmatched_hits=0
    sess=requests.Session()
    sess.headers["User-Agent"]="island-chapter1-license-audit/1.0 (metadata-only academic reproducibility)"
    # Group by source; query exact IDs in modest batches. A PMC id may be
    # indexed in a MED record, so query PMCID and check the returned PMCID.
    groups=C.defaultdict(list)
    for key in sorted(keys):groups[key[0]].append(key)
    for s,group in sorted(groups.items()):
        for i in range(0,len(group),35):
            chunk=group[i:i+35]
            clauses=[f'PMCID:"{k}"' if s=="PMC" else f'(SRC:{s} AND EXT_ID:"{k}")' for _,k in chunk]
            params={"query":" OR ".join(clauses),"format":"json","resultType":"core","pageSize":1000,"synonym":"false"}
            response=None
            for attempt in range(3):
                try:
                    r=sess.get(ROOT,params=params,timeout=(12,40))
                    r.raise_for_status()
                    data=r.json()
                    hits=data.get("resultList",{}).get("result",[])
                    if int(data.get("hitCount",-1))!=len(hits):
                        raise ValueError("truncated or malformed search response")
                    response=hits;requests_ok+=1;break
                except (requests.RequestException,ValueError,TypeError):
                    if attempt<2:time.sleep(2**attempt)
            if response is None:
                failures+=1
                for k in chunk:records[k]={"status":"request_failed"}
            else:
                for rec in response:
                    if not any(exact_match(rec,k) for k in chunk):unmatched_hits+=1
                for k in chunk:
                    hits=[rec for rec in response if exact_match(rec,k)]
                    if not hits:records[k]={"status":"not_returned"};continue
                    # Retain only minimal public bibliographic/rights metadata.
                    keep=("source","id","pmcid","doi","license","isOpenAccess","inPMC","firstPublicationDate","pubYear","versionNumber")
                    projected=[{f:rec.get(f) for f in keep} for rec in hits]
                    if len({json.dumps(rec,sort_keys=True) for rec in projected})!=1:
                        duplicate_conflicts+=1;records[k]={"status":"conflicting_metadata"};continue
                    rec=projected[0]
                    records[k]={"status":classify(rec,k),"record":rec}
            time.sleep(.35)
    assert set(records)==set(keys)
    return records,{"successful_batches":requests_ok,"failed_batches":failures,"conflicting_matches":duplicate_conflicts,"unexpected_returned_records":unmatched_hits}
def load_inputs(folder):
    z={}
    for k,(name,h) in PINS.items():
        p=folder/name
        if sha(p)!=h:raise ValueError("archive digest mismatch: "+k)
        z[k]=zipfile.ZipFile(p)
    l=read(z["provenance"],"SOURCE_LINEAGE_DETAILS.csv")
    tokens=set(l.loc[l.source_family.eq(FAMILY),"source_lineage"])
    keyed={t:exact_id(t) for t in tokens}
    if len(tokens)!=2024 or len(set(keyed.values()))!=2024 or dict(C.Counter(k[0] for k in keyed.values()))!=SRC_COUNTS:
        raise ValueError("historical identifier roster drift")
    ledger=read(z["ledger"],"species_axis_coverage.csv.gz")
    rights=read(z["rights"],"CELL_RELEASE_RIGHTS.csv.gz")
    recovery=read(z["provenance"],"VALIDATED_LOW_PROVENANCE_RECOVERY.csv")
    cols=["accepted_species","axis"]
    if any(x.duplicated(cols).any() for x in [ledger,rights,recovery]):raise ValueError("duplicate cell keys")
    r=rights[rights.blocked_source_families.map(lambda x:FAMILY in x.split("|"))].copy()
    x=r.merge(ledger[cols+["source_lineages"]],on=cols,validate="one_to_one").merge(recovery[cols+["support_source_lineages"]],on=cols,how="left",validate="one_to_one").fillna("")
    assert len(x)==18826 and x.release_status.eq("rights_review_required").all()
    cells=[]
    for row in x.itertuples():
        direct=set(row.source_lineages.split("|"))
        upstream=set(json.loads(row.support_source_lineages)) if row.support_source_lineages else set()
        found=(direct|upstream)&tokens
        unknown={t for t in (direct|upstream) if "europepmc.org" in t and t not in tokens}
        if not found or unknown:raise ValueError("cell has unresolved Europe PMC lineage scope")
        cells.append((frozenset(keyed[t] for t in found),row.axis,row.quality,row.blocked_source_families==FAMILY,row.has_validated_low_lineage=="True"))
    assert sum(c[3] for c in cells)==883
    assert sum(c[3] and c[1]=="reproductive_assurance" for c in cells)==781
    assert sum(c[4] for c in cells)==8431
    return set(keyed.values()),cells
def analyze(keys,cells,records,request_summary):
    stat=C.Counter(v["status"] for v in records.values())
    by_source={s:dict(C.Counter(records[k]["status"] for k in keys if k[0]==s)) for s in sorted(SRC_COUNTS)}
    candidates={k for k,v in records.items() if v["status"]=="permissive_metadata_candidate"}
    touched=[c for c in cells if c[0]&candidates]
    all_epmc=[c for c in cells if c[0]<=candidates]
    sole=[c for c in all_epmc if c[3]]
    def summarize(rows):
        return {"cells":len(rows),"by_axis":dict(C.Counter(c[1] for c in rows)),"by_quality":dict(C.Counter(c[2] for c in rows)),"derived_cells":sum(c[4] for c in rows)}
    return {
      "timestamp_utc":datetime.now(timezone.utc).isoformat(),
      "input_archive_sha256":{k:h for k,(_,h) in PINS.items()},
      "historical_article_ids":len(keys),"historical_ids_by_source":SRC_COUNTS,
      "affected_cells":summarize(cells),"sole_family_blocked_cells":summarize([c for c in cells if c[3]]),
      "records_by_metadata_status":dict(stat),"source_by_status":by_source,
      "candidate_normalized_license_labels":dict(C.Counter(normalize_license(records[k]["record"].get("license"))[0] for k in candidates)),
      "candidate_explicit_version_count":sum(bool(normalize_license(records[k]["record"].get("license"))[1]) for k in candidates),
      "candidate_touched_cells":summarize(touched),
      "all_required_epmc_ids_are_metadata_candidates":summarize(all_epmc),
      "hypothetical_sole_family_cell_ceiling":summarize(sole),
      "request_summary":request_summary,
      "policy_changed":False,"new_cells_authorized":0,"public_release":False,
      "interpretation":"Current article metadata only. Exact historical source ID match is not proof of the exact acquisition version or applicable source-data rights. CC-BY/CC0/CC-BY-SA candidates require version, license statement, third-party content, attribution and actual evidence scope review. For derived cells every upstream lineage is required; other original source-family blockers and semantic gates remain in force. No permission is granted by this audit.",
    }
def self_test():
    assert normalize_license("CC BY-NC-SA")==("cc-by-nc-sa","")
    assert normalize_license("https://creativecommons.org/licenses/by/4.0/")==("cc-by","4.0")
    assert normalize_license("CC BY-ND 4.0")[0]=="cc-by-nd"
    assert exact_id("url:https://europepmc.org/article/MED/123")==("MED","123")
    r={"source":"MED","id":"123","pmcid":"PMC456","doi":"10.1/x","license":"cc by","isOpenAccess":"Y"}
    assert exact_match(r,("MED","123")) and not exact_match(r,("AGR","123"))
    assert classify(r,("MED","123"))=="permissive_metadata_candidate"
    assert classify(dict(r,license="cc by nc"),("MED","123"))=="nc_or_nd_metadata"
    assert classify(dict(r,isOpenAccess="N"),("MED","123"))=="permissive_label_scope_unverified"
    print("EPMC_METADATA_SELF_TEST_PASS")
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source-dir",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--private-metadata",type=Path,help="Optional local custodian receipt. Never enable in public CI.")
    p.add_argument("--roster-only",action="store_true")
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:self_test();raise SystemExit
    keys,cells=load_inputs(a.source_dir)
    if a.roster_only:
        print("EPMC_ROSTER_PASS "+json.dumps({"ids":len(keys),"source_counts":SRC_COUNTS,"affected_cells":len(cells),"sole_blocked":sum(c[3] for c in cells),"derived":sum(c[4] for c in cells)}))
        raise SystemExit
    records,req=get_metadata(keys)
    result=analyze(keys,cells,records,req)
    a.out.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
    if a.private_metadata:
        a.private_metadata.write_text(json.dumps([{"historical_source":s,"historical_id":k,**r} for (s,k),r in sorted(records.items())],indent=2,sort_keys=True)+"\n")
    print("EPMC_LICENSE_METADATA_RESULT "+json.dumps(result,sort_keys=True,allow_nan=False),flush=True)
    if not req["successful_batches"]:raise SystemExit("No metadata batches succeeded; no article licenses verified")
