#!/usr/bin/env python3
"""Prioritize exact cell-rights blockers, without granting rights or altering data.

Only aggregate counts and preselected provider names leave the temporary runner.
A hypothetical clearance releases a cell only if *every* frozen blocker is cleared.
Neither an unresolved provenance token nor a hypothetical clearance is permission.
"""
from __future__ import annotations
import argparse, collections as c, csv, gzip, hashlib, io, json, zipfile
from pathlib import Path
from urllib.parse import urlsplit

RIGHTS_SHA = "6c5f2577e85e95e841b8659427e27c72ce631f59a3da5a335f54f505de5ea32a"
PROVENANCE_SHA = "5d3e5ea80e1cf9757cc7b3bbb64db2ca67f0b8ecce7f7f52668472511f07ad36"
AXES = ("floral_structural_complexity", "flower_colour", "reproductive_assurance")
QUALITIES = ("high", "medium", "low")
SOURCES = (
    "domain:efloras.org", "domain:gobotany.nativeplanttrust.org",
    "provider_treatment:flora_malesiana", "database:pladias",
    "ken-fern-useful-plants:species", "domain:gbif.org",
    "provider_treatment:flora_of_china", "database:baseflor",
    "provider:efloras", "provider_treatment:flora_of_panama",
    "database:floraweb_biolflor", "domain:nzpcn.org.nz", "grassbase:url",
    "provider_treatment:nybg_floraneotropica", "provider_treatment:nybg_memoirs",
    "unresolved:other_doi", "unresolved:citation_hash", "domain:europepmc.org",
    "domain:plantnet.rbgsyd.nsw.gov.au", "domain:nparks.gov.sg",
    "domain:plants.ces.ncsu.edu", "database:ecoflora",
    "wave42-secondary-reproductive-databases:shared-redistribution-guard",
)
def label(s):
    # Explicit allowlist prevents arbitrary lineage tokens, species names and URLs
    # from entering the anonymous public receipt.
    return s if s in SOURCES else "withheld_family_"+hashlib.sha256(s.encode()).hexdigest()[:16]
def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest()
def check(p,h):
    if sha(p)!=h:raise ValueError("source ZIP SHA mismatch")
def counts(rows):
    sets=c.Counter(); byaxis=c.Counter(); byquality=c.Counter(); status=c.Counter(); incidence=c.Counter(); singles=c.Counter(); byaxisstatus=c.Counter(); seen=set()
    for r in rows:
        key=(r["accepted_species"],r["axis"])
        assert key not in seen;seen.add(key)
        axis,quality,st=r["axis"],r["quality"],r["release_status"]
        assert axis in AXES and quality in QUALITIES
        status[st]+=1;byaxisstatus[axis,st]+=1
        if st!="rights_review_required":continue
        b=frozenset(x for x in r["blocked_source_families"].split("|") if x)
        assert b and len(b)==int(r["blocked_source_family_count"])
        sets[b]+=1;byaxis[b,axis]+=1;byquality[b,quality]+=1
        for s in b:incidence[s]+=1
        if len(b)==1:singles[next(iter(b))]+=1
    return sets,byaxis,byquality,status,incidence,singles,byaxisstatus
def greedy(patterns,n):
    selected=set();steps=[]
    for _ in range(n):
        marginal=c.Counter()
        for b,count in patterns.items():
            left=b-selected
            if len(left)==1:marginal[next(iter(left))]+=count
        if not marginal:break
        k,gain=min(marginal.items(),key=lambda x:(-x[1],x[0]))
        selected.add(k);steps.append((k,gain))
    return steps
def analyze(rights,provenance):
    check(rights,RIGHTS_SHA);check(provenance,PROVENANCE_SHA)
    with zipfile.ZipFile(rights) as z:
        with gzip.GzipFile(fileobj=z.open("CELL_RELEASE_RIGHTS.csv.gz")) as g:
            sets,byaxis,byquality,status,inc,single,bas=counts(csv.DictReader(io.TextIOWrapper(g,encoding="utf-8")))
        frozen=json.loads(z.read("FINAL_RELEASE_GATE.json"))
    assert dict(status)==frozen["release_status_counts"] and sum(status.values())==222688
    assert len(inc)==5032 and sum(single.values())==63174 and len(sets)==7898
    assert status["redistributable"]==46274 and status["semantic_review_required"]==1
    repro=c.Counter({b:n for (b,a),n in byaxis.items() if a=="reproductive_assurance"})
    gre=greedy(sets,20);gr=greedy(repro,8)
    def summarize(group):
        s=set(group);n=sum(v for b,v in sets.items() if b<=s)
        return {"families":len(s),"additional_cells":n,"potential_total":46274+n,
                "touched_cells":sum(v for b,v in sets.items() if b&s),
                "additional_by_axis":{a:sum(v for (b,k),v in byaxis.items() if k==a and b<=s) for a in AXES},
                "additional_by_quality":{q:sum(v for (b,k),v in byquality.items() if k==q and b<=s) for q in QUALITIES}}
    def ranked(d,n):return [k for k,_ in sorted(d.items(),key=lambda x:(-x[1],x[0]))[:n]]
    scenarios={
      "efloras_three_labels":summarize(["domain:efloras.org","provider:efloras","provider_treatment:flora_of_china"]),
      "european_four_databases":summarize(["database:pladias","database:floraweb_biolflor","database:baseflor","database:ecoflora"]),
      "four_item_license_candidates":summarize(["provider_treatment:flora_malesiana","database:baseflor","domain:gbif.org","grassbase:url"]),
      "top10_incidence":summarize(ranked(inc,10)),
      "top10_singleton":summarize(ranked(single,10)),
      "all_axes_greedy10":summarize([k for k,_ in gre[:10]]),
      "all_axes_greedy20":summarize([k for k,_ in gre]),
      "reproductive_greedy4":summarize([k for k,_ in gr[:4]]),
      "reproductive_greedy8":summarize([k for k,_ in gr]),
    }
    family_stats=[{"family":s,"affected":inc[s],"sole_blocker":single[s],
      "sole_blocker_reproductive":repro[frozenset([s])],
      "sole_blocker_high_medium":sum(byquality[frozenset([s]),q] for q in ("high","medium"))}
      for s in SOURCES]
    # Exact provenance-source scope diagnostic: use stored historical maps only.
    with zipfile.ZipFile(provenance) as z:
        details=c.Counter();gbif_paths=c.Counter();fm_packages=set();fm_n=0
        with io.TextIOWrapper(z.open("SOURCE_LINEAGE_DETAILS.csv"),encoding="utf-8") as f:
            for r in csv.DictReader(f):
                s=r["source_family"]
                if s in SOURCES:details[s]+=1
                if s=="domain:gbif.org":
                    token=r["source_lineage"].removeprefix("url:")
                    p=urlsplit(token).path.strip("/").split("/")
                    p=[x for x in p if x not in ("en","ja","fr","es","de")]
                    gbif_paths["species_page" if p and p[0]=="species" else "other"]+=1
        with io.TextIOWrapper(z.open("FLORML_LINEAGE_FAMILY_MAP.csv"),encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["source_family"]=="provider_treatment:flora_malesiana":
                    fm_n+=1;fm_packages.add(r["source_url"])
    assert gbif_paths["species_page"]==1867 and gbif_paths["other"]==0
    assert fm_n==3393 and len(fm_packages)==1
    # Only the public provider-package URL (no cell/lineage rows) is retained in prose.
    result={"source_rights_zip_sha256":RIGHTS_SHA,"source_provenance_zip_sha256":PROVENANCE_SHA,
      "status_counts":dict(status),"blocker_family_count":len(inc),"unique_blocker_sets":len(sets),
      "cells_with_one_blocker":sum(single.values()),"reproductive_cells_with_one_blocker":sum(n for b,n in repro.items() if len(b)==1),
      "axis_baseline":{a:{"resolved":sum(v for (j,k),v in bas.items() if j==a),"authorized":bas[a,"redistributable"]} for a in AXES},
      "selected_family_stats":family_stats,"scenarios":scenarios,
      "all_axes_greedy_steps":[{"family":label(k),"marginal_cells":n} for k,n in gre],
      "reproductive_greedy_steps":[{"family":label(k),"marginal_reproductive_cells":n} for k,n in gr],
      "source_scope":{"GBIF_historical_lineage_rows":dict(gbif_paths),"Flora_Malesiana_mapped_lineages":fm_n,
        "Flora_Malesiana_distinct_WFO_packages":len(fm_packages),"historical_source_detail_counts":dict(details)},
      "decision":"Counterfactual priorities only. All blockers must be legally resolved; unknown citation families are not license-granting providers. No promotion, policy change, cell release, new analysis, or editor contact."}
    return result
def self_test():
    p=c.Counter({frozenset(["a"]):5,frozenset(["b"]):2,frozenset(["a","b"]):9,frozenset(["c","d"]):30})
    assert greedy(p,2)==[("a",5),("b",11)]
    assert label("arbitrary:species or URL").startswith("withheld_family_")
    print("RIGHTS_PRIORITY_SELF_TEST_PASS")
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--rights",type=Path)
    p.add_argument("--provenance",type=Path)
    p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:self_test()
    else:
        r=analyze(a.rights,a.provenance)
        a.out.write_text(json.dumps(r,indent=2,sort_keys=True,allow_nan=False)+"\n")
        print("RIGHTS_PRIORITY_RESULT "+json.dumps(r,sort_keys=True,allow_nan=False))
