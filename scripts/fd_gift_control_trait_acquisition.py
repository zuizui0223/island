#!/usr/bin/env python3
"""Protected GIFT 3.2 direct fruit-dryness / growth-form acquisition and eligibility.
No animal/trait association or species-level output. Frozen GIFT4505 roster unchanged.
"""
from __future__ import annotations
import argparse, collections, concurrent.futures, hashlib, io, json, math, pathlib, re, sys, time, urllib.parse, urllib.request
import pandas as pd
import numpy as np
from fd_original_tropical_pollinator_bridge import read_pin

ENDPOINT="https://gift.uni-goettingen.de/api/extended/index3.2.php"
VERSIONS="https://gift.uni-goettingen.de/api/index.php?query=versions"
TARGETS={"3.18.1":{"dry","fleshy"},"1.2.1":{"herb","shrub","tree","other"}}
# Scope fixed before receiving new fruit records. No fruit_type->dryness inference.
def txt(x):return "" if x is None else " ".join(str(x).split())
def flag(x):return txt(x).lower() in {"1","true","yes"}
def sha(b):return hashlib.sha256(b).hexdigest()
def url(**kwargs):return ENDPOINT+"?"+urllib.parse.urlencode(kwargs)
def fetch(u,cap=100_000_000):
    for attempt in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(u,headers={"User-Agent":"GIFT island research source audit"}),timeout=90) as r:
                b=r.read(cap+1)
            if len(b)>cap:raise ValueError("source_response_over_cap")
            j=json.loads(b)
            if not isinstance(j,list) or any(not isinstance(x,dict) for x in j):raise ValueError("unexpected_api_schema")
            return j,sha(b)
        except Exception:
            if attempt==2:raise
            time.sleep(2**attempt)
def norm(s):
    v=txt(s).split()
    if len(v)!=2 or not re.fullmatch(r"[A-Za-z-]+",v[0]) or not re.fullmatch(r"[A-Za-z-]+",v[1]):return None
    if v[1].lower() in {"sp","spp","cf","aff"}:return None
    return " ".join(v).casefold()
def pairs_for(meta,refs,trait):
    pairs=set()
    restricted=set()
    for row in meta:
        if txt(row.get("bias")) not in {"","0"}:continue
        ref=txt(row.get("ref_ID"))
        if not ref or not any(str(k).startswith("trait") and txt(v)==trait for k,v in row.items()):continue
        r=refs.get(ref,{})
        if txt(r.get("restricted"))!="0" or not txt(r.get("ref_long")):restricted.add(ref);continue
        pairs.add(ref)
    return sorted(pairs,key=lambda s:(len(s),s)),len(restricted)
def parse_record(row,trait,ref,refs,species):
    if txt(row.get("ref_ID"))!=ref or txt(row.get("trait_ID"))!=trait:return None,"scope_mismatch"
    if flag(row.get("derived")) or flag(row.get("bias_deriv")):return None,"derived"
    if txt(refs.get(ref,{}).get("restricted"))!="0" or not txt(refs.get(ref,{}).get("ref_long")):return None,"restricted_or_no_citation"
    if not flag(row.get("matched")) or not flag(row.get("resolved")):return None,"unresolved"
    if any(flag(row.get(k)) for k in ("cf_genus","cf_species","aff_species")) or txt(row.get("subtaxon")):return None,"uncertain_or_infraspecific"
    work=txt(row.get("work_ID"));name=norm(species.get(work))
    if not name:return None,"missing_work_name"
    state=txt(row.get("trait_value")).lower()
    if state not in TARGETS[trait]:return None,"non_single_state"
    rid=txt(row.get("trait_derived_ID"))
    if not rid:return None,"missing_record_id"
    return (name,state,ref,rid,work),None
def consolidate(rows,ambiguous_names):
    # Repeated same-state records provide provenance, not repeated species.
    d=collections.defaultdict(set);refs=collections.defaultdict(set);ids=collections.defaultdict(set)
    for name,state,ref,rid,work in rows:
        d[name].add(state);refs[name].add(ref);ids[name].add(work)
    ambiguous={n for n,v in ids.items() if len(v)!=1}|set(ambiguous_names)
    conflict={n for n,v in d.items() if len(v)!=1}
    states={n:next(iter(v)) for n,v in d.items() if n not in ambiguous|conflict}
    return states,{"candidate_names":len(d),"conflicting_state_names":len(conflict),"ambiguous_work_names":len(set(d)&ambiguous),
                   "unique_unambiguous_names":len(states),"names_with_multiple_references":sum(len(refs[n])>=2 for n in states)}
def coverage(flora,geo,states,scope):
    f=flora if scope=="all_recorded" else flora[flora.origin_status.eq("native")]
    f=f[f.name.notna()].copy();f["state"]=f.name.map(states)
    base=f.groupby("island_id").accepted_species.nunique().rename("roster")
    scored=f[f.state.notna()].groupby("island_id").agg(n=("accepted_species","nunique"),states=("state","nunique"))
    d=geo[["island_id","analysis_regime"]].merge(base.reset_index(),on="island_id",how="inner",validate="one_to_one").merge(scored.reset_index(),on="island_id",how="left",validate="one_to_one")
    d["n"]=d.n.fillna(0).astype(int);d["states"]=d.states.fillna(0).astype(int)
    d["fraction"]=d.n/d.roster
    out={}
    for reg,z in [("all",d),*d.groupby("analysis_regime",dropna=False)]:
        out[str(reg)]={"islands_with_roster":len(z),"islands_ge1":int(z.n.ge(1).sum()),"islands_ge5":int(z.n.ge(5).sum()),
                      "islands_ge30":int(z.n.ge(30).sum()),"islands_ge5_and_ge50pct":int((z.n.ge(5)&z.fraction.ge(.5)).sum()),
                      "islands_ge5_multistate":int((z.n.ge(5)&z.states.ge(2)).sum()),
                      "median_scored_among_ge5":float(z.loc[z.n.ge(5),"n"].median()) if z.n.ge(5).any() else None}
    return out
def self_test():
    ref={"7":{"restricted":"0","ref_long":"source"}}
    sp={"1":"Abc def"}
    row=dict(ref_ID="7",trait_ID="3.18.1",derived="0",bias_deriv="0",matched="1",resolved="1",work_ID="1",trait_derived_ID="10",trait_value="fleshy")
    r,e=parse_record(row,"3.18.1","7",ref,sp);assert r and e is None
    assert parse_record(dict(row,trait_value="dry/fleshy"),"3.18.1","7",ref,sp)[1]=="non_single_state"
    assert parse_record(dict(row,derived="1"),"3.18.1","7",ref,sp)[1]=="derived"
    assert consolidate([r,("abc def","dry","8","11","1")],set())[0]=={}
    assert consolidate([r,r],set())[0]=={"abc def":"fleshy"}
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path)
    p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--trait",required=True,choices=tuple(TARGETS))
    a=p.parse_args();self_test();a.out.parent.mkdir(parents=True,exist_ok=True)
    f=pd.read_csv(io.BytesIO(read_pin(a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert len(geo)==8264 and geo.island_id.is_unique
    f["name"]=f.accepted_species.map(norm)
    # Never use fuzzy names or synthesize missing fruit morphology.
    ambiguous=set(f[f.name.notna()].groupby("name").accepted_species.nunique().loc[lambda x:x>1].index)
    out={"status":"source_audit_only","trait_id":a.trait,"gift_version":"3.2","botanical_parent":4505,"geography_parent":8264,
         "scope":"public raw non-derived records; direct candidates not adjudicated to high/medium evidence",
         "source_pins":{},"coverage":None}
    def checkpoint():
        a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    checkpoint()
    try:
        versions,h=fetch(VERSIONS);out["source_pins"]["versions"]=h
        if "3.2" not in {txt(x.get("version")) for x in versions}:raise ValueError("version_not_available")
        for query in ("traits_meta","references","reference_traits"):
            records,h=fetch(url(query=query));out["source_pins"][query]=h
            if query=="traits_meta":meta=records
            elif query=="references":refs={txt(x.get("ref_ID")):x for x in records if txt(x.get("ref_ID"))}
            else:ref_traits=records
        target=[x for x in meta if txt(x.get("Lvl3"))==a.trait]
        assert len(target)==1
        out["trait_metadata"]={k:target[0].get(k) for k in ("Lvl3","Trait2","Units","count")}
        refs_ok,nrestricted=pairs_for(ref_traits,refs,a.trait)
        out["reference_pairs"]={"unrestricted_with_citation":len(refs_ok),"restricted_or_missing":nrestricted}
        if not refs_ok:raise ValueError("no_eligible_public_references")
        checkpoint()
        # Read species dictionary in stable versioned pages; no raw rows leave this job.
        species={};bad_ids=set()
        for start in range(0,600_000,100_000):
            records,h=fetch(url(query="species",startat=start,limit=100_000))
            out["source_pins"]["species_"+str(start)]=h
            for row in records:
                k=txt(row.get("work_ID"));n=txt(row.get("work_species"))
                if not k or not n:continue
                if k in species and species[k]!=n:bad_ids.add(k)
                species[k]=n
            if len(records)<100_000:break
        for k in bad_ids:species.pop(k,None)
        out["species_dictionary"]={"distinct_work_ids":len(species),"conflicting_ids":len(bad_ids)}
        print("GIFT_CONTROL_METADATA",json.dumps({"trait_id":a.trait,"reference_pairs":out["reference_pairs"],"species_dictionary":out["species_dictionary"]}),flush=True)
        # Limited parallelism, process every eligible reference; failed sources stay explicit.
        def one(ref):
            try:
                rows,h=fetch(url(query="traits_raw",traitid=a.trait,deriv=0,biasderiv=0,refid=ref))
                return ref,rows,h,None
            except Exception as e:return ref,None,None,type(e).__name__
        good=[];reasons=collections.Counter();nraw=0;failed=collections.Counter();hashes=[];succeeded=0
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            for ref,rows,h,err in executor.map(one,refs_ok):
                if err:failed[err]+=1;continue
                succeeded+=1;nraw+=len(rows);hashes.append((ref,h))
                for row in rows:
                    rec,err=parse_record(row,a.trait,ref,refs,species)
                    if err:reasons[err]+=1
                    else:good.append(rec)
                if succeeded%25==0:
                    print("GIFT_CONTROL_PROGRESS",json.dumps({"trait_id":a.trait,"references_read":succeeded,"references_total":len(refs_ok),"raw_rows":nraw}),flush=True)
        out["acquisition"]={"references_succeeded":succeeded,"references_failed":sum(failed.values()),"failure_types":dict(failed),
                            "raw_rows":nraw,"direct_candidate_rows":len(good),"exclusion_reasons":dict(reasons),
                            "sorted_reference_response_hash_digest":sha(json.dumps(sorted(hashes),separators=(",",":")).encode())}
        states,checks=consolidate(good,ambiguous)
        # Unambiguous source names also must match the frozen accepted-name universe.
        matched_names=set(f.name.dropna())&set(states)
        out["reconciliation"]=checks|{"exact_matched_original_names":len(matched_names),"accepted_candidate_state_counts":dict(collections.Counter(states[n] for n in matched_names))}
        out["coverage"]={scope:coverage(f,geo,states,scope) for scope in ("all_recorded","known_native")}
        out["status"]="candidate_coverage_complete" if succeeded==len(refs_ok) else "candidate_coverage_partial_sources"
        out["primary_e2_status"]="not_run_unreviewed_evidence_and_animal_input"
    except Exception as e:
        out["status"]="acquisition_failed";out["error_type"]=type(e).__name__;out["error_code"]=str(e) if str(e) in {"source_response_over_cap","unexpected_api_schema","version_not_available","no_eligible_public_references"} else "see_source_pipeline"
    checkpoint()
    print("GIFT_CONTROL_AUDIT",json.dumps(out,ensure_ascii=False,allow_nan=False),flush=True)
if __name__=="__main__":main()
