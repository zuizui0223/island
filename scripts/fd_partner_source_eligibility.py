#!/usr/bin/env python3
"""Eligibility only: original GIFT4505 flower/fruit/life-form coverage and Marshall Bombus lists.
No effect fit and no protected identities exported. Contract frozen in fd commit ab0ee26f.
"""
import argparse, hashlib, io, json, pathlib, re, urllib.request, zipfile
import numpy as np
import pandas as pd
from fd_original_tropical_pollinator_bridge import read_pin, get, pair_states

BEE_MD5="372c22e2c24e203600a29d7f874a907d"
MODEL_SHA="a7ddae4cdb41453313d20811371cf64bd38b42b788fc89ae96c284d8d667ddd7"
SPECIES_SHA="95330e8408cf3f3816f65789c8bb374633b64342657e10051ea7e053818ac8e8"
GIFT_META="https://gift.uni-goettingen.de/api/extended/index3.2.php?query=traits_meta"

def digest(b):return hashlib.sha256(b).hexdigest()

def decode(b):
    for enc in ("utf-8-sig","cp1252","latin-1"):
        try:return b.decode(enc)
        except UnicodeDecodeError:pass
    raise ValueError("unknown source encoding")

def fruit_audit(t,flora):
    x=t[t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    x=x[x.normalized_value.str.strip().ne("")&~x.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    counts=x.groupby("trait_name").accepted_species.nunique().astype(int).to_dict()
    candidate={k:v for k,v in counts.items() if re.search(r"fruit|dispers|growth.form|life.form",k,re.I)}
    out={"ledger_direct_distinct_species_by_trait":counts,"candidate_direct_traits":candidate,
         "fruit_primary_status":"no_certified_direct_fleshy_fruit_endpoint" if not candidate else "requires_ontology_and_source_review",
         "warning":"Trait names or dispersal mode alone do not certify fruit morphology. No effect was fitted."}
    try:
        raw=get(GIFT_META,5_000_000);j=json.loads(raw)
        if isinstance(j,dict):
            poss=[v for v in j.values() if isinstance(v,list)]
            if len(poss)!=1:raise ValueError("unrecognized metadata schema")
            j=poss[0]
        if not isinstance(j,list) or not j or not isinstance(j[0],dict):raise ValueError("unrecognized metadata records")
        # Only public trait metadata fields; no species/provider records.
        fields=("Lvl1","Category","Lvl2","Trait1","Lvl3","Trait2","Units","type","comment","count","trait_ID")
        matched=[]
        for row in j:
            if any(re.search(r"fruit|dispers|growth.form|life.form",str(row.get(k,"")),re.I) for k in ("Trait1","Trait2")):
                matched.append({k:row[k] for k in fields if k in row})
        out["gift_3_2_metadata"]={"status":"read","sha256":digest(raw),"records":len(j),"candidates":matched,
          "warning":"Published metadata counts, not direct/eligible island coverage. Need raw un-derived un-restricted records and fixed ontology."}
    except Exception as e:
        out["gift_3_2_metadata"]={"status":"unavailable","error_type":type(e).__name__}
    return out

def fetch_bee():
    rec=json.loads(get("https://zenodo.org/api/records/23222636",2_000_000))
    z=next(v for v in rec["files"] if v["key"]=="Data.zip")
    assert z["size"]==262554015 and z["checksum"]=="md5:"+BEE_MD5
    raw=get(z["links"].get("content") or z["links"]["self"],300_000_000)
    assert len(raw)==z["size"] and hashlib.md5(raw).hexdigest()==BEE_MD5
    with zipfile.ZipFile(io.BytesIO(raw)) as f:
        m=f.read("Data/model_table.csv");s=f.read("Data/islandspecieslist.csv")
    assert digest(m)==MODEL_SHA and digest(s)==SPECIES_SHA
    m=pd.read_csv(io.StringIO(decode(m)),dtype=str,keep_default_na=False)
    s=pd.read_csv(io.StringIO(decode(s)),dtype=str,keep_default_na=False)
    assert len(m)==305 and m.OBJECTID_1.is_unique and set(s.OBJECTID_1)<=set(m.OBJECTID_1)
    return m,s

def count_bees(m,s):
    for c in ("NonNative","Uncertain"):
        s[c]=s[c].str.strip()
    s["name"]=s.Corrected.str.replace("_"," ",regex=False).str.strip().str.replace(r"\s+"," ",regex=True)
    # Fail closed for checklist qualification if any native label/name is unresolved.
    s["native"]=s.NonNative.eq("0")
    s["certain"]=s.Uncertain.eq("0")
    s["valid_binomial"]=s.name.str.match(r"^[A-Z][A-Za-z-]+ [a-z][a-z-]+(?: |$)") & ~s.name.str.contains(r"\b(?:sp|spp|cf|aff)\b",regex=True)
    s["bombus"]=s.name.str.split().str[0].eq("Bombus")
    rows=[]
    for island,v in s.groupby("OBJECTID_1"):
        nat=v[v.native];certain=nat[nat.certain&nat.valid_binomial]
        allnames=set(nat.name);names=set(certain.name);bom=set(certain.loc[certain.bombus,"name"])
        rows.append(dict(bee_id=island,native_all=len(allnames),certain_native=len(names),bombus=len(bom),
            other=len(names-bom),uncertain_bombus=int((nat.bombus&~(nat.certain&nat.valid_binomial)).sum()),
            invalid_native=int((~nat.valid_binomial).sum()),unknown_native_flag=int((~v.NonNative.isin(["0","1"])).sum()),
            unknown_uncertain_flag=int((~nat.Uncertain.isin(["0","1"])).sum())))
    x=pd.DataFrame(rows)
    m=m.rename(columns={"OBJECTID_1":"bee_id","Archipelago":"group","FinalRichNat":"source_native","Cov":"coverage","AreaModel":"area"})
    d=m[["bee_id","group","source_native","coverage","area"]].merge(x,on="bee_id",how="left",validate="one_to_one")
    for c in ("source_native","coverage","area"):d[c]=pd.to_numeric(d[c],errors="coerce")
    # No rows in a checklist is not a qualified nonrecord. Source coverage scale is explicitly audited.
    d["reconciled"]=d.native_all.notna()&d.source_native.eq(d.native_all)
    d["qualified"]=d.reconciled&d.coverage.between(.90,1)&d.area.ge(1)&d.certain_native.gt(0)&d.uncertain_bombus.eq(0)&d.invalid_native.eq(0)&d.unknown_native_flag.eq(0)&d.unknown_uncertain_flag.eq(0)&d.group.str.strip().ne("")&~d.group.str.lower().isin(["na","nan","none","unknown","unresolved"])
    return d,s

def support(d):
    counts=d.groupby("group").bombus.agg(lambda v:int(v.gt(0).nunique()))
    return {"islands":len(d),"groups":int(d.group.nunique()),"positive":int(d.bombus.gt(0).sum()),
        "nonrecord":int(d.bombus.eq(0).sum()),"informative_presence_groups":int(counts.gt(1).sum()),
        "gate_200_islands_30_groups":bool(len(d)>=200 and counts.gt(1).sum()>=30)}

def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",type=pathlib.Path,required=True);p.add_argument("--geography",type=pathlib.Path,required=True);p.add_argument("--bridge",type=pathlib.Path,required=True);p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args()
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    bridge=json.loads(read_pin(a.bridge,"bridge"))
    g=pd.DataFrame(bridge["geometry"],columns=bridge["public_bridge_columns"]).rename(columns={"gshhg_geometry_hash_id":"island_id"})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert len(geo)==8264 and geo.island_id.is_unique and len(g)==234 and g.island_id.is_unique and g.bee_id.is_unique
    assert (g[["bee_overlap","gshhg_overlap"]]>=.9-1e-8).all().all()
    a.out.parent.mkdir(parents=True,exist_ok=True)
    out={"status":"source_eligibility_only","contract_fd_commit":"ab0ee26f99c854ee7c8e73c3dc686a77a5a80f57","botanical_parent":4505,"geography_parent":8264,
         "fruit_and_growth":fruit_audit(t,f)}
    print("FRUIT_GROWTH_AUDIT",json.dumps(out["fruit_and_growth"],ensure_ascii=True),flush=True)
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    try:
        m,s=fetch_bee();d,s=count_bees(m,s)
    except Exception as e:
        out["bombus"]={"status":"source_unavailable","error_type":type(e).__name__}
        print("BOMBUS_SOURCE_UNAVAILABLE",type(e).__name__,flush=True)
        a.out.write_text(json.dumps(out,indent=2,allow_nan=False));return
    x=d.merge(g,on="bee_id",validate="one_to_one").merge(geo[["island_id","analysis_regime"]],on="island_id",validate="one_to_one")
    # Original same-species pair states and scored support only. Never inspect trait values by animal exposure.
    pair=pair_states(t);roster=f[f.accepted_species.isin(pair.index)]
    n=roster.groupby("island_id").accepted_species.nunique().rename("scored_n")
    x=x.join(n,on="island_id")
    x["scored_n"]=x.scored_n.fillna(0)
    north=x[x.analysis_regime.eq("northern_midlatitude")].copy()
    # Original regime spelling must be supported before deciding a model is impossible.
    regimes=x.analysis_regime.value_counts(dropna=False).to_dict()
    if "northern_midlatitude" not in regimes:
        # archived geography may label the northern regime 'north_midlatitude'; report literal labels only
        out["regime_warning"]="verify exact original north label; no automatic alternate regime chosen"
    bx=s[s.native&s.certain&s.valid_binomial&s.bombus]
    record_counts={"source_model_rows":len(d),"species_rows":len(s),"native_bombus_distinct_names":int(bx.name.nunique()),
        "native_bombus_recorded_model_islands":int(d.bombus.gt(0).sum()),"native_bombus_uncertain_islands":int(d.uncertain_bombus.gt(0).sum()),
        "model_reconciled_islands":int(d.reconciled.sum()),"model_qualified_islands":int(d.qualified.sum()),
        "coverage_nonmissing":int(d.coverage.notna().sum()),"coverage_outside_0_1":int((~d.coverage.between(0,1)&d.coverage.notna()).sum()),
        "coverage_range":[float(d.coverage.min()),float(d.coverage.max())] if d.coverage.notna().any() else None,
        "native_flag_values":s.NonNative.value_counts().to_dict(),"uncertain_flag_values":s.Uncertain.value_counts().to_dict(),
        "geometry_matches":len(x),"geometry_by_regime":{str(k):int(v) for k,v in regimes.items()}}
    out["bombus"]={"status":"eligibility_only","source":record_counts,"north_matched_with_any_checklist":support(north[north.bombus.notna()]),
        "north_qualified":support(north[north.qualified]),"north_qualified_pair_ge5":support(north[north.qualified&north.scored_n.ge(5)]),
        "all_region_qualified":support(x[x.qualified]),"warning":"Nonrecord is checklist qualified, not verified absence; no E1 or E2 fitted. No plant traits compared by Bombus exposure."}
    print("BOMBUS_ELIGIBILITY",json.dumps(out["bombus"],allow_nan=False),flush=True)
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False))

if __name__=="__main__":main()
