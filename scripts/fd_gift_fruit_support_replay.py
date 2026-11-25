#!/usr/bin/env python3
"""Replay the frozen public direct fruit-dryness source; audit genus and documentation support.
No animal trait association; no export of plant names, island IDs, or source rows.
"""
import argparse,collections,concurrent.futures,hashlib,io,json,pathlib
import pandas as pd
from fd_original_tropical_pollinator_bridge import read_pin,accepted_genus
from fd_gift_control_trait_acquisition import fetch,url,txt,norm,pairs_for,parse_record,consolidate,sha,coverage
TRAIT="3.18.1"
SOURCE_DIGEST="bb55f76c47172fd52a74bb2a44dd78aef7b712d102ccd3f56523488a6b7093f0"
META_PIN="ee20bafd237be38b4cb97d65219f65d1f02ecc5fca09867c9cdfc1640fb30af4"
PAIR_PIN="e09622191683a5cb51f90684e0d6b29db88e093cf6e93e7fd4862dadb844967c"
SPECIES_PINS=["482b979ee66bbe0f9ae3d2a115dc09d54219ca364951d771ca2b0bdf25adc3fb","e2c437f4c9af81f7a6078a464e6cbb5ad985ffeb392837612ca7beb346348191","95cb3c407b18e32d972cd3be3d7771df3c5b33207c8477e54b2e16a99e64d57c","7ab61c4dd24d73533b0aebad29e962f7a67fb081c60f43ff9c9712f6156df51d"]
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path)
    p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path)
    a=p.parse_args();a.out.parent.mkdir(parents=True,exist_ok=True)
    f=pd.read_csv(io.BytesIO(read_pin(a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any() and len(geo)==8264
    f["name"]=f.accepted_species.map(norm);f["genus"]=f.accepted_species.map(accepted_genus)
    amb=set(f[f.name.notna()].groupby("name").accepted_species.nunique().loc[lambda x:x>1].index)
    meta,h=fetch(url(query="traits_meta"));assert h==META_PIN
    rt,h=fetch(url(query="reference_traits"));assert h==PAIR_PIN
    ref_rows,refhash=fetch(url(query="references"))
    refs={txt(x.get("ref_ID")):x for x in ref_rows if txt(x.get("ref_ID"))}
    rp,excluded=pairs_for(rt,refs,TRAIT);assert len(rp)==3
    species={};bad=set()
    for i,expected in enumerate(SPECIES_PINS):
        rows,h=fetch(url(query="species",startat=i*100000,limit=100000));assert h==expected
        for r in rows:
            k=txt(r.get("work_ID"));n=txt(r.get("work_species"))
            if not k or not n:continue
            if k in species and species[k]!=n:bad.add(k)
            species[k]=n
    assert not bad and len(species)==372098
    good=[];hashes=[]
    for ref in rp:
        rows,h=fetch(url(query="traits_raw",traitid=TRAIT,deriv=0,biasderiv=0,refid=ref))
        hashes.append((ref,h))
        for row in rows:
            r,e=parse_record(row,TRAIT,ref,refs,species)
            if r:good.append(r)
    digest=sha(json.dumps(sorted(hashes),separators=(",",":")).encode())
    assert digest==SOURCE_DIGEST, "fruit raw responses changed; do not silently re-freeze"
    states,checks=consolidate(good,amb)
    matched=set(f.name.dropna())&set(states);assert len(matched)==2183
    assert collections.Counter(states[k] for k in matched)=={"dry":993,"fleshy":1190}
    refnames={ref:set(r[0] for r in good if r[2]==ref)&matched for ref in rp}
    rmask=collections.defaultdict(set)
    for ref,names in refnames.items():
        for name in names:rmask[name].add(ref)
    # Reuse frozen public animal bridge only for geographic support (no diet or fruit fraction regression).
    b=json.loads(read_pin(a.bridge,"bridge"))
    g=pd.DataFrame(b["geometry"],columns=b["public_bridge_columns"]).rename(columns={"gshhg_geometry_hash_id":"island_id"})
    birds=pd.DataFrame(b["bird_candidate"])[["bee_id","group"]]
    join=g[["island_id","bee_id"]].merge(birds,on="bee_id",validate="one_to_one")
    assert len(join)<=52
    f=f.merge(geo[["island_id","analysis_regime"]],on="island_id",how="inner",validate="many_to_one")
    f["state"]=f.name.map(states)
    counts={}
    for scope,d in [("all_recorded",f),("known_native",f[f.origin_status.eq("native")])]:
        counts[scope]={}
        for region,z in [("all",d),*d.groupby("analysis_regime")]:
            pl=z[z.state.notna()&z.genus.notna()].copy()
            tab=pl.drop_duplicates("accepted_species")
            n=pl.groupby("island_id").accepted_species.nunique().rename("scored")
            vary=set(tab.groupby("genus").state.nunique().loc[lambda x:x>1].index)
            nvar=pl[pl.genus.isin(vary)].groupby("island_id").accepted_species.nunique()
            orig=z.groupby("island_id").accepted_species.nunique()
            n=n.reindex(orig.index,fill_value=0);nvar=nvar.reindex(orig.index,fill_value=0)
            eligible=n.ge(5);m=join[join.island_id.isin(n[eligible].index)]
            gr=m.groupby("group").size()
            refs_diag=[]
            for idx,ref in enumerate(rp,1):
                names=refnames[ref];only={k for k in matched if rmask[k]=={ref}}
                nz=pl[~pl.name.isin(only)].groupby("island_id").accepted_species.nunique().reindex(orig.index,fill_value=0)
                refs_diag.append({"reference_ordinal":idx,"recorded_candidate_species_in_region":int(tab.name.isin(names).sum()),"exclusive_species_in_region":int(tab.name.isin(only).sum()),"islands_ge5_after_reference_omission":int(nz.ge(5).sum())})
            counts[scope][str(region)]={"species":len(tab),"genera":int(tab.genus.nunique()),"variable_genera":len(vary),"species_in_variable_genera":int(tab.genus.isin(vary).sum()),
                "islands_ge5":int(eligible.sum()),"islands_ge5_with_variable_genus_carrier":int((eligible&nvar.gt(0)).sum()),
                "median_coverage_ge5":float((n[eligible]/orig[eligible]).median()) if eligible.any() else None,
                "bird_geographic_candidate_islands_ge5":len(m),"bird_candidate_groups":len(gr),"bird_candidate_groups_ge2_upper_bound":int(gr.ge(2).sum()),"reference_omission":refs_diag}
    out={"status":"replayed_direct_candidates_support_only","original_parent":4505,"geography_parent":8264,"source_run":45530336990,"fruit_digest":digest,"reference_metadata_sha256":refhash,
        "matched_names":len(matched),"counts":counts,"warning":"Variable genus without evidence-grade stratification is an optimistic movable-carrier upper bound; current candidate records not yet graded. Bird candidate geography provisional. No effect or animal diet comparison."}
    a.out.write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")
    print("FRUIT_SUPPORT_AUDIT",json.dumps(out,allow_nan=False),flush=True)
if __name__=="__main__":main()
