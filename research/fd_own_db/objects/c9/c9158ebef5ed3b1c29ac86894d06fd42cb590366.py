"""Family-local diet completeness, without imputing missing diets in other families.

Family assignment is an exact, conservative MDD nomenclatural classification.
The source diet family's label must agree when a donor exists. If a mapped
taxon has unresolved family, that island cannot certify focal-family closure.
Species-level diet transfer still uses the existing exact/singleton rule.
"""
from collections import Counter
import numpy as np
from bat_family_context import FAMILIES,CONTROLS,OUTCOMES
from bat_mdd_synonym_audit import col,index_names

ALLOWED={"available","name_combination","nomen_novum","justified_emendation","as_emended"}
def taxonomy(checklist,syn,map_names,aug):
    # Strict status allowlist for this NEW family assignment; do not revise
    # the previously archived species-level sensitivity in place.
    s=col(syn,"MDD_nomenclature_status")
    keep=syn[s].map(lambda x: bool(x.strip()) and set(t.strip().casefold() for t in x.split("|"))<=ALLOWED)
    idx,_=index_names(checklist,syn[keep])
    order=col(checklist,"order","MDD_order")
    c=checklist[checklist[order].str.casefold().eq("chiroptera")]
    i=col(c,"id","MDD_species_id","MDD_id")
    fam=col(c,"family","MDD_family")
    family_by_id=dict(zip(c[i],c[fam].str.strip().str.casefold()))
    assert len(family_by_id)==len(c) and all(family_by_id.values())
    result={};audit=Counter()
    for n in sorted(map_names):
        f={family_by_id[j] for j in idx.get(n,set())}
        if len(f)!=1:
            audit["unresolved_family"]+=1
            continue
        name=next(iter(f))
        if n in aug.index and str(aug.loc[n,"family"]).strip().casefold()!=name:
            audit["diet_family_conflict"]+=1
            continue
        result[n]=name
        audit["assigned"]+=1
        if name in FAMILIES:audit[name]+=1
    assert len(result)+audit["unresolved_family"]+audit["diet_family_conflict"]==len(map_names)
    return result,dict(audit)

def allocate(n):
    # [nonnectar, nectar, diet unknown]; separately count unknown family.
    return {f:np.zeros((n,3),dtype=int) for f in FAMILIES}|{"unresolved":np.zeros(n,dtype=int)}

def add(a,f,k,idx):
    if f is None:
        a["unresolved"][idx]+=1
    elif f in FAMILIES:
        a[f][idx,2 if k==3 else int(k>0)]+=1

def attach(frame,a,prefix):
    frame[prefix+"_unknown_family"]=a["unresolved"]
    for f in FAMILIES:
        for j,t in enumerate(("other","nectar","unknown")):
            frame[prefix+"_"+f+"_"+t]=a[f][:,j]
    summed=a["unresolved"]+sum(a[f].sum(1) for f in FAMILIES)
    assert (summed<=frame.total.to_numpy()).all()

def panels(fl,tab,geo,estimate):
    result={}
    for scope,v in fl.items():
        d=v.merge(tab[.5],on="island_id",validate="one_to_one").merge(geo,on="island_id",validate="one_to_one")
        d=d[d.analysis_regime.eq("tropical")&d.range_area.ge(20)&d.n_species.ge(5)&d.distance_to_continent_km.ge(0)].copy()
        d["block"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
        d["log_area"]=np.log2(d.range_area);d["log_distance"]=np.log2(1+d.distance_to_continent_km)
        d["log_n"]=np.log2(d.n_species);d["abs_lat"]=d.island_latitude.abs()
        for prefix in ("local_exact","local_mdd"):
            for family in FAMILIES:
                s=d.copy()
                s["focal_nectar"]=s[prefix+"_"+family+"_nectar"]
                s["focal_other"]=s[prefix+"_"+family+"_other"]
                s["focal_unknown"]=s[prefix+"_"+family+"_unknown"]
                s["focal_total"]=s.focal_nectar+s.focal_other+s.focal_unknown
                s=s[s.focal_total.gt(0)&s[prefix+"_unknown_family"].eq(0)].copy()
                s["log_family_nectar"]=np.log1p(s.focal_nectar)
                s["log_family_other"]=np.log1p(s.focal_other)
                s["log_nonfamily"]=np.log1p(s.total-s.focal_total)
                s=s.dropna(subset=["island_latitude","island_longitude",*CONTROLS,*OUTCOMES]).copy()
                s=s[np.isfinite(s[CONTROLS+list(OUTCOMES)].to_numpy(float)).all(1)]
                full=s[s.focal_unknown.eq(0)].copy()
                # Same strict taxonomy in both fits; only all-bat diet
                # completeness differs. Original exact-match closure inferred
                # from preserved unknown count, not the augmented diet count.
                all_complete=full[full[prefix+"_global_unknown"].eq(0)]
                fits={}
                for tag,q in (("all_bat_diets_complete",all_complete),("focal_family_diets_complete",full)):
                    fits[tag]={"cohort":len(q),"focal_nectar_positive":int(q.focal_nectar.gt(0).sum()),
                               "focal_other_positive":int(q.focal_other.gt(0).sum()),
                               "fit":estimate(q,["log_family_nectar","log_family_other"],CONTROLS)}
                result[scope+"_"+prefix+"_"+family]={"taxonomically_resolved_family_positive":len(s),
                        "focal_diet_complete":len(full),"newly_eligible_due_to_family_local_closure":len(full)-len(all_complete),
                        "fits":fits}
    return {"panels":result,"warning":"Post-pilot strict family taxonomy, all map names assigned to a family and focal family diet complete. Other families may have missing diet; only total nonfocal species count used. No diet imputation. 30-block/200 informative-island gate unchanged; no evolutionary or pollination claim."}

def self_test():
    import pandas as pd
    a=allocate(3);add(a,"pteropodidae",1,np.array([0,1]));add(a,"pteropodidae",3,np.array([1]));add(a,None,3,np.array([2]))
    d=pd.DataFrame({"total":[1,2,1]});attach(d,a,"test")
    assert d.test_pteropodidae_nectar.tolist()==[1,1,0]
    assert d.test_pteropodidae_unknown.tolist()==[0,1,0]
    assert d.test_unknown_family.tolist()==[0,0,1]
    print("BAT_FAMILY_LOCAL_CLOSURE_TEST_PASS",flush=True)
if __name__=="__main__":self_test()
