"""Two independent nectar-bearing bat family contexts, post-pilot exploration.

Counts are from the same original range overlap and EltonTraits diet labels.
The family label is the diet source's taxonomy (including reciprocal singleton
donors only in the separately reported sensitivity). Family presence is not a
realm assignment; introduced/accidental occurrences cannot be inferred.
"""
import numpy as np

FAMILIES=("pteropodidae","phyllostomidae")
CONTROLS=["log_area","log_distance","log_n","abs_lat",*[f"climate_pc{i}" for i in range(1,5)],"log_nonfamily"]
OUTCOMES=("white","white_open_brush","red_pink","tubular")

def allocate(n):
    return {f:np.zeros((n,2),dtype=int) for f in FAMILIES}

def add_counts(counts,family,diet_class,indices):
    f=str(family).strip().casefold()
    if f in counts and diet_class in (0,1,2):
        counts[f][indices,int(diet_class>0)]+=1

def attach_counts(frame,counts):
    for f,a in counts.items():
        frame[f+"_other"]=a[:,0]
        frame[f+"_nectar"]=a[:,1]
        assert (frame[f+"_other"]+frame[f+"_nectar"]<=frame.total).all()
    assert sum(frame[f+"_nectar"] for f in FAMILIES).le(frame.nectar).all()
    assert sum(frame[f+"_other"] for f in FAMILIES).le(frame.other).all()

def clade_panels(fl,tab,geo,estimate):
    ans={}
    # Fixed spatial support and min-name cohort, not a search over cutoffs.
    for scope,v in fl.items():
        d=v.merge(tab[.5],on="island_id",validate="one_to_one").merge(geo,on="island_id",validate="one_to_one")
        d=d[d.analysis_regime.eq("tropical")&d.complete&d.range_area.ge(20)&d.n_species.ge(5)&d.distance_to_continent_km.ge(0)].copy()
        d["block"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
        d["log_area"]=np.log2(d.range_area);d["log_distance"]=np.log2(1+d.distance_to_continent_km)
        d["log_n"]=np.log2(d.n_species);d["abs_lat"]=d.island_latitude.abs()
        for family in FAMILIES:
            s=d.copy()
            count=s[family+"_nectar"]+s[family+"_other"]
            # Family-positive only: compare nectar and non-nectar species in
            # the same family, not continent-specific taxonomic absence.
            s=s[count.gt(0)].copy();count=count.loc[s.index]
            s["log_family_nectar"]=np.log1p(s[family+"_nectar"])
            s["log_family_other"]=np.log1p(s[family+"_other"])
            s["log_nonfamily"]=np.log1p(s.total-count)
            s=s.dropna(subset=["island_latitude","island_longitude",*CONTROLS,*OUTCOMES]).copy()
            s=s[np.isfinite(s[CONTROLS+list(OUTCOMES)].to_numpy(float)).all(1)]
            # Count as observed only within the fully diet-mapped cohort.
            assert (s[family+"_nectar"]+s[family+"_other"]<=s.total).all()
            key=scope+"_"+family
            ans[key]={"cohort":int(len(s)),"family_nectar_positive":int(s[family+"_nectar"].gt(0).sum()),
                      "family_other_positive":int(s[family+"_other"].gt(0).sum()),
                      "nonfamily_nectar_positive":int((s.nectar-s[family+"_nectar"]).gt(0).sum()),
                      "fit":estimate(s,["log_family_nectar","log_family_other"],CONTROLS)}
    return {"panels":ans,
      "warning":"Exploratory family-conditioned comparison using coarse floral states and diet-source taxonomy. Two family contexts are not independent evolutionary tests; no realm equivalence, measured visitation, species-level trait matching, or plant-genus adjustment. Same 30-block/200-informative-island gate retained."}

def self_test():
    import pandas as pd
    a=allocate(3)
    add_counts(a,"Pteropodidae",1,np.array([0,1]))
    add_counts(a,"Pteropodidae",0,np.array([1]))
    add_counts(a,"Phyllostomidae",2,np.array([2]))
    add_counts(a,"Phyllostomidae",3,np.array([0]))
    d=pd.DataFrame({"total":[2,2,1],"nectar":[1,1,1],"other":[0,1,0]})
    attach_counts(d,a)
    assert d.pteropodidae_nectar.tolist()==[1,1,0]
    assert d.pteropodidae_other.tolist()==[0,1,0]
    assert d.phyllostomidae_nectar.tolist()==[0,0,1]
    print("BAT_FAMILY_COUNT_TEST_PASS",flush=True)

if __name__=="__main__":self_test()
