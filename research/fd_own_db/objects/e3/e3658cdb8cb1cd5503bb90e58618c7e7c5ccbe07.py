"""Diagnose why the outcome-blind synonym-expanded cohort changes the effect.
No change to original estimand or support gate. Not a causal decomposition.
"""
import numpy as np
import pandas as pd

def diagnostic(f,old,new,geo,estimate):
    v=f["direct"]
    controls=["log_area","log_distance","log_n","abs_lat",*[f"climate_pc{i}" for i in range(1,5)]]
    outcomes=["white","white_open_brush","red_pink","tubular"]
    features=["log_nectar","log_other"]
    result={}
    for t in (.5,.1):
        a=old[t].set_index("island_id");b=new[t].set_index("island_id").loc[a.index]
        assert (a.total==b.total).all() and b.loc[a.complete,"complete"].all()
        assert (a.loc[a.complete,["nectar","other"]]==b.loc[a.complete,["nectar","other"]]).all().all()
        merged=b.reset_index().merge(v,on="island_id",validate="one_to_one").merge(geo,on="island_id",validate="one_to_one")
        d=merged[merged.analysis_regime.eq("tropical")&merged.complete&merged.n_species.ge(5)&merged.range_area.ge(20)&merged.distance_to_continent_km.ge(0)].copy()
        d["old_complete"]=d.island_id.isin(a.index[a.complete])
        d["log_area"]=np.log2(d.range_area);d["log_distance"]=np.log2(1+d.distance_to_continent_km);d["log_n"]=np.log2(d.n_species);d["abs_lat"]=d.island_latitude.abs()
        d["log_nectar"]=np.log1p(d.nectar);d["log_other"]=np.log1p(d.other)
        d=d.dropna(subset=controls+features+outcomes+["island_latitude","island_longitude"]).copy()
        d=d[np.isfinite(d[controls+features+outcomes].to_numpy(float)).all(1)].copy()
        d["block"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
        original=d[d.old_complete].copy()
        old_blocks=set(original.block)
        informative=set(k for k,z in original.groupby("block") if len(z)>=2 and z.log_nectar.nunique()>1)
        new_only=d[~d.old_complete]
        q={
          "original":original,
          "expanded_in_original_informative_blocks":d[d.block.isin(informative)],
          "expanded_in_all_original_blocks":d[d.block.isin(old_blocks)],
          "expanded_all":d,
          "new_islands_only":new_only}
        ans={k:estimate(z,features,controls) for k,z in q.items()}
        def desc(z):
            return {"islands":len(z),"blocks":int(z.block.nunique()),
              "white_mean":float(z.white.mean()) if len(z) else None,
              "nectar_median":float(z.nectar.median()) if len(z) else None,
              "other_median":float(z.other.median()) if len(z) else None,
              "abs_lat_median":float(z.abs_lat.median()) if len(z) else None,
              "area_median":float(z.range_area.median()) if len(z) else None}
        result[str(t)]={"original_covariate_complete":desc(original),
          "added_covariate_complete":desc(new_only),
          "added_in_original_informative_blocks":int(new_only.block.isin(informative).sum()),
          "added_in_other_original_blocks":int((new_only.block.isin(old_blocks)&~new_only.block.isin(informative)).sum()),
          "added_in_new_blocks":int((~new_only.block.isin(old_blocks)).sum()),
          "fits":ans}
    return {"panels":result,"warning":"Sequential sample composition diagnostic. Models are refitted with changing group demeaning and weights; increments are not additive causal contributions. Original complete islands retain identical diet counts. No post-hoc relaxed support gate or independent replication."}
