#!/usr/bin/env python3
"""Public historical colour × form syndrome-proxy coherence pilot.
No mainland plant pool, full-core phenotype, genus structure, or actual guild is inferred.
"""
import argparse,hashlib,io,json,math,pathlib
import numpy as np
import pandas as pd
from plant_syndrome_core import CORE,TEMPLATES,encode,components
from tropical_pollinator_source_preflight import get
from plant_syndrome_q2_original import GEO_SHA,SEED

BASE="https://raw.githubusercontent.com/zuizui0223/island/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/results/"
PINS={"direct":"6c834cb58105993b3bda10b6816198a6023a25ff4ce407eebccee5b93d5d1c7f",
      "direct_single_state":"ebb615091b971606a6626a488c44322340e6f0526c330b86b88ffee287eb9bc2"}
def sha(b):return hashlib.sha256(b).hexdigest()

def pilot_fit(d,draws=1999):
    cols=["log_distance","log_area","log_n","abs_lat",*[f"climate_pc{i}" for i in range(1,5)]]
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+["block","coherence"]).copy()
    groups=[]
    for _,g in d.groupby("block",sort=True):
        if len(g)<2 or g.log_distance.nunique()<2:continue
        x=g[cols].to_numpy(float);y=g.coherence.to_numpy(float)
        groups.append(((x-x.mean(0))/math.sqrt(len(x)),(y-y.mean())/math.sqrt(len(y)),len(y)))
    k=len(groups);n=sum(g[2] for g in groups)
    out={"eligible":len(d),"informative_groups":k,"informative_islands":n}
    def fit(w):
        j=np.flatnonzero(w)
        if not len(j) or sum(groups[i][2]-1 for i in j)-len(cols)<4:return None
        X=np.vstack([groups[i][0]*math.sqrt(w[i]) for i in j])
        Y=np.concatenate([groups[i][1]*math.sqrt(w[i]) for i in j])
        if np.linalg.matrix_rank(X)<len(cols):return None
        return float(np.linalg.lstsq(X,Y,rcond=None)[0][0])
    b=fit(np.ones(k,int))
    if b is None:return out|{"status":"not_estimable"}
    rng=np.random.default_rng(SEED)
    vals=[]
    for _ in range(draws):
        q=fit(np.bincount(rng.integers(0,k,k),minlength=k))
        if q is not None:vals.append(q)
    leave=[]
    for i in range(k):
        w=np.ones(k,int);w[i]=0;q=fit(w)
        if q is not None:leave.append(q)
    support=k>=30 and n>=200 and len(vals)>=.9*draws and len(leave)==k
    return out|{"status":"exploratory_estimable" if support else "insufficient_support",
        "slope_per_log2_distance":b,"bootstrap_interval":np.quantile(vals,[.025,.975]).tolist() if len(vals)>=.9*draws else None,
        "leave_group_range":[min(leave),max(leave)] if leave else None,
        "valid_bootstrap":len(vals),"valid_leave_group":len(leave)}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args()
    raw=get(BASE+"floral_space_20261106/all/profile_coordinates.csv",maxbytes=2_000_000)
    ref=pd.read_csv(io.BytesIO(raw));assert len(ref)==496 and ref.profile_id.is_unique
    graw=get(BASE+"geography_20260924/corrected_geography_covariates.csv",maxbytes=8_000_000)
    assert sha(graw)==GEO_SHA
    geo=pd.read_csv(io.BytesIO(graw),dtype={"island_id":str});assert len(geo)==8264 and geo.island_id.is_unique
    results={}
    for scope,pin in PINS.items():
        raw_c=get(BASE+"floral_space_20261106/"+scope+"/island_profile_carriers.csv.gz",maxbytes=2_000_000)
        assert sha(raw_c)==pin
        c=pd.read_csv(io.BytesIO(raw_c),compression="gzip",dtype={"island_id":str})
        assert not c.duplicated(["island_id","profile_id"]).any() and c.n_species.gt(0).all()
        pids=set(c.profile_id);r=ref[ref.profile_id.isin(pids)].copy().reset_index(drop=True)
        assert len(r)==len(pids)
        if scope=="direct_single_state":
            assert all(not r[v].str.contains("|",regex=False).any() for v in CORE[:2])
        idx=dict(zip(r.profile_id,r.index))
        X,pairs,weights,labels,_,_=encode(r,CORE[:2],multistate=scope=="direct")
        rows=[]
        for i,v in c.groupby("island_id",sort=True):
            n=int(v.n_species.sum())
            if n<10:continue
            ix=np.array([idx[k] for k in v.profile_id])
            q=components(X[ix],pairs,counts=v.n_species.to_numpy())
            rows.append({"island_id":i,"n":n,"coherence":float(weights@q["total"])})
        d=pd.DataFrame(rows).merge(geo,on="island_id",validate="one_to_one")
        d=d[d.area_km2.gt(0)&d.distance_to_continent_km.ge(0)].dropna(subset=["island_latitude","island_longitude","analysis_regime"]).copy()
        d["log_area"]=np.log2(d.area_km2);d["log_distance"]=np.log2(1+d.distance_to_continent_km)
        d["log_n"]=np.log2(d.n);d["abs_lat"]=d.island_latitude.abs()
        d["block"]=d.analysis_regime.astype(str)+"_"+np.floor(d.island_latitude/10).astype(int).astype(str)+"_"+np.floor(d.island_longitude/10).astype(int).astype(str)
        out={"source_carrier_islands":int(c.island_id.nunique()),"geographic_matched_n10":len(d),
             "island_equal_mean_coherence":float(d.coherence.mean()),"all":pilot_fit(d),"regions":{}}
        for k,v in d.groupby("analysis_regime",sort=True):
            out["regions"][str(k)]={"islands":len(v),"island_equal_mean_coherence":float(v.coherence.mean()),"fit":pilot_fit(v)}
        results[scope]=out
        print("PUBLIC_PLANT_Q2",scope,json.dumps(out,allow_nan=False),flush=True)
    out={"panels":results,"source_pins":PINS|{"geography":GEO_SHA},
         "primary_pilot":"direct_single_state","warning":"Only historical public colour/form profile name counts, not original 4-axis Q2 or source-relative cohort. Direct multistate is state-set membership, not measured co-expression; no plant genus, source flora, or pollination data. Isolation slope is exploratory and not proof of coherence breakdown."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
if __name__=="__main__":main()
