"""Post-pilot missing-diet sensitivity. Never relabel unknown species as observed.
Keep a common island/block cohort across three explicitly hypothetical codings.
No positive mechanism claim and no formal identification bounds from extremes.
"""
import numpy as np
import pandas as pd

FEATURES=("known","unknown_as_other","unknown_as_nectar")
OUTCOMES=("white","white_open_brush","red_pink","tubular")
SEED=20261125
CONTROLS=["log_area","log_distance","log_n","abs_lat",*[f"climate_pc{i}" for i in range(1,5)],"log_unknown"]

def fit_common(d,draws=1999):
    # Unknown count is retained as a distinct covariate in every scenario.
    # Counts under each extreme are coherent across all islands.
    base={"complete_covariate_islands":len(d),"islands_with_unknown":int(d.unknown.gt(0).sum())}
    if d.empty:return base|{"status":"not_estimable"}
    assert (d[["nectar","other","unknown"]]>=0).all().all()
    x0=d[CONTROLS].to_numpy(float)
    f=np.stack([np.log1p(d[["nectar","other"]].to_numpy(float)),
                np.log1p(np.column_stack([d.nectar,d.other+d.unknown])),
                np.log1p(np.column_stack([d.nectar+d.unknown,d.other]))],axis=0)
    yy=d[list(OUTCOMES)].to_numpy(float)
    blocks=[]
    for _,v in d.groupby("block",sort=True):
        ix=v.index.to_numpy()
        if len(ix)<2 or not all(np.ptp(f[j,ix,0])>0 for j in range(3)):continue
        xx=np.concatenate([f[:,ix,:],np.broadcast_to(x0[ix],(3,len(ix),x0.shape[1]))],axis=2)
        y=yy[ix]
        xx=(xx-xx.mean(1,keepdims=True))/np.sqrt(len(ix));y=(y-y.mean(0))/np.sqrt(len(ix))
        blocks.append((xx,y,len(ix)))
    k=len(blocks);n=sum(v[2] for v in blocks)
    base|={"common_informative_blocks":k,"common_informative_islands":n}
    def solve(w):
        active=np.flatnonzero(w);p=2+len(CONTROLS)
        if not len(active) or sum(blocks[i][2] for i in active)-len(active)-p<4:return None
        X=np.concatenate([blocks[i][0]*np.sqrt(w[i]) for i in active],axis=1)
        Y=np.vstack([blocks[i][1]*np.sqrt(w[i]) for i in active])
        if any(np.linalg.matrix_rank(x)<p for x in X):return None
        B=np.stack([np.linalg.lstsq(x,Y,rcond=None)[0] for x in X])
        return B[:,0]-B[:,1]
    obs=solve(np.ones(k,dtype=int))
    if obs is None:return base|{"status":"not_estimable"}
    rng=np.random.default_rng(SEED)
    boot=[]
    for _ in range(draws):
        z=solve(np.bincount(rng.integers(0,k,k),minlength=k))
        if z is not None:boot.append(z)
    leave=[]
    for i in range(k):
        w=np.ones(k,dtype=int);w[i]=0;z=solve(w)
        if z is not None:leave.append(z)
    gate=k>=30 and n>=200 and len(boot)>=.9*draws and len(leave)==k
    ans={}
    for j,o in enumerate(OUTCOMES):
        ans[o]={}
        for a,s in enumerate(FEATURES):
            lv=np.asarray(leave)[:,a,j] if leave else []
            ans[o][s]={"contrast":float(obs[a,j]),
                "paired_group_bootstrap_interval":np.quantile(np.asarray(boot)[:,a,j],[.025,.975]).tolist() if len(boot)>=.9*draws else None,
                "leave_block_range":[float(min(lv)),float(max(lv))] if len(lv) else None}
        ans[o]["all_scenarios_same_positive_point_sign"]=bool((obs[:,j]>0).all())
        ans[o]["all_scenarios_positive_bootstrap_lower"]=bool(len(boot)>=.9*draws and (np.quantile(np.asarray(boot)[:,:,j],.025,axis=0)>0).all())
    return base|{"status":"exploratory_estimable" if gate else "insufficient_support",
                 "valid_bootstrap":len(boot),"valid_leave_block":len(leave),"endpoints":ans}

def missingness_panels(f,tab,geo):
    result={}
    # Keep primary spatial resolution and carrier floor fixed; this is a
    # focused missingness audit, not another threshold search.
    for scope,v in f.items():
        for overlap,c in tab.items():
            if overlap!=.5:continue
            d=v.merge(c,on="island_id",validate="one_to_one").merge(geo,on="island_id",validate="one_to_one")
            d=d[d.analysis_regime.eq("tropical")&d.total.gt(0)&d.range_area.gt(0)&d.distance_to_continent_km.ge(0)].copy()
            d["block"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str))
            d["log_area"]=np.log2(d.range_area);d["log_distance"]=np.log2(1+d.distance_to_continent_km)
            d["log_n"]=np.log2(d.n_species);d["abs_lat"]=d.island_latitude.abs();d["log_unknown"]=np.log1p(d.unknown)
            for amin in (20,100):
                for nmin in (5,):
                    q=d[d.range_area.ge(amin)&d.n_species.ge(nmin)].copy()
                    for coverage in (.5,.8,.95):
                        a=q[q.match_share.ge(coverage)].copy()
                        b=a.dropna(subset=CONTROLS+list(OUTCOMES)+["nectar","other","unknown","island_latitude","island_longitude"]).reset_index(drop=True)
                        b=b[np.isfinite(b[CONTROLS+list(OUTCOMES)].to_numpy(float)).all(1)].reset_index(drop=True)
                        name=f"{scope}_cover{overlap}_area{amin}_n{nmin}_diet{coverage}"
                        result[name]={"before_coverage":len(q),"after_coverage":len(a),
                          "original_complete":int(q.complete.sum()),
                          "fit":fit_common(b)}
    return {"panels":result,"warning":"Post-pilot sensitivity. Known counts do not represent complete diet. Extreme assignments are hypothetical, not coefficient identification bounds; intermediate species-specific assignments can produce values outside their range. Common cohort and paired bootstrap across all scenarios; no independent replication or pollination claim."}

def self_test():
    rng=np.random.default_rng(SEED);n=300
    d=pd.DataFrame({"block":np.repeat(np.arange(30),10),"nectar":rng.integers(0,8,n),"other":rng.integers(0,10,n),"unknown":rng.integers(0,5,n)})
    for c in CONTROLS:d[c]=rng.normal(size=n)
    d["log_unknown"]=np.log1p(d.unknown)
    for o in OUTCOMES:d[o]=.5*np.log1p(d.nectar)-.1*np.log1p(d.other)+.01*rng.normal(size=n)
    r=fit_common(d,draws=99)
    assert r["status"]=="exploratory_estimable" and abs(r["endpoints"]["white"]["known"]["contrast"]-.6)<.01
    d["nectar"]=0
    assert fit_common(d,draws=99)["status"]=="not_estimable"
    print("BAT_UNKNOWN_DIET_SENSITIVITY_TEST_PASS",flush=True)
if __name__=="__main__":self_test()
