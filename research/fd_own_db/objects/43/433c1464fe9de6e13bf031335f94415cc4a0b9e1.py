#!/usr/bin/env python3
"""Frozen exploratory Bombus-specific test, using existing geometry-certified public flower bridge.
No original protected GIFT individual records; no checklist absence outside the source cohort.
"""
import argparse, hashlib, io, json, pathlib, time, zipfile
import numpy as np
import pandas as pd
from tropical_public_flora_guild_association import read_artifact, get_public_flora
from tropical_pollinator_source_preflight import get
from tropical_bee_bird_crosswalk import BEE_MD5, MODEL_SHA, SPECIES_SHA, decode_source

OUTCOMES=("bilabiate_or_papilionaceous","open_radial","tubular_only","rao_q_without_replacement")
CONTROLS=("log_area","log_distance","log_n","bee_coverage","joint_coverage","AMT","AP")
SEED=20261125

def source():
    meta=json.loads(get("https://zenodo.org/api/records/23222636"))
    f=next(x for x in meta["files"] if x["key"]=="Data.zip")
    assert f["size"]==262554015 and f["checksum"]=="md5:"+BEE_MD5
    for trial in range(3):
        try:
            raw=get(f["links"].get("content") or f["links"]["self"],maxbytes=300_000_000)
            assert len(raw)==f["size"] and hashlib.md5(raw).hexdigest()==BEE_MD5
            break
        except Exception:
            if trial==2: raise
            time.sleep(5)
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        mb=z.read("Data/model_table.csv"); ob=z.read("Data/islandspecieslist.csv")
    assert hashlib.sha256(mb).hexdigest()==MODEL_SHA
    assert hashlib.sha256(ob).hexdigest()==SPECIES_SHA
    model=pd.read_csv(io.StringIO(decode_source(mb)),dtype={"OBJECTID_1":str})
    occ=pd.read_csv(io.StringIO(decode_source(ob)),dtype=str,keep_default_na=False)
    assert len(model)==305 and model.OBJECTID_1.is_unique
    assert {"OBJECTID_1","Corrected","NonNative","Uncertain"}<=set(occ)
    occ["name"]=occ.Corrected.str.strip().str.replace(r"\s+"," ",regex=True)
    occ=occ[occ.name.ne("") & occ.OBJECTID_1.isin(model.OBJECTID_1)].copy()
    occ["bombus"]=occ.name.str.split().str[0].str.casefold().eq("bombus")
    result=model.rename(columns={"OBJECTID_1":"bee_id","Archipelago":"group","FinalRichNat":"model_native","AreaModel":"area","DistModel":"distance","Cov":"bee_coverage"})
    for tag,sub in (("native",occ[occ.NonNative.str.strip().eq("0")]),
                    ("certain",occ[occ.NonNative.str.strip().eq("0")&occ.Uncertain.str.strip().eq("0")])):
        # Distinct name units, never rows; retain source non-Bombus controls.
        sub=sub.drop_duplicates(["OBJECTID_1","name"])
        sums=sub.groupby("OBJECTID_1").bombus.agg(["sum","count"])
        result[tag+"_bombus"]=result.bee_id.map(sums["sum"]).fillna(0).astype(int)
        result[tag+"_other"]=result.bee_id.map(sums["count"]-sums["sum"]).fillna(0).astype(int)
    result["reconciled"]=(result.native_bombus+result.native_other)==pd.to_numeric(result.model_native)
    assert int(result.reconciled.sum())==304
    summary={"model_islands":len(result),"reconciled":int(result.reconciled.sum()),
             "native_name_pairs":int((result.native_bombus+result.native_other).sum()),
             "bombus_positive_source_units":int(result.native_bombus.gt(0).sum()),
             "bombus_native_name_pairs":int(result.native_bombus.sum()),
             "bombus_certain_name_pairs":int(result.certain_bombus.sum())}
    assert summary["native_name_pairs"]==11064
    cols=["bee_id","group","model_native","area","distance","bee_coverage","AMT","AP","reconciled","native_bombus","native_other","certain_bombus","certain_other"]
    return result[cols],summary

def estimate(d,controls,draws=1999):
    cols=["log_bombus","log_other",*controls]
    d=d.dropna(subset=["group",*cols,*OUTCOMES]).copy()
    blocks=[]
    for _,g in d.groupby("group",sort=True):
        # A group must have within-group variation in Bombus, not merely total bees.
        if len(g)<2 or g.log_bombus.nunique()<2: continue
        x=g[cols].to_numpy(float); y=g[list(OUTCOMES)].to_numpy(float)
        if not (np.isfinite(x).all() and np.isfinite(y).all()): continue
        blocks.append(((x-x.mean(0))/np.sqrt(len(x)),(y-y.mean(0))/np.sqrt(len(x)),len(x)))
    k=len(blocks); n=sum(b[2] for b in blocks)
    base={"eligible_islands":len(d),"informative_groups":k,"informative_islands":n,
          "bombus_positive_eligible":int(d.log_bombus.gt(0).sum()),
          "controls":list(controls)}
    def solve(w):
        active=np.flatnonzero(w)
        if not len(active) or sum(blocks[i][2] for i in active)-len(active)-len(cols)<4:return None
        x=np.vstack([blocks[i][0]*np.sqrt(w[i]) for i in active])
        y=np.vstack([blocks[i][1]*np.sqrt(w[i]) for i in active])
        if np.linalg.matrix_rank(x)<len(cols):return None
        b=np.linalg.lstsq(x,y,rcond=None)[0]
        return b[0]-b[1]
    v=solve(np.ones(k,dtype=int))
    if v is None:return base|{"status":"not_estimable"}
    rng=np.random.default_rng(SEED); boot=[]; leave=[]
    for _ in range(draws):
        w=np.bincount(rng.integers(0,k,k),minlength=k)
        z=solve(w)
        if z is not None:boot.append(z)
    for i in range(k):
        w=np.ones(k,dtype=int);w[i]=0
        z=solve(w)
        if z is not None:leave.append(z)
    out={}
    for i,e in enumerate(OUTCOMES):
        ci=np.quantile(np.asarray(boot)[:,i],[.025,.975]).tolist() if len(boot)>=.9*draws else None
        lv=np.asarray(leave)[:,i].tolist() if leave else []
        out[e]={"contrast_bombus_minus_other":float(v[i]),"group_bootstrap_interval":ci,
                "leave_group_range":[min(lv),max(lv)] if lv else None}
    gate=k>=10 and n>=30 and len(boot)>=.9*draws and len(leave)==k
    p=out[OUTCOMES[0]]
    supported=gate and p["group_bootstrap_interval"][0]>0 and p["leave_group_range"][0]>0
    return base|{"status":"exploratory_support" if supported else "exploratory_not_supported" if gate else "insufficient_support",
                 "bootstrap_success":len(boot),"leave_group_success":len(leave),"endpoints":out}

def panels(fl,bees,draws=1999):
    out={}
    for scope,frame in fl.items():
        base=frame.merge(bees,on="bee_id",validate="one_to_one")
        # Exclude unreconciled model/source island and zero native checklist rows;
        # zero Bombus is meaningful only inside a positive, reconciled bee list.
        base=base[base.reconciled & (base.native_bombus+base.native_other).gt(0)].copy()
        for c in ("area","distance","n_names","bee_coverage","joint_coverage","AMT","AP"):
            base[c]=pd.to_numeric(base[c],errors="coerce")
        base=base[base.area.gt(0)&base.distance.ge(0)&base.n_names.gt(0)].copy()
        base["log_area"]=np.log2(base.area);base["log_distance"]=np.log2(1+base.distance);base["log_n"]=np.log2(base.n_names)
        # Use all exact matched islands for primary; region splits predeclared, never choose by outcomes.
        regions={"all":base,**{str(k):v for k,v in base.groupby("analysis_regime",sort=True)}}
        for rg,region in regions.items():
            for nmin in (5,50):
                for quality in ("native","certain"):
                    d=region[region.n_names.ge(nmin)].copy()
                    d=d[(d[quality+"_bombus"]+d[quality+"_other"]).gt(0)].copy()
                    d["log_bombus"]=np.log1p(d[quality+"_bombus"])
                    d["log_other"]=np.log1p(d[quality+"_other"])
                    key=f"{scope}_{rg}_n{nmin}_{quality}"
                    out[key]=estimate(d,CONTROLS,draws)
    return out

def self_test():
    rng=np.random.default_rng(SEED)
    n=96
    d=pd.DataFrame({"group":np.repeat(np.arange(12),8),
                    "log_bombus":np.tile(np.arange(8)/7,12),
                    "log_other":rng.normal(size=n)})
    for c in CONTROLS:d[c]=rng.normal(size=n)
    for o in OUTCOMES:d[o]=.5*d.log_bombus-.1*d.log_other+rng.normal(0,.001,n)
    v=estimate(d,CONTROLS,draws=199)
    assert v["status"]=="exploratory_support" and v["informative_groups"]==12
    assert abs(v["endpoints"][OUTCOMES[0]]["contrast_bombus_minus_other"]-.6)<.01
    d["log_bombus"]=0.
    assert estimate(d,CONTROLS,draws=199)["status"]=="not_estimable"
    print("SYNTHETIC_ESTIMATOR_TEST_PASS",flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path);p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test()
        return
    if a.out is None:p.error("--out is required for source run")
    a.out.mkdir(parents=True,exist_ok=True)
    fl,profile_sha=get_public_flora(read_artifact())
    bees,counts=source()
    print("BOMBUS_SOURCE",json.dumps(counts),flush=True)
    result=panels(fl,bees)
    primary="direct_all_n5_native"
    r={"source":counts,"profile_sha256":profile_sha,"primary_panel":primary,
       "primary_endpoint":OUTCOMES[0],"panels":result,
       "warning":"Exploratory, source-positive bee checklist cohort; Bombus genus membership not effective pollination. Bilabiate/papilionaceous historical name-unit state-set membership is not actual symmetry or tube depth. No plant genus-conditioned test or historical shift. Other endpoints/regions/scopes are sensitivity only."}
    (a.out/"bombus_public_flora.json").write_text(json.dumps(r,indent=2,allow_nan=False))
    for k,v in result.items():print("BOMBUS_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
