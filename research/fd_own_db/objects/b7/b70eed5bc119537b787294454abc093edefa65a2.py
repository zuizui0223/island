#!/usr/bin/env python3
"""Public aggregate-only floral phenotype bridge, source-matched exploratory diagnostics.
No private species data, no assignment of floral morphology to an effective pollinator.
"""
import argparse, collections, hashlib, io, json, math, os, pathlib, subprocess, zipfile
import numpy as np
import pandas as pd
from tropical_pollinator_source_preflight import get
from tropical_three_guild_candidate_audit import fetch
from tropical_bee_bird_crosswalk import decode_source,BEE_MD5,MODEL_SHA

ART=(14903203565,"b661f1cd414a63bdc4450b8c7b1691530fd23afdf8f66c405129c189587c98a6")
BASE="https://raw.githubusercontent.com/zuizui0223/island/d4a5dafab5f0737ea6c1c1c0263fc20decf6035c/results/floral_space_20261106/"
SCOPES=("direct","direct_single_state")
ENDPOINTS=("tubular_only","tubular_or_salverform","extended_access_shape","open_radial","bilabiate_or_papilionaceous","red_pink","red_pink_and_extended_shape","rao_q_without_replacement")
MASKS={
"tubular_only":{"tubular"},
"tubular_or_salverform":{"tubular","salverform"},
"extended_access_shape":{"tubular","salverform","funnel_trumpet","spurred"},
"open_radial":{"open_radial"},
"bilabiate_or_papilionaceous":{"bilabiate","papilionaceous"}}
def sha(b):return hashlib.sha256(b).hexdigest()
def read_artifact():
    aid,digest=ART
    b=subprocess.check_output(["curl","--fail","--silent","--show-error","--location","--max-time","90","--max-filesize","2000000","-H","Authorization: Bearer "+os.environ["GH_TOKEN"],"-H","Accept: application/vnd.github+json",f"https://api.github.com/repos/zuizui0223/fd/actions/artifacts/{aid}/zip"])
    assert sha(b)==digest
    out={}
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        for scope in SCOPES:
            name=scope+"_bee_public_flower_bridge.csv"
            p=[n for n in z.namelist() if pathlib.PurePosixPath(n).name==name]
            assert len(p)==1
            out[scope]=pd.read_csv(z.open(p[0]),dtype={"bee_id":str})
    assert len(out["direct"])==231 and len(out["direct_single_state"])==230
    return out
def get_public_flora(bridge):
    ref=get(BASE+"all/profile_coordinates.csv",maxbytes=2_000_000)
    profiles=pd.read_csv(io.BytesIO(ref))
    assert profiles.profile_id.is_unique
    forms=profiles.floral_form.str.split("|",regex=False).map(set)
    colors=profiles.flower_primary_color.str.split("|",regex=False).map(set)
    for name,mask in MASKS.items():profiles[name]=forms.map(lambda s:bool(s&mask))
    profiles["red_pink"]=colors.map(lambda s:"red_pink" in s)
    profiles["red_pink_and_extended_shape"]=profiles["red_pink"]&profiles["extended_access_shape"]
    out={}
    pins={"direct":"6c834cb58105993b3bda10b6816198a6023a25ff4ce407eebccee5b93d5d1c7f","direct_single_state":"ebb615091b971606a6626a488c44322340e6f0526c330b86b88ffee287eb9bc2"}
    for scope in SCOPES:
        raw=get(BASE+scope+"/island_profile_carriers.csv.gz",maxbytes=2_000_000)
        assert sha(raw)==pins[scope]
        f=pd.read_csv(io.BytesIO(raw),compression="gzip")
        f=f[f.island_id.isin(bridge[scope].island_id)].copy()
        assert not f.duplicated(["island_id","profile_id"]).any()
        tab=f.merge(profiles[["profile_id",*ENDPOINTS[:-1]]],on="profile_id",validate="many_to_one")
        for k in ENDPOINTS[:-1]:tab[k]=tab[k].astype(int)*tab.n_species
        counts=tab.groupby("island_id")[["n_species",*ENDPOINTS[:-1]]].sum()
        m=bridge[scope].merge(counts.reset_index(),on="island_id",validate="one_to_one")
        assert len(m)==len(bridge[scope]) and (m.n_names==m.n_species).all()
        for k in ENDPOINTS[:-1]:m[k]/=m.n_species
        out[scope]=m
    return out,sha(ref)
def get_bees():
    # Signed archive source; individual native species not exported or needed.
    meta=json.loads(get("https://zenodo.org/api/records/23222636"))
    f=next(f for f in meta["files"] if f["key"]=="Data.zip")
    assert f["size"]==262554015 and f["checksum"]=="md5:"+BEE_MD5
    raw=get(f["links"].get("content") or f["links"]["self"],maxbytes=300_000_000)
    assert hashlib.md5(raw).hexdigest()==BEE_MD5
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        b=z.read("Data/model_table.csv");assert sha(b)==MODEL_SHA
    d=pd.read_csv(io.StringIO(decode_source(b)),dtype={"OBJECTID_1":str}).rename(columns={"OBJECTID_1":"bee_id","Archipelago":"group","FinalRichNat":"native_bees","AreaModel":"area","DistModel":"distance","Cov":"bee_coverage","n_records":"bee_records","TypeModel":"geo_origin"})
    assert len(d)==305 and d.bee_id.is_unique
    for c in ("native_bees","area","distance","bee_coverage","bee_records","AMT","AP"):d[c]=pd.to_numeric(d[c],errors="coerce")
    return d
def design(d,features,controls,weights="group",resample=None):
    """Within-archipelago estimates; rank and residual degrees-of-freedom gates."""
    predictors=features+controls
    ds=[];ys=[];support=[]
    for group,v in d.groupby("group"):
        if len(v)<2:continue
        if resample is not None and not resample.get(group,0):continue
        x=v[predictors].to_numpy(float)
        y=v[list(ENDPOINTS)].to_numpy(float)
        if not np.isfinite(x).all() or not np.isfinite(y).all():continue
        if np.max(np.ptp(x[:,:len(features)],axis=0))<1e-12:continue
        m=1 if resample is None else resample[group]
        wt=math.sqrt(m/(len(v) if weights=="group" else 1))
        ds.append((x-x.mean(0))*wt);ys.append((y-y.mean(0))*wt);support.append(group)
    if not ds:return None
    X=np.vstack(ds);Y=np.vstack(ys);dof=len(X)-len(ds)-X.shape[1]
    if dof<4 or np.linalg.matrix_rank(X)<X.shape[1]:return None
    beta=np.linalg.lstsq(X,Y,rcond=None)[0]
    return {"beta":beta[:len(features)].T,"groups":support,"dof":dof,"n_informative":len(X)}
def fit(d,features,controls,draws=499):
    d=d.dropna(subset=features+controls+list(ENDPOINTS)).copy()
    o={"n":len(d),"archipelagos":int(d.group.nunique()),"features":features,"controls":controls,"weighting":"equal group"}
    v=design(d,features,controls)
    if v is None:return o|{"status":"insufficient_rank_or_support"}
    groups=v["groups"];o.update(status="exploratory",informative_groups=len(groups),n_informative=v["n_informative"],dof=v["dof"],prespecified_minimum_support=(len(groups)>=10 and v["n_informative"]>=30))
    # Cache each demeaned group design. The bootstrap is algebraically the
    # identical weighted least squares on replicated groups, without redoing
    # pandas groupby and centering hundreds of times.
    p=len(features)+len(controls)
    grams=[];cross=[];sizes=[];blocks=[]
    for g in groups:
        q=d[d.group==g]
        x=q[features+controls].to_numpy(float)
        y=q[list(ENDPOINTS)].to_numpy(float)
        x=(x-x.mean(0))/math.sqrt(len(q));y=(y-y.mean(0))/math.sqrt(len(q))
        grams.append(x.T@x);cross.append(x.T@y);sizes.append(len(q));blocks.append((x,y))
    grams=np.asarray(grams);cross=np.asarray(cross);sizes=np.asarray(sizes)
    def fast(m):
        # Rank via the original stacked design rather than squared condition number.
        active=np.flatnonzero(m)
        if not len(active) or int(sizes[active].sum())-len(active)-p<4:return None
        x=np.vstack([blocks[i][0]*math.sqrt(m[i]) for i in active])
        if np.linalg.matrix_rank(x)<p:return None
        y=np.vstack([blocks[i][1]*math.sqrt(m[i]) for i in active])
        return np.linalg.lstsq(x,y,rcond=None)[0][:len(features)].T
    same=fast(np.ones(len(groups),dtype=int))
    assert same is not None and np.allclose(same,v["beta"],rtol=1e-8,atol=1e-10)
    rng=np.random.default_rng(20261124);boot=[];drop=[]
    for i in range(draws):
        m=np.bincount(rng.integers(0,len(groups),len(groups)),minlength=len(groups))
        z=fast(m)
        if z is not None:boot.append(z)
    for g in range(len(groups)):
        m=np.ones(len(groups),dtype=int);m[g]=0
        z=fast(m)
        if z is not None:drop.append(z)
    rows={}
    for i,y in enumerate(ENDPOINTS):
        rows[y]={}
        for j,x in enumerate(features):
            a=float(v["beta"][i,j])
            b=[float(z[i,j]) for z in boot]
            c=[float(z[i,j]) for z in drop]
            rows[y][x]={"slope":a,"resample_q025_q975":np.quantile(b,[.025,.975]).tolist() if len(b)>=50 else None,"bootstrap_draws":len(b),"leave_one_group_min_max":[min(c),max(c)] if c else None,"same_sign_drops":int(sum(np.sign(z)==np.sign(a) for z in c)),"drop_count":len(c)}
    return o|{"endpoints":rows}
def prep(d):
    d=d.copy();d["log_area"]=np.log2(d.area);d["log_distance"]=np.log2(1+d.distance);d["log_n_names"]=np.log2(d.n_names)
    d["log_bees"]=np.log1p(d.native_bees);d["log_bee_records"]=np.log1p(d.bee_records)
    d["joint_coverage"]=pd.to_numeric(d.joint_coverage,errors="coerce")
    return d
def summarize_panels(fl,bees,bird,but):
    res={}; counts={}
    for scope,v in fl.items():
        base=v.merge(bees[["bee_id","group","native_bees","area","distance","bee_coverage","bee_records","geo_origin","AMT","AP"]],on="bee_id",validate="one_to_one")
        for nmin in (5,50):
            # Separate guild-pair estimands, no zero-imputation and no cherry-picked source masks.
            b=prep(base[(base.analysis_regime=="tropical")&(base.n_names>=nmin)&(base.area>0)&(base.distance>=0)&(base.native_bees>=0)&(base.bee_records>=0)])
            tag=f"{scope}_n{nmin}"
            counts[tag]={"bee":len(b)}
            cov=["log_area","log_distance","log_n_names"]
            res[tag+"_bee_raw"]=fit(b,["log_bees"],cov)
            res[tag+"_bee_coverage"]=fit(b,["log_bees"],cov+["bee_coverage","joint_coverage"])
            res[tag+"_bee_climate"]=fit(b,["log_bees"],cov+["AMT","AP","bee_coverage","joint_coverage"])
            # Bird provisional identity, same bee ID and separate count definitions. Negative-control nonnectar birds.
            bird_cols=["bee_id","nectar_gt0","nectar_ge20","nectar_ge50","all_birds","match_type","area_ratio"]
            x=b.merge(bird[bird_cols],on="bee_id",validate="one_to_one")
            x["log_nonnectar"]=np.log1p(x.all_birds-x.nectar_gt0)
            x["log_bird"]=np.log1p(x.nectar_gt0)
            x["log_bird20"]=np.log1p(x.nectar_ge20);x["log_bird50"]=np.log1p(x.nectar_ge50)
            counts[tag]["bird_candidates"]=len(x)
            for label,z in [("all",x),("tight",x[x.area_ratio.between(.9,1.1)&x.match_type.eq("normalized_name")])]:
                for f in ("log_bird","log_bird20","log_bird50","log_nonnectar"):
                    # Incremental bird coefficient conditional on bee count, identical covariates.
                    res[tag+"_bird_"+label+"_"+f]=fit(z,["log_bees",f],cov)
                res[tag+"_bird_"+label+"_paired_control"]=fit(z,["log_bees","log_bird","log_nonnectar"],cov)
            for mask in ("all","precision_le1km","no_gbifid_with_reference"):
                col=mask+"_records";ncol=mask+"_recorded_species"
                y=b.merge(but[["bee_id",col,ncol]],on="bee_id",validate="one_to_one")
                for minrec in (10,50):
                    q=y[(y[col]>=minrec)&(y[ncol]>0)].copy()
                    q["log_butterfly"]=np.log1p(q[ncol]);q["log_butterfly_records"]=np.log1p(q[col])
                    k=f"{tag}_butterfly_{mask}_{minrec}"
                    counts[k]={"positive_units":len(q)}
                    res[k+"_raw"]=fit(q,["log_bees","log_butterfly"],cov)
                    res[k+"_effort"]=fit(q,["log_bees","log_butterfly"],cov+["log_butterfly_records"])
    return res,counts
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    br=read_artifact();fl,profile_sha=get_public_flora(br)
    bird=fetch("bee_bird");but=fetch("butterfly");bees=get_bees()
    res,counts=summarize_panels(fl,bees,bird,but)
    r={"artifact_pin":ART,"profile_sha256":profile_sha,"counts":counts,"panels":res,"warning":"Public historical PR321 colour/form category state-set aggregates only, not original 4505 radial symmetry x deep tube. Shape indicators mean documented category membership, NOT tube length, realised guild or a morph-level colour/form observation. Bird matches provisional; butterfly count positive records only. All tests exploratory, many dependent panels, percentile resampling not calibrated inference, no causality or guild-specific selection."}
    (a.out/"public_flower_guild_association.json").write_text(json.dumps(r,indent=2,allow_nan=False))
    # Short log of primary rows and all matched cohort counts. Full cells in artifact.
    primary={k:v for k,v in res.items() if k.startswith("direct_n50_") and ("_bee_" in k or "_bird_all_" in k or "_butterfly_" in k)}
    print("FLORAL_GUILD_COUNTS",json.dumps(counts),flush=True)
    for k,v in primary.items():
        short={x:y for x,y in v.items() if x!="endpoints"}
        short["endpoints"]={x:y for x,y in v.get("endpoints",{}).items() if x in ("extended_access_shape","tubular_only","red_pink_and_extended_shape","rao_q_without_replacement")}
        print("FLORAL_GUILD_PANEL",k,json.dumps(short),flush=True)
if __name__=="__main__":main()
