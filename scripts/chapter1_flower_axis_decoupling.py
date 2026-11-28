#!/usr/bin/env python3
"""Post-hoc own-DB H1 extension: paired floral-axis isolation contrasts.
No official GIFT island or mainland flora, no animal data. Anonymous aggregates only.
"""
import argparse, hashlib, json, math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import t

LOCK={"flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
      "ledger":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
      "geography":"08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7"}
REGIONS=("northern_midlatitude","northern_high_latitude","tropical","southern_extratropical")
DIST="log_distance_to_continent_km"
COV=("log_island_area_km2","climate_pc1","climate_pc2","climate_pc3","climate_pc4")
TRAITS=("floral_symmetry","tube_depth_class")
METRICS=("radial","shallow_open","radial_deep","bilateral_deep","radial_shallow","bilateral_shallow","axis_difference","conditional_difference")
SCOPES=("all_recorded","literal_native","high_only")
ENDPOINTS=("axis_difference","conditional_difference")
SEED=20261125
DRAWS=499

def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def trait_pairs(ledger, high_only=False):
    z=ledger.loc[ledger.trait_name.isin(TRAITS)&ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(["high"] if high_only else ["high","medium"])&
        ledger.normalized_value.notna(),["accepted_species","trait_name","normalized_value"]].copy()
    z["normalized_value"]=z.normalized_value.str.strip()
    z=z[(z.trait_name.eq(TRAITS[0])&z.normalized_value.isin(["actinomorphic","zygomorphic"]))|
        (z.trait_name.eq(TRAITS[1])&z.normalized_value.isin(["absent_or_open","shallow","deep"]))]
    if z.duplicated(["accepted_species","trait_name"]).any():raise ValueError("duplicate retained species-trait evidence")
    v=z.pivot(index="accepted_species",columns="trait_name",values="normalized_value").dropna(subset=list(TRAITS))
    # Genus token is a grouping proxy, not a phylogenetic tree.
    genus=v.index.to_series().str.extract(r"^([A-Z][A-Za-z-]+)\s+[a-z][A-Za-z-]+",expand=False)
    valid=genus.notna();v=v.loc[valid].copy()
    v["genus"]=genus.loc[valid];v["A"]=v[TRAITS[0]].eq("actinomorphic").astype(int)
    v["O"]=v[TRAITS[1]].isin(["absent_or_open","shallow"]).astype(int)
    return v

def cells(A,O):
    return np.column_stack([A*(1-O),(1-A)*(1-O),A*O,(1-A)*O]).astype(float)

def rates(c):
    n=c.sum(1);deep=c[:,0]+c[:,1];shallow=c[:,2]+c[:,3]
    p=np.divide(c,n[:,None],out=np.full_like(c,np.nan),where=n[:,None]>0)
    a=p[:,0]+p[:,2];o=p[:,2]+p[:,3]
    cond=np.divide(c[:,0],deep,out=np.full(len(c),np.nan),where=deep>0)-np.divide(c[:,2],shallow,out=np.full(len(c),np.nan),where=shallow>0)
    return np.column_stack([a,o,p,a-o,cond])

def design(d, regional=False):
    cols=[DIST,*COV,"log_carriers","coverage"]
    X=d[cols].to_numpy(float)
    # Raw H1 distance transform, same unit across regions; scale other columns
    # solely for numerical stability. Constant nuisance terms fail closed.
    sd=X.std(0)
    if np.any(sd<=1e-12):return None,None
    X=(X-X.mean(0))/sd;X[:,0]*=sd[0]
    X=np.column_stack([np.ones(len(X)),X])
    if regional:
        z=d.analysis_regime.eq("tropical").to_numpy(float)
        X=np.column_stack([X,X*z[:,None]])
        j=X.shape[1]//2+1
    else:j=1
    return X,j

def estimate(X,Y,blocks,j):
    n,k=X.shape;g=len(set(blocks))
    if n<50 or g<10 or n<=k or g<=k or np.linalg.matrix_rank(X)<k:return None
    inv=np.linalg.inv(X.T@X);B=inv@X.T@Y;E=Y-X@B
    # A single paired response is the direct slope difference, not a
    # difference of two separate significance tests.
    infl=(X@inv[:,j])[:,None]*E
    sums=np.stack([infl[blocks==u].sum(0) for u in sorted(set(blocks))])
    se=np.sqrt(np.sum(sums*sums,axis=0)*g/(g-1)*(n-1)/(n-k))
    b=B[j];crit=t.ppf(.975,g-1)
    out={}
    for i in range(Y.shape[1]):
        if se[i]<=1e-14:out[str(i)]={"beta":float(b[i]),"status":"degenerate_se","p":1.}
        else:out[str(i)]={"beta":float(b[i]),"se":float(se[i]),"interval":[float(b[i]-crit*se[i]),float(b[i]+crit*se[i])],"p":float(2*t.sf(abs(b[i]/se[i]),g-1)),"status":"fit"}
    return out

def bootstrap(X,j,M,C,gcode,blocks,metric,seed):
    # Coherent positive genus weights shared across all occurrences + common
    # spatial-block resample. Fixed original cohort and covariates.
    # Sensitivity distribution, NOT a calibrated confidence interval.
    rng=np.random.default_rng(seed)
    _,bi=np.unique(blocks,return_inverse=True);ng=int(gcode.max())+1;nb=int(bi.max())+1
    draws=[];geo_only=[]
    obs=rates(np.asarray(M@C))[:,metric]
    for _ in range(DRAWS):
        gw=rng.exponential(1,ng);w=gw[gcode]
        r=rates(np.asarray(M@(C*w[:,None])))[:,metric]
        bw=np.bincount(rng.integers(nb,size=nb),minlength=nb)[bi]
        x=X*np.sqrt(bw[:,None])
        if np.linalg.matrix_rank(x)<X.shape[1] or np.any(~np.isfinite(r)):continue
        draw=np.linalg.lstsq(x,np.column_stack([r,obs])*np.sqrt(bw[:,None]),rcond=None)[0][j]
        draws.append(float(draw[0]));geo_only.append(float(draw[1]))
    def report(x):
        return {"valid_draws":len(x),"range_025_975":np.quantile(x,[.025,.975]).tolist() if len(x)>=.9*DRAWS else None}
    return {"joint_genus_geography":report(draws),"geography_only":report(geo_only),"not_calibrated_ci":True}

def holm(rows):
    order=sorted(range(len(rows)),key=lambda i:rows[i].get("p",1.))
    last=0.
    for j,i in enumerate(order):
        last=max(last,min(1.,rows[i].get("p",1.)*(len(rows)-j)))
        rows[i]["holm_p"]=last

def make_frame(flora,v,geo,scope):
    f=flora if scope!="literal_native" else flora.loc[flora.origin_status.eq("native")]
    f=f[["island_id","accepted_species"]]
    species=v.index.to_list();idx=pd.Series(np.arange(len(species)),index=species)
    j=f.merge(idx.rename("si"),left_on="accepted_species",right_index=True,validate="many_to_one")
    islandids=sorted(j.island_id.unique());ix=pd.Series(np.arange(len(islandids)),index=islandids)
    M=sparse.csr_matrix((np.ones(len(j)),(j.island_id.map(ix).to_numpy(),j.si.to_numpy())),shape=(len(islandids),len(species)))
    if M.nnz!=len(j) or np.any(M.data!=1):raise ValueError("duplicate incidence")
    C=cells(v.A.to_numpy(),v.O.to_numpy());counts=np.asarray(M@C);n=counts.sum(1)
    d=pd.DataFrame({"island_id":islandids,"mi":np.arange(len(islandids)),"carrier_n":n,
        "deep_n":counts[:,:2].sum(1),"shallow_n":counts[:,2:].sum(1)})
    d["coverage"]=d.carrier_n/d.island_id.map(f.groupby("island_id").size())
    d["log_carriers"]=np.log(d.carrier_n)
    d=d.merge(geo[["island_id","analysis_regime","spatial_block",DIST,*COV]],on="island_id",validate="one_to_one")
    for k in (DIST,*COV):d[k]=pd.to_numeric(d[k],errors="coerce")
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=[DIST,*COV,"spatial_block","coverage","analysis_regime"])
    d=d.loc[d.analysis_regime.isin(REGIONS)&d.spatial_block.astype(str).str.strip().ne("")]
    gcode=pd.factorize(v.genus,sort=True)[0]
    return d,M,C,gcode

def fit_panel(d,M,C,gcode,scope,endpoint,region,nmin,regional=False):
    metric=METRICS.index(endpoint);d=d.loc[d.carrier_n.ge(nmin)].copy()
    if endpoint=="conditional_difference":d=d.loc[d.deep_n.ge(5)&d.shallow_n.ge(5)]
    support_window=None
    if regional:
        d=d.loc[d.analysis_regime.isin(["tropical","northern_midlatitude"])]
        ss=[d.loc[d.analysis_regime.eq(k),DIST] for k in ["tropical","northern_midlatitude"]]
        if min(map(len,ss))<50:return {"scope":scope,"endpoint":endpoint,"region":region,"nmin":nmin,"n":len(d),"status":"insufficient_common_support","p":1.}
        lo=max(x.quantile(.05) for x in ss);hi=min(x.quantile(.95) for x in ss)
        support_window=[float(lo),float(hi)]
        d=d.loc[d[DIST].between(lo,hi)]
    else:d=d.loc[d.analysis_regime.eq(region)]
    d=d.sort_values("mi")
    b=d.spatial_block.to_numpy(str);n=len(d);g=len(set(b))
    out={"scope":scope,"endpoint":endpoint,"region":region,"nmin":nmin,"n":n,"blocks":g,"p":1.,"status":"insufficient_support","distance_common_support":support_window}
    if regional:
        rr={k:{"n":int((d.analysis_regime==k).sum()),"blocks":d.loc[d.analysis_regime.eq(k),"spatial_block"].nunique()} for k in ["tropical","northern_midlatitude"]}
        out["regional_support"]=rr
        if any(v["n"]<50 or v["blocks"]<10 for v in rr.values()):return out
    if n<50 or g<10:return out
    X,j=design(d,regional)
    if X is None:return out|{"status":"constant_predictor"}
    mat=M[d.mi.to_numpy(int)]
    Y=rates(np.asarray(mat@C))
    # All seven unconditional metrics share exactly this roster. Conditional
    # metric appears only where both depth subsets are adequately represented.
    use=list(range(7))+([7] if endpoint=="conditional_difference" else [])
    y=Y[:,use]
    if np.any(~np.isfinite(y)):raise ValueError("nonfinite supported response")
    fit=estimate(X,y,b,j)
    if fit is None:return out|{"status":"rank_or_df_failure"}
    metrics={METRICS[k]:fit[str(i)] for i,k in enumerate(use)}
    assert abs(metrics["axis_difference"]["beta"]-(metrics["radial"]["beta"]-metrics["shallow_open"]["beta"]))<1e-10
    assert abs(sum(metrics[k]["beta"] for k in METRICS[2:6]))<1e-10
    out.update(metrics[endpoint]);out["metrics"]=metrics
    out["mean_carriers"]=float(d.carrier_n.mean());out["median_coverage"]=float(d.coverage.median())
    out["distinct_scored_species"]=int(np.count_nonzero(np.asarray(mat.sum(0)).ravel()))
    # Pre-fixed sensitivity for every estimable primary and regional test.
    if nmin==10 and out["status"]=="fit":
        out["weight_sensitivity"]=bootstrap(X,j,mat,C,gcode,b,metric,SEED)
    return out

def run(flora,ledger,geo):
    assert not flora.duplicated(["island_id","accepted_species"]).any()
    assert geo.island_id.is_unique
    results=[];supports=[]
    for scope in SCOPES:
        v=trait_pairs(ledger,scope=="high_only")
        d,M,C,gcode=make_frame(flora,v,geo,scope)
        supports.append({"scope":scope,"trait_pair_species":len(v),"islands_matched_geography":len(d),"islands_n10":int(d.carrier_n.ge(10).sum()),"islands_two_depth_groups_n5":int((d.deep_n.ge(5)&d.shallow_n.ge(5)).sum())})
        for nmin in [10,30]:
            local=[];cross=[]
            for endpoint in ENDPOINTS:
                for region in REGIONS:local.append(fit_panel(d,M,C,gcode,scope,endpoint,region,nmin))
                cross.append(fit_panel(d,M,C,gcode,scope,endpoint,"tropical_minus_north_mid",nmin,True))
            holm(local);holm(cross);results.extend(local+cross)
    assert len(results)==60
    return {"owner":"island","sources":LOCK,"supports":supports,"panels":results,
        "primary_family":"all_recorded n>=10: eight regional paired contrasts, Holm including unsupported",
        "secondary_family":"two tropical-minus-northern-midlatitude contrasts on common 5–95% H1 log-distance support, separate Holm per scope/threshold",
        "warning":"Post-hoc outcome-motivated cross-sectional assemblage analysis. Direct High/Medium single-state H1-compatible pair carriers. Equal-island OLS share contrast vs raw H1 log distance, area, climate, coverage and carrier count, CR1 spatial cluster t. Coherent genus multiplier+geography sensitivity not calibrated CI. No named pollinator, measured access, within-species change, or independent replication."}

def self_test():
    a=np.array([1,0,1,0]);o=np.array([0,0,1,1])
    c=cells(a,o);v=rates(c.sum(0,keepdims=True))[0]
    assert np.allclose(v,[.5,.5,.25,.25,.25,.25,0,0])
    x=rates(np.array([[3.,1.,1.,3.]]))[0]
    assert abs(x[7]-.5)<1e-12
    assert np.isclose(x[6],x[2]-x[5])
    rng=np.random.default_rng(10);n=240
    xx=rng.normal(size=n);zz=rng.normal(size=n);X=np.column_stack([np.ones(n),xx,zz]);b=np.repeat(np.arange(30),8)
    y=np.column_stack([.2+.03*xx+rng.normal(0,.1,n),.3-.02*xx+rng.normal(0,.1,n)])
    y=np.column_stack([y,y[:,0]-y[:,1]]);r=estimate(X,y,b,1)
    assert abs(r["2"]["beta"]-(r["0"]["beta"]-r["1"]["beta"]))<1e-12
    rows=[{"p":.01},{"p":.03},{"p":1.}];holm(rows);assert np.allclose([v["holm_p"] for v in rows],[.03,.06,1])
    print("FLORAL_DECOUPLING_SELF_TEST_PASS")

def main():
    p=argparse.ArgumentParser()
    for k in ("flora","ledger","geography","out"):p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    for k,h in LOCK.items():
        if digest(getattr(a,k))!=h:raise ValueError("source hash mismatch "+k)
    flora=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    ledger=pd.read_csv(a.ledger,dtype=str)
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    assert len(flora)==1039757 and flora.island_id.nunique()==4505 and len(geo)==8264
    z=run(flora,ledger,geo)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("FLORAL_DECOUPLING_SUPPORT "+json.dumps(z["supports"],allow_nan=False),flush=True)
    for v in z["panels"]:print("FLORAL_DECOUPLING_PANEL "+json.dumps(v,allow_nan=False),flush=True)
    print("FLORAL_DECOUPLING_COMPLETE 60",flush=True)

if __name__=="__main__":main()
