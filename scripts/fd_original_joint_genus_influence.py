"""Protected original GIFT: exact finite-sample joint Rao partition.
Descriptive island metrics, geographic association conditional on covariates.
No island/species identifiers exported; no evolutionary or pollinator inference.
"""
import argparse, io, json, pathlib
import numpy as np
import pandas as pd
from fd_original_tropical_pollinator_bridge import read_pin, PIN, accepted_genus

PAIRS=[("floral_form","floral_symmetry")]
SEED=20261125

def partition(a,b,g):
    n=len(a)
    if n<2: raise ValueError("need two species")
    _,a=np.unique(a,return_inverse=True);_,b=np.unique(b,return_inverse=True)
    A=a.max()+1;B=b.max()+1
    joint=np.bincount(a*B+b,minlength=A*B).reshape(A,B)
    ca=joint.sum(1);cb=joint.sum(0);den=n*(n-1)
    collision=lambda c: np.sum(c*(c-1))
    q=1-collision(joint)/den
    qa=1-collision(ca)/den;qb=1-collision(cb)/den
    q0=1-(1-qa)*(1-qb)
    # Expected ordered matching pairs under independent within-stratum
    # permutation of B, conditional on both local marginal tables.
    expected=0.; prev=np.zeros((A,B))
    for k in np.unique(g):
        m=g==k;h=int(m.sum())
        x=np.bincount(a[m],minlength=A);y=np.bincount(b[m],minlength=B)
        if h>1: expected+=collision(x)*collision(y)/(h*(h-1))
        z=np.outer(x,y)/h
        expected+=2*np.sum(prev*z);prev+=z
    qg=1-expected/den
    return np.array([q,qa,qb,q0,qg-q0,q-qg,q-q0])

NAMES=[]

def cell_metrics(h,x,y,cells):
    n=len(h)
    obs=h.groupby([x,y]).size().to_dict()
    ax=h[x].value_counts().to_dict();by=h[y].value_counts().to_dict()
    expected={}
    for _,g in h.groupby("stratum"):
        for a,na in g[x].value_counts().items():
            for b,nb in g[y].value_counts().items():
                expected[a,b]=expected.get((a,b),0)+na*nb/len(g)
    out=[]
    for a,b in cells:
        o=obs.get((a,b),0)/n; e=expected.get((a,b),0)/n
        out.extend([o,ax.get(a,0)*by.get(b,0)/n**2,e,o-e])
    z=np.array(out).reshape(-1,4)
    assert np.allclose(z.sum(0),[1,1,1,0])
    return z.ravel()

def self_test():
    import itertools
    a=np.array([0,0,1,1,0,1]);b=np.array([0,1,0,1,1,1]);g=np.array([0,0,0,0,1,1])
    v=partition(a,b,g)
    vals=[]
    for x in itertools.permutations(b[:4]):
        for y in itertools.permutations(b[4:]):
            vals.append(partition(a,np.array(x+y),g)[0])
    assert np.isclose(np.mean(vals),v[3]+v[4])
    assert np.isclose(v[0],sum(v[3:6]))
    assert np.isclose(partition(a,b,np.arange(len(a)))[5],0)
    assert np.isclose(partition(a,b,np.zeros(len(a)))[4],0)
    h=pd.DataFrame({"a":[0,0,1,1],"b":[0,0,1,1],"stratum":["g"]*4})
    z=cell_metrics(h,"a","b",[(0,0),(0,1),(1,0),(1,1)]).reshape(-1,4)
    assert np.allclose(z[:,3],[.25,-.25,-.25,.25])
    print("JOINT_CELL_SELF_TEST_PASS",flush=True)

def fit(d):
    cols=["distance","area","log_n","coverage","abs_lat","climate_pc1","climate_pc2","climate_pc3","climate_pc4"]
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+NAMES+["block"])
    xs=[];ys=[]
    for _,v in d.groupby("block"):
        if len(v)<2 or v.distance.nunique()<2:continue
        x=v[cols].to_numpy(float);y=v[NAMES].to_numpy(float)
        xs.append((x-x.mean(0))/np.sqrt(len(v)));ys.append((y-y.mean(0))/np.sqrt(len(v)))
    K=len(xs);n=sum(map(len,xs));out={"complete_islands":len(d),"informative_blocks":K,"informative_islands":n}
    def calc(w):
        ids=np.flatnonzero(w)
        if not len(ids) or sum(len(xs[i])-1 for i in ids)-len(cols)<4:return None
        X=np.vstack([xs[i]*np.sqrt(w[i]) for i in ids]);Y=np.vstack([ys[i]*np.sqrt(w[i]) for i in ids])
        if np.linalg.matrix_rank(X)<len(cols):return None
        return np.linalg.lstsq(X,Y,rcond=None)[0][0]
    p=calc(np.ones(K))
    if p is None:return out|{"status":"not_estimable"}
    rng=np.random.default_rng(SEED);boot=[];leave=[]
    for _ in range(999):
        z=calc(np.bincount(rng.integers(K,size=K),minlength=K))
        if z is not None:boot.append(z)
    for i in range(K):
        w=np.ones(K);w[i]=0;z=calc(w)
        if z is not None:leave.append(z)
    ci=np.quantile(boot,[.025,.975],axis=0).T if len(boot)>=900 else [None]*len(p)
    lr=np.stack([np.min(leave,0),np.max(leave,0)],axis=1) if leave else [None]*len(p)
    return out|{"status":"descriptive_supported" if K>=30 and n>=200 and len(boot)>=900 and len(leave)==K else "insufficient_support","valid_bootstrap":len(boot),"valid_leave":len(leave),"slopes":{k:{"estimate":float(p[i]),"interval":None if ci[i] is None else list(ci[i]),"leave_range":None if lr[i] is None else list(lr[i])} for i,k in enumerate(NAMES)}}

def influence(d,j):
    cols=["distance","area","log_n","coverage","abs_lat","climate_pc1","climate_pc2","climate_pc3","climate_pc4"]
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=cols+["block"]).copy()
    d=d.groupby("block").filter(lambda v:len(v)>=2 and v.distance.nunique()>1).reset_index(drop=True)
    if len(d)<15:return {"status":"not_estimable"}
    ix={v:i for i,v in enumerate(d.island_id)}
    x=d[cols].to_numpy(float);X=np.zeros_like(x);scale=np.zeros(len(d))
    for ids in d.groupby("block").indices.values():
        scale[ids]=1/np.sqrt(len(ids));X[ids]=(x[ids]-x[ids].mean(0))*scale[ids,None]
    if np.linalg.matrix_rank(X)<len(cols):return {"status":"not_estimable"}
    w=np.linalg.pinv(X)[0]*scale
    for ids in d.groupby("block").indices.values():w[ids]-=w[ids].mean()
    targets=[("bilabiate","zygomorphic"),("tubular","zygomorphic"),("funnel_trumpet","zygomorphic")]
    q=j[j.island_id.isin(ix)].copy()
    q["a"]=q.floral_form.eq("bilabiate").astype(int)
    q["b"]=q.floral_symmetry.eq("zygomorphic").astype(int)
    result={}
    for a,b in targets:
        q["a"]=q.floral_form.eq(a).astype(int);q["b"]=q.floral_symmetry.eq(b).astype(int);q["ab"]=q.a*q.b
        z=q.groupby(["island_id","genus","grade"]).agg(n=("a","size"),a=("a","sum"),b=("b","sum"),ab=("ab","sum")).reset_index()
        z["res"]=z.ab-z.a*z.b/z.n
        z["i"]=z.island_id.map(ix)
        z["contribution"]=z.res*w[z.i.to_numpy()]/d.n.to_numpy()[z.i.to_numpy()]
        c=z.groupby("genus").contribution.sum()
        point=float(c.sum());ref=float(w@d[a+"__"+b+"__residual"].to_numpy())
        assert np.isclose(point,ref,atol=1e-10)
        active=c[c.abs()>1e-14];ordered=active.abs().sort_values(ascending=False)
        # Count independently variable strata, not merely species-rich islands.
        informative=z[(z.a>0)&(z.a<z.n)&(z.b>0)&(z.b<z.n)]
        local_support={"variable_strata":len(informative),"variable_islands":int(informative.island_id.nunique()),"variable_genera":int(informative.genus.nunique()),"variable_blocks":int(d.iloc[informative.i.unique()].block.nunique())}
        total=float(ordered.sum())
        drops={}
        # Anonymous concentration by geographic block, same original design.
        block_c=z.assign(block=d.block.to_numpy()[z.i.to_numpy()]).groupby("block").contribution.sum()
        block_abs=block_c.abs().sort_values(ascending=False)
        block_share=float(block_abs.iloc[0]/block_abs.sum()) if block_abs.sum()>0 else None
        for k in [1,3,5,10]:
            names=ordered.index[:k];drops[str(k)]=float(point-c.loc[names].sum())
        # Fixed original denominator and design: remove contribution, not taxa from roster.
        result[a+"__"+b]={"slope":point,"local_support":local_support,"top_block_absolute_share":block_share,"nonzero_contributing_genera":len(active),"absolute_contribution_sum":total,"top1_absolute_share":float(ordered.iloc[0]/total) if total else None,"top5_absolute_share":float(ordered.iloc[:5].sum()/total) if total else None,"fixed_design_drop_top_k":drops,"single_genus_drop_range":[float((point-active).min()),float((point-active).max())] if len(active) else None,"islands_with_nonzero_local_residual":int(z.loc[z.res.abs()>1e-14,"island_id"].nunique())}
    return {"status":"descriptive","islands":len(d),"blocks":int(d.block.nunique()),"targets":result}

def main():
    p=argparse.ArgumentParser();p.add_argument("--self-test",action="store_true");p.add_argument("--original-root",type=pathlib.Path);p.add_argument("--geography",type=pathlib.Path);p.add_argument("--out",type=pathlib.Path);a=p.parse_args()
    if a.self_test:self_test();return
    r=a.original_root
    f=pd.read_csv(io.BytesIO(read_pin(r/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(read_pin(r/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")),compression="gzip",dtype=str)
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert f.island_id.nunique()==4505 and not f.duplicated(["island_id","accepted_species"]).any()
    assert geo.island_id.nunique()==8264 and geo.island_id.is_unique
    t=t[t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    t["normalized_value"]=t.normalized_value.str.strip()
    t=t[t.normalized_value.ne("")&~t.normalized_value.str.contains("|",regex=False)&~t.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    t=t[t.trait_name.isin(set(sum((list(x) for x in PAIRS),[])))]
    assert not t.duplicated(["accepted_species","trait_name"]).any()
    v=t.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    q=t.pivot(index="accepted_species",columns="trait_name",values="evidence_quality")
    panels={}
    global NAMES
    cats=sorted(set(v["floral_form"].dropna()))
    syms=sorted(set(v["floral_symmetry"].dropna()))
    cells=[(a,b) for a in cats for b in syms]
    NAMES=[f"{a}__{b}__{k}" for a,b in cells for k in ("observed","margin","genus_grade","residual")]
    for x,y in PAIRS:
        z=v[[x,y]].dropna().copy()
        z["grade"]=q.loc[z.index,x]+"_"+q.loc[z.index,y]
        z["genus"]=[accepted_genus(s) for s in z.index];z=z.dropna()
        z["stratum"]=z.genus+"_"+z.grade
        for scope in ["all_recorded","original_native","high_quality_only"]:
            ff=f[f.origin_status.eq("native")] if scope=="original_native" else f
            zz=z[z.grade.eq("high_high")] if scope=="high_quality_only" else z
            totals=ff.groupby("island_id").accepted_species.nunique()
            j=ff.merge(zz,left_on="accepted_species",right_index=True,validate="many_to_one")
            rows=[]
            for i,h in j.groupby("island_id"):
                n=len(h)
                if n<5:continue
                m=cell_metrics(h,x,y,cells)
                rows.append(dict(island_id=i,n=n,coverage=n/totals[i],**dict(zip(NAMES,m))))
            d=pd.DataFrame(rows,columns=["island_id","n","coverage",*NAMES]).merge(geo,on="island_id",validate="one_to_one")
            d["distance"]=np.log2(1+d.distance_to_continent_km.where(d.distance_to_continent_km.ge(0)))
            # Use pinned geographic area column; fail closed if schema differs.
            area_col=next((c for c in ["island_area_km2","area_km2","area"] if c in d),None)
            if area_col is None:raise ValueError("Missing pinned area column")
            d["area"]=np.log2(d[area_col].where(d[area_col].gt(0)))
            d["log_n"]=np.log2(d.n);d["abs_lat"]=d.island_latitude.abs()
            d["block"]=(np.floor(d.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(d.island_longitude/10).astype("Int64").astype(str)).where(d.island_latitude.notna()&d.island_longitude.notna())
            for region in ["all","tropical","northern_extratropical","southern_extratropical"]:
                # Region labels are audited separately; no silent relabeling.
                if region=="all":u=d
                elif region=="tropical":u=d[d.analysis_regime.eq("tropical")]
                elif region=="northern_extratropical":u=d[d.analysis_regime.ne("tropical")&d.island_latitude.gt(23.5)]
                else:u=d[d.analysis_regime.ne("tropical")&d.island_latitude.lt(-23.5)]
                for nmin in [5,30]:
                    w=u[u.n.ge(nmin)]
                    key="__".join([x,y,scope,region,str(nmin)])
                    panels[key]=influence(w,j)
    out={"botanical_parent":4505,"geographic_parent":8264,"pins":PIN,"panels":panels,"warning":"Exploratory geographic association. Fixed-design genus contribution deletion, not full roster refit; no species-coherent inference. No causal or confirmatory inference. Native/high and n thresholds nested. Region definitions use tropical source label and latitude thresholds."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False,indent=2))
    for k,v in panels.items():print("JOINT_INFLUENCE_PANEL",k,json.dumps(v,allow_nan=False),flush=True)
    print("JOINT_INFLUENCE_COMPLETE",json.dumps({"panels":len(panels),"supported":sum(x["status"]=="descriptive_supported" for x in panels.values())}),flush=True)
if __name__=="__main__":main()
