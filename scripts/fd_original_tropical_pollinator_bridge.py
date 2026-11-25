#!/usr/bin/env python3
"""Protected original GIFT4505 radial/bilateral-deep flower endpoints vs independent animal records.
The public geography bridge is external. All plant rows, island/plant/genus identifiers and
matched membership stay inside protected source runner. Only aggregate coefficients and counts exit.
"""
import argparse,collections,hashlib,io,json,math,pathlib,urllib.request,zipfile
import numpy as np
import pandas as pd
from scipy import sparse

PIN={"flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
"traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
"assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
"geography":"08218617b5880153986012e79817de65a74dd4e44b5a65e8a537dab8cb2cd0a7",
"bridge":"467f111ca14be773bc7dbe6cd76880c443d5f8d6a81c0dad55ca6b8799302776"}
BEE_MD5="372c22e2c24e203600a29d7f874a907d"
MODEL_SHA="a7ddae4cdb41453313d20811371cf64bd38b42b788fc89ae96c284d8d667ddd7"
OUTCOMES=("radial_deep","bilateral_deep","joint_Q")
TRAITS=("floral_symmetry","tube_depth_class")
def sha(b):return hashlib.sha256(b).hexdigest()
def read_pin(p,key):
    b=p.read_bytes()
    if sha(b)!=PIN[key]:raise ValueError("Source hash mismatch "+key)
    return b
def get(url,cap):
    req=urllib.request.Request(url,headers={"User-Agent":"fd-protected-analysis/1.0"})
    with urllib.request.urlopen(req,timeout=120) as r:
        b=r.read(cap+1)
    if len(b)>cap:raise ValueError("download cap")
    return b
def bee_model():
    rec=json.loads(get("https://zenodo.org/api/records/23222636",2_000_000))
    f=next(f for f in rec["files"] if f["key"]=="Data.zip")
    assert f["size"]==262554015 and f["checksum"]=="md5:"+BEE_MD5
    b=get(f["links"].get("content") or f["links"]["self"],300_000_000)
    assert hashlib.md5(b).hexdigest()==BEE_MD5
    with zipfile.ZipFile(io.BytesIO(b)) as z:m=z.read("Data/model_table.csv")
    assert sha(m)==MODEL_SHA
    try:s=m.decode("utf-8-sig")
    except UnicodeDecodeError:s=m.decode("cp1252")
    d=pd.read_csv(io.StringIO(s),dtype={"OBJECTID_1":str}).rename(columns={"OBJECTID_1":"bee_id","Archipelago":"group","FinalRichNat":"bee","AreaModel":"area","DistModel":"distance","Cov":"bee_cov"})
    assert len(d)==305 and d.bee_id.is_unique
    return d[["bee_id","group","bee","area","distance","bee_cov"]]
def public_animals(p):
    a=json.loads(read_pin(p,"bridge"))
    x=pd.DataFrame(a["geometry"],columns=a["public_bridge_columns"]).rename(columns={"gshhg_geometry_hash_id":"island_id"})
    assert len(x)==234 and x.bee_id.is_unique and x.island_id.is_unique
    # Rounded public overlap values; originally strict >=.9 in the source geometry runner.
    assert (x[["bee_overlap","gshhg_overlap"]]>=.9-1e-8).all().all()
    d=x.merge(bee_model(),on="bee_id",validate="one_to_one")
    b=pd.DataFrame(a["bird_candidate"]).drop(columns=["group"])
    t=pd.DataFrame(a["butterfly_record_opportunity"])
    assert len(b)==52 and len(t)==98
    d=d.merge(b,on="bee_id",how="left",validate="one_to_one").merge(t,on="bee_id",how="left",validate="one_to_one")
    for c in ("bee","area","distance","bee_cov","nectar_gt0","nectar_ge50","all_birds"):
        d[c]=pd.to_numeric(d[c],errors="coerce")
    assert d.bee.notna().all() and d.bee.ge(0).all()
    d["log_bee"]=np.log1p(d.bee)
    d["log_area"]=np.log2(d.area);d["log_distance"]=np.log2(1+d.distance)
    d["log_nectar"]=np.log1p(d.nectar_gt0);d["log_nectar50"]=np.log1p(d.nectar_ge50)
    d["log_nonnectar"]=np.log1p(d.all_birds-d.nectar_gt0)
    for mask in ("all","precision_le1km","no_gbifid_with_reference"):
        d[mask+"_records"]=pd.to_numeric(d[mask+"_records"],errors="coerce")
        d[mask+"_recorded_species"]=pd.to_numeric(d[mask+"_recorded_species"],errors="coerce")
    return d
def accepted_genus(s):
    p=s.strip().split()
    return p[0] if len(p)>=2 and p[0][0].isupper() and p[0].replace("-","").isalpha() and p[1].lower() not in {"sp.","spp.","cf.","aff.","?"} else None
def pair_states(t):
    z=t[t.trait_name.isin(TRAITS)&t.evidence_scope.eq("species_direct")&t.evidence_quality.isin(["high","medium"])&t.normalized_value.notna()].copy()
    z["normalized_value"]=z.normalized_value.str.strip()
    z=z[z.normalized_value.ne("")&~z.normalized_value.str.contains("|",regex=False)&~z.normalized_value.str.lower().isin(["none","nan","null","unknown","unresolved"])]
    assert not z.duplicated(["accepted_species","trait_name"]).any() and z.accepted_species.notna().all()
    a=z.pivot(index="accepted_species",columns="trait_name",values="normalized_value").dropna(subset=list(TRAITS))
    q=z.pivot(index="accepted_species",columns="trait_name",values="evidence_quality").reindex(a.index)
    a["grade"]=np.where(q[list(TRAITS)].eq("high").all(axis=1),"high","medium")
    a["genus"]=[accepted_genus(s) for s in a.index]
    a=a[a.genus.notna()]
    a["pair"]=list(zip(a[TRAITS[0]],a[TRAITS[1]]))
    assert a.index.is_unique
    return a
def fit_block(d,xcols):
    idx=[];xs=[];mult=[];groups=[]
    for g,v in d.groupby("group"):
        if len(v)<2:continue
        x=v[xcols].to_numpy(float)
        if np.max(np.ptp(x,axis=0))<1e-12:continue
        x=x-x.mean(0)
        idx.append(v.index.to_numpy());xs.append(x/math.sqrt(len(v)));mult.append(1/math.sqrt(len(v)));groups.append(g)
    if not xs:return None
    X=np.vstack(xs);p=X.shape[1]
    dof=len(X)-len(xs)-p
    if dof<4 or np.linalg.matrix_rank(X)<p:return None
    inverse=np.linalg.pinv(X)
    B=np.zeros((p,len(d)));pos=0
    for ids,w in zip(idx,mult):
        b=inverse[:,pos:pos+len(ids)]*w;B[:,ids]=b-b.mean(1)[:,None];pos+=len(ids)
    return B,groups,dof,len(X)
def null_p(draw,obs):
    return float(min(1.,2*min((1+np.sum(draw<=obs))/(len(draw)+1),(1+np.sum(draw>=obs))/(len(draw)+1))))
def analyze(d,pair,rosters,features,controls,nperm,seed):
    xcols=features+controls
    d=d.replace([np.inf,-np.inf],np.nan).dropna(subset=xcols+["group"]).reset_index(drop=True)
    result={"island_count":len(d),"archipelago_count":int(d.group.nunique()),"features":features,"controls":controls}
    if len(d)<8:return result|{"status":"insufficient"}
    fit=fit_block(d,xcols)
    if fit is None:return result|{"status":"insufficient_rank"}
    B,groups,dof,ninfo=fit
    result|={"status":"exploratory","informative_groups":len(groups),"n_informative":ninfo,"dof":dof,"primary_support":len(groups)>=10 and ninfo>=30}
    rs=[rosters[i] for i in d.island_id]
    spp=sorted(set.union(*rs));index={s:i for i,s in enumerate(spp)}
    ii=[];jj=[]
    for i,s in enumerate(rs):
        ii.extend([i]*len(s));jj.extend(index[z] for z in s)
    n=np.array(list(map(len,rs)))
    A=sparse.csr_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(d),len(spp)))
    A=sparse.diags(1/n)@A
    tab=pair.loc[spp]
    cats=sorted(set(tab.pair));catidx={c:i for i,c in enumerate(cats)}
    code=np.array([catidx[c] for c in tab.pair])
    target=np.array([[c==("actinomorphic","deep"),c==("zygomorphic","deep")] for c in cats],float)
    eye=np.eye(len(cats))
    def y_for(c):
        p=A@eye[c]
        return np.column_stack([p@target,(1-np.square(p).sum(1))*n/(n-1)])
    Y=y_for(code);obs=(B@Y)[:len(features)]
    result["distinct_scored_species"]=len(spp);result["median_scored_species"]=float(np.median(n))
    result["outcome_means"]={o:float(Y[:,i].mean()) for i,o in enumerate(OUTCOMES)}
    rng=np.random.default_rng(seed)
    tests={}
    for label,fields in [("genus",["genus"]),("genus_quality",["genus","grade"])]:
        groups0=[np.array(v,dtype=int) for v in tab.reset_index().groupby(fields,sort=True).indices.values() if len(v)>1]
        groups0=[g for g in groups0 if len(np.unique(code[g]))>1]
        draws=np.empty((nperm,len(features),len(OUTCOMES)))
        for k in range(nperm):
            c=code.copy()
            for ix in groups0:c[ix]=rng.permutation(c[ix])
            draws[k]=(B@y_for(c))[:len(features)]
        tests[label]={"draws":draws,"movable_species":sum(map(len,groups0))}
    out={}
    for j,x in enumerate(features):
        out[x]={}
        for i,y in enumerate(OUTCOMES):
            row={"coefficient":float(obs[j,i])}
            for label,t in tests.items():
                v=t["draws"][:,j,i]
                row[label]={"movable_species":t["movable_species"],"null_mean":float(v.mean()),"null_q025_q975":np.quantile(v,[.025,.975]).tolist(),"observed_minus_null_mean":float(obs[j,i]-v.mean()),"p_two_sided":null_p(v,obs[j,i])}
            out[x][y]=row
    # Leave-one-archipelago diagnostic on same original observed trait rows.
    for g in groups:
        keep=d.group.ne(g).to_numpy()
        z=fit_block(d.loc[keep].reset_index(drop=True),xcols)
        if z is None:continue
        beta=(z[0]@Y[keep])[:len(features)]
        for j,x in enumerate(features):
            for i,y in enumerate(OUTCOMES):
                out[x][y].setdefault("leave_one_group_values",[]).append(float(beta[j,i]))
    for x in out:
        for row in out[x].values():
            v=row.pop("leave_one_group_values",[])
            row["leave_one_group_min_max"]=[min(v),max(v)] if v else None
            row["leave_one_group_n"]=len(v)
    result["outcomes"]=out
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--original-root",required=True,type=pathlib.Path);p.add_argument("--geography",required=True,type=pathlib.Path);p.add_argument("--bridge",required=True,type=pathlib.Path);p.add_argument("--out",required=True,type=pathlib.Path);p.add_argument("--permutations",type=int,default=199);a=p.parse_args()
    root=a.original_root
    fb=read_pin(root/"fixed/canonical/input/chapter1_status_flora.csv.gz","flora")
    tb=read_pin(root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz","traits")
    read_pin(root/"fixed/source/result/island_source_assignments.csv.gz","assign")
    flora=pd.read_csv(io.BytesIO(fb),compression="gzip",dtype=str,usecols=["island_id","accepted_species","origin_status"])
    t=pd.read_csv(io.BytesIO(tb),compression="gzip",dtype=str,usecols=["accepted_species","trait_name","normalized_value","evidence_scope","evidence_quality"])
    geo=pd.read_csv(io.BytesIO(read_pin(a.geography,"geography")),dtype={"island_id":str})
    assert flora.island_id.nunique()==4505 and not flora.duplicated(["island_id","accepted_species"]).any()
    assert geo.island_id.nunique()==8264 and geo.island_id.is_unique and "analysis_regime" in geo
    animals=public_animals(a.bridge).merge(geo[["island_id","analysis_regime"]],on="island_id",validate="one_to_one")
    pair=pair_states(t);panels={};support={}
    tropical=set(geo.loc[geo.analysis_regime.eq("tropical"),"island_id"])
    source_tropical=flora[flora.island_id.isin(tropical)]
    assert source_tropical.island_id.nunique()==1558
    animals=animals[animals.analysis_regime.eq("tropical")].copy()
    for scope in ("all_recorded","original_native","high_quality_only"):
        f=source_tropical if scope!="original_native" else source_tropical[source_tropical.origin_status.eq("native")]
        total=f.groupby("island_id").accepted_species.nunique()
        pa=pair if scope!="high_quality_only" else pair[pair.grade.eq("high")]
        selected=f[f.accepted_species.isin(pa.index)&f.island_id.isin(animals.island_id)]
        roster=selected.groupby("island_id").accepted_species.agg(set).to_dict()
        d=animals[animals.island_id.isin(roster)].copy()
        d["n"]=d.island_id.map(lambda i:len(roster[i]))
        d["log_n"]=np.log2(d.n);d["direct_coverage"]=d.n/d.island_id.map(total)
        support[scope]={"original_tropical_flora_islands":int(f.island_id.nunique()),"matched_pair_islands":len(d),"n5":int(d.n.ge(5).sum()),"n30":int(d.n.ge(30).sum())}
        for threshold in ((5,30) if scope!="high_quality_only" else (5,)):
            q=d[d.n.ge(threshold)].copy()
            cv=["log_area","log_distance","log_n"]
            cases=[("bee",q,["log_bee"],cv),("bee_coverage",q,["log_bee"],cv+["bee_cov","direct_coverage"]),
                   ("bird",q,["log_bee","log_nectar"],cv),
                   ("bird_control",q,["log_bee","log_nectar","log_nonnectar"],cv)]
            for mask in ("all","precision_le1km","no_gbifid_with_reference"):
                b=q[q[mask+"_records"].ge(10)&q[mask+"_recorded_species"].gt(0)].copy()
                b["log_butterfly"]=np.log1p(b[mask+"_recorded_species"])
                b["log_butterfly_records"]=np.log1p(b[mask+"_records"])
                cases.append(("butterfly_"+mask,b,["log_bee","log_butterfly"],cv+["log_butterfly_records"]))
            for label,b,fs,cs in cases:
                key=f"{scope}_n{threshold}_{label}"
                panels[key]=analyze(b,pa,roster,fs,cs,a.permutations,20261124)
    # Adjust within the complete reported genus-quality family; deliberately exploratory.
    ps=[]
    for v in panels.values():
        for x in v.get("outcomes",{}).values():
            for row in x.values():
                ps.append(row["genus_quality"])
    order=sorted(ps,key=lambda r:r["p_two_sided"]);last=1.
    for i in range(len(order)-1,-1,-1):
        last=min(last,order[i]["p_two_sided"]*len(order)/(i+1));order[i]["BH_q_exploratory"]=float(last)
    out={"status":"ORIGINAL_GIFT4505_EXTERNAL_ANIMAL_BRIDGE_EXPLORATORY","original_parent_botanical_islands":4505,"original_geography_islands":8264,"original_tropical_botanical_islands":1558,"public_geometry_bridge":234,"source_hashes":PIN,"support":support,"panels":panels,"permutations":a.permutations,"genus_quality_test_count":len(ps),"warning":"Original genuine species Direct joint pair states, not local phenotype variation or effective pollination. Fixed incidence genotype-free taxonomy shuffles conditional on matched carrier pool; no historical adaptation or causality. Provisional bird geography and positive-only butterfly records. No individual botanical or geographic identifiers exported."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,indent=2,allow_nan=False))
    print("ORIGINAL_PAIR_SUPPORT",json.dumps(support),flush=True)
    for k,v in panels.items():print("ORIGINAL_PAIR_PANEL",k,json.dumps(v),flush=True)
if __name__=="__main__":main()
