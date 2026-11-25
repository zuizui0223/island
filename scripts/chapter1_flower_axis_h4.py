#!/usr/bin/env python3
"""H4 post-hoc exact-species floral-axis interaction, outcome-blind support gate.

Use frozen corrected H3 effect-row identity and original GloPL accepted names.
No new experimental effect, no causal interpretation. Only anonymous summaries.
"""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t
from chapter1_flower_axis_decoupling import LOCK,digest,trait_pairs,holm
from chapter1_flower_axis_h2 import AUDIT_SHA,build_scores

GLOPL_SHA="264ca8c2237f126b3c209fcba039e19b291324ddf6e7a96992dc8235a55984de"
H3_SHA="f375b0c394cb2de8df991a6ef6107fde33a8cef909a379520f1c820ecb465d1c"
MEAS=["PL_Effect_Size_Type1","PL_Effect_Size_Type2","Constant_added","Level_of_Supplementation"]
KEY=["study_key","site_key","analysis_regime","log1p_distance_to_major_continent_km"]
GROUP=KEY+["species_key","radial","deep",*MEAS]
PRIMARY=["global","tropical"]
MIN_SPECIES=10
MIN_PUBS=5

def key(v):
    if pd.isna(v):return ""
    return " ".join(str(v).replace("_"," ").split()).casefold()

def support(z):
    cats=[]
    for a in (0,1):
        for d in (0,1):
            x=z.loc[z.radial.eq(a)&z.deep.eq(d)]
            cats.append({"radial":a,"deep":d,"species":x.species_key.nunique(),"publications":x.study_key.nunique(),"sites":x.site_key.nunique(),"effect_rows":len(x)})
    n=z.species_key.nunique();p=z.study_key.nunique()
    return {"species":n,"publications":p,"rows":len(z),"four_cells":cats,
      "pass":bool(n>=50 and p>=30 and all(x["species"]>=MIN_SPECIES and x["publications"]>=MIN_PUBS for x in cats))}

def columns(z,conditioned):
    # Fixed contrast, no trait dropping. Nuisance dummies may be redundant
    # in a restricted cohort; deterministic rank check uses no effect values.
    a=z.radial.to_numpy(float);d=z.deep.to_numpy(float)
    X=[np.ones(len(z)),a,d,a*d]
    names=["intercept","radial","deep","radial_x_deep"]
    if conditioned:
        X.append(z.repro.to_numpy(float));names.append("selfing_core")
    if np.linalg.matrix_rank(np.column_stack(X))!=len(X):return None,None,None
    dist=z.log1p_distance_to_major_continent_km.to_numpy(float)
    cand=[(dist-dist.mean(),"log_distance")]
    for c in ["analysis_regime",*MEAS]:
        lev=sorted(z[c].astype(str).unique())
        cand.extend((z[c].astype(str).eq(x).to_numpy(float),c+"="+x) for x in lev[1:])
    skipped=[]
    for v,n in cand:
        if np.linalg.matrix_rank(np.column_stack([*X,v]))>len(X):
            X.append(v);names.append(n)
        else:skipped.append(n.split("=")[0])
    return np.column_stack(X),names,sorted(set(skipped))

def fit(z,conditioned):
    g=GROUP+(["repro"] if conditioned else [])
    cells=z.groupby(g,dropna=False,as_index=False).PL_Effect_Size.mean()
    counts=cells.groupby("study_key").study_key.transform("size").to_numpy(float)
    w=1/counts
    X,names,skip=columns(cells,conditioned)
    if X is None:return {"status":"nonidentifiable_trait_design","p":1.}
    n,k=X.shape;b=cells.study_key.to_numpy(str);groups=np.unique(b);G=len(groups)
    if n<=k+2 or G<=k+1:return {"status":"insufficient_rank_or_df","p":1.}
    y=cells.PL_Effect_Size.to_numpy(float)
    inv=np.linalg.inv(X.T@(w[:,None]*X));beta=inv@(X.T@(w*y));res=y-X@beta
    score=np.stack([X[b==u].T@(w[b==u]*res[b==u]) for u in groups])
    V=G/(G-1)*(n-1)/(n-k)*inv@score.T@score@inv
    contrasts={
      "radial_advantage_shallow":np.eye(k)[1],
      "radial_advantage_deep":np.eye(k)[1]+np.eye(k)[3],
      "difference_in_radial_advantage_deep_minus_shallow":np.eye(k)[3],
    }
    out={}
    for lab,c in contrasts.items():
        v=float(c@beta);se=float(np.sqrt(max(0,c@V@c)))
        if se<=1e-12:out[lab]={"beta":v,"p":1.,"status":"degenerate_se"};continue
        q=t.ppf(.975,G-1)
        out[lab]={"beta":v,"se":se,"interval":[v-q*se,v+q*se],"p":float(2*t.sf(abs(v/se),G-1)),"status":"fit"}
    main=out["difference_in_radial_advantage_deep_minus_shallow"]
    return {"status":main["status"],"p":main["p"],"cells":n,"publications":G,"species":cells.species_key.nunique(),"nuisance_rank_omissions":skip,"contrasts":out}

def self_test():
    z=pd.DataFrame({"radial":[0,0,1,1],"deep":[0,1,0,1],"species_key":["a","b","c","d"],"study_key":["u"]*4,"site_key":["v"]*4})
    assert not support(z)["pass"]
    assert [x["species"] for x in support(z)["four_cells"]]==[1]*4
    assert key(" Abies_alba ")=="abies alba"
    print("H4_JOINT_SUPPORT_TEST_PASS")

def main():
    p=argparse.ArgumentParser()
    for k in ("ledger","state-audit","h3","glopl","out"):p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true");a=p.parse_args()
    if a.self_test:self_test();return
    for path,h in [(a.ledger,LOCK["ledger"]),(a.state_audit,AUDIT_SHA),(a.h3,H3_SHA),(a.glopl,GLOPL_SHA)]:
        if digest(path)!=h:raise ValueError("frozen source SHA mismatch")
    ledger=pd.read_csv(a.ledger,dtype=str);v=trait_pairs(ledger)
    v=v.assign(species_key=v.index.map(key),radial=v.A,deep=1-v.O)
    assert v.species_key.is_unique
    w=build_scores(pd.read_csv(a.state_audit,dtype=str))
    w=w.assign(species_key=w.index.map(key),repro=w.min2)
    assert w.species_key.is_unique
    v=v.merge(w[["species_key","repro"]],on="species_key",how="left",validate="one_to_one")
    # Source identity and all support computed before reading effect column.
    raw=pd.read_csv(a.glopl,encoding="latin-1",dtype=str,usecols=["Species_accepted_names"])
    r=pd.read_csv(a.h3,usecols=["row_id",*KEY,*MEAS],dtype={x:str for x in [*KEY[:3],*MEAS]})
    assert len(r)==2969 and r.row_id.is_unique and r.row_id.min()>=0 and r.row_id.max()<len(raw)
    r["species_key"]=raw.iloc[r.row_id.to_numpy(int)].Species_accepted_names.map(key).to_numpy()
    for c in MEAS:r[c]=r[c].fillna("").str.strip()
    z=r.merge(v[["species_key","radial","deep","repro"]],on="species_key",validate="many_to_one")
    cohorts={};rows=[]
    for region in PRIMARY:
        q=z if region=="global" else z.loc[z.analysis_regime.eq("tropical")]
        for conditioned in (False,True):
            x=q.dropna(subset=["repro"]) if conditioned else q
            s=support(x);cohorts[region,conditioned]=x
            rows.append({"region":region,"conditioned":conditioned,"support":s,"status":"insufficient_support","p":1.})
    print("H4_JOINT_PREFLIGHT "+json.dumps(rows,allow_nan=False),flush=True)
    if any(x["support"]["pass"] for x in rows):
        h=pd.read_csv(a.h3,usecols=["row_id","PL_Effect_Size"])
        orig=pd.read_csv(a.glopl,encoding="latin-1",usecols=["PL_Effect_Size"])
        check=pd.to_numeric(orig.iloc[h.row_id.to_numpy(int)].PL_Effect_Size,errors="coerce")
        assert np.allclose(h.PL_Effect_Size.to_numpy(float),check.to_numpy(float),rtol=0,atol=1e-12)
        for row in rows:
            if row["support"]["pass"]:
                x=cohorts[row["region"],row["conditioned"]].merge(h,on="row_id",validate="one_to_one")
                row.update(fit(x,row["conditioned"]))
    holm(rows)
    result={"owner":"island","sources":{"own_ledger":LOCK["ledger"],"own_state_audit":AUDIT_SHA,"corrected_H3_rows":H3_SHA,"original_GloPL":GLOPL_SHA},"panels":rows,
        "primary":"tropical unconditioned interaction; global and min2-conditioned same-source diagnostics; Holm over 4 fixed slots",
        "warning":"Post-hoc current PL association on exact normalized accepted species only. Minimum 10 distinct species and 5 publications per morphology cell, overall 50 species and 30 publications. Publication equal weighting/finite cluster t; repeated taxa may remain dependent. H3/H4 original experiments reused, no independent replication, causal mediation, or measured visitor mechanics. No outcome fit if support fails."}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for row in rows:print("H4_JOINT_PANEL "+json.dumps(row,allow_nan=False),flush=True)
    print("H4_JOINT_COMPLETE 4",flush=True)

if __name__=="__main__":main()
