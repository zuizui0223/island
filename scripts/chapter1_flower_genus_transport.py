#!/usr/bin/env python3
"""Held-out genus profile transport of PR354's tropical conditional contrast."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from chapter1_flower_axis_decoupling import LOCK, REGIONS, digest, trait_pairs, make_frame, design, estimate, rates

MODES = ("external", "other_tropical")
SCOPES = (("all_recorded", 10), ("all_recorded", 30), ("high_only", 10), ("literal_native", 10))
MODELS = ("blind", "entry", "loading")
ALPHA = 5.
DRAWS = 499
SEED = 20261128

def train(C, genus, donor, embargo):
    allowed = np.setdiff1d(np.unique(donor), np.unique(embargo))
    assert not np.intersect1d(allowed, embargo).size
    if len(allowed) < 100:
        return None
    z = C[allowed].sum(0)
    if z[:2].sum() == 0 or z[2:].sum() == 0:
        return None
    prior = (z + .5) / (len(allowed) + 2.)
    ng = int(genus.max()) + 1
    count = np.bincount(genus[allowed], minlength=ng)
    sums = np.zeros((ng, 4))
    np.add.at(sums, genus[allowed], C[allowed])
    probs = (sums + ALPHA * prior) / (count[:, None] + ALPHA)
    assert np.allclose(probs.sum(1), 1)
    return prior, probs, count, allowed

def predict(prior, probs, count, target_genus):
    g, n = np.unique(target_genus, return_counts=True)
    entry = probs[g].mean(0)
    loading = np.average(probs[g], weights=n, axis=0)
    p = np.stack([prior, entry, loading])
    assert np.all(p > 0) and np.allclose(p.sum(1), 1)
    return p, float(np.mean(count[target_genus] == 0)), float(np.mean(count[g] == 0))

def contrast(p):
    return p[..., 0] / p[..., :2].sum(-1) - p[..., 2] / p[..., 2:].sum(-1)

def block_gain(values, blocks):
    by = np.stack([values[blocks == b].mean(0) for b in np.unique(blocks)])
    rng = np.random.default_rng(SEED)
    boot = by[rng.integers(len(by), size=(DRAWS, len(by)))].mean(1)
    return {"mean": by.mean(0).tolist(), "range_025_975": np.quantile(boot, [.025,.975], axis=0).T.tolist(),
            "blocks": len(by), "not_calibrated_ci": True}

def gain(x):
    return np.stack([x[...,0]-x[...,2], x[...,0]-x[...,1], x[...,1]-x[...,2]], axis=-1)

def donor_rosters(f, v, geo, scope):
    # Embargo across ALL recorded status rows, even when the target is native.
    index = pd.Series(np.arange(len(v)), index=v.index, name="si")
    gg = geo[["island_id", "analysis_regime", "spatial_block"]].dropna()
    gg = gg.loc[gg.analysis_regime.isin(REGIONS) & gg.spatial_block.astype(str).str.strip().ne("")]
    z = f.merge(index, left_on="accepted_species", right_index=True, validate="many_to_one")
    z = z.merge(gg, on="island_id", validate="many_to_one")
    embargo = {str(b): np.unique(x.si.to_numpy(int)) for b,x in z.groupby("spatial_block")}
    if scope == "literal_native":
        z = z.loc[z.origin_status.eq("native")]
    rosters = {mode: {str(b): np.unique(x.si.to_numpy(int)) for b,x in
                     z.loc[z.analysis_regime.ne("tropical") if mode=="external" else z.analysis_regime.eq("tropical")].groupby("spatial_block")}
               for mode in MODES}
    return embargo, rosters

def panel(f, ledger, geo, parent, scope, nmin):
    v = trait_pairs(ledger, scope=="high_only")
    d, M, C, genus = make_frame(f, v, geo, scope)
    d = d.loc[d.analysis_regime.eq("tropical") & d.carrier_n.ge(nmin) &
              d.deep_n.ge(5) & d.shallow_n.ge(5)].sort_values("mi")
    old = next(x for x in parent["panels"] if x["scope"]==scope and x["nmin"]==nmin and
               x["region"]=="tropical" and x["endpoint"]=="conditional_difference")
    blocks = d.spatial_block.to_numpy(str)
    n = len(d)
    if n != old["n"] or len(set(blocks)) != old["blocks"]:
        raise ValueError("parent mask mismatch")
    base = {"scope":scope, "nmin":nmin, "islands":n,"blocks":len(set(blocks))}
    if n < 50 or len(set(blocks)) < 10:
        return [base | {"mode":m,"status":"insufficient_parent_support"} for m in MODES]
    X,j = design(d)
    mm = M[d.mi.to_numpy(int)].tocsr()
    q = np.asarray(mm @ C) / d.carrier_n.to_numpy()[:,None]
    obs = contrast(q)
    original = estimate(X, obs[:,None], blocks,j)
    if original is None or not np.isclose(original["0"]["beta"], old["beta"], atol=1e-10,rtol=0):
        raise ValueError("parent slope mismatch")
    ix = [mm.indices[mm.indptr[i]:mm.indptr[i+1]] for i in range(n)]
    embargo, rosters = donor_rosters(f,v,geo,scope)
    panels=[]
    for mode in MODES:
        ps=np.full((n,3,4), np.nan); unseen=np.full((n,2),np.nan)
        support=[]; failed=0
        for b in np.unique(blocks):
            test=np.flatnonzero(blocks==b)
            ban=embargo[b]
            target=np.unique(np.concatenate([ix[i] for i in test]))
            assert np.all(np.isin(target,ban))
            donor_groups=[np.setdiff1d(s,ban) for k,s in rosters[mode].items() if k!=b]
            donor_groups=[s for s in donor_groups if len(s)]
            if len(donor_groups)<10:
                failed+=1; continue
            train_result=train(C,genus,np.unique(np.concatenate(donor_groups)),ban)
            if train_result is None:
                failed+=1; continue
            prior, profiles, counts, allowed=train_result
            assert not np.intersect1d(allowed,target).size
            support.append([len(allowed),len(donor_groups)])
            for i in test:
                ps[i],unseen[i,0],unseen[i,1]=predict(prior,profiles,counts,genus[ix[i]])
        mask=np.isfinite(ps).all(axis=(1,2))
        out=base|{"mode":mode,"scored_islands":int(mask.sum()),"failed_folds":failed,
                  "status":"fit" if mask.all() else "incomplete_fold_support"}
        if support:
            out["donor_species_min_median_max"]=np.quantile(np.asarray(support)[:,0],[0,.5,1]).tolist()
            out["donor_blocks_min"]=int(np.asarray(support)[:,1].min())
        if not mask.all():
            panels.append(out);continue
        out["unknown_genus_weight_mean"]=dict(zip(("loading","entry"),unseen.mean(0).tolist()))
        pred=contrast(ps)
        y=np.column_stack([obs,pred, pred[:,1]-pred[:,0],pred[:,2]-pred[:,1],obs-pred[:,2]])
        labels=("observed","blind","entry","loading","entry_increment","loading_increment","residual")
        assert np.allclose(y[:,0],y[:,1]+y[:,4]+y[:,5]+y[:,6])
        fit=estimate(X,y,blocks,j)
        if fit is None:
            raise ValueError("same design fit failed")
        slopes={k:float(fit[str(i)]["beta"]) for i,k in enumerate(labels)}
        assert np.isclose(slopes["observed"],sum(slopes[k] for k in ("blind","entry_increment","loading_increment","residual")),atol=1e-10)
        out["slopes"]=slopes
        sq=np.sum((ps-q[:,None,:])**2,axis=2)
        log=-np.sum(q[:,None,:]*np.log(ps),axis=2)
        ce=(pred-obs[:,None])**2
        out["block_equal_gains"]={k:block_gain(gain(a),blocks) for k,a in
                                    (("category_squared_error",sq),("log_loss",log),("conditional_squared_error",ce))}
        # Proper Brier and log-score gains per actual species, averaged across its
        # held-out occurrences, then species- or genus-equal. No new significance.
        tall=np.zeros((len(v),6)); occurrences=np.zeros(len(v))
        for i, ids in enumerate(ix):
            loss=np.sum((ps[i][None,:,:]-C[ids][:,None,:])**2,axis=2)
            ll=-(C[ids]@np.log(ps[i]).T)
            gg=np.column_stack([gain(loss),gain(ll)])
            tall[ids]+=gg; occurrences[ids]+=1
        test=occurrences>0
        per=tall[test]/occurrences[test,None]; tg=genus[test]
        bygen=np.stack([per[tg==g].mean(0) for g in np.unique(tg)])
        out["target_taxon_equal_gain"]={"species":int(test.sum()),"genera":int(len(bygen)),
                                       "species_equal":per.mean(0).tolist(),
                                       "genus_equal":bygen.mean(0).tolist()}
        # Proper score identity: squared category distance differs from expected
        # Brier only by observed irreducible impurity, shared across models.
        assert np.allclose(gain(np.sum(ps**2,axis=2)-2*np.sum(ps*q[:,None,:],axis=2)+1),gain(sq))
        out["conditional_rmse_by_model"]=dict(zip(MODELS,np.sqrt(ce.mean(0)).tolist()))
        out["category_mse_by_model"]=dict(zip(MODELS,sq.mean(0).tolist()))
        panels.append(out)
    return panels

def self_test():
    C=np.eye(4)[np.tile([0,0,1,2,3],30)]
    g=np.repeat([0,1,2],50)
    t=train(C,g,np.arange(150),np.arange(10))
    assert t is not None and len(t[3])==140 and min(t[3])>=10
    p,u,e=predict(*t[:3],np.array([0,0,1]))
    assert np.allclose(p.sum(1),1)
    assert np.allclose(gain(np.array([3.,2.,1.])),[2,1,1])
    # Nonuniform genus representation changes mixture; unknown genus uses prior.
    q=np.array([[.7,.1,.1,.1],[.1,.1,.1,.7],[.25,.25,.25,.25]])
    p,u,e=predict(q[2],q,np.array([5,5,0]),np.array([0,0,1,2]))
    assert np.isclose(u,.25) and np.isclose(e,1/3)
    assert not np.allclose(p[1],p[2])
    assert np.isclose(contrast(np.array([.375,.125,.125,.375])),.5)
    print("FLORAL_GENUS_TRANSPORT_SELF_TEST_PASS")

def main():
    p=argparse.ArgumentParser()
    for k in ("flora","ledger","geography","parent","out"):
        p.add_argument("--"+k,type=Path)
    p.add_argument("--self-test",action="store_true")
    a=p.parse_args()
    if a.self_test:
        self_test();return
    for k,h in LOCK.items():
        if digest(getattr(a,k))!=h:raise ValueError("source hash mismatch "+k)
    parent=json.loads(a.parent.read_text())
    assert parent["run"]==45579499069
    f=pd.read_csv(a.flora,dtype=str,usecols=["island_id","accepted_species","origin_status"])
    ledger=pd.read_csv(a.ledger,dtype=str)
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    assert len(f)==1039757 and f.island_id.nunique()==4505 and len(geo)==8264
    assert not f.duplicated(["island_id","accepted_species"]).any() and geo.island_id.is_unique
    rows=[]
    for scope,nmin in SCOPES:
        for r in panel(f,ledger,geo,parent,scope,nmin):
            rows.append(r)
            print("FLORAL_GENUS_TRANSPORT_PANEL "+json.dumps(r,allow_nan=False),flush=True)
    assert len(rows)==8
    result={"owner":"island","parent_run":45579499069,"source_hashes":LOCK,
            "primary":"all_recorded n10 external; loading versus blind category squared-error gain",
            "gain_order":["loading_vs_blind","entry_vs_blind","loading_vs_entry"],
            "target_taxon_gain_order":["Brier loading/blind","Brier entry/blind","Brier loading/entry","log loading/blind","log entry/blind","log loading/entry"],
            "panels":rows,
            "warning":"Post-outcome known-genus conditional prediction with full accepted paired-species block embargo. Observed target genus and inherited eligibility supplied; not prospective establishment. Bootstrap with frozen folds descriptive. No first-source grouping or independent biological replication. No raw exports."}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("FLORAL_GENUS_TRANSPORT_COMPLETE 8",flush=True)

if __name__=="__main__":main()
