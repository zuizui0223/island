"""Audit genus-heldout GloBI documented guild shares one whole family at a time.

Input names and providers stay private. Five outcomes are GLOBAL reference-share
components, NOT island-specific pollinator occurrence, pollen transfer or fitness.
Family bootstrap is one-family-one-vote; within-family morphology permutations
are posthoc and do not remove source-provider or subfamily dependence.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

GUILDS=("bombus","non_bombus_bees","lepidoptera","flower_visiting_birds","diptera")

def _project(q):
    u=np.sort(q,axis=1)[:,::-1]
    cs=np.cumsum(u,axis=1)-1.
    k=np.arange(1,q.shape[1]+1)
    rho=(u-cs/k>0).sum(axis=1)
    theta=cs[np.arange(len(q)),rho-1]/rho
    return np.maximum(q-theta[:,None],0.)

def audit_guilds(features:pd.DataFrame, families:pd.Series, documented:np.ndarray, *,
                 source_fields:tuple[str,...], floral_field:str="form_open",
                 ridge_alpha:float=100., bootstrap:int=4999,
                 shuffles:int=999, seed:int=53161)->pd.DataFrame:
    """Require fully observed and unambiguous botanical source evidence."""
    if len(features)<50 or not features.index.is_unique or not features.index.equals(families.index):
        raise ValueError("Missing or misaligned unique genus units")
    if families.isna().any() or families.astype(str).str.contains("AMBIG:|unknown|unresolved",case=False,regex=True).any():
        raise ValueError("Unresolved or conflicting source family")
    fam=families.astype(str).to_numpy()
    groups=[np.flatnonzero(fam==f) for f in sorted(np.unique(fam))]
    if len(groups)<20 or bootstrap<199 or shuffles<99 or not (np.isfinite(ridge_alpha) and ridge_alpha>0):
        raise ValueError("Insufficient spatial/taxon support or invalid simulation")
    columns=[*source_fields,floral_field]
    if len(columns)!=len(set(columns)) or not set(columns).issubset(features):
        raise ValueError("Missing or duplicate source/evidence morphology predictors")
    X=features[columns].to_numpy(float)
    Y=np.asarray(documented,float)
    if not np.isfinite(X).all():
        raise ValueError("Unknown floral form or source effort must not be imputed")
    if Y.shape!=(len(features),5) or not np.isfinite(Y).all() or (Y<0).any() or not np.allclose(Y.sum(axis=1),1):
        raise ValueError("Invalid five GloBI documented proportions")
    base=X[:,:-1];form=X[:,-1]
    folds=[(np.setdiff1d(np.arange(len(X)),te),te) for te in groups]
    def prediction(A):
        p=np.empty_like(Y)
        for tr,te in folds:
            train=A[tr];test=A[te]
            mu=train.mean(axis=0);sd=train.std(axis=0);sd[sd<1e-12]=1
            train=(train-mu)/sd;test=(test-mu)/sd
            intercept=Y[tr].mean(axis=0)
            beta=np.linalg.solve(train.T@train+ridge_alpha*np.eye(train.shape[1]),
                                 train.T@(Y[tr]-intercept))
            p[te]=intercept+test@beta
        return _project(p)
    err0=(prediction(base)-Y)**2
    err1=(prediction(X)-Y)**2
    residual=err0-err1
    movable=[ix for ix in groups if len(ix)>=2 and np.unique(form[ix]).size>1]
    rng=np.random.default_rng(seed)
    null=np.empty((shuffles,5))
    for t in range(shuffles):
        f=form.copy()
        for ix in movable:
            f[ix]=rng.permutation(form[ix])
        null[t]=(err0-(prediction(np.column_stack((base,f)))-Y)**2).mean(axis=0)
    rows=[]
    for j,name in enumerate((*GUILDS,"five_channel_mean")):
        delta=residual[:,j] if j<5 else residual.mean(axis=1)
        gen_err0=err0[:,j] if j<5 else err0.mean(axis=1)
        gen_err1=err1[:,j] if j<5 else err1.mean(axis=1)
        fammean=np.array([delta[ix].mean() for ix in groups])
        boot=fammean[rng.integers(len(fammean),size=(bootstrap,len(fammean)))].mean(axis=1)
        ctrl=null[:,j] if j<5 else null.mean(axis=1)
        rows.append(dict(channel=name,n_genera=len(X),n_families=len(groups),
            n_morphology_permutable_families=len(movable),
            n_morphology_permutable_genera=sum(map(len,movable)),
            genus_equal_baseline_MSE=float(gen_err0.mean()),
            genus_equal_model_MSE=float(gen_err1.mean()),
            genus_equal_MSE_gain=float(delta.mean()),
            family_equal_MSE_gain=float(fammean.mean()),
            family_equal_CI025=float(np.quantile(boot,.025)),
            family_equal_CI975=float(np.quantile(boot,.975)),
            n_positive_families=int((fammean>0).sum()),
            within_family_null_median=float(np.median(ctrl)),
            within_family_null_025=float(np.quantile(ctrl,.025)),
            within_family_null_975=float(np.quantile(ctrl,.975)),
            posthoc_one_sided_tail=float((1+(ctrl>=delta.mean()).sum())/(shuffles+1)),
            provider_independent=False,local_pollination_measured=False))
    return pd.DataFrame.from_records(rows)
