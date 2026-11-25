"""Source-region reuse audit for fd functional traits.

Original island x source IDs and individual predictions stay private.
Returned grouped summaries contain no identities or inferential p-values.
Connected components of overlapping GIFT source labels are conservative
stress tests, NOT demonstrated independent ancestral histories.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

SCORE_FIELDS = ("region","scope","trait","ledger","island_id","spatial_block","seed","brier_advantage")
SOURCE_FIELDS = ("island_id","source_mode","source_rank","entity_ID")
COLUMNS = (
    "region","scope","trait","ledger","n_islands","n_spatial_blocks",
    "n_source_signatures","top_two_source_signature_islands","top_two_fraction",
    "block_equal_gain","signature_equal_gain","signature_bootstrap_lower",
    "signature_bootstrap_upper","fraction_signatures_positive",
    "remaining_islands_drop_top_two","block_equal_gain_drop_top_two",
    "components_any_donor_overlap","largest_component_any_overlap",
    "components_three_of_five_overlap","largest_component_three_of_five_overlap",
    "largest_donor_component_islands","component_equal_gain",
    "fraction_donor_components_positive",
)

def _component_labels(signatures: list[tuple[str, ...]], threshold: int):
    n=len(signatures)
    if n==0:
        return 0,np.array([],dtype=int)
    s=[set(x) for x in signatures]
    edges=[(i,i) for i in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if len(s[i]&s[j])>=threshold:
                edges.extend(((i,j),(j,i)))
    arr=np.asarray(edges,dtype=int)
    adj=csr_matrix((np.ones(len(arr),dtype=np.int8),(arr[:,0],arr[:,1])),shape=(n,n))
    nc,labels=connected_components(adj,directed=False)
    return int(nc),labels

def _components(signatures: list[tuple[str, ...]], threshold: int):
    nc,labels=_component_labels(signatures,threshold)
    return nc,(int(np.bincount(labels).max()) if nc else 0)

def audit_source_overlap(scores: pd.DataFrame, assignments: pd.DataFrame, *,
    control_seed: int=20261110, bootstrap_draws: int=4999,
    min_islands: int=40, min_blocks: int=10, min_signatures: int=5
) -> pd.DataFrame:
    """One preselected placebo seed only; suppress thin cells.

    'signature bootstrap' resamples the *identical* five-region source
    signatures, which share donor regions and thus are NOT independent units.
    Its interval is a descriptive sensitivity, not confirmatory uncertainty.
    """
    if min_islands<40 or min_blocks<10 or min_signatures<5 or bootstrap_draws<199:
        raise ValueError("Privacy/support gate too small")
    if not set(SCORE_FIELDS).issubset(scores):
        raise ValueError("Score fields missing")
    if not set(SOURCE_FIELDS).issubset(assignments):
        raise ValueError("Source fields missing")
    x=scores[list(SCORE_FIELDS)].copy()
    if x.empty:
        raise ValueError("Empty scores")
    for field in ("region","scope","trait","ledger","island_id","spatial_block"):
        if x[field].isna().any() or x[field].astype(str).str.strip().eq("").any():
            raise ValueError("Missing score metadata")
        x[field]=x[field].astype(str).str.strip()
    for field in ("seed","brier_advantage"):
        x[field]=pd.to_numeric(x[field],errors="coerce")
        if not np.isfinite(x[field].to_numpy(dtype=float)).all():
            raise ValueError("Invalid score or seed")
    x=x[x.seed.eq(control_seed)].copy()
    if x.empty:
        raise ValueError("Predeclared placebo seed absent")
    if x.duplicated(["region","scope","trait","ledger","island_id"]).any():
        raise ValueError("Duplicate island-level scores")
    if (x.groupby("island_id").region.nunique().gt(1).any()
        or x.groupby("island_id").spatial_block.nunique().gt(1).any()):
        raise ValueError("Inconsistent focal geography")
    y=assignments[list(SOURCE_FIELDS)].copy()
    for field in ("island_id","source_mode","entity_ID"):
        if y[field].isna().any() or y[field].astype(str).str.strip().eq("").any():
            raise ValueError("Missing source metadata")
        y[field]=y[field].astype(str).str.strip()
    y["source_rank"]=pd.to_numeric(y.source_rank,errors="coerce")
    if y.source_rank.isna().any():
        raise ValueError("Missing source rank")
    y=y[y.source_mode.eq("geo_k5")&y.source_rank.between(1,5)].copy()
    if y.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate rank")
    if y.duplicated(["island_id","entity_ID"]).any():
        raise ValueError("Duplicate source donor")
    sig=y.groupby("island_id").entity_ID.agg(lambda v:tuple(sorted(set(v))))
    if not sig.map(len).eq(5).all():
        raise ValueError("Every source signature must have exactly five donors")
    x=x.join(sig.rename("signature"),on="island_id")
    if x.signature.isna().any():
        raise ValueError("Missing source signature")
    rng=np.random.default_rng(control_seed)
    out=[]
    for (region,scope,trait,ledger),g in x.groupby(
        ["region","scope","trait","ledger"],sort=True
    ):
        ni=g.island_id.nunique(); nb=g.spatial_block.nunique()
        avg=g.groupby("signature").brier_advantage.mean()
        ns=len(avg)
        if ni<min_islands or nb<min_blocks or ns<min_signatures:
            continue
        sizes=g.groupby("signature").size().sort_values(ascending=False)
        top=list(sizes.index[:2])
        rem=g[g.signature.map(lambda t:t not in top)]
        val=avg.to_numpy()
        reps=val[rng.integers(0,len(val),size=(bootstrap_draws,len(val)))].mean(axis=1)
        signatures=list(avg.index)
        ncomp,labels=_component_labels(signatures,threshold=1)
        ncomp3,largest3=_components(signatures,threshold=3)
        groupids=dict(zip(signatures,labels))
        withcomp=g.assign(component=g.signature.map(groupids))
        groupavg=withcomp.groupby("component").brier_advantage.mean()
        out.append(dict(
            region=region,scope=scope,trait=trait,ledger=ledger,
            n_islands=int(ni),n_spatial_blocks=int(nb),n_source_signatures=int(ns),
            top_two_source_signature_islands=int(sizes.iloc[:2].sum()),
            top_two_fraction=float(sizes.iloc[:2].sum()/ni),
            block_equal_gain=float(g.groupby("spatial_block").brier_advantage.mean().mean()),
            signature_equal_gain=float(avg.mean()),
            signature_bootstrap_lower=float(np.quantile(reps,.025)),
            signature_bootstrap_upper=float(np.quantile(reps,.975)),
            fraction_signatures_positive=float((avg>0).mean()),
            remaining_islands_drop_top_two=int(len(rem)),
            block_equal_gain_drop_top_two=float(rem.groupby("spatial_block").brier_advantage.mean().mean()) if len(rem) else np.nan,
            components_any_donor_overlap=int(ncomp),
            largest_component_any_overlap=int(np.bincount(labels).max()),
            components_three_of_five_overlap=int(ncomp3),
            largest_component_three_of_five_overlap=int(largest3),
            largest_donor_component_islands=int(withcomp.groupby("component").size().max()),
            component_equal_gain=float(groupavg.mean()),
            fraction_donor_components_positive=float((groupavg>0).mean()),
        ))
    return pd.DataFrame(out,columns=COLUMNS)
