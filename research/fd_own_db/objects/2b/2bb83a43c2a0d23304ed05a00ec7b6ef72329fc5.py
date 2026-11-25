"""Paired, same-island contrast of source-genus predictive gains.

Inputs must be predictions fitted outside focal accepted species, archipelago
and phenotype provider. This is a scorer, not a model fitter or causal test.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

REQUIRED=("island_id","archipelago","spatial_block","first_source",
          "trait_provider","sector","n_target","n_positive",
          "p_pooled","p_genus","target_species_heldout",
          "archipelago_heldout","source_provider_heldout")
SECTORS=("actinomorphic_deep","tubular_bell_deep")

def paired_source_genus_gain(frame:pd.DataFrame, *, draws:int=999,
                             min_islands:int=30,min_groups:int=10)->dict:
    """Compare per-island mean Brier gains on the EXACT common island set.

    Equal-weight groups are archipelagos; no island-pair pseudo-replication.
    If group support is too small, no estimate is reported.
    """
    missing=set(REQUIRED)-set(frame)
    if missing:raise ValueError(f"missing fields {sorted(missing)}")
    if frame.empty or draws<199 or min_islands<1 or min_groups<2:
        raise ValueError("invalid cohort or thresholds")
    if frame.island_id.isna().any() or frame.duplicated(["island_id","sector"]).any():
        raise ValueError("island-sector duplicate or missing island")
    if not frame.sector.isin(SECTORS).all():
        raise ValueError("unknown flower architecture")
    for col in ("archipelago","spatial_block","first_source","trait_provider"):
        if frame[col].isna().any() or frame[col].astype(str).str.strip().eq("").any():
            raise ValueError("missing provenance group")
    for col in ("target_species_heldout","archipelago_heldout","source_provider_heldout"):
        if not frame[col].map(lambda v:isinstance(v,(bool,np.bool_)) and bool(v)).all():
            raise ValueError(f"unverified independent holdout: {col}")
    n=pd.to_numeric(frame.n_target,errors="coerce").to_numpy(float)
    k=pd.to_numeric(frame.n_positive,errors="coerce").to_numpy(float)
    if not np.isfinite(n+k).all() or np.any(n<1) or np.any(n!=np.floor(n)) or np.any(k<0) or np.any(k>n) or np.any(k!=np.floor(k)):
        raise ValueError("invalid target phenotype counts")
    d=frame.copy()
    for c in ("p_pooled","p_genus"):
        p=pd.to_numeric(d[c],errors="coerce").to_numpy(float)
        if not np.isfinite(p).all() or np.any((p<0)|(p>1)):
            raise ValueError("invalid heldout probability")
        d[c+"_brier"]=(k*(1-p)**2+(n-k)*p**2)/n
    d["gain"]=d.p_pooled_brier-d.p_genus_brier
    ids=d.groupby("island_id").sector.nunique()
    common=ids.index[ids.eq(2)]
    paired=d.loc[d.island_id.isin(common)].copy()
    # Same island has the same source/archipelago/provider for both endpoints.
    for col in ("archipelago","spatial_block","first_source","trait_provider"):
        if paired.groupby("island_id")[col].nunique().gt(1).any():
            raise ValueError(f"non-comparable endpoint provenance: {col}")
    wide=paired.pivot(index="island_id",columns="sector",values="gain")
    delta=(wide[SECTORS[0]]-wide[SECTORS[1]]).rename("delta")
    info=paired.drop_duplicates("island_id").set_index("island_id")
    z=info[["archipelago","spatial_block","first_source","trait_provider"]].join(delta)
    out={"n_original_islands":int(d.island_id.nunique()),
         "n_common_islands":int(len(z)),
         "n_archipelagos":int(z.archipelago.nunique()),
         "n_spatial_blocks":int(z.spatial_block.nunique()),
         "n_first_sources":int(z.first_source.nunique()),
         "n_providers":int(z.trait_provider.nunique())}
    if len(z)<min_islands or min(out["n_archipelagos"],out["n_spatial_blocks"],out["n_first_sources"])<min_groups or out["n_providers"]<2:
        out["status"]="NOT_IDENTIFIABLE_PAIRED_INDEPENDENT_SUPPORT"
        return out
    v=z.groupby("archipelago",sort=True).delta.mean().to_numpy(float)
    rng=np.random.default_rng(20261117)
    bs=v[rng.integers(len(v),size=(draws,len(v)))].mean(axis=1)
    out.update(status="DESCRIPTIVE_PAIRED_OUT_OF_FOLD_CONTRAST",
               delta_gain_radial_minus_tubular=float(v.mean()),
               bootstrap95=[float(x) for x in np.quantile(bs,[.025,.975])],
               note="Group-bootstrap conditional interval, not a causal or preregistered test")
    return out
