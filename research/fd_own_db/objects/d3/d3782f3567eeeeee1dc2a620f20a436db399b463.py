"""Score independently held-out floral predictions; never infer guild absence."""
from __future__ import annotations
import numpy as np
import pandas as pd

COLUMNS = ("island_id","archipelago","spatial_block","first_source",
           "trait_provider","n_target","n_positive","p_baseline",
           "p_genus","p_guild","guild_supported",
           "target_species_heldout","archipelago_heldout","source_provider_heldout")

def evaluate(frame: pd.DataFrame, *, draws: int=999, min_islands: int=30,
             min_groups: int=10) -> dict:
    """Evaluate precomputed, out-of-archipelago predictions only.

    Caller is responsible for target-species, source-provider, and archipelago
    holdout during model training. This evaluates predictions, not causal effects.
    """
    if not set(COLUMNS).issubset(frame):
        raise ValueError("missing columns")
    if frame.empty or frame.island_id.isna().any() or frame.island_id.duplicated().any():
        raise ValueError("one row per unique island required")
    for c in ("archipelago","spatial_block","first_source","trait_provider"):
        if frame[c].isna().any() or frame[c].astype(str).str.strip().eq("").any():
            raise ValueError("missing independent group")
    if not frame.guild_supported.map(lambda x:isinstance(x,(bool,np.bool_))).all():
        raise ValueError("guild_supported must be an explicit boolean")
    for c in ("target_species_heldout","archipelago_heldout","source_provider_heldout"):
        if not frame[c].map(lambda x:isinstance(x,(bool,np.bool_)) and bool(x)).all():
            raise ValueError(f"out-of-fold provenance not proven: {c}")
    if draws<199 or min_islands<1 or min_groups<2:
        raise ValueError("invalid support thresholds")
    n=pd.to_numeric(frame.n_target,errors="coerce").to_numpy(float)
    k=pd.to_numeric(frame.n_positive,errors="coerce").to_numpy(float)
    if not np.isfinite(n+k).all() or np.any(n<1) or np.any(n!=np.floor(n)) or np.any(k<0) or np.any(k>n) or np.any(k!=np.floor(k)):
        raise ValueError("invalid species counts")
    d=frame.copy()
    for c in ("p_baseline","p_genus","p_guild"):
        p=pd.to_numeric(frame[c],errors="coerce").to_numpy(float)
        # No guild measurement means no guild model, not a zero effect.
        eligible=d.guild_supported.to_numpy(dtype=bool) if c=="p_guild" else np.ones(len(d),bool)
        if not np.isfinite(p[eligible]).all() or np.any(p[eligible]<0) or np.any(p[eligible]>1):
            raise ValueError("invalid probability on supported model rows")
        if c=="p_guild" and np.isfinite(p[~eligible]).any():
            raise ValueError("unsupported guild rows must not contain predictions")
        d[c+"_brier"]=(k*(1-p)**2+(n-k)*p**2)/n
    d["gain_genus"]=d.p_baseline_brier-d.p_genus_brier
    d["gain_guild"]=d.p_genus_brier-d.p_guild_brier
    def summarize(sub, guild):
        out={"n_islands":len(sub),"n_archipelagos":sub.archipelago.nunique(),
             "n_blocks":sub.spatial_block.nunique(),
             "n_sources":sub.first_source.nunique(),
             "n_providers":sub.trait_provider.nunique()}
        if len(sub)<min_islands or min(out["n_archipelagos"],out["n_blocks"],out["n_sources"])<min_groups or (guild and out["n_providers"]<2):
            out["status"]="NOT_IDENTIFIABLE_INDEPENDENT_SUPPORT"
            return out
        v=sub.groupby("archipelago",sort=True)[["gain_genus","gain_guild"]].mean().to_numpy()
        avg=v.mean(axis=0)
        rng=np.random.default_rng(20261117+int(guild))
        boot=v[rng.integers(len(v),size=(draws,len(v)))].mean(axis=1)
        for j,c in enumerate(("gain_genus","gain_guild") if guild else ("gain_genus",)):
            out[c]=float(avg[j])
            out[c+"_bootstrap95"]=[float(x) for x in np.quantile(boot[:,j],[.025,.975])]
        out["status"]="DESCRIPTIVE_HELDOUT_PREDICTIVE_COMPARISON"
        return out
    return {"n_input_islands":len(d),"source_genus":summarize(d,False),
            "guild_subset":summarize(d.loc[d.guild_supported],True),
            "not_causal":True,"no_guild_absence_inferred":True}
