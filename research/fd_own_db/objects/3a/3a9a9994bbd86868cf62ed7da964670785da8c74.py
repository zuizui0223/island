"""Exploratory ecosystem-scale island floral FD versus experimental pollen limitation.

Island co-location does NOT mean the experiment measured one of the directly
typed plant species. No temporal pollinator decline, visit identity,
flower phenotype evolution or mediation is identified by this model.
Outcome is GloPL *hand pollen* response, not insect counts.

One island is one prediction/association unit. Publication×island cells are
first averaged, then each island gets equal weight. All uncertainties are
conditional spatial-block resamples, not independent evolutionary histories.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

PRIMARY=("source_delta_marginal_rao","source_delta_joint_rao")
SECONDARY=("source_delta_open_radial_share",)
NEEDED=(
    "island_id","analysis_regime","spatial_block",
    "log_distance_to_continent_km","log_island_area_km2",
    *PRIMARY,*SECONDARY,
    "pollen_limitation_island_mean","n_independent_GloPL_studies"
)


def _residual_beta(group:pd.DataFrame, endpoint:str) -> float:
    """OLS partial slope of a z-scored floral Q delta, equal island weight."""
    x=group.copy()
    z=pd.to_numeric(x[endpoint],errors="coerce").to_numpy(float)
    sd=float(z.std(ddof=0))
    if sd<1e-10 or not np.isfinite(sd):
        raise ValueError("No independent island variation in functional FD")
    z=(z-z.mean())/sd
    y=x.pollen_limitation_island_mean.to_numpy(float)
    d=x[["log_distance_to_continent_km","log_island_area_km2"]].to_numpy(float)
    d=(d-d.mean(axis=0))/np.maximum(d.std(axis=0),1e-10)
    regime=pd.get_dummies(x.analysis_regime.astype(str),dtype=float)
    cols=sorted(regime.columns)
    regime=regime.reindex(columns=cols[1:]).to_numpy(float)
    X=np.column_stack((np.ones(len(x)),z,d,regime))
    if np.linalg.matrix_rank(X)<X.shape[1]:
        raise ValueError("Original climate/regime source data not full rank")
    return float(np.linalg.lstsq(X,y,rcond=None)[0][1])


def audit_crosssectional_pollen_fd(
    island_fd:pd.DataFrame,
    pollen:pd.DataFrame,
    *,
    seed:int=20261113, block_draws:int=1999,
    min_islands:int=25, min_blocks:int=15,
) -> dict:
    """Pollen input must be original publication×island mean cells.

    Input columns: island_id,study_key,PL_Effect_Size_mean. Study record counts
    do not turn repeated rows or multiple study measurements into new islands.
    """
    if min_islands!=25 or min_blocks!=15 or block_draws<999:
        raise ValueError("Frozen ecological association support gate changed")
    if not set(NEEDED[:5]+PRIMARY+SECONDARY).issubset(island_fd):
        raise ValueError("Island source relative floral FD or geography missing")
    if not {"island_id","study_key","PL_Effect_Size_mean"}.issubset(pollen):
        raise ValueError("Independent experimental publication×island mean missing")
    if island_fd.island_id.duplicated().any():
        raise ValueError("Island floral FD must have one row per actual island")
    if pollen.duplicated(["island_id","study_key"]).any():
        raise ValueError("Multiple experiment cells from same original publication/island")
    if pollen[["island_id","study_key"]].isna().any().any():
        raise ValueError("Publication/island identity incomplete")
    p=pollen.copy()
    p["PL_Effect_Size_mean"]=pd.to_numeric(p.PL_Effect_Size_mean,errors="coerce")
    p=p.loc[np.isfinite(p.PL_Effect_Size_mean)].copy()
    ag=p.groupby("island_id",as_index=False).agg(
        pollen_limitation_island_mean=("PL_Effect_Size_mean","mean"),
        n_independent_GloPL_studies=("study_key","nunique"))
    work=island_fd.merge(ag,on="island_id",how="inner",validate="one_to_one")
    for key in NEEDED[3:]:
        work[key]=pd.to_numeric(work[key],errors="coerce")
    work=work.replace([np.inf,-np.inf],np.nan).dropna(subset=list(NEEDED))
    counts=work.groupby("analysis_regime").island_id.nunique().to_dict()
    nblocks=work.spatial_block.nunique()
    npublications=p.loc[p.island_id.isin(work.island_id),"study_key"].nunique()
    report={
        "status":"ISLAND_CROSSSECTIONAL_EXPERIMENTAL_PL_VS_FLORAL_FD",
        "n_islands_with_source_Q_and_finite_publication_effects":int(len(work)),
        "n_spatial_blocks":int(nblocks),
        "n_independent_original_publications":int(npublications),
        "n_original_publication_island_cells":int(len(p.loc[p.island_id.isin(work.island_id)])),
        "regions_n_islands":{k:int(v) for k,v in counts.items()},
        "same_flower_species_as_experiment_is_not_assumed":True,
        "temporal_decline_or_effective_visitor_identity_tested":False,
        "ecological_mediation_or_adaptation_identified":False,
        "original_source_taxa_sites_island_ids_exported":False,
    }
    if (len(work)<min_islands or nblocks<min_blocks or npublications<20
        or counts.get("northern_midlatitude",0)<10 or counts.get("tropical",0)<10):
        report.update(status="ECOLOGICAL_ASSOCIATION_NOT_IDENTIFIABLE",
            endpoints=[])
        return report
    results=[]
    for j,key in enumerate((*PRIMARY,*SECONDARY)):
        try:
            observed=_residual_beta(work,key)
        except ValueError:
            results.append(dict(endpoint=key,status="NONESTIMABLE_RANK_OR_CONSTANT"))
            continue
        blocks=sorted(work.spatial_block.unique())
        rng=np.random.default_rng(seed+j)
        estimates=[]
        for _ in range(block_draws):
            chosen=rng.choice(blocks,size=len(blocks),replace=True)
            d=pd.concat([work.loc[work.spatial_block.eq(b)] for b in chosen],
                        ignore_index=True)
            try:
                estimates.append(_residual_beta(d,key))
            except ValueError:
                pass
        leave=[]
        for b in blocks:
            other=work.loc[~work.spatial_block.eq(b)]
            try:
                leave.append(_residual_beta(other,key))
            except ValueError:
                pass
        if len(estimates)<.8*block_draws:
            results.append(dict(endpoint=key,
              status="BOOTSTRAP_IDENTIFIABILITY_FAILED",
              original_sd_scale_slope=observed,
              n_valid_block_draws=len(estimates)))
            continue
        lo,hi=np.quantile(estimates,[.025,.975])
        results.append(dict(endpoint=key,
          status="POSTHOC_CONDITIONAL_ASSOCIATION_NOT_CAUSAL",
          independent_primary_endpoint=key in PRIMARY,
          original_sd_scale_slope=observed,
          spatial_block_bootstrap_025=float(lo),
          spatial_block_bootstrap_975=float(hi),
          n_block_bootstrap_replicates_used=len(estimates),
          n_leave_one_block_results=len(leave),
          fraction_leave_block_out_same_slope_sign=(
            float(np.mean(np.sign(leave)==np.sign(observed))) if leave else None)))
    report["endpoints"]=results
    return report
