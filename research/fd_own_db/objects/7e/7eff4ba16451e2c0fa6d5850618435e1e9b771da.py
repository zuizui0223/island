"""Joint geology×LGM source-overlap/positivity audit on ALL 4505 GIFT floras.

The independently Roeble/Zizka classified geological origin of a GIFT island
does NOT equal its separately Weigelt-2013 mainland glacial connection.
One physical island may be oceanic but historically LGM-connected according
to an imperfect crosswalk; this is audited, never silently relabelled.

For a source-relative TRUE same-plant pair morphology:
- report the real original island/support/variance counts in all 2×2 cells;
- require eligible flora >=30, >=10 10° blocks and >=10 first GIFT source
  groups in BOTH cells before estimating ONE within-stratum contrast;
- require the FOUR cells to meet these gates before an equal-strata adjusted
  geological/LGM effect or interaction can be identified.
- bootstrap the same UNION of physical/source groups across all four cells,
  not islands as independent observations.

No causal history, native status, local pollinator or adaptive origin claim.
"""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd

from .geology_genus_flower_pair import (
    AXES,FOCAL_PAIRS,MODES,METRICS,MIN_CLASS_ISLANDS,
    MIN_INDEPENDENT_BLOCKS,MIN_INDEPENDENT_SOURCES,
    MIN_DIRECT_ISLAND,MIN_DIRECT_SOURCE,_contrast,
)

REGIONS=("all","northern_midlatitude","tropical")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
GEO=("continental","oceanic")
LGM=("lgm_connected","lgm_unconnected")
CELLS=tuple((g,l) for g in GEO for l in LGM)
COLUMNS=[f"{m}__delta_{t}" for m in MODES for t in METRICS]


def _qualified(q:pd.DataFrame)->dict:
    return dict(
        n_islands=int(q.island_id.nunique()),
        n_geographic_blocks=int(q.spatial_block.nunique()),
        n_first_mainland_sources=int(q.first_source.nunique()),
        n_within_genus_variable_flora_islands=int(
            q.island__n_variable_genera.gt(0).sum()),
        n_same_informative_genus_island_and_both_sources=int(
            q.n_shared_informative_genera_both_source_modes.gt(0).sum()),
    )


def _passes(x:dict)->bool:
    return (x["n_islands"]>=MIN_CLASS_ISLANDS and
            x["n_geographic_blocks"]>=MIN_INDEPENDENT_BLOCKS and
            x["n_first_mainland_sources"]>=MIN_INDEPENDENT_SOURCES)


def _four_cell_adjusted(q:pd.DataFrame,*,cluster:str,seed:str,
                        draws:int)->dict:
    """Two-way equal-strata contrast on FOUR documented cross-class cells."""
    g=q.dropna(subset=[cluster]).groupby(
        [cluster,"geological_origin","LGM_connectivity"],sort=True)[COLUMNS].mean()
    clusters=sorted(q[cluster].dropna().unique())
    X=np.full((len(clusters),4,len(COLUMNS)),np.nan)
    for i,k in enumerate(clusters):
        for j,(geol,l) in enumerate(CELLS):
            if (k,geol,l) in g.index:
                X[i,j,:]=g.loc[(k,geol,l)].to_numpy(float)
    if not np.isfinite(X).any():
        return dict(status="NOT_IDENTIFIABLE_JOINT_CLASS_SOURCE_GROUPS")
    means=np.nanmean(X,axis=0)
    if not np.isfinite(means).all():
        return dict(status="NOT_IDENTIFIABLE_JOINT_CLASS_SOURCE_GROUPS")
    # index cell order: Cont-connected, Cont-unconnected,
    # oceanic-connected, oceanic-unconnected.
    def contrasts(m):
        c_c,c_u,o_c,o_u=m
        origin=((o_c-c_c)+(o_u-c_u))/2
        lgm=((c_u-c_c)+(o_u-o_c))/2
        interaction=(o_u-c_u)-(o_c-c_c)
        return np.stack((origin,lgm,interaction),axis=0)
    estimates=contrasts(means)
    rng=np.random.default_rng(int.from_bytes(
        sha256(seed.encode()).digest()[:8],"big"))
    idx=rng.integers(len(clusters),size=(draws,len(clusters)))
    with np.errstate(invalid="ignore"):
        sampled=X[idx]
        cnt=np.isfinite(sampled).sum(axis=1)
        valid=(cnt[:,:,0]>0).all(axis=1)
    if valid.sum()<draws*.99:
        return dict(status="NOT_IDENTIFIABLE_FOUR_CELL_CLUSTER_BOOTSTRAP",
                    n_clusters=int(len(clusters)),
                    n_valid_cluster_bootstraps=int(valid.sum()))
    boot=np.nanmean(sampled[valid],axis=1)
    effects=contrasts(boot.transpose(1,0,2)).transpose(1,0,2)
    names=("oceanic_minus_continental_equal_LGM_strata",
           "LGM_unconnected_minus_connected_equal_geology_strata",
           "geology_by_LGM_interaction")
    result=dict(status="DESCRIPTIVE_FULL_FOUR_CELL_GEOL_LGM_OVERLAP",
                n_union_independent_clusters=int(len(clusters)),
                n_valid_cluster_bootstraps=int(valid.sum()))
    for ei,ename in enumerate(names):
        for ci,col in enumerate(COLUMNS):
            result[ename+"__"+col]=float(estimates[ei,ci])
            result[ename+"__"+col+"_bootstrap95"]=[
                float(v) for v in np.quantile(effects[:,ei,ci],[.025,.975])]
    return result


def audit_geology_lgm_common_support(frame:pd.DataFrame, *,
                                     population:int=4505,
                                     draws:int=999)->dict:
    need={"island_id","regime","scope","sector","geological_origin",
          "LGM_connectivity","spatial_block","first_source",
          "island_total","island__n_direct","island__n_variable_genera",
          "n_shared_informative_genera_both_source_modes",
          *COLUMNS}
    for m in MODES:need|={f"{m}__n_direct",f"{m}__total"}
    if missing:=need-set(frame):
        raise ValueError("missing independent original GIFT source cross-history fields "+str(sorted(missing)))
    if (frame.island_id.nunique()!=population or
        len(frame)!=population*3*11 or
        frame.duplicated(["island_id","scope","sector"]).any()):
        raise ValueError("original complete 4505-island species Direct pair denominator contracted")
    if not frame.geological_origin.fillna("unresolved").isin(
        ("oceanic","continental","mixed","unresolved")).all():
        raise ValueError("independent geology categorization corrupted")
    if not frame.LGM_connectivity.fillna("unresolved").isin(
        ("lgm_connected","lgm_unconnected","unresolved")).all():
        raise ValueError("independent glacial GMMC source label corrupted")
    if draws<199:raise ValueError("fewer than 199 original source bootstraps")
    panels=[]
    for region in REGIONS:
      for scope in SCOPES:
       for sector in FOCAL_PAIRS:
        subset=frame.loc[frame.scope.eq(scope)&frame.sector.eq(sector)].copy()
        if region!="all":subset=subset.loc[subset.regime.eq(region)]
        q=subset.loc[subset.island__n_direct.ge(MIN_DIRECT_ISLAND)]
        for mode in MODES:
            q=q.loc[q[f"{mode}__n_direct"].ge(MIN_DIRECT_SOURCE)]
        q=q.dropna(subset=COLUMNS)
        both=q.loc[q.geological_origin.isin(GEO)&q.LGM_connectivity.isin(LGM)].copy()
        group_support={}
        for geol,l in CELLS:
            cell=both.loc[both.geological_origin.eq(geol)&
                          both.LGM_connectivity.eq(l)]
            group_support[f"{geol}__{l}"]=_qualified(cell)
        same_cohort_counts=dict(
            n_region_original_global_islands=int(subset.island_id.nunique()),
            n_original_pair_source_comparable_before_two_history_labels=int(q.island_id.nunique()),
            n_both_geological_and_LGM_history_qualified=int(both.island_id.nunique()),
            n_geologically_classified_but_LGM_unknown=int(q.loc[
                q.geological_origin.isin(GEO)&~q.LGM_connectivity.isin(LGM),
                "island_id"].nunique()),
            n_LGM_classified_but_geology_unknown=int(q.loc[
                ~q.geological_origin.isin(GEO)&q.LGM_connectivity.isin(LGM),
                "island_id"].nunique()),
            n_same_variable_genus_on_island_and_both_sources=int(
                both.n_shared_informative_genera_both_source_modes.gt(0).sum()),
        )
        outcomes={}
        stratifications=(
            ("geological_origin","LGM_connectivity",LGM,GEO),
            ("LGM_connectivity","geological_origin",GEO,LGM),
        )
        for target_axis,stratum_axis,stratum_values,categories in stratifications:
            for v in stratum_values:
                sub=both.loc[both[stratum_axis].eq(v)]
                cross=[_qualified(sub.loc[sub[target_axis].eq(cat)])
                       for cat in categories]
                name=f"{target_axis}__within__{stratum_axis}__{v}"
                if not all(map(_passes,cross)):
                    outcomes[name]=dict(
                        status="NOT_IDENTIFIABLE_COMMON_CLASS_SOURCE_OR_BLOCK_SUPPORT",
                        categories=list(categories),
                        category_support=cross,
                    )
                    continue
                res=dict(status="DESCRIPTIVE_HISTORY_CONDITIONAL_BOTANICAL_SOURCE_FLORA",
                         categories=list(categories),
                         category_support=cross)
                for grouping in ("spatial_block","first_source"):
                    stat=_contrast(sub,axis=target_axis,categories=categories,
                        unit=grouping,columns=COLUMNS,
                        seed=f"{region}|{scope}|{sector}|{name}|{grouping}",
                        draws=draws)
                    res[grouping]=stat
                outcomes[name]=res
        joint_ok=all(_passes(group_support[f"{geol}__{l}"])
                     for geol,l in CELLS)
        joint=dict(status=("DESCRIPTIVE_FOUR_CLASS_SUPPORT_BEFORE_BOOTSTRAP"
                            if joint_ok else "NOT_IDENTIFIABLE_ALL_FOUR_HISTORY_CELLS"),
                   n_four_eligible_cells=sum(_passes(group_support[f"{geol}__{l}"])
                                             for geol,l in CELLS))
        if joint_ok:
            for grouping in ("spatial_block","first_source"):
                joint[grouping]=_four_cell_adjusted(
                    both,cluster=grouping,
                    seed=f"joint|{region}|{scope}|{sector}|{grouping}",
                    draws=draws)
        panels.append(dict(
            region=region,scope=scope,sector=sector,
            **same_cohort_counts,
            original_four_joint_history_cell_support=group_support,
            within_one_history_class_contrasts=outcomes,
            simultaneous_origin_and_LGM_equal_strata_adjustment=joint,
        ))
    return dict(
        status="GIFT4505_INDEPENDENT_GEOLOGY_X_LGM_FOUR_CELL_IDENTIFIABILITY",
        n_original_global_botanical_floras=population,
        n_original_full_flora_scope_pair_rows=len(frame),
        n_reporting_region_scope_flower_panels=len(panels),
        n_history_stratified_panels_supported=sum(
            x["status"]=="DESCRIPTIVE_HISTORY_CONDITIONAL_BOTANICAL_SOURCE_FLORA"
            for p in panels for x in p["within_one_history_class_contrasts"].values()),
        n_both_history_axes_jointly_supported=sum(
            p["simultaneous_origin_and_LGM_equal_strata_adjustment"]["status"]==
            "DESCRIPTIVE_FOUR_CLASS_SUPPORT_BEFORE_BOOTSTRAP"
            for p in panels),
        panels=panels,
        limits=dict(
            geological_origin_Roeble_Zizka_independent_not_distance_proxy=True,
            glacial_connection_Weigelt_distinct_from_geological_origin=True,
            all_four_cells_need_independent_source_and_block_support=True,
            conditional_history_contrasts_observational_not_causal=True,
            source_species_identity_not_historical_colonization=True,
            original_native_status_nonrandomly_documented=True,
            no_flower_pollinator_syndrome_or_guild_shift_observed=True,
            no_source_species_archipelago_provider_outer_holdout=True,
            all_original_4505_floras_preserved=True,
        ))
