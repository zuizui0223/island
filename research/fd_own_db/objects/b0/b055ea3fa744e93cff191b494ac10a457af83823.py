"""Source-verified 4,505-GIFT-island floral trait-pair genus covariance by origin.

Geological island class (independently matched Roeble/Zizka) and last-glacial
mainland connectivity (Weigelt GMMC) are DIFFERENT variables and NEVER
substituted. Source-relative true SAME-species floral access composition
and exact within/between botanical genus pairing are compared descriptively
among observational island cohorts, with whole physical 10-degree block and
first-candidate-mainland-source resampling.

No source/genus component is historical colonization, actual effective
pollination or local island phenotypic evolution. Missing Direct traits, island
origin, source provenance or clusters cannot be filled as zero.
"""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd

FOCAL_PAIRS=("actinomorphic_deep","zygomorphic_deep","tubular_bell_deep")
AXES={
    "geological_origin":("continental","oceanic"),
    "LGM_connectivity":("lgm_connected","lgm_unconnected"),
}
REGIMES=("all","northern_midlatitude","tropical")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
MODES=("geo5","climate10")
METRICS=("joint","coupling","between_genus_coupling","within_genus_coupling")
MIN_DIRECT_ISLAND=5
MIN_DIRECT_SOURCE=10
MIN_CLASS_ISLANDS=30
MIN_INDEPENDENT_BLOCKS=10
MIN_INDEPENDENT_SOURCES=10


def original_geology_join(original:pd.DataFrame,crosswalk:pd.DataFrame,
                          *,population:int=4505)->pd.DataFrame:
    """Retain every original island; no distance proxy supplies geology."""
    needed={"island_id","geological_origin","LGM_connectivity",
            "area_km2","distance_to_continent_km"}
    if not needed<=set(crosswalk):
        raise ValueError("original independently verified geology/GMMC columns missing")
    if crosswalk.island_id.duplicated().any():
        raise ValueError("duplicate original GSHHG physical island geology")
    if original.island_id.nunique()!=population:
        raise ValueError("all original GIFT accepted floral island floras required")
    if "geological_origin" in original or "LGM_connectivity" in original:
        raise ValueError("geological origin not independently joined")
    if not crosswalk.geological_origin.fillna("unresolved").isin(
        (*AXES["geological_origin"],"mixed","unresolved")).all():
        raise ValueError("unknown original independent geology category")
    if not crosswalk.LGM_connectivity.fillna("unresolved").isin(
        (*AXES["LGM_connectivity"],"unresolved")).all():
        raise ValueError("unknown original GMMC mainland connectivity")
    include=["island_id","geological_origin","LGM_connectivity",
             "area_km2","distance_to_continent_km"]
    out=original.merge(crosswalk[include],on="island_id",how="left",
                       validate="many_to_one",indicator=True)
    if len(out)!=len(original) or out.island_id.nunique()!=population:
        raise AssertionError("original global botanical island population narrowed")
    # Missing original physical geography stays unresolved: never invent class.
    out["geological_origin"]=out.geological_origin.fillna("unresolved")
    out["LGM_connectivity"]=out.LGM_connectivity.fillna("unresolved")
    if "mixed" in AXES["geological_origin"]:
        raise AssertionError("mixed cannot be used as oceanic or continental")
    return out.drop(columns="_merge")


def _group_means(frame:pd.DataFrame, class_col:str,unit:str, cols:list[str]):
    # One vote per real 10-degree physical block within each class, or
    # per true first GIFT mainland donor-region identifier. A region/group
    # may contain several islands; no pseudo-independent island pairs.
    grouped=frame.dropna(subset=[unit]).groupby([unit,class_col],sort=True)[cols].mean()
    keys=sorted(set(grouped.index.get_level_values(0)))
    outcome={}
    for label in frame[class_col].unique():
        arr=np.full((len(keys),len(cols)),np.nan,dtype=float)
        for i,key in enumerate(keys):
            if (key,label) in grouped.index:
                arr[i]=grouped.loc[(key,label)].to_numpy(dtype=float)
        outcome[label]=arr
    return outcome


def _contrast(frame:pd.DataFrame,*,axis:str,categories:tuple[str,str],
              unit:str,columns:list[str],seed:str,draws:int)->dict:
    arrays=_group_means(frame,axis,unit,columns)
    a,b=categories
    arr_a=arrays[a];arr_b=arrays[b]
    na=int(np.isfinite(arr_a[:,0]).sum())
    nb=int(np.isfinite(arr_b[:,0]).sum())
    out=dict(unit=unit,n_clusters_category_a=na,n_clusters_category_b=nb,
             n_union_clusters=len(arr_a))
    if min(na,nb)<MIN_INDEPENDENT_BLOCKS:
        out["status"]="NOT_IDENTIFIABLE_ORIGIN_CLUSTER_SUPPORT"
        return out
    m_a=np.nanmean(arr_a,axis=0)
    m_b=np.nanmean(arr_b,axis=0)
    obs=m_b-m_a
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    draw_ids=rng.integers(len(arr_a),size=(draws,len(arr_a)))
    # Joint resampling of the UNION of exact original source clusters
    # preserves correlations where one source block contains both labels.
    A=arr_a[draw_ids]
    B=arr_b[draw_ids]
    valid_a=np.isfinite(A[:,:,0]).sum(axis=1)>0
    valid_b=np.isfinite(B[:,:,0]).sum(axis=1)>0
    valid=valid_a&valid_b
    if int(valid.sum())<int(draws*.99):
        out.update(status="NOT_IDENTIFIABLE_ORIGIN_BOOTSTRAP_COVERAGE",
                   n_bootstrap_valid=int(valid.sum()))
        return out
    diff=np.nanmean(B[valid],axis=1)-np.nanmean(A[valid],axis=1)
    for j,col in enumerate(columns):
        out[col+"__category_a_mean"]=float(m_a[j])
        out[col+"__category_b_mean"]=float(m_b[j])
        out[col+"__category_b_minus_a"]=float(obs[j])
        out[col+"__bootstrap95"]=[
            float(v) for v in np.quantile(diff[:,j],[.025,.975])]
    for mode in MODES:
        pref=lambda metric:f"{mode}__delta_{metric}"
        lhs=out[pref("coupling")+"__category_b_minus_a"]
        rhs=(out[pref("between_genus_coupling")+"__category_b_minus_a"]+
             out[pref("within_genus_coupling")+"__category_b_minus_a"])
        if abs(lhs-rhs)>1e-10:
            raise AssertionError("exact origin-class flower genus partition lost")
    out.update(status="DESCRIPTIVE_GIFT_ORIGIN_SOURCE_RELATIVE_FLORAL_PAIRING",
               n_bootstrap_valid=int(valid.sum()))
    return out


def original_pair_by_geology(frame:pd.DataFrame,*,population:int=4505,
                             draws:int=999)->dict:
    """All 4,505 island×3 scope×11 pair rows; 3 biological FOCAL pairs."""
    needed={"island_id","regime","scope","sector","spatial_block","first_source",
            "island__n_direct","island__n_variable_genera","island_total",
            *AXES}
    for m in MODES:
        needed|={f"{m}__n_direct",f"{m}__total"}
        needed|={f"{m}__delta_{metric}" for metric in METRICS}
    if missing:=needed-set(frame):
        raise ValueError(f"original source/geo phenotype fields missing: {sorted(missing)}")
    if (frame.island_id.nunique()!=population or
        len(frame)!=population*3*11 or
        frame.duplicated(["island_id","scope","sector"]).any()):
        raise ValueError("original full global flora/pair provenance lost")
    if draws<199:
        raise ValueError("insufficient whole-source-cluster bootstrap")
    panels=[]
    columns=[f"{mode}__delta_{m}" for mode in MODES for m in METRICS]
    for axis,categories in AXES.items():
      for region in REGIMES:
       for scope in SCOPES:
        for sector in FOCAL_PAIRS:
            data=frame.loc[frame.scope.eq(scope)&frame.sector.eq(sector)]
            if region!="all":data=data.loc[data.regime.eq(region)]
            all_original=int(data.island_id.nunique())
            eligible=data.loc[data[axis].isin(categories)].copy()
            q=eligible.loc[eligible.island__n_direct.ge(MIN_DIRECT_ISLAND)]
            for mode in MODES:
                q=q.loc[q[f"{mode}__n_direct"].ge(MIN_DIRECT_SOURCE)]
            q=q.dropna(subset=columns)
            # Do not compute a categorical origin contrast if the smallest
            # class lacks both independent source and physical-block support.
            per_category={}
            for label in categories:
                sub=q.loc[q[axis].eq(label)]
                n_blocks=int(sub.spatial_block.nunique())
                n_first=int(sub.first_source.nunique())
                per_category[label]=dict(
                    n_islands=int(len(sub)),
                    n_spatial_blocks=n_blocks,n_first_sources=n_first,
                    n_islands_with_some_real_variable_genus=int(
                        sub.island__n_variable_genera.gt(0).sum()))
            if min(per_category[x]["n_islands"] for x in categories)<MIN_CLASS_ISLANDS or (
                min(per_category[x]["n_spatial_blocks"] for x in categories)<MIN_INDEPENDENT_BLOCKS or
                min(per_category[x]["n_first_sources"] for x in categories)<MIN_INDEPENDENT_SOURCES):
                status="NOT_IDENTIFIABLE_ORIGIN_CLASS_SOURCE_SUPPORT"
            else:
                status="DESCRIPTIVE_ORIGIN_CLASS_BOTH_SOURCE_MODES_SUPPORTED"
            entry=dict(axis=axis,categories=list(categories),
                       region=region,scope=scope,sector=sector,
                       n_original_botanical_islands_in_region=all_original,
                       n_independently_classified_island_floras=int(len(eligible)),
                       n_source_pair_eligible_with_both_modes=int(len(q)),
                       class_support=per_category,status=status)
            for grouping in ("spatial_block","first_source"):
                if status.startswith("NOT_IDENTIFIABLE"):
                    entry[grouping]=dict(status=status)
                else:
                    entry[grouping]=_contrast(
                        q,axis=axis,categories=categories,unit=grouping,
                        columns=columns,seed=f"{axis}|{region}|{scope}|{sector}|{grouping}",
                        draws=draws)
            panels.append(entry)
    return dict(
        status="GIFT_4505_GEOLOGICAL_ORIGIN_AND_LGM_SEPARATE_SOURCE_GENUS_PAIRING",
        original_global_botanical_island_floras=population,
        n_original_island_status_true_pair_rows=len(frame),
        n_independent_GIFT_islands_with_oceanic_or_continental_geology=int(
            frame.loc[frame.scope.eq("all_recorded")&
                      frame.sector.eq(FOCAL_PAIRS[0])&
                      frame.geological_origin.isin(AXES["geological_origin"]),"island_id"].nunique()),
        n_original_islands_with_LGM_mainland_connectivity=int(
            frame.loc[frame.scope.eq("all_recorded")&
                      frame.sector.eq(FOCAL_PAIRS[0])&
                      frame.LGM_connectivity.isin(AXES["LGM_connectivity"]),"island_id"].nunique()),
        n_report_panels=len(panels),panels=panels,
        identification_limits=dict(
            geological_oceanic_continental_independent_publication_only=True,
            LGM_connectivity_never_substituted_for_geological_origin=True,
            unresolved_geology_retained_in_4505_population_but_not_classes=True,
            source_candidate_not_historical_colonizer=True,
            genus_pairing_partition_exact_not_causes=True,
            sampling_origin_and_direct_coverage_nonrandom=True,
            source_block_and_first_source_both_require_support=True,
            no_pollinator_guild_classification_or_effective_pollination=True,
            no_island_adaptation_identified=True,
            original_4505_botanical_parent_never_shrunk=True,
        ))
