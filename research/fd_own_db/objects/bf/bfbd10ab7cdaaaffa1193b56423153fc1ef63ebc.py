"""Outcome-blind genus concentration controls for NHM numeric congener prediction.

Input contains PRIVATE genus-resolved island-scope-source sufficient statistics,
NEVER original accepted species names in public output. Concentration rankings
use counts of original Direct-documented focal plants ONLY, not prediction
accuracy or bird/butterfly guild outcomes. Sensitivity analyses are posthoc
exploratory and retain the full 4,505-island population denominator.
"""
from __future__ import annotations
from hashlib import sha256
import numpy as np
import pandas as pd

REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
TIERS=("any_valid_NHM_numeric_tube_length",
       "two_source_corpora_consistent_tube_length")
MODES=("geo5","climate10")
MIN_TARGET=5
MIN_GENERA=2
MIN_DONORS=10
MIN_ISLANDS=30
MIN_BLOCKS=10


def score_genus_sufficient_stats(group:pd.DataFrame)->dict:
    """Recompute two source baselines and Brier scores on same target species.

    A source genus g supplies only (numeric donor count, total source genus
    checklist size, mean probability from held-g-out NHM-cm calibrator).
    The target supplies only its original n observed ordinal depth labels and
    n positive 'deep'. No direct trait target enters prediction.
    """
    n=int(group.n_target.sum())
    ng=int(len(group))
    d={"n_target":n,"n_genera":ng,
       "largest_genus_carrier_fraction":float(group.n_target.max()/n) if n else None,
       "effective_genus_number":float(n*n/(group.n_target.astype(float)**2).sum())
       if n else None}
    if ng==0 or n==0:return {**d,"status":"NOT_IDENTIFIABLE_DONOR_SUPPORT"}
    for mode in MODES:
        donors=int(group[mode+"__n_numeric"].sum())
        source_total=int(group[mode+"__n_original_source_genus"].sum())
        d[mode+"__n_numeric_donors"]=donors
        d[mode+"__n_source_species_in_supported_genera"]=source_total
        d[mode+"__n_genera_with_two_numeric_donors"]=int(
            group[mode+"__n_numeric"].ge(2).sum()
        )
    if n<MIN_TARGET or ng<MIN_GENERA or min(
        d[mode+"__n_numeric_donors"] for mode in MODES
    )<MIN_DONORS:
        return {**d,"status":"NOT_IDENTIFIABLE_DONOR_SUPPORT"}
    ngs=group.n_target.to_numpy(float)
    ys=group.n_target_deep.to_numpy(float)
    d["island_documentation_fraction"]=float(n/np.sum(
        group.n_target_original_genus.to_numpy(float)
    ))
    for mode in MODES:
        dn=group[mode+"__n_numeric"].to_numpy(float)
        p=group[mode+"__genus_pdeep"].to_numpy(float)
        pooled=float(np.sum(dn*p)/np.sum(dn))
        pooled_genus_equal=float(np.mean(p))
        exact=lambda pred: ys - 2.*pred*ys + ngs*pred**2
        e_genus=exact(p)
        e_pooled=exact(pooled)
        e_pool_equal=exact(pooled_genus_equal)
        d[mode+"__brier_gain_species_equal"]=float(
            np.sum(e_pooled-e_genus)/n
        )
        d[mode+"__brier_gain_genus_equal"]=float(
            np.mean((e_pooled-e_genus)/ngs)
        )
        d[mode+"__brier_gain_donor_genus_equal_pooled"]=float(
            np.sum(e_pool_equal-e_genus)/n
        )
        denom=d[mode+"__n_source_species_in_supported_genera"]
        d[mode+"__source_documentation_fraction"]=float(
            d[mode+"__n_numeric_donors"]/denom
        )
    d["documentation_10pp_both"]=bool(all(
        abs(d["island_documentation_fraction"]-
            d[mode+"__source_documentation_fraction"])<=.10
        for mode in MODES
    ))
    d["source_2plus_genus_both"]=bool(all(
        d[mode+"__n_genera_with_two_numeric_donors"]>=2 for mode in MODES
    ))
    d["status"]="ELIGIBLE_PAIRED_SOURCE_GENUS_CONCENTRATION"
    return d


def _cluster_summary(rows:pd.DataFrame,*,group:str,seed:str,bootstrap:int)->dict:
    grouped=rows.dropna(subset=[group]).groupby(group,sort=True)
    r={"n_original_floras":int(len(rows)),"n_geography_or_source_groups":int(grouped.ngroups)}
    if len(rows)<MIN_ISLANDS or grouped.ngroups<MIN_BLOCKS:
        r["status"]="NOT_IDENTIFIABLE_GROUP_SUPPORT"
        return r
    cols=[mode+"__"+key for mode in MODES for key in (
        "brier_gain_species_equal","brier_gain_genus_equal",
        "brier_gain_donor_genus_equal_pooled",
    )]
    vals=grouped[cols].mean().to_numpy(float)
    mean=vals.mean(0)
    rng=np.random.default_rng(int.from_bytes(
        sha256(seed.encode()).digest()[:8],"big"))
    boot=vals[rng.integers(len(vals),size=(bootstrap,len(vals)))].mean(1)
    for i,name in enumerate(cols):
        r[name+"__block_mean"]=float(mean[i])
        r[name+"__bootstrap95"]=[
            float(x) for x in np.quantile(boot[:,i],[.025,.975])
        ]
    r["status"]="DESCRIPTIVE_GENUS_CONCENTRATION_SENSITIVITY"
    return r


def audit_genus_concentration(
    private:pd.DataFrame,roster:pd.DataFrame,*,population=4505,bootstrap=999
)->dict:
    """Private genus/plant labels remain in source CI; only aggregates leave it."""
    need={"island_id","regime","scope","tier","genus",
          "spatial_block","first_source","n_target",
          "n_target_deep","n_target_original_genus"}
    for mode in MODES:
        need|={mode+"__n_numeric",mode+"__n_original_source_genus",
               mode+"__genus_pdeep"}
    if m:=need-set(private):
        raise ValueError(f"missing protected genus-score field {sorted(m)}")
    if private.duplicated(["island_id","scope","tier","genus"]).any():
        raise ValueError("duplicate original island scope tier genus")
    if not private.regime.isin(REGIMES).all() or (
        not private.scope.isin(SCOPES).all()) or (
        not private.tier.isin(TIERS).all()):
        raise ValueError("original flora region/nativeness/NHM source tier drift")
    if roster.island_id.nunique()!=population or roster.island_id.duplicated().any():
        raise ValueError("original complete 4505 island roster contracted")
    if not private.island_id.isin(roster.island_id).all():
        raise ValueError("source genus rows outside full original botanical islands")
    for col in ("n_target","n_target_deep","n_target_original_genus"):
        vals=pd.to_numeric(private[col],errors="raise")
        if (vals<0).any() or (vals%1!=0).any():
            raise ValueError(f"invalid botanical genus target counts {col}")
    if ((private.n_target_deep>private.n_target)|
        (private.n_target>private.n_target_original_genus)).any():
        raise ValueError("known original target Direct depth cardinalities inconsistent")
    for mode in MODES:
        n=private[mode+"__n_numeric"]
        tot=private[mode+"__n_original_source_genus"]
        p=private[mode+"__genus_pdeep"]
        if (n<=0).any() or (tot<n).any() or p.isna().any() or (
            ~p.between(0,1)).any():
            raise ValueError("source numeric genus donor record invalid")
    if bootstrap<199:raise ValueError("insufficient original geographical uncertainty draws")

    result=[]
    for region in REGIMES:
        for scope in SCOPES:
            for tier in TIERS:
                part=private.loc[
                    private.regime.eq(region)&private.scope.eq(scope)&
                    private.tier.eq(tier)
                ].copy()
                grouped={i:q for i,q in part.groupby("island_id",sort=True)}
                primary=[]
                for island,q in grouped.items():
                    scored=score_genus_sufficient_stats(q)
                    if scored["status"]!="ELIGIBLE_PAIRED_SOURCE_GENUS_CONCENTRATION":
                        continue
                    row={k:scored[k] for k in scored}
                    row.update(
                        island_id=island,regime=region,scope=scope,tier=tier,
                        spatial_block=q.spatial_block.iloc[0],
                        first_source=q.first_source.iloc[0],
                    )
                    primary.append(row)
                base=pd.DataFrame(primary)
                genus_rank=part.loc[
                    part.island_id.isin(set(base.island_id) if len(base) else set())
                ].groupby("genus",sort=True).n_target.sum().sort_values(
                    ascending=False,kind="mergesort"
                )
                top=genus_rank.index.tolist()
                baseline_count=int(sum(genus_rank.values))
                base_frac=float(genus_rank.head(5).sum()/baseline_count) if baseline_count else None
                out={
                    "region":region,"scope":scope,"tier":tier,
                    "n_original_region_floras":int(roster.regime.eq(region).sum()),
                    "n_original_islands_source_numeric_paired_support":int(len(base)),
                    "n_distinct_original_target_genera_in_eligible_islands":int(len(top)),
                    "top_five_original_genera_fraction_of_target_direct_occurrences":base_frac,
                    "n_original_islands_Direct_balance_both_10pp":int(
                        base.documentation_10pp_both.sum()) if len(base) else 0,
                    "n_original_islands_source_two_plus_genera_both":int(
                        base.source_2plus_genus_both.sum()) if len(base) else 0,
                    "n_original_islands_effective_three_plus_genera":int(
                        base.effective_genus_number.ge(3).sum()) if len(base) else 0,
                    "n_original_islands_largest_genus_at_most_half":int(
                        base.largest_genus_carrier_fraction.le(.5).sum()) if len(base) else 0,
                    "largest_genus_fraction_median":float(
                        base.largest_genus_carrier_fraction.median()) if len(base) else None,
                    "effective_genus_number_median":float(
                        base.effective_genus_number.median()) if len(base) else None,
                    "original_taxon_genus_ranks_never_exported":True,
                }
                sets={"all_original_eligible":base}
                if len(base):
                    sets.update({
                        "documentation_balance10":base.loc[base.documentation_10pp_both],
                        "source_2plus_genera":base.loc[base.source_2plus_genus_both],
                        "effective_genus_at_least3":base.loc[base.effective_genus_number.ge(3)],
                        "largest_genus_at_most_half":base.loc[base.largest_genus_carrier_fraction.le(.5)],
                    })
                else:
                    for k in (
                        "documentation_balance10","source_2plus_genera",
                        "effective_genus_at_least3","largest_genus_at_most_half"
                    ):
                        sets[k]=base
                # Top global genera are selected from ORIGINAL carrier counts
                # only, without inspecting response traits, predictions or scores.
                baseids=set(base.island_id) if len(base) else set()
                for k in (1,5):
                    forbidden=set(top[:k])
                    survivors=[]
                    for island,q in grouped.items():
                        if island not in baseids:
                            continue
                        stats=score_genus_sufficient_stats(
                            q.loc[~q.genus.isin(forbidden)]
                        )
                        if stats["status"]!="ELIGIBLE_PAIRED_SOURCE_GENUS_CONCENTRATION":
                            continue
                        stats.update(
                            island_id=island,spatial_block=q.spatial_block.iloc[0],
                            first_source=q.first_source.iloc[0],
                        )
                        survivors.append(stats)
                    sets[f"remove_top_{k}_original_carrier_genera"]=pd.DataFrame(survivors)
                # Same-island leave largest genus: rank only number of focal
                # documented plants, independent of phenotype or gain.
                survivors=[]
                baseids=set(base.island_id) if len(base) else set()
                for island,q in grouped.items():
                    if island not in baseids:continue
                    winner=q.sort_values(
                        ["n_target","genus"],ascending=[False,True],
                        kind="mergesort"
                    ).genus.iloc[0]
                    sub=q.loc[q.genus.ne(winner)]
                    stats=score_genus_sufficient_stats(sub)
                    if stats["status"]!="ELIGIBLE_PAIRED_SOURCE_GENUS_CONCENTRATION":
                        continue
                    stats.update(
                        island_id=island,spatial_block=q.spatial_block.iloc[0],
                        first_source=q.first_source.iloc[0],
                    )
                    survivors.append(stats)
                sets["remove_largest_within_each_island"]=pd.DataFrame(survivors)
                for label,frame in sets.items():
                    for grouping in ("spatial_block","first_source"):
                        if frame.empty:
                            out[f"{label}__{grouping}"]={
                                "n_original_floras":0,
                                "n_geography_or_source_groups":0,
                                "status":"NOT_IDENTIFIABLE_GROUP_SUPPORT",
                            }
                        else:
                            out[f"{label}__{grouping}"]=_cluster_summary(
                                frame,group=grouping,
                                seed=f"{region}|{scope}|{tier}|{label}|{grouping}",
                                bootstrap=bootstrap,
                            )
                result.append(out)
    return {
        "status":"ORIGINAL_GIFT_4505_EXTERNAL_NUMERIC_TUBE_GENUS_CONCENTRATION_AUDIT",
        "original_global_botanical_floras":population,
        "n_original_region_scope_corpus_panels":len(result),
        "panels":result,
        "study_limits":{
            "original_target_species_genera_private_not_exported":True,
            "genus_ranked_by_target_documented_species_count_not_phenotype":True,
            "source_pool_recalculated_after_genus_removal":True,
            "genus_equal_not_phylogenetically_independent":True,
            "source_NHM_numeric_and_original_Direct_may_share_publications":True,
            "no_unobserved_flower_morphology_imputed":True,
            "source_origin_gift_checklist_not_actual_seed_arrival":True,
            "original_4505_floras_retained":True,
            "exploratory_posthoc_dominant_genus_sensitivity":True,
            "no_pollinator_identity_or_adaptation_inferred":True,
        },
    }
