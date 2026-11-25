"""Same-island cohort inference audit for 4,505 GIFT floras × qualified GBIF guilds.

Rationale: the previous raw Bombus–floral colour×form result used 760 vs23
islands, while the source-relative comparison used only 50 vs16 DIFFERENT
islands. This explicitly computes raw and source-relative category contrasts
on the SAME selected islands, then conditions on physical 10-degree blocks
and independent first GIFT candidate-region signatures where feasible.

No ecological treatment is randomized: sign swaps and intervals are
descriptive sensitivity diagnostics, never a causal or historical-shift test.
All 4,505 original floras remain in initial coverage denominators; source-
matched cases are a transparent additional *analysis opportunity* stratum.
"""
from __future__ import annotations

from hashlib import sha256

import numpy as np
import pandas as pd

from .global_partial_pollinator_guild_bridge import (
    normalize_guild_evidence,CHANNELS
)
from .global_partial_trait_flora_fd import TRAITS,PAIR_TRAITS

ENDPOINTS=tuple(TRAITS)+tuple("__".join(z) for z in PAIR_TRAITS)
GUILDS=("bombus","non_bombus_bees","flower_visiting_birds")
PAIR_METRICS=(
    ("observed_finite5","source_relative_finite5"),
    ("observed_genus_equal_finite5","source_relative_genus_equal_finite5"),
    ("observed_two_trait_pairing","source_relative_pairing"),
)
MIN_EACH=10
MIN_BLOCKS=5
MIN_SOURCES=5
# Fixed source-qualified eight-slot family is exploratory post-hoc, not confirmatory.
PRIMARY_SPECS=(
  ("bombus","flower_primary_color__floral_form","observed_finite5","source_relative_finite5"),
  ("bombus","flower_primary_color__floral_form","observed_two_trait_pairing","source_relative_pairing"),
  ("bombus","floral_form","observed_finite5","source_relative_finite5"),
  ("bombus","floral_form","observed_genus_equal_finite5","source_relative_genus_equal_finite5"),
  ("flower_visiting_birds","flower_primary_color__floral_form","observed_finite5","source_relative_finite5"),
  ("flower_visiting_birds","flower_primary_color__floral_form","observed_two_trait_pairing","source_relative_pairing"),
  ("flower_visiting_birds","floral_symmetry","observed_finite5","source_relative_finite5"),
  ("flower_visiting_birds","floral_symmetry","observed_genus_equal_finite5","source_relative_genus_equal_finite5"),
)


def _seed(*parts:object)->int:
    return int.from_bytes(sha256("|".join(map(str,parts)).encode()).digest()[:8],"big")


def group_status_contrast(frame:pd.DataFrame,metric:str,cluster:str)->dict|None:
    src=frame.loc[frame.channel_state.eq("retained"),[cluster,metric]].groupby(
        cluster)[metric].mean().dropna()
    no=frame.loc[frame.channel_state.eq("disrupted"),[cluster,metric]].groupby(
        cluster)[metric].mean().dropna()
    if min(len(src),len(no))<5:
        return None
    return dict(mean=float(src.mean()-no.mean()),
                n_detected_clusters=int(len(src)),
                n_qualified_nondetected_clusters=int(len(no)))


def matched_block_delta(frame:pd.DataFrame,metric:str,cluster:str,*,
                        bootstrap:int,seed:int)->dict:
    """Within-cluster observed-minus-adequate-nondetection differences.

    This controls the fixed region/block (or GIFT first-source) aggregate
    in the same island cohort. No causal assignment is implied.
    """
    temp=frame[[cluster,"channel_state",metric]].dropna(subset=[metric,cluster])
    sums=temp.groupby([cluster,"channel_state"])[metric].mean().unstack()
    if not {"retained","disrupted"}.issubset(sums):
        return dict(n_mixed_clusters=0,status="NO_WITHIN_CLUSTER_MIXING")
    d=(sums["retained"]-sums["disrupted"]).dropna().to_numpy(float)
    n=len(d)
    if n<MIN_BLOCKS:
        return dict(n_mixed_clusters=n,status="MIXED_CLUSTER_POWER_TOO_LOW")
    rng=np.random.default_rng(seed)
    boot=d[rng.integers(n,size=(bootstrap,n))].mean(axis=1)
    effect=float(d.mean())
    # Block signs are NOT randomized geology or pollinator treatment.
    # Enumeration is a conservative, transparent *diagnostic* only.
    if n<=16:
        patterns=np.arange(1<<n,dtype=np.uint32)[:,None]
        bit=np.arange(n,dtype=np.uint32)[None,:]
        sign=np.where(((patterns>>bit)&1)==1,1.,-1.)
        null=np.abs(sign@d/n)
        p=float((null>=abs(effect)-1e-12).mean())
        mechanism="EXACT_FULL_TWO_SIDED_CLUSTER_SIGN_DIAGNOSTIC"
        simulations=int(1<<n)
    else:
        total=9999
        sign=rng.choice(np.array([-1.,1.]),size=(total,n))
        null=np.abs(sign@d/n)
        p=float((1+np.count_nonzero(null>=abs(effect)-1e-12))/(total+1))
        mechanism="NUMERICAL_TWO_SIDED_CLUSTER_SIGN_DIAGNOSTIC"
        simulations=total
    leave=(np.sum(d)-d)/(n-1)
    return dict(status="DESCRIPTIVE_NOT_CAUSAL",
                n_mixed_clusters=n,mean=effect,
                bootstrap95=[float(x) for x in np.quantile(boot,[.025,.975])],
                two_sided_cluster_sign_diagnostic_p=p,
                sign_method=mechanism,n_sign_permutations=simulations,
                leave_one_mixed_cluster_range=[
                  float(leave.min()),float(leave.max())])


def _support(a:pd.DataFrame)->dict:
    rr=a.loc[a.channel_state.eq("retained")]
    dd=a.loc[a.channel_state.eq("disrupted")]
    b=[int(z.spatial_block.nunique()) for z in (rr,dd)]
    s=[int(z.source_signature.nunique()) for z in (rr,dd)]
    eligible=(min(len(rr),len(dd))>=MIN_EACH and
              min(b)>=MIN_BLOCKS and min(s)>=MIN_SOURCES)
    return dict(n_GBIF_detected=int(len(rr)),
                n_qualified_nondetection=int(len(dd)),
                n_block_by_state=b,n_first_GIFT_source_by_state=s,
                sufficient_independent_support=bool(eligible))


def _is_real_value(a:pd.Series)->pd.Series:
    # All original effect-size calculations produce float.NaN for
    # inapplicable metrics; never turn these into true biological zero.
    return pd.to_numeric(a,errors="coerce").notna()


def _bh(p:list[float])->list[float]:
    n=len(p)
    order=np.argsort(p)
    out=np.zeros(n)
    prev=1.
    for rank in range(n-1,-1,-1):
        k=int(order[rank])
        prev=min(prev,p[k]*n/(rank+1))
        out[k]=prev
    return list(map(float,out))


def source_cohort_and_guild_robustness(
    private_scores:pd.DataFrame, original_channels:pd.DataFrame,*,
    seed:int=20261114,draws:int=999
)->dict:
    if draws<999:
        raise ValueError("Original gross and source-conditional diagnostic needs >=999 samples")
    required={"island_id","regime","geological_origin","original_scope",
              "outcome","spatial_block","source_signature",
              "n_recorded_direct_species","observed_direct_coverage",
              *[x for a in PAIR_METRICS for x in a]}
    if not required.issubset(private_scores):
        raise ValueError("Missing original source FD paired metrics")
    raw=private_scores
    n_global=int(raw.island_id.nunique())
    if n_global<4000 or len(raw)!=n_global*2*len(ENDPOINTS):
        raise ValueError("Global original island flora population was narrowed")
    if raw.duplicated(["island_id","original_scope","outcome"]).any():
        raise ValueError("Original island species trait duplicated")
    states=normalize_guild_evidence(original_channels,set(raw.island_id))
    states=states.loc[states.channel_id.isin(GUILDS)&
                      states.channel_state.isin(("retained","disrupted"))]
    out=[]
    # Do not generate a gigantic 585k copy, since this audit explicitly
    # requires only a handful of biologically interpretable fixed metrics.
    for region in ["GLOBAL_ALL","tropical","northern_midlatitude",
                   "northern_high_latitude","southern_extratropical"]:
        f=raw if region=="GLOBAL_ALL" else raw.loc[raw.regime.eq(region)]
        if f.empty:continue
        for scope in ("all_recorded","original_native"):
            u=f.loc[f.original_scope.eq(scope)]
            for guild in GUILDS:
                channels=states.loc[states.channel_id.eq(guild),
                                    ["island_id","channel_state","source_state"]]
                if channels.empty:continue
                base=u.merge(channels,on="island_id",validate="many_to_one",
                             how="inner")
                # Both statuses must have independent candidate-source
                # opportunity, else detecting only one may be biased.
                base=base.loc[base.source_state.eq("available")]
                for trait in ENDPOINTS:
                    group=base.loc[base.outcome.eq(trait)]
                    for raw_metric,relative_metric in PAIR_METRICS:
                        full=group.loc[_is_real_value(group[raw_metric])]
                        # The same island membership, channel state, source
                        # ranking and valid trait subtype for both quantities.
                        same=full.loc[_is_real_value(full[relative_metric])]
                        count_full=_support(full)
                        count_same=_support(same)
                        row=dict(
                            region=region,original_source_flora_status=scope,
                            GBIF_guild=guild,trait_or_pair=trait,
                            original_raw_metric=raw_metric,
                            original_GIFT_source_relative_metric=relative_metric,
                            all_4505_original_floras_remain_main_denominator=n_global,
                            full_raw=count_full,exact_same_island_source=count_same,
                            n_original_all_trait_scope_islands=int(
                                len(u.loc[u.outcome.eq(trait)])),
                        )
                        if count_full["sufficient_independent_support"]:
                            row["original_raw_full_group_block_difference"]=group_status_contrast(
                                full,raw_metric,"spatial_block")
                            # IMPORTANT: previous global contrasts average
                            # detected/adequately-nondetected BLOCKS SEPARATELY.
                            # They do NOT ensure any block is observed in
                            # BOTH states. Test mixed block overlap BEFORE
                            # attributing the difference to a GBIF guild.
                            row["physical_block_within_full_raw"]=matched_block_delta(
                                full,raw_metric,"spatial_block",bootstrap=draws,
                                seed=_seed(seed,region,scope,guild,trait,
                                           raw_metric,"full_block"))
                            row["GIFT_first_source_within_full_raw"]=matched_block_delta(
                                full,raw_metric,"source_signature",bootstrap=draws,
                                seed=_seed(seed,region,scope,guild,trait,
                                           raw_metric,"full_source"))
                        if count_same["sufficient_independent_support"]:
                            for label,col in (
                                ("physical_block","spatial_block"),
                                ("GIFT_first_nearest_region","source_signature")):
                                row[label+"_same_raw"]=group_status_contrast(
                                    same,raw_metric,col)
                                row[label+"_same_source_relative"]=group_status_contrast(
                                    same,relative_metric,col)
                                row[label+"_within_raw"]=matched_block_delta(
                                    same,raw_metric,col,bootstrap=draws,seed=_seed(
                                      seed,region,scope,guild,trait,raw_metric,"raw",col))
                                row[label+"_within_source_relative"]=matched_block_delta(
                                    same,relative_metric,col,bootstrap=draws,seed=_seed(
                                      seed,region,scope,guild,trait,raw_metric,"source",col))
                            # Compare recorded trait rates on the EXACT SAME
                            # islands; difference is observation opportunity.
                            cov=same.groupby("channel_state")[
                                "observed_direct_coverage"].mean().to_dict()
                            row["same_island_direct_plant_trait_record_fraction_by_guild_status"]=dict(
                                GBIF_detected=float(cov.get("retained",float("nan"))),
                                qualified_nondetection=float(cov.get("disrupted",float("nan"))))
                        out.append(row)
    # Transparent post-hoc 8-slot contrast family, not confirmatory a priori:
    # source-relative diagnostic requires >=5 genuinely MIXED physical blocks.
    panel=[]
    for guild,trait,raw_metric,source_metric in PRIMARY_SPECS:
        row=next((x for x in out if x["region"]=="GLOBAL_ALL"
             and x["original_source_flora_status"]=="all_recorded"
             and x["GBIF_guild"]==guild and x["trait_or_pair"]==trait
             and x["original_raw_metric"]==raw_metric and
             x["original_GIFT_source_relative_metric"]==source_metric),None)
        if row is None:raise ValueError("Original predeclared posthoc diagnostic missing")
        k=row.get("physical_block_within_source_relative",{})
        stat=(k.get("two_sided_cluster_sign_diagnostic_p")
              if k.get("status")=="DESCRIPTIVE_NOT_CAUSAL" else None)
        panel.append(dict(
            guild=guild,trait_or_pair=trait,metric_pair=[raw_metric,source_metric],
            n_source_same_detected=row["exact_same_island_source"]["n_GBIF_detected"],
            n_source_same_qualified_nondetected=row["exact_same_island_source"][
                "n_qualified_nondetection"],
            n_physical_blocks_with_both_statuses=int(k.get("n_mixed_clusters",0)),
            source_relative_within_block_sign_diagnostic_p=stat,
            testable_within_physical_blocks=stat is not None,
        ))
    ps=[x["source_relative_within_block_sign_diagnostic_p"]
        if x["testable_within_physical_blocks"] else 1.
        for x in panel]
    for t,q in zip(panel,_bh(ps)):
        # p=1 for untestable registered slot: no silent dropping
        t["BH_q_with_untestable_slots_conservatively_p1"]=float(q) if (
            t["testable_within_physical_blocks"]) else None
    return dict(
        status="GLOBAL_4505_ORIGINAL_GUILD_RAW_VS_SOURCE_SAME_COHORT_CONFOUNDING_AUDIT",
        n_original_individual_botanical_islands=n_global,
        n_source_cohort_contrast_panels=len(out),
        n_source_same_island_panels_with_support=sum(
            p["exact_same_island_source"]["sufficient_independent_support"] for p in out),
        original_all_cohort_and_same_GIFT_source_analyses=out,
        source_fixed_posthoc_primary_eight_slot_family=panel,
        n_original_primary_eight_slot_testable_with_mixed_blocks=sum(
            p["testable_within_physical_blocks"] for p in panel),
        source_limits=dict(
            exploratory_posthoc_eight_slots_not_hypotheses_pre_registered_before_previous_results=True,
            exact_same_islands_raw_and_source_reference_compared=True,
            block_cluster_must_contain_BOTH_GBIF_statuses=True,
            previous_separately_clustered_comparisons_do_not_prove_any_within_block_overlap=True,
            full_original_raw_cohort_and_source_comparable_subcohort_both_receive_mixed_block_diagnostic=True,
            GIFT_first_nearest_region_must_contain_BOTH_GBIF_statuses=True,
            original_source_GBIF_qualified_nondetection_not_historical_loss=True,
            original_native_flora_status_not_imputed=True,
            source_species_flower_phenotype_not_island_local=True,
            source_guild_GBIF_taxon_presence_not_plant_visitor=True,
            uncertainty_not_adjusted_for_detection_effort_or_phylogenetic_relatedness=True,
            sign_flip_not_a_randomized_or_causal_p_value=True,
            untestable_primary_slots_conservatively_counted_as_p_equals_1_for_BH=True,
            original_global_4505_flora_denominator_not_cut_to_50_islands=True,
            restricted_island_taxon_GBIF_source_group_IDs_stay_private=True,
        )
    )
