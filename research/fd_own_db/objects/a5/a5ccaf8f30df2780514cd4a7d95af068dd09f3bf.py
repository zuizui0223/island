"""Geologically classified island endemic vs native-nonendemic source-congener trait discordance.

Unlike Direct flower-triple complete cases, use ONE original high/medium
species_direct single-state trait at a time. The outcome is difference
between documented named plant species and OTHER documented source-region
congeneric species, NOT within-species change, confirmed ancestral state,
colonization failure, reproductive fitness, or pollinator-mediated adaptation.

Preserve original native vs endemic source evidence; source lists are
candidate GIFT mainland botanical checklists, not known evolutionary donors.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import numpy as np
import pandas as pd

from .geological_endemic_trait_opportunity import TRAITS, ORIGINS
from .regional_native_provenance_volume import cluster_equal_bootstrap

PRIMARY_TRAITS=("flower_primary_color","floral_symmetry")
MIN_PAIRS_PER_ISLAND=3
MIN_ISLANDS_PER_ORIGIN=10
MIN_BLOCKS_PER_ORIGIN=5
MIN_SOURCES_PER_ORIGIN=5


def _seed(*keys:object)->int:
    return int.from_bytes(sha256(":".join(map(str,keys)).encode()).digest()[:8],"big")


def _genus(name:str)->str:
    parts=name.strip().split()
    return parts[0] if len(parts)>=2 else ""


def strict_direct_trait_states(ledger:pd.DataFrame)->dict[str,dict[str,str]]:
    fields={"accepted_species","trait_name","normalized_value",
            "evidence_scope","evidence_quality"}
    if not fields.issubset(ledger):
        raise ValueError("Original direct evidence fields missing")
    z=ledger.loc[
        ledger.trait_name.isin(TRAITS)&
        ledger.evidence_scope.eq("species_direct")&
        ledger.evidence_quality.isin(("high","medium"))&
        ledger.normalized_value.notna()&
        ~ledger.normalized_value.astype(str).str.contains("|",regex=False)
    ].copy()
    if z.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("Ambiguous source Direct plant×trait state")
    return {t:z.loc[z.trait_name.eq(t)].set_index(
        "accepted_species").normalized_value.to_dict()
        for t in TRAITS}


def plant_source_discordance(
    names:set[str], source_states:dict[str,str], focal_states:dict[str,str]
)->list[dict]:
    """One datum per actually measured source-listed ACCEPTED focal plant.

    Exact discordance against a uniformly sampled distinct GIFT source
    congener (expected value). This avoids favorable post-hoc matching on
    the flower phenotype, but donor reference depth still influences it.
    """
    donor_genus=defaultdict(list)
    for key,state in source_states.items():
        g=_genus(key)
        if g:
            donor_genus[g].append(state)
    out=[]
    for taxon in sorted(names):
        if taxon not in focal_states:
            continue
        genus=_genus(taxon)
        if not genus or genus not in donor_genus:
            continue
        x=focal_states[taxon]
        v=donor_genus[genus]
        out.append(dict(
            accepted_species=taxon,
            source_genus_congener_count=len(v),
            discordance=float(np.mean([x!=s for s in v])),
            _private_original_state=x,
        ))
    return out


def _status_contrast_within_identical_genus(
        endemic:list[dict],native_nonendemic:list[dict]
)->tuple[int,int,int,float]:
    """Pre-donor trait-information audit for common source-genus contrasts.

    The original species-direct phenotype exists once per accepted plant
    species. Under same-genus matching, a zero contrast can be a deterministic
    consequence of zero within-genus floral-category variation; count that
    inability to falsify as support loss, not biological equivalence.

    Each shared genus is equally weighted, regardless of number of island
    species; original source conspecifics are excluded upstream.
    """
    e=defaultdict(list)
    n=defaultdict(list)
    for row in endemic:
        e[_genus(row["accepted_species"])].append(row)
    for row in native_nonendemic:
        n[_genus(row["accepted_species"])].append(row)
    common=sorted((set(e)&set(n))-{""})
    n_status_variable=0
    n_nonzero_source_discordance=0
    paired=[]
    for genus in common:
        e_states={p["_private_original_state"] for p in e[genus]}
        n_states={p["_private_original_state"] for p in n[genus]}
        # A genus with only one total state is STRICTLY uninformative:
        # source-relative discordance must be the same in both status groups.
        if len(e_states|n_states)>1:
            n_status_variable+=1
        local_delta=(np.mean([p["discordance"] for p in e[genus]])-
                     np.mean([p["discordance"] for p in n[genus]]))
        if abs(local_delta)>1e-12:
            n_nonzero_source_discordance+=1
        paired.append(float(local_delta))
    diff=float(np.mean(paired)) if len(common)>=2 else float("nan")
    return len(common),n_status_variable,n_nonzero_source_discordance,diff


def geological_endemic_state_discordance(
    direct:pd.DataFrame,raw:pd.DataFrame,geology:pd.DataFrame,
    candidate_assignments:pd.DataFrame,mainland:pd.DataFrame,*,
    seed:int=20261114,boots:int=1999,
)->dict:
    """Private island-level outcomes, public REGION×GEOLOGY×TRAIT aggregates."""
    if boots<999:
        raise ValueError("Conditional original source bootstrap support weakened")
    expected={"island_id","accepted_species","origin_status",
              "endemic_status","floristic_status"}
    if not expected.issubset(raw) or raw.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original island native/endemic source ledger invalid")
    if not raw.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original native-origin category invalid")
    if not raw.endemic_status.isin(("endemic","nonendemic","unresolved")).all():
        raise ValueError("Original independent endemism category invalid")
    if ((raw.floristic_status.eq("endemic") &
        (~raw.origin_status.eq("native") | ~raw.endemic_status.eq("endemic"))).any() or
        (raw.floristic_status.eq("native_nonendemic") &
        (~raw.origin_status.eq("native") | ~raw.endemic_status.eq("nonendemic"))).any()):
        raise ValueError("Endemic source status contradicts original GIFT native origin")
    if not {"island_id","geological_origin","analysis_regime","spatial_block"}.issubset(geology):
        raise ValueError("Independent Roeble/Zizka geological data missing")
    if geology.island_id.duplicated().any():
        raise ValueError("Original geology island duplicates")
    if not {"island_id","source_mode","source_rank","entity_ID"}.issubset(candidate_assignments):
        raise ValueError("Original nearest GIFT mainland candidate assignments missing")
    if not {"entity_ID","work_species"}.issubset(mainland):
        raise ValueError("Original source donor plant roster missing")
    states=strict_direct_trait_states(direct)
    typed=geology.loc[geology.geological_origin.isin(ORIGINS),[
        "island_id","geological_origin","analysis_regime","spatial_block"]]
    if len(typed)<4:
        raise ValueError("No independent geology types")
    observed=raw.loc[raw.island_id.isin(typed.island_id)].merge(
        typed,on="island_id",validate="many_to_one",how="inner")
    assignments=candidate_assignments.loc[
        candidate_assignments.source_mode.eq("geo_k5")].copy()
    assignments["source_rank"]=pd.to_numeric(assignments.source_rank,errors="raise")
    assignments=assignments.loc[assignments.source_rank.between(1,5)]
    if assignments.duplicated(["island_id","source_rank"]).any():
        raise ValueError("GIFT candidate source rank duplicates")
    near=assignments.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    donor_by_source=mainland.drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()

    island_rows=[]
    for island,frame in observed.groupby("island_id",sort=True):
        source_ids=near.get(str(island),())
        if len(source_ids)!=5 or len(set(source_ids))!=5:
            continue # source unavailable, NOT ecological absence
        every_island_taxon=set(frame.accepted_species)
        donor_names=set().union(*(donor_by_source.get(k,set()) for k in source_ids))-every_island_taxon
        status={
            "endemic":set(frame.loc[frame.floristic_status.eq("endemic"),"accepted_species"]),
            "native_nonendemic":set(frame.loc[
                frame.floristic_status.eq("native_nonendemic"),"accepted_species"]),
        }
        for trait in TRAITS:
            ds=states[trait]
            donors={name:ds[name] for name in donor_names if name in ds}
            bystatus={}
            for role in ("endemic","native_nonendemic"):
                data=plant_source_discordance(status[role],donors,ds)
                bystatus[role]=data
            a,b=bystatus["endemic"],bystatus["native_nonendemic"]
            row=dict(
                island_id=str(island),regime=frame.analysis_regime.iloc[0],
                geological_origin=frame.geological_origin.iloc[0],
                spatial_block=frame.spatial_block.iloc[0],
                source_signature=source_ids[0],
                trait=trait,
                n_endemic_direct_source_congener_rows=len(a),
                n_native_nonendemic_direct_source_congener_rows=len(b),
                endemic_mean_discordance=float(np.mean([q["discordance"] for q in a]))
                    if a else np.nan,
                native_nonendemic_mean_discordance=float(
                    np.mean([q["discordance"] for q in b])) if b else np.nan,
                endemic_median_source_congener_count=float(
                    np.median([q["source_genus_congener_count"] for q in a]))
                    if a else np.nan,
                native_nonendemic_median_source_congener_count=float(
                    np.median([q["source_genus_congener_count"] for q in b]))
                    if b else np.nan,
            )
            (common_genus_count,informative_genus_count,
             nonzero_genus_count,common_genus_diff)=_status_contrast_within_identical_genus(a,b)
            row["n_same_island_same_source_genus_endemic_native_groups"]=common_genus_count
            row["n_same_island_genus_with_status_group_floral_state_variation"]=informative_genus_count
            row["n_same_island_genus_with_nonzero_endemic_native_source_discordance"]=nonzero_genus_count
            row["same_island_same_genus_endemic_minus_native_discordance"]=common_genus_diff
            row["same_island_endemic_minus_native_discordance"]=(
                row["endemic_mean_discordance"]-
                row["native_nonendemic_mean_discordance"]
                if min(len(a),len(b))>=MIN_PAIRS_PER_ISLAND else np.nan)
            # The exact same focal accepted species may appear on several
            # islands, so each physical-island average is the analysis row.
            island_rows.append(row)

    all_data=pd.DataFrame(island_rows)
    if all_data.empty:
        raise ValueError("No original island/trait/source congener coverage")
    summaries=[]
    for (region,origin,trait),g in all_data.groupby(
        ["regime","geological_origin","trait"],sort=True):
        paired=g.loc[g.same_island_endemic_minus_native_discordance.notna()]
        matched=g.loc[g.same_island_same_genus_endemic_minus_native_discordance.notna()]
        support=dict(
            regime=region,geological_origin=origin,trait=trait,
            n_original_islands_with_5_source_regions=len(g),
            n_islands_with_any_endemic_source_congener=int(
                g.n_endemic_direct_source_congener_rows.gt(0).sum()),
            n_islands_with_any_native_nonendemic_source_congener=int(
                g.n_native_nonendemic_direct_source_congener_rows.gt(0).sum()),
            n_endemic_island_taxon_source_congener_rows=int(
                g.n_endemic_direct_source_congener_rows.sum()),
            n_native_nonendemic_island_taxon_source_congener_rows=int(
                g.n_native_nonendemic_direct_source_congener_rows.sum()),
            n_islands_with_both_ge3_congener_comparisons=int(len(paired)),
            n_islands_with_at_least_two_same_genus_endemic_nonendemic_pairs=int(len(matched)),
            n_blocks_with_same_genus_pairs=int(matched.spatial_block.nunique()),
            n_sources_with_same_genus_pairs=int(matched.source_signature.nunique()),
            n_common_endemic_nonendemic_genus_units_in_this_group=int(
                g.n_same_island_same_source_genus_endemic_native_groups.sum()),
            n_status_variable_endemic_nonendemic_genus_island_units=int(
                g.n_same_island_genus_with_status_group_floral_state_variation.sum()),
            n_genus_island_units_with_nonzero_source_discordance=int(
                g.n_same_island_genus_with_nonzero_endemic_native_source_discordance.sum()),
            n_islands_with_any_within_genus_status_flower_variation=int(
                g.n_same_island_genus_with_status_group_floral_state_variation.gt(0).sum()),
            n_independent_blocks_with_both=int(paired.spatial_block.nunique()),
            n_source_signatures_with_both=int(paired.source_signature.nunique()),
            n_recorded_genus_opportunity_less_than_3_endemic_sites=int(
                g.n_endemic_direct_source_congener_rows.lt(3).sum()),
        )
        okay=(len(paired)>=MIN_ISLANDS_PER_ORIGIN and
              paired.spatial_block.nunique()>=MIN_BLOCKS_PER_ORIGIN and
              paired.source_signature.nunique()>=MIN_SOURCES_PER_ORIGIN)
        support["paired_status_comparison_reportable"]=bool(okay)
        same_genus_ok=(len(matched)>=MIN_ISLANDS_PER_ORIGIN and
                       matched.spatial_block.nunique()>=MIN_BLOCKS_PER_ORIGIN and
                       matched.source_signature.nunique()>=MIN_SOURCES_PER_ORIGIN)
        support["same_genus_endemic_native_comparison_reportable"]=bool(same_genus_ok)
        if same_genus_ok:
            for label,group in (("block","spatial_block"),("source","source_signature")):
                support[label+"_matched_genus_endemic_minus_native_discordance"]=cluster_equal_bootstrap(
                    matched,"same_island_same_genus_endemic_minus_native_discordance",group,
                    draws=boots,seed=_seed(seed,"within_genus",region,origin,trait,label))
        if okay:
            for label,group in (("block","spatial_block"),("source","source_signature")):
                support[label+"_endemic_minus_native_discordance"]=cluster_equal_bootstrap(
                    paired,"same_island_endemic_minus_native_discordance",group,
                    draws=boots,seed=_seed(seed,"paired",region,origin,trait,label))
        summaries.append(support)
    results=[]
    for (region,trait),g in all_data.groupby(["regime","trait"],sort=True):
        left=g.loc[g.geological_origin.eq("continental")&
                   g.same_island_endemic_minus_native_discordance.notna()]
        right=g.loc[g.geological_origin.eq("oceanic")&
                    g.same_island_endemic_minus_native_discordance.notna()]
        okay=all(len(q)>=MIN_ISLANDS_PER_ORIGIN and
                 q.spatial_block.nunique()>=MIN_BLOCKS_PER_ORIGIN and
                 q.source_signature.nunique()>=MIN_SOURCES_PER_ORIGIN
                 for q in (left,right))
        within_left=g.loc[g.geological_origin.eq("continental")&
                    g.same_island_same_genus_endemic_minus_native_discordance.notna()]
        within_right=g.loc[g.geological_origin.eq("oceanic")&
                    g.same_island_same_genus_endemic_minus_native_discordance.notna()]
        row=dict(regime=region,trait=trait,primary_trait=bool(trait in PRIMARY_TRAITS),
            n_oceanic_same_genus_paired_islands=int(len(within_right)),
            n_continental_same_genus_paired_islands=int(len(within_left)),
            n_oceanic_paired_islands=int(len(right)),
            n_continental_paired_islands=int(len(left)),
            n_oceanic_blocks=int(right.spatial_block.nunique()),
            n_continental_blocks=int(left.spatial_block.nunique()),
            origin_comparison_identifiable_in_supported_cohort=bool(okay))
        if okay:
            for label,cluster in (("block","spatial_block"),("source","source_signature")):
                a=left.groupby(cluster).same_island_endemic_minus_native_discordance.mean().to_numpy(float)
                b=right.groupby(cluster).same_island_endemic_minus_native_discordance.mean().to_numpy(float)
                rng=np.random.default_rng(_seed(seed,"origin",region,trait,label))
                resampled_b=b[rng.integers(len(b),size=(boots,len(b)))].mean(axis=1)
                resampled_a=a[rng.integers(len(a),size=(boots,len(a)))].mean(axis=1)
                row[label+"_oceanic_minus_continental_same_island_discordance"]=dict(
                    mean=float(b.mean()-a.mean()),
                    ci95=[float(v) for v in np.quantile(resampled_b-resampled_a,[.025,.975])])
        results.append(row)
    return dict(
        status="ORIGINAL_GEOLOGICAL_ENDEMIC_SINGLE_TRAIT_CONGENER_DISCORDANCE",
        n_geology_typed_islands=int(len(typed)),
        n_original_islands_with_5_source_region_trait_eligibility=int(all_data.island_id.nunique()),
        primary_single_trait_hypotheses=list(PRIMARY_TRAITS),
        source_provenance=dict(
            original_endemic_status_separate_from_native=True,
            species_direct_high_medium_unambiguous_only=True,
            source_genus_congeners_selected_without_trait_outcome_matching=True,
            plant_state_once_per_accepted_species_not_island_population=True,
            donor_state_frequency_discordance_is_not_adaptation=True,
            same_genus_zero_contrast_power_audited_with_original_focal_state_variation=True,
            original_source_GIFT_nearest_five_not_ancestral=True,
            regional_clustering_is_not_phylogenetic_replication=True,
            geological_origin_is_not_LGM_connectivity=True,
            conditional_intervals_not_adjusted_for_multiple_exploratory_traits=True,
        ),
        endemic_vs_native_same_island_congener_support=summaries,
        geological_origin_difference_of_endemic_status_effect=results,
        no_raw_species_island_genus_or_source_id_exported=True,
    )
