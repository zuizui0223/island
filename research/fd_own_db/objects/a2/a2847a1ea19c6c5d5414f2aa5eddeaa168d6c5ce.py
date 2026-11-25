"""Exactly partition GLOBAL rare-profile *metric coverage* by genus turnover.

For each source-qualified island and its five-nearest GIFT candidate native
mainland flora, remove ALL accepted island plant species from donor source.
Hold 30 original Direct plant names per island/source via independent rarefaction.

The union-of-Hamming-balls categorical reference profile coverage has an
exact two-player SHAPLEY partition into:
  (A) genera represented on BOTH island and donor source ("shared")
  (B) genera found in only one recorded source flora ("unshared").
For reference bin R, v(S)=fraction of original observed R profiles within
one categorical mismatch of any species in group S. Partition:
phi_shared=0.5[v(shared)+v(shared+unshared)-v(unshared)]
phi_unshared=0.5[v(unshared)+v(shared+unshared)-v(shared)].
Thus contributions sum EXACTLY to total volume per island/source, including
profile overlap, without using an invalid additive Rao Q proxy.

An independent fixed-*five distinct matched genera* diagnostic draws one
different-source species per genus per island/source group, so genus labels
and genus counts are identical without requiring 15 or 30 matched species.
This is a selected, sparse, n=5 cohort and NOT within-species adaptation or
a universal source-genus-conditional null. A fixed-8-genus sensitivity
measures identifiability only if enough unique genera exist.

Reference profiles are original Direct high/medium 2,234 species/100 actual
categorical triples, rarity classes from original Direct species frequency.
No identity-level source, genus, island or taxon rows are returned.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .discrete_functional_volume import TRAITS, _resample_group
from .rare_profile_decomposition import RarityReference, genus, _seedstr

KINDS = ("rare", "common")
SCOPES=("all_recorded","source_native")
CONTROLS=("shared","unshared","total")
MATCHED_FIVE=5
MATCHED_EIGHT=8


def _source_genus(names: set[str]) -> set[str]:
    return {g for n in names if (g := genus(n))}


def _ball_mask(reference:RarityReference, names:list[str],
               triples:pd.DataFrame, *,
               original_species_mask:dict[str,int]|None=None)->int:
    mask=0
    if original_species_mask is not None:
        # Pre-compile source original trait state -> radius-1 integer bitmap
        # ONCE at the start of a real 928-island analysis; it is exact and
        # avoids high-cost pandas lookups for ~100k rarefied network draws.
        for species in names:
            if species not in original_species_mask:
                raise ValueError("Original Direct species missing from reference")
            mask |= original_species_mask[species]
    elif len(names):
        for p in set(map(tuple,triples.loc[names,list(TRAITS)].to_numpy(dtype=str))):
            if p not in reference.volume.masks:
                raise ValueError("Original direct phenotype not in global reference")
            mask |= reference.volume.masks[p][1]
    return mask


def _bin_value(ref:RarityReference,mask:int,bin_name:str)->float:
    return ((mask & ref.bins[bin_name]).bit_count()
            /ref.n_bins[bin_name])


def two_group_shapley(ref:RarityReference,shared_mask:int,
                      exclusive_mask:int)->dict[str,float]:
    """Exact total=shared+unshared for every original reference rarity bin."""
    everything=shared_mask|exclusive_mask
    values={}
    for rarity in KINDS:
        s=_bin_value(ref,shared_mask,rarity)
        u=_bin_value(ref,exclusive_mask,rarity)
        total=_bin_value(ref,everything,rarity)
        # Shapley over genus-shared and genus-exclusive species groups;
        # overlap is assigned 50/50 and cannot be counted twice.
        shared=.5*(s+total-u)
        unshared=.5*(u+total-s)
        if abs(shared+unshared-total)>1e-12:
            raise AssertionError("Categorical coverage Shapley attribution not exact")
        values[f"{rarity}_shared"]=float(shared)
        values[f"{rarity}_unshared"]=float(unshared)
        values[f"{rarity}_total"]=float(total)
    for kind in CONTROLS:
        values[f"gap_{kind}"]=values[f"rare_{kind}"]-values[f"common_{kind}"]
    if abs(values["gap_shared"]+values["gap_unshared"]-values["gap_total"])>1e-12:
        raise AssertionError("Rare-common exact Shapley decomposition failed")
    return values


def _genus_groups(names:set[str])->dict[str,tuple[str,...]]:
    groups=defaultdict(list)
    for name in sorted(names):
        g=genus(name)
        if g:
            groups[g].append(name)
    return {g:tuple(n) for g,n in groups.items()}


def _sample_components(
    reference:RarityReference,triples:pd.DataFrame,
    local_names:set[str],source_names:set[str],*,
    n:int=30,repeats:int=40,seed:int=20261113,
    original_species_mask:dict[str,int]|None=None
)->dict[str,float]:
    if n!=30 or repeats<20:
        raise ValueError("Frozen original macro sample/rarefaction is fixed")
    if len(local_names)<30 or len(source_names)<30:
        raise ValueError("Inadequate original Direct species on island or GIFT source")
    if local_names.intersection(source_names):
        raise ValueError("Source contains exact focal-island taxon-name leakage")
    shared_genus=_source_genus(local_names)&_source_genus(source_names)
    rng=np.random.default_rng(seed)
    rows=[]
    for _ in range(repeats):
        i=rng.choice(sorted(local_names),size=n,replace=False).tolist()
        m=rng.choice(sorted(source_names),size=n,replace=False).tolist()
        distinct_ig=set(map(genus,i))
        distinct_mg=set(map(genus,m))
        # Group belonging is frozen from FULL SOURCE FLORA before sample;
        # do NOT flip a genus to exclusive because n=30 misses it.
        input_groups=[]
        for pool in (i,m):
            shared=[t for t in pool if genus(t) in shared_genus]
            exclusive=[t for t in pool if genus(t) not in shared_genus]
            sh=_ball_mask(reference,shared,triples,
                            original_species_mask=original_species_mask)
            un=_ball_mask(reference,exclusive,triples,
                            original_species_mask=original_species_mask)
            input_groups.append(two_group_shapley(reference,sh,un))
        diff={k:input_groups[0][k]-input_groups[1][k]
              for k in input_groups[0]}
        diff.update(n_sampled_island_shared_genus=int(sum(genus(x) in shared_genus for x in i)),
                    n_sampled_source_shared_genus=int(sum(genus(x) in shared_genus for x in m)))
        rows.append(diff)
    scores=pd.DataFrame(rows).mean(numeric_only=True).to_dict()
    scores["n_shared_source_genera_in_original_flora"]=len(shared_genus)
    scores["n_original_local_genera"]=len(_source_genus(local_names))
    scores["n_original_source_genera"]=len(_source_genus(source_names))
    return {k:float(v) for k,v in scores.items()}


def _match_unique_genera(
    ref:RarityReference,triple:pd.DataFrame,
    local_names:set[str],donor_names:set[str],*,
    n:int=5,repeats:int=40,seed:int=20261113,
    original_species_mask:dict[str,int]|None=None
)->dict[str,float]:
    if n not in (MATCHED_FIVE,MATCHED_EIGHT) or repeats<20:
        raise ValueError("Frozen distinct genus source-match size/rarefaction changed")
    a=_genus_groups(local_names)
    b=_genus_groups(donor_names-local_names)
    shared=sorted(set(a)&set(b))
    if len(shared)<n:
        raise ValueError("Original source has too few DISTINCT shared genera")
    rng=np.random.default_rng(seed)
    values=[]
    for _ in range(repeats):
        selected=sorted(rng.choice(shared,size=n,replace=False).tolist())
        aa=[str(rng.choice(a[g])) for g in selected]
        bb=[str(rng.choice(b[g])) for g in selected]
        if len(set(aa))!=n or len(set(bb))!=n or set(aa)&set(bb):
            raise ValueError("Genus matched source taxon names not independent")
        p=_ball_mask(ref,aa,triple,
                     original_species_mask=original_species_mask)
        q=_ball_mask(ref,bb,triple,
                     original_species_mask=original_species_mask)
        values.append((
            _bin_value(ref,p,"rare")-_bin_value(ref,p,"common")
            -_bin_value(ref,q,"rare")+_bin_value(ref,q,"common"),
            _bin_value(ref,p,"rare")-_bin_value(ref,q,"rare"),
            _bin_value(ref,p,"common")-_bin_value(ref,q,"common")
        ))
    x=np.asarray(values,float).mean(axis=0)
    return dict(matched_genus_gap=float(x[0]),
                matched_genus_rare_coverage=float(x[1]),
                matched_genus_common_coverage=float(x[2]),
                n_global_source_shared_genera=len(shared))


def original_genus_shapley_macro(
    triple:pd.DataFrame,flora:pd.DataFrame,geo:pd.DataFrame,
    sources:pd.DataFrame,mainland:pd.DataFrame,*,
    repeats:int=40,reference_draws:int=399,
    bootstraps:int=1999,seed:int=20261113,
)->dict:
    if repeats<20 or reference_draws<199 or bootstraps<999:
        raise ValueError("Original source rare floral support contract cannot be weakened")
    for frame,need in (
      (flora,{"island_id","accepted_species","origin_status"}),
      (geo,{"island_id","analysis_regime","spatial_block"}),
      (sources,{"island_id","source_mode","source_rank","entity_ID"}),
      (mainland,{"entity_ID","work_species"})):
        if not need.issubset(frame):
            raise ValueError("Original source genus decomposition required fields absent")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original accepted species island ledger duplicated")
    if not flora.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Original native status invalid")
    original=flora.loc[flora.accepted_species.isin(triple.index)].copy()
    eligible=set(original.groupby("island_id").accepted_species.nunique().pipe(
        lambda x:x[x>=30]).index)
    if len(eligible)!=928:
        raise ValueError("Original 928 island source floral cohort not found")
    original=original.loc[original.island_id.isin(eligible)]
    geog=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates()
    if geog.island_id.duplicated().any() or geog.isna().any().any():
        raise ValueError("Original geography ambiguous")
    if not eligible.issubset(set(geog.island_id)):
        raise ValueError("One of the source original islands has no physical block")
    bygeo=geog.set_index("island_id")
    candidate=sources.loc[
       sources.island_id.isin(eligible)&
       sources.source_mode.eq("geo_k5")&
       sources.source_rank.between(1,5),
       ["island_id","source_rank","entity_ID"]]
    if candidate.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original source rank duplicated")
    donors=candidate.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(donors.get(i,()))!=5 or len(set(donors[i]))!=5 for i in eligible):
        raise ValueError("Source 5-nearest candidate native GIFT floras incomplete")
    gift=mainland.loc[
       mainland.work_species.isin(triple.index),
       ["entity_ID","work_species"]].drop_duplicates().groupby(
       "entity_ID").work_species.agg(set).to_dict()
    byisland={k:x for k,x in original.groupby("island_id")}
    ref=RarityReference.from_original(triple,null_draws=reference_draws,seed=seed)
    # Compiled exact per-species category-neighborhood lookup: each bit means
    # one of 100 actually observed original global reference profiles.
    species_mask={
        str(species):ref.volume.masks[tuple(row)][1]
        for species,row in zip(triple.index,
             triple[list(TRAITS)].to_numpy(dtype=str))
    }
    if len(species_mask)!=len(triple):
        raise ValueError("Direct global species IDs not unique")
    all_rows=[]
    genus_rows=[]
    for island in sorted(eligible):
        group=byisland[island]
        all_original_names=set(group.accepted_species)
        raw_source=set().union(*(gift.get(j,set()) for j in donors[island]))
        source=(raw_source&set(triple.index))-all_original_names
        if len(source)<30:
            continue
        for scope in SCOPES:
            members=(all_original_names if scope=="all_recorded" else
                     set(group.loc[group.origin_status.eq("native"),"accepted_species"]))
            if len(members)<30:
                continue
            metadata=dict(
                island_id=island,scope=scope,
                regime=str(bygeo.at[island,"analysis_regime"]),
                spatial_block=str(bygeo.at[island,"spatial_block"]),
                source_signature=donors[island][0],
                n_local_complete_Direct_species=len(members),
                n_source_complete_Direct_species=len(source))
            stats=_sample_components(ref,triple,members,source,
                n=30,repeats=repeats,seed=_seedstr(seed,island,scope,"shapley"),
                original_species_mask=species_mask)
            entry={**metadata,**stats}
            all_rows.append(entry)
            genus_count=len(_source_genus(members)&_source_genus(source))
            support={**metadata,"n_shared_original_genera":genus_count,
                "five_distinct_genera_support":genus_count>=MATCHED_FIVE,
                "eight_distinct_genera_support":genus_count>=MATCHED_EIGHT}
            for k in (MATCHED_FIVE,MATCHED_EIGHT):
                if genus_count>=k:
                    matched=_match_unique_genera(ref,triple,members,source,
                        n=k,repeats=repeats,seed=_seedstr(seed,island,scope,"g",k),
                        original_species_mask=species_mask)
                    for metric,v in matched.items():
                        support[f"g{k}_{metric}"]=v
            genus_rows.append(support)
    p=pd.DataFrame(all_rows)
    m=pd.DataFrame(genus_rows)
    regions=[]
    for (scope,regime),q in p.groupby(["scope","regime"]):
        if len(q)<10 or q.spatial_block.nunique()<5:
            continue
        rec=dict(scope=scope,regime=regime,n_original_eligible_islands=len(q),
            n_physical_blocks=q.spatial_block.nunique(),
            n_first_GIFT_source_signatures=q.source_signature.nunique(),
            median_shared_genera=int(q.n_shared_source_genera_in_original_flora.median()),
            median_local_genera=int(q.n_original_local_genera.median()),
            median_source_genera=int(q.n_original_source_genera.median()),
            n_original_source_disjoint_plant_names_ge_30=int(len(q)))
        for bin_ in KINDS:
            total=(q[f"{bin_}_shared"]+q[f"{bin_}_unshared"])-q[f"{bin_}_total"]
            if not np.allclose(total.to_numpy(),0,atol=1e-11):
                raise ValueError("Original Shapley and recorded profile support not additive")
        for j,key in enumerate((
            "gap_total","gap_shared","gap_unshared",
            "rare_total","rare_shared","rare_unshared",
            "common_total","common_shared","common_unshared")):
            rec["median_"+key]=float(q[key].median())
            for typ,cluster in (("block","spatial_block"),
                                ("source","source_signature")):
                a,lo,hi=_resample_group(q,key,cluster,seed+41*j+
                    (0 if typ=="block" else 1),draws=bootstraps)
                rec[typ+"_mean_"+key]=a
                rec[typ+"_ci025_"+key]=lo
                rec[typ+"_ci975_"+key]=hi
        regions.append(rec)
    matched=[]
    for (scope,regime),q in m.groupby(["scope","regime"]):
        for k in (MATCHED_FIVE,MATCHED_EIGHT):
            name=f"{'five' if k==5 else 'eight'}_distinct_genera_support"
            chosen=q.loc[q[name]].copy()
            nblocks=int(chosen.spatial_block.nunique())
            nsources=int(chosen.source_signature.nunique())
            accepted=(len(chosen)>=10 and nblocks>=5 and nsources>=5)
            rec=dict(scope=scope,regime=regime,
                matched_distinct_genera_per_island=k,
                n_eligible_source_disjoint_flora=len(q),
                n_islands_matching_genus_support=len(chosen),
                n_spatial_blocks=nblocks,n_first_GIFT_sources=nsources,
                reportable_under_10islands_5blocks_5sources=bool(accepted))
            if accepted:
                for j,metric in enumerate(("matched_genus_gap",
                      "matched_genus_rare_coverage","matched_genus_common_coverage")):
                    key=f"g{k}_{metric}"
                    rec["median_"+key]=float(chosen[key].median())
                    for typ,cluster in (("block","spatial_block"),
                                        ("source","source_signature")):
                        a,lo,hi=_resample_group(chosen,key,cluster,
                            seed+151*j+(0 if typ=="block" else 1),draws=bootstraps)
                        rec[typ+"_mean_"+key]=a
                        rec[typ+"_ci025_"+key]=lo
                        rec[typ+"_ci975_"+key]=hi
            matched.append(rec)
    return dict(
       status="ORIGINAL_ISLAND_RARE_PROFILE_GENUS_SHAPLEY_ACCOUNTING",
       n_original_eligible_Direct_islands=928,
       n_original_reference_Direct_species=len(triple),
       n_original_reference_triple_profiles=len(ref.volume.profiles),
       n_original_island_scope_rows_decomposed=len(p),
       reference_rarity_bins=ref.n_bins,
       all_source_disjoint_flora_genus_shapley_regions=regions,
       restricted_distinct_genus_source_pairing_panels=matched,
       method={
         "equal_sample_n_per_original_island_and_GIFT_source":30,
         "original_direct_rarefactions":repeats,
         "strict_same_plant_name_source_exclusion":True,
         "global_category_reference_all_original_Direct":True,
         "coverage_radius_exactly_one_categorical_trait_mismatch":True,
         "shapley_2_player_shared_vs_source_exclusive_original_genera":True,
         "exact_additive_shapley_coverage_proven":True,
         "shared_genus_decision_uses_original_complete_flora_not_sampled_traits":True,
         "five_unique_genus_one_island_species_one_source_species_control":True,
         "eight_unique_genus_support_sensitivity":True,
         "sampled_genus_matched_control_is_a_selected_subset":True,
         "block_and_nearest_GIFT_source_bootstrap_draws":bootstraps,
       },
       source_genera_are_observed_assemblage_not_verified_phylogeny=True,
       no_genetically_informed_within_species_adaptation_test=True,
       no_actual_pollinator_service_or_historical_decline_tested=True,
       no_original_island_taxon_species_genus_source_identity_released=True,
    )
