"""Strict same-genus species-match sensitivity for rare categorical floral volume.

This is a source/lineage-controlled **15-species restricted cohort** study,
not a full 928-island matched-genus null. It matches the EXACT number of plant
species per genus between each island and the five-nearest original native
GIFT candidate mainland regions; no matching on flower phenotype, outcome or
treatment. Every source plant taxon has a different accepted name from every
recorded local island plant (including unresolved-origin names).

Because the 928-island global source cannot generally support full n30
same-genus controls in the tropics, this runs a predeclared separate n15
estimand, reports attrition, and never extrapolates to the other islands.
"""
from __future__ import annotations

from collections import Counter,defaultdict
import numpy as np
import pandas as pd

from .rare_profile_decomposition import RarityReference,genus,_seedstr
from .discrete_functional_volume import TRAITS,_resample_group

N_MATCH=15
DRAWS=40
SUPPORT_MIN_ISLANDS=10
SUPPORT_MIN_BLOCKS=5
SUPPORT_MIN_SOURCE_SIGNATURES=5
METRICS=("rare_gap","k_residual_gap","frequency_weighted_k_residual_gap",
         "rare_coverage","common_coverage")


def _group_species(names:set[str])->dict[str,tuple[str,...]]:
    by=defaultdict(list)
    for x in sorted(names):
        key=genus(x)
        if key:
            by[key].append(x)
    return {k:tuple(v) for k,v in by.items()}


def _one_sample_score(
    reference:RarityReference,species_names:list[str],triple:pd.DataFrame,
)->np.ndarray:
    profiles=triple.loc[species_names,list(TRAITS)].to_numpy(dtype=str)
    coverage,k=reference.measure(profiles)
    raw=float(coverage[0]-coverage[2])
    expected=reference.expected_uniform[k]
    weighted=reference.expected_speciesfreq[k]
    return np.asarray([
        raw,
        raw-float(expected[0]-expected[2]),
        raw-float(weighted[0]-weighted[2]),
        float(coverage[0]),float(coverage[2])
    ],dtype=float)


def _genus_slots(
    island_species:set[str],source_species:set[str]
)->tuple[dict[str,tuple[str,...]],dict[str,tuple[str,...]],list[str]]:
    local=_group_species(island_species)
    mainland=_group_species(source_species-island_species)
    slots=[]
    for g in sorted(local):
        count=min(len(local[g]),len(mainland.get(g,())))
        slots.extend([g]*count)
    return local,mainland,slots


def _simulate(
    reference:RarityReference,triple:pd.DataFrame,
    island_names:set[str],source_names:set[str],*,n:int,repeats:int,seed:int
)->dict:
    """Return only local result (private). Source donor ≠ any island taxon."""
    if n!=15 or repeats<20:
        raise ValueError("Original genus-conditioned analysis floor changed")
    local,src,slots=_genus_slots(island_names,source_names)
    if len(slots)<n:
        raise ValueError("Insufficient outcome-blind same-genus source controls")
    if len(island_names)<n or len(source_names)<n:
        raise ValueError("Insufficient full nonmatched source controls")
    rand=np.random.default_rng(seed)
    matched=[]
    ordinary=[]
    for draw in range(repeats):
        genus_idx=rand.choice(len(slots),size=n,replace=False)
        chosen=Counter(slots[int(i)] for i in genus_idx)
        local_names=[]
        donor_names=[]
        for g in sorted(chosen):
            m=chosen[g]
            local_names.extend(rand.choice(local[g],size=m,replace=False).tolist())
            donor_names.extend(rand.choice(src[g],size=m,replace=False).tolist())
        if len(set(local_names))!=n or len(set(donor_names))!=n:
            raise ValueError("Nonindependent same-plant draw or genus mismatch")
        if Counter(map(genus,local_names))!=Counter(map(genus,donor_names)):
            raise ValueError("Source and island genus composition differs")
        matched.append(
            _one_sample_score(reference,local_names,triple)-
            _one_sample_score(reference,donor_names,triple)
        )
        # Unrestricted Direct island/source same N=15 to identify the
        # selected-cohort effect of *requiring shared source genera*.
        island_draw=rand.choice(sorted(island_names),size=n,replace=False).tolist()
        source_draw=rand.choice(sorted(source_names),size=n,replace=False).tolist()
        ordinary.append(
            _one_sample_score(reference,island_draw,triple)-
            _one_sample_score(reference,source_draw,triple)
        )
    return dict(
        matched=np.asarray(matched,float).mean(axis=0),
        unrestricted=np.asarray(ordinary,float).mean(axis=0),
        n_genus_source_opportunity_slots=len(slots),
        n_shared_source_genera=len({i for i in slots}),
        n_original_local_species=len(island_names),
        n_independent_source_species=len(source_names),
    )


def original_genus_matched_rare_reanalysis(
    triple:pd.DataFrame,original_flora:pd.DataFrame,
    geo:pd.DataFrame,assignments:pd.DataFrame,mainland:pd.DataFrame,*,
    repeats:int=40,reference_draws:int=399,
    bootstrap:int=1999,seed:int=20261113,
)->dict:
    if repeats<20 or reference_draws<199 or bootstrap<999:
        raise ValueError("Original genus-comparable 15 species reference design weakened")
    for t,cols in (
      (original_flora,("island_id","accepted_species","origin_status")),
      (geo,("island_id","analysis_regime","spatial_block")),
      (assignments,("island_id","source_mode","source_rank","entity_ID")),
      (mainland,("entity_ID","work_species"))
    ):
        if not set(cols).issubset(t):
            raise ValueError("Original source taxonomic/geographic contract mismatch")
    original=original_flora.loc[
        original_flora.accepted_species.isin(triple.index)
    ].copy()
    if original.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original accepted plant duplicated on island")
    if not original.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Unrecognized original native status")
    full=original.groupby("island_id").accepted_species.nunique()
    eligible=set(full.loc[full>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original 928 macro flora population drift")
    original=original.loc[original.island_id.isin(eligible)]
    geo=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates()
    if geo.island_id.duplicated().any() or not eligible.issubset(set(geo.island_id)):
        raise ValueError("Original blocks do not cover exactly eligible flora")
    loc=geo.set_index("island_id")
    assignment=assignments.loc[
       assignments.island_id.isin(eligible)&
       assignments.source_mode.eq("geo_k5")&
       assignments.source_rank.between(1,5),
       ["island_id","source_rank","entity_ID"]]
    if assignment.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate candidate GIFT source rank")
    near=assignment.sort_values(["island_id","source_rank"]).groupby(
         "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(near.get(i,()))!=5 or len(set(near[i]))!=5 for i in eligible):
        raise ValueError("Original five source entities missing")
    gifts=mainland.loc[
        mainland.work_species.isin(triple.index),
        ["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    reference=RarityReference.from_original(triple,
        null_draws=reference_draws,seed=seed)
    grouped={i:f for i,f in original.groupby("island_id")}
    rows=[]
    opportunity=[]
    for island in sorted(eligible):
        nativepool=set(grouped[island].loc[
            grouped[island].origin_status.eq("native"),"accepted_species"])
        allpool=set(grouped[island].accepted_species)
        candidate=set().union(*(gifts.get(k,set()) for k in near[island]))
        candidate=candidate.intersection(triple.index)-allpool
        for status,names in (("all_recorded",allpool),("source_native",nativepool)):
            local,donors,slots=_genus_slots(names,candidate)
            row=dict(scope=status,regime=str(loc.at[island,"analysis_regime"]),
               island_id=island,spatial_block=str(loc.at[island,"spatial_block"]),
               source_signature=near[island][0],
               n_same_genus_slots=len(slots),n_shared_genera=len(set(slots)),
               n_local_names=len(names),n_source_disjoint_names=len(candidate),
               n_matchable_15=len(slots)>=N_MATCH,
               n_matchable_30=len(slots)>=30)
            opportunity.append(row)
            if len(slots)<N_MATCH:
                continue
            score=_simulate(reference,triple,names,candidate,
                n=N_MATCH,repeats=repeats,seed=_seedstr(seed,island,status))
            for k,measure in enumerate(METRICS):
                row["matched15_delta_"+measure]=float(score["matched"][k])
                row["unrestricted15_delta_"+measure]=float(score["unrestricted"][k])
                row["genus_matching_vs_unrestricted15_"+measure]=float(
                    score["matched"][k]-score["unrestricted"][k])
            rows.append(row)
    potential=pd.DataFrame(opportunity)
    private=pd.DataFrame(rows)
    output=[]
    support=[]
    for (scope,regime),q in potential.groupby(["scope","regime"]):
        matching=q.loc[q.n_matchable_15]
        candidate_count=int(len(matching))
        blocks=int(matching.spatial_block.nunique())
        signatures=int(matching.source_signature.nunique())
        eligibility=(candidate_count>=SUPPORT_MIN_ISLANDS and
             blocks>=SUPPORT_MIN_BLOCKS and
             signatures>=SUPPORT_MIN_SOURCE_SIGNATURES)
        report=dict(scope=scope,regime=regime,
           n_original_direct_30_species_islands=int(len(q)),
           n_with_complete_same_genus_match_15_islands=candidate_count,
           n_with_complete_same_genus_match_30_islands=int(q.n_matchable_30.sum()),
           n_original_physical_blocks_match15=blocks,
           n_first_GIFT_source_signatures_match15=signatures,
           median_max_same_genus_source_slots=float(q.n_same_genus_slots.median()),
           eligible_to_report_conditional_genus15_comparison=bool(eligibility))
        support.append(report)
        if not eligibility:
            continue
        gp=private.loc[
           private.scope.eq(scope)&private.regime.eq(regime)]
        if len(gp)!=candidate_count:
            raise ValueError("Genus matched source execution lost eligible island")
        rec=dict(report)
        for j,metric in enumerate(METRICS):
            for m,ctype in enumerate(("matched15_delta_",
                  "unrestricted15_delta_","genus_matching_vs_unrestricted15_")):
                key=ctype+metric
                rec["median_"+key]=float(gp[key].median())
                for weight,group in (("block","spatial_block"),
                                     ("source","source_signature")):
                    mean,lo,hi=_resample_group(gp,key,group,
                          seed+100*m+20*j+(weight=="source"),draws=bootstrap)
                    rec[weight+"_equal_mean_"+key]=mean
                    rec[weight+"_ci025_"+key]=lo
                    rec[weight+"_ci975_"+key]=hi
        output.append(rec)
    return {
        "status":"GENUS_MATCHED_SOURCE_DISJOINT_CATEGORICAL_RARITY_15_SPECIES",
        "n_original_complete_Direct_reference_species":len(triple),
        "n_original_928_eligible_islands":len(eligible),
        "same_genus_source_support_by_scope_region":support,
        "same_genus_source_matched15_effect_by_scope_region":output,
        "n_genus_matched_source_effect_region_panels_reported":len(output),
        "design":{
            "n_species_exact_per_island_and_mainland":N_MATCH,
            "n_fixed_direct_morphology_rarefactions":repeats,
            "n_speciesfreq_and_uniform_reference_k_null_draws":reference_draws,
            "n_conditional_geography_source_bootstraps":bootstrap,
            "require_at_least_10_physical_islands_and_5_blocks_and_5_sources":True,
            "identical_genus_counts_per_island_source_draw":True,
            "original_exact_local_plant_names_removed_from_all_candidate_sources":True,
            "source_identity_not_verified_as_ancestor":True,
            "sample_selection_restricted_to_source_genus_exchangeable_flora":True,
            "reported_only_if_source_genus_support_gates_pass":True,
            "not_global_full_928_genus_conditional_null":True,
        },
        "phylogeny_provenance_observation_effort_not_controlled":True,
        "actual_pollinator_loss_or_evolution_not_measured":True,
        "protected_taxon_species_island_genus_or_GIFT_entity_never_exported":True,
    }
