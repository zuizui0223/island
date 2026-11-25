"""Island versus GIFT donor rarity breadth with EQUAL genus richness and species N.

Unlike exact SAME-genera controls, this isolates the effect of different
genus IDENTITIES while both samples contain exactly 15 DISTINCT genera and
exactly one different accepted Direct plant species from each genus.

Other estimands:
 - fixed n=30 source draw can have unequal genus richness;
 - fixed n=5 identical genera tests different species within overlapping
   recorded genera;
 - this separate fixed n=15 distinct genera permits unlike genus IDs.
All three are source-dependent descriptive trait-composition contrasts.
No plant species phenotype varies by location in the source database.
"""
from __future__ import annotations

from collections import defaultdict
import numpy as np
import pandas as pd

from .discrete_functional_volume import TRAITS, _resample_group
from .rare_profile_decomposition import RarityReference,genus,_seedstr

PRIMARY_G=15
SENSITIVITY_G=10
OUTCOMES=("rare","common","rare_minus_common",
          "k_residual_rare_minus_common",
          "speciesfreq_k_residual_rare_minus_common")


def _by_genus(accepted_names:set[str])->dict[str,tuple[str,...]]:
    d=defaultdict(list)
    for n in sorted(accepted_names):
        g=genus(n)
        if g:
            d[g].append(n)
    return {g:tuple(v) for g,v in d.items()}


def _sample_one_per_genus(by_genus:dict[str,tuple[str,...]],*,n:int,
                          rng:np.random.Generator)->list[str]:
    if n not in (PRIMARY_G,SENSITIVITY_G):
        raise ValueError("Original distinct genus count contract changed")
    if len(by_genus)<n:
        raise ValueError("Too few Direct species-supported genera")
    chosen=rng.choice(sorted(by_genus),size=n,replace=False)
    species=[str(rng.choice(by_genus[str(g)])) for g in chosen]
    if len(set(species))!=n or len({genus(name) for name in species})!=n:
        raise AssertionError("Independent genus-unit draw not one taxon/genus")
    return species


def _one_score(reference:RarityReference,triple:pd.DataFrame,
               names:list[str], *,
               source_profile_bitmap:dict[str,tuple[tuple[str,str,str],int]]|None=None)->np.ndarray:
    if source_profile_bitmap is None:
        p=triple.loc[names,list(TRAITS)].to_numpy(dtype=str)
        coverage,k=reference.measure(p)
    else:
        seen=set()
        covered=0
        for name in names:
            profile,bitmap=source_profile_bitmap[name]
            seen.add(profile)
            covered |= bitmap
        k=len(seen)
        coverage=np.asarray([
            (covered & reference.bins[c]).bit_count()/reference.n_bins[c]
            for c in ("rare","middle","common")
        ],dtype=float)
    u=reference.expected_uniform[k]
    w=reference.expected_speciesfreq[k]
    return np.asarray([float(coverage[0]),float(coverage[2]),
        float(coverage[0]-coverage[2]),
        float(coverage[0]-coverage[2]-u[0]+u[2]),
        float(coverage[0]-coverage[2]-w[0]+w[2])])


def _one_island(
    reference:RarityReference,triple:pd.DataFrame,
    local:set[str],source:set[str],*,n:int=15,
    draws:int=40,seed:int=20261113,
    source_profile_bitmap:dict[str,tuple[tuple[str,str,str],int]]|None=None
)->dict:
    if draws<20:
        raise ValueError("Insufficient original reference genus-unit thinnings")
    if local&source:
        raise ValueError("Original local accepted plant must be disjoint from GIFT source")
    a=_by_genus(local); b=_by_genus(source)
    if len(a)<n or len(b)<n:
        raise ValueError("Original genus-richness matching not feasible")
    rng=np.random.default_rng(seed)
    results=[]
    for _ in range(draws):
        island_names=_sample_one_per_genus(a,n=n,rng=rng)
        donor_names=_sample_one_per_genus(b,n=n,rng=rng)
        results.append(_one_score(reference,triple,island_names,
                           source_profile_bitmap=source_profile_bitmap)-
                       _one_score(reference,triple,donor_names,
                           source_profile_bitmap=source_profile_bitmap))
    z=np.asarray(results,float).mean(axis=0)
    return dict(zip(OUTCOMES,map(float,z)))


def original_genus_richness_normalization(
    triple:pd.DataFrame,flora:pd.DataFrame,
    geo:pd.DataFrame,sources:pd.DataFrame,mainland:pd.DataFrame,*,
    repeats:int=40,reference_draws:int=399,
    bootstrap:int=1999,seed:int=20261113,
)->dict:
    if repeats<20 or reference_draws<199 or bootstrap<999:
        raise ValueError("Original source genus richness control support weakened")
    for tab,cols in (
        (flora,{"island_id","accepted_species","origin_status"}),
        (geo,{"island_id","analysis_regime","spatial_block"}),
        (sources,{"island_id","source_mode","source_rank","entity_ID"}),
        (mainland,{"entity_ID","work_species"})):
        if not cols.issubset(tab):
            raise ValueError("Missing original GIFT flora/genus geographic source columns")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate species/island origin data")
    if not flora.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Original source nativeness unsupported")
    original=flora.loc[flora.accepted_species.isin(triple.index)]
    n=original.groupby("island_id").accepted_species.nunique()
    eligible=set(n[n>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original 928 Direct island sampling cohort drifted")
    original=original.loc[original.island_id.isin(eligible)].copy()
    geo=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates()
    if geo.island_id.duplicated().any() or not eligible.issubset(set(geo.island_id)):
        raise ValueError("Original region/spatial-block source invalid")
    loc=geo.set_index("island_id")
    assignments=sources.loc[sources.island_id.isin(eligible) &
        sources.source_mode.eq("geo_k5") & sources.source_rank.between(1,5),
        ["island_id","source_rank","entity_ID"]].copy()
    if assignments.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Conflicting original source rank within island")
    donors=assignments.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(donors.get(i,()))!=5 or len(set(donors[i]))!=5 for i in eligible):
        raise ValueError("Original five-nearest GIFT lists not preserved")
    gift=(mainland.loc[mainland.work_species.isin(triple.index),
        ["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict())
    byisland={i:g for i,g in original.groupby("island_id")}
    ref=RarityReference.from_original(triple,null_draws=reference_draws,seed=seed)
    source_profile_bitmap={}
    for name,row in zip(triple.index,triple[list(TRAITS)].to_numpy(dtype=str)):
        profile=tuple(row)
        source_profile_bitmap[str(name)]=(profile,ref.volume.masks[profile][1])
    if len(source_profile_bitmap)!=len(triple):
        raise ValueError("Original Direct plant species identifiers duplicated")
    support=[]
    rows=[]
    for i in sorted(eligible):
        group=byisland[i]
        all_local=set(group.accepted_species)
        all_source=set().union(*(gift.get(k,set()) for k in donors[i]))
        all_source=(all_source & set(triple.index))-all_local
        donor_g=_by_genus(all_source)
        for scope,names in (
            ("all_recorded",all_local),
            ("source_native",set(group.loc[group.origin_status.eq("native"),"accepted_species"]))
        ):
            if len(names)<30:
                continue
            local_g=_by_genus(names)
            meta=dict(
                island_id=i,scope=scope,regime=str(loc.at[i,"analysis_regime"]),
                spatial_block=str(loc.at[i,"spatial_block"]),
                source_signature=donors[i][0],
                n_local_genus=len(local_g),
                n_source_genus=len(donor_g),
                n_local_species=len(names),
                n_donor_species=len(all_source))
            for n_genus in (PRIMARY_G,SENSITIVITY_G):
                record=dict(meta,n_genus=n_genus,
                  has_genus_richness_support=(
                    len(local_g)>=n_genus and len(donor_g)>=n_genus))
                support.append(record)
                if not record["has_genus_richness_support"]:
                    continue
                z=_one_island(ref,triple,names,all_source,n=n_genus,
                    draws=repeats,seed=_seedstr(seed,i,scope,n_genus),
                    source_profile_bitmap=source_profile_bitmap)
                rows.append(dict(record,**z))
    p=pd.DataFrame(support)
    actual=pd.DataFrame(rows)
    out=[]
    for (scope,regime,n_genus),group in p.groupby(["scope","regime","n_genus"]):
        eligible_group=group.loc[group.has_genus_richness_support]
        rec=dict(scope=scope,regime=regime,n_distinct_genera_per_side=n_genus,
            n_original_30_species_direct_flora_islands=len(group),
            n_supported_islands=len(eligible_group),
            n_spatial_blocks=eligible_group.spatial_block.nunique(),
            n_first_nearest_GIFT_source_signatures=eligible_group.source_signature.nunique(),
            median_original_island_genus=int(group.n_local_genus.median()),
            median_original_source_genus=int(group.n_source_genus.median()))
        sufficient=(len(eligible_group)>=10 and
            eligible_group.spatial_block.nunique()>=5 and
            eligible_group.source_signature.nunique()>=5)
        rec["regional_comparison_reportable"]=bool(sufficient)
        if sufficient:
            sub=actual.loc[
                actual.scope.eq(scope)&actual.regime.eq(regime)&
                actual.n_genus.eq(n_genus)]
            if len(sub)!=len(eligible_group):
                raise ValueError("Genus richness controlled site omitted from source")
            for j,metric in enumerate(OUTCOMES):
                rec["median_"+metric]=float(sub[metric].median())
                for weight,cluster in (("block","spatial_block"),
                                        ("source","source_signature")):
                    mean,lo,hi=_resample_group(sub,metric,cluster,
                        seed+100*n_genus+12*j+(weight=="source"),
                        draws=bootstrap)
                    rec[weight+"_mean_"+metric]=mean
                    rec[weight+"_ci025_"+metric]=lo
                    rec[weight+"_ci975_"+metric]=hi
        out.append(rec)
    return dict(
        status="ORIGINAL_GENUS_RICHNESS_STANDARDIZED_RARE_CATEGORY_SUPPORT",
        original_complete_Direct_species=len(triple),
        n_original_eligible_islands=len(eligible),
        original_global_profile_count=len(ref.volume.profiles),
        equal_genus_richness_panels=out,
        n_protected_source_genus_unit_rows=len(actual),
        design={
            "one_different_accepted_plant_species_per_genus":True,
            "same_genus_number_but_genus_IDENTITY_NOT_matched":True,
            "primary_15_distinct_genera_each_side":True,
            "sensitivity_10_distinct_genera_each_side":True,
            "all_local_accepted_species_names_excluded_from_candidate_sources":True,
            "reference_original_100_profile_global_Direct_rarity":True,
            "fixed_sample_rarefactions":repeats,
            "fixed_k_null_reference_uniform_and_species_frequency":reference_draws,
            "conditionally_equal_source_and_spatial_block_bootstraps":bootstrap,
        },
        observed_genus_species_source_not_true_phylogenetic_ancestry=True,
        original_species_triple_state_not_population_specific=True,
        no_actual_pollinator_loss_or_fitness_inferred=True,
        no_original_source_genus_taxon_island_or_GIFT_name_exported=True,
    )
