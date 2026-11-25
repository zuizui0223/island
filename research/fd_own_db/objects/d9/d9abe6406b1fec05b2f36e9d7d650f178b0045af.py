"""Is the endemic island-vs-mainland floral discordance different by GEOLOGY in shared genera?

This is an outcome-blind, accepted-name-disjoint and GIFT-source-qualified
cross-origin taxonomic matching audit, not a demonstration of adaptation.
Geological origin is independently Roeble/Zizka classified, NOT LGM connection.
The sample unit for a genus-level contrast is the DISTINCT botanical genus,
not the duplicated island×plant record. Distinct species and geographic/source
group counts are all audited before evaluating directional contrasts.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from itertools import product

import numpy as np
import pandas as pd

from .geological_endemic_trait_opportunity import TRAITS
from .geological_endemic_trait_divergence import strict_direct_trait_states, plant_source_discordance

GEOLOGIES=("oceanic","continental")
PRIMARY=("flower_primary_color","floral_symmetry")
SAMPLES=1999
MIN_SHARED_GENERA=5
MIN_GEO_BLOCKS=5
MIN_GIFT_SOURCES=5


def _seed(*parts):
    return int.from_bytes(sha256(":".join(map(str,parts)).encode()).digest()[:8],"big")


def _genus(name):
    a=str(name).split()
    return a[0] if len(a)>=2 else ""


def _exact_signed_genus_sensitivity(differences:np.ndarray)->dict:
    """Exact 2^N paired-genus sign sensitivity, noncausal observational null.

    The original 10-degree geographic/source grouping remains potentially
    shared, so this is a descriptive stress test, not an island-formation
    treatment-randomization p-value.
    """
    a=np.asarray(differences,dtype=float)
    if a.ndim!=1 or len(a)<5 or len(a)>16 or not np.isfinite(a).all():
        raise ValueError("Signed genus audit requires 5..16 finite matched genera")
    n=len(a)
    indicators=np.arange(1<<n,dtype=np.uint32)[:,None]
    bits=np.arange(n,dtype=np.uint32)[None,:]
    signed=np.where(((indicators >> bits)&1)==1,1.,-1.)
    null=np.abs(signed@a/n)
    observed=abs(float(np.mean(a)))
    loo=(np.sum(a)-a)/(n-1)
    return dict(
        exact_two_sided_genus_label_sign_flip_p=float(
            np.mean(null >= observed-1e-12)),
        n_exact_sign_reassignments=int(1<<n),
        leave_one_genus_out_effect_range=[
            float(np.min(loo)),float(np.max(loo))],
        leave_one_genus_out_same_direction_fraction=float(
            np.mean(np.sign(loo)==np.sign(np.mean(a)))),
        n_genus_pairs_with_nonzero_difference=int(np.count_nonzero(
            np.abs(a)>1e-12)),
    )


def calculate_endemic_shared_genus_geology(
    ledger:pd.DataFrame,
    original_status:pd.DataFrame,
    geography:pd.DataFrame,
    candidate_sources:pd.DataFrame,
    native_mainland:pd.DataFrame,
    *,seed:int=20261114,draws:int=SAMPLES,
)->dict:
    """Analysis of observed SINGLE-trait variation, not a trait-profile expansion.

    All original accepted endemic plant names stay private. Counts are
    independent taxon-unit evidence gates; no candidate mainland taxon is
    declared the historical ancestor.
    """
    if draws<999:
        raise ValueError("Insufficient cluster bootstrap draws")
    required={"island_id","accepted_species","origin_status",
              "endemic_status","floristic_status"}
    if not required.issubset(original_status) or original_status.duplicated(
        ["island_id","accepted_species"]).any():
        raise ValueError("Original plant endemicity source nonunique/incomplete")
    orig=original_status
    if not orig.origin_status.isin(("native","introduced","unresolved")).all() or not (
        orig.endemic_status.isin(("endemic","nonendemic","unresolved")).all()):
        raise ValueError("Original native/endemic status labels changed")
    if ((orig.floristic_status.eq("endemic") &
          (~orig.origin_status.eq("native") | ~orig.endemic_status.eq("endemic"))).any()):
        raise ValueError("Original endemic plant native origin not proved")
    geo_cols={"island_id","geological_origin","analysis_regime","spatial_block"}
    if not geo_cols.issubset(geography) or geography.island_id.duplicated().any():
        raise ValueError("Independent geological geography invalid")
    typed=geography.loc[geography.geological_origin.isin(GEOLOGIES)]
    if len(typed)<4:
        raise ValueError("No source verified geology")
    known=orig.loc[orig.island_id.isin(typed.island_id)&
                   orig.floristic_status.eq("endemic")].merge(
                       typed[["island_id","geological_origin","analysis_regime","spatial_block"]],
                       on="island_id",how="inner",validate="many_to_one")
    if known.empty:
        raise ValueError("No native+endemic original labelled species on geological islands")
    states=strict_direct_trait_states(ledger)
    req={"island_id","entity_ID","source_mode","source_rank"}
    if not req.issubset(candidate_sources):
        raise ValueError("Original five-nearest donor evidence missing")
    nearby=candidate_sources.loc[
        candidate_sources.source_mode.eq("geo_k5")].copy()
    nearby["source_rank"]=pd.to_numeric(nearby.source_rank,errors="raise")
    nearby=nearby.loc[nearby.source_rank.between(1,5)]
    if nearby.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Repeated source rank on island")
    byid=nearby.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if not {"entity_ID","work_species"}.issubset(native_mainland):
        raise ValueError("Original GIFT source flora identity missing")
    source_pool=native_mainland[["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    all_by_island=orig.loc[orig.island_id.isin(typed.island_id)].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    private=[]
    for island,group in known.groupby("island_id",sort=True):
        five=byid.get(str(island),())
        if len(five)!=5 or len(set(five))!=5:
            continue
        source=set().union(*(source_pool.get(x,set()) for x in five))
        source-=all_by_island[str(island)]
        endemic=set(group.accepted_species)
        for trait in TRAITS:
            st=states[trait]
            local=endemic & set(st)
            if not local:
                continue
            donor_states={name:st[name] for name in source if name in st}
            if not donor_states:
                continue
            matched=plant_source_discordance(local,donor_states,st)
            for row in matched:
                name=row["accepted_species"]
                private.append(dict(
                    region=group.analysis_regime.iloc[0],
                    geology=group.geological_origin.iloc[0],
                    trait=trait,
                    island_id=str(island),
                    block=group.spatial_block.iloc[0],
                    first_source=five[0],
                    genus=_genus(name),
                    taxon=name,
                    state=st[name],
                    discordance=row["discordance"],
                    n_source_congeners=row["source_genus_congener_count"],
                ))
    data=pd.DataFrame(private)
    if data.empty:
        raise ValueError("No original source-Direct endemic congener pairs")
    if data[["island_id","trait","taxon"]].duplicated().any():
        raise ValueError("Original focal endemic trait repeated within island")
    strata=[]
    for (regime,geology,trait),g in data.groupby(
        ["region","geology","trait"],sort=True
    ):
        strata.append(dict(
            region=regime,geological_origin=geology,trait=trait,
            original_island_taxon_rows=int(len(g)),
            n_original_individual_islands=int(g.island_id.nunique()),
            n_original_distinct_endemic_plant_species=int(g.taxon.nunique()),
            n_original_distinct_endemic_botanical_genera=int(g.genus.nunique()),
            n_original_independent_10degree_blocks=int(g.block.nunique()),
            n_original_first_GIFT_donor_region_signatures=int(g.first_source.nunique()),
            median_original_source_congeners=float(g.n_source_congeners.median()),
            ratio_original_island_taxon_rows_to_distinct_endemic_species=float(
                len(g)/g.taxon.nunique()),
        ))
    contrasts=[]
    for (regime,trait),g in data.groupby(["region","trait"],sort=True):
        o=g.loc[g.geology.eq("oceanic")].copy()
        c=g.loc[g.geology.eq("continental")].copy()
        # Species recorded on BOTH geological types in the same stratum
        # cannot establish distinct endemic plant lineages across origins.
        both=set(o.taxon)&set(c.taxon)
        o=o.loc[~o.taxon.isin(both)]
        c=c.loc[~c.taxon.isin(both)]
        og=set(o.genus)
        cg=set(c.genus)
        common=sorted((og&cg)-{""})
        # On each geology, first collapse repeated island records of the
        # SAME accepted species; then equal weight distinct species within
        # each genus; lastly one weight per matched genus.
        estimates={}
        for grp,label in ((o,"oceanic"),(c,"continental")):
            genus_vals={}
            genus_states={}
            for genus,sub in grp.groupby("genus"):
                if genus not in common:
                    continue
                t=sub.groupby("taxon").discordance.mean()
                genus_vals[genus]=float(t.mean())
                genus_states[genus]=set(sub.state)
            estimates[label]=(genus_vals,genus_states)
        common=[g for g in common if
                g in estimates["oceanic"][0] and g in estimates["continental"][0]]
        differences=np.array([
            estimates["oceanic"][0][genus]-estimates["continental"][0][genus]
            for genus in common],dtype=float)
        n_variable=sum(
            len(estimates["oceanic"][1][genus]|
                estimates["continental"][1][genus])>1 for genus in common)
        so=o.loc[o.genus.isin(common)]
        sc=c.loc[c.genus.isin(common)]
        rec=dict(
            region=regime,trait=trait,
            primary_trait=bool(trait in PRIMARY),
            n_endemic_species_shared_by_both_geologies_excluded=len(both),
            n_distinct_oceanic_plant_species_source_qualified=int(o.taxon.nunique()),
            n_distinct_continental_plant_species_source_qualified=int(c.taxon.nunique()),
            n_distinct_oceanic_source_qualified_genera=int(len(og)),
            n_distinct_continental_source_qualified_genera=int(len(cg)),
            n_same_genus_shared_across_geological_types=len(common),
            n_shared_genera_with_different_original_endemic_trait_categories=int(n_variable),
            n_oceanic_islands_covered_by_shared_genera=int(so.island_id.nunique()),
            n_continental_islands_covered_by_shared_genera=int(sc.island_id.nunique()),
            n_oceanic_blocks_covered_by_shared_genera=int(so.block.nunique()),
            n_continental_blocks_covered_by_shared_genera=int(sc.block.nunique()),
            n_oceanic_first_source_groups_covered_by_shared_genera=int(so.first_source.nunique()),
            n_continental_first_source_groups_covered_by_shared_genera=int(sc.first_source.nunique()),
        )
        support=(len(common)>=MIN_SHARED_GENERA and n_variable>=2 and
                 min(so.block.nunique(),sc.block.nunique())>=MIN_GEO_BLOCKS and
                 min(so.first_source.nunique(),sc.first_source.nunique())>=MIN_GIFT_SOURCES)
        rec["qualified_to_compare_source_endemic_discordance_within_shared_genera"]=bool(support)
        if support:
            rng=np.random.default_rng(_seed(seed,"shared_genus",regime,trait))
            draw=differences[rng.integers(len(differences),size=(draws,len(differences)))].mean(axis=1)
            rec["paired_genus_oceanic_minus_continental_discordance"]=dict(
                mean=float(differences.mean()),
                ci95=[float(x) for x in np.quantile(draw,[.025,.975])],
                n_original_independent_genus_units=len(differences),
                n_independent_shared_state_variable_genera=n_variable,
            )
            if len(differences)<=16:
                rec["descriptive_exact_genus_sign_and_loo_sensitivity"]=(
                    _exact_signed_genus_sensitivity(differences))
        contrasts.append(rec)
    return dict(
        status="ORIGINAL_GEOLOGICAL_ENDEMIC_SHARED_GENUS_SOURCE_PAIRING_AUDIT",
        original_geological_islands_with_valid_labels=int(len(typed)),
        n_original_distinct_endemic_source_qualified_trait_island_plant_rows=len(data),
        primary_traits=list(PRIMARY),
        original_endemic_species_and_independent_genus_opportunity=strata,
        original_shared_genus_geological_comparisons=contrasts,
        method=dict(
            accepted_plant_species_not_island_records_are_taxon_units=True,
            repeated_same_accepted_species_deduplicated_before_genus_means=True,
            same_accepted_species_in_both_geological_classes_excluded=True,
            genus_equal_bootstrap_not_unjustified_island_independence=True,
            per_trait_not_all_three_flower_traits_required=True,
            source_assignment_GIFT_geo_k5_excludes_all_original_island_taxa=True,
            correct_original_floristic_endemic_and_native_source_status=True,
            original_focal_trait_variation_required_for_effect_reporting=True,
            source_congeners_are_not_true_ancestors=True,
            original_plant_state_fixed_per_species_not_local_population=True,
            no_pollinator_shift_or_evolved_adaptation_identified=True,
            no_privileged_botanical_taxa_or_island_IDs_exported=True,
        ),
    )
