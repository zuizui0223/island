"""Source-only NHM floral COLOUR×FORM joint-data feasibility against all 4505 GIFT island floras.

Use all five pinned separate, independently released NHM text corpora. A donor
must be an original named GIFT native mainland region-list species with BOTH
single-state whole-flower primary colour and floral shape on the same exact
accepted species, and must be absent from EVERY original 4505 island flora.
BHL, eFloras and Ecoflora are text documentation sources, NOT independent
mainland field phenotype populations, ancestral source or pollen effectiveness.
Never fabricate missing colour/form cross-products or allow island test-name
reuse in any source region. No rights-protected taxon/region/island IDs exported.
"""
from __future__ import annotations
import argparse
import json
import hashlib
from collections import defaultdict
from pathlib import Path
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_nhm_original_provider_audit import (
    PROVIDERS,BLOB_PINS,read_provider_file,consensus_provider_profiles,
)
from fd_nhm_external_flora_preflight import (
    nhm_git_blob_sha,nhm_external_categories
)

MIN_DIRECT=5


def cross_corpus_joint(by_provider:dict):
    """Return strong same-corpus and weaker different-corpus trait opportunity."""
    if not by_provider:raise ValueError("No independently pinned corpus")
    primary,primary_dual,primary_conflicts,_=consensus_provider_profiles(
        by_provider,axis="primary")
    form,form_dual,form_conflicts,_=consensus_provider_profiles(
        by_provider,axis="form")
    combined={name:(primary[name],form[name])
              for name in set(primary)&set(form)}
    raw_same=defaultdict(set)
    two_both=defaultdict(set)
    for provider,axes in by_provider.items():
        direct_same=set(axes["primary"])&set(axes["form"])
        for name in direct_same:
            raw_same[name].add((axes["primary"][name],axes["form"][name]))
            two_both[name].add(provider)
    same={name:next(iter(pairs)) for name,pairs in raw_same.items()
          if len(pairs)==1 and name in combined
          and next(iter(pairs))==combined[name]}
    twice={name:same[name] for name,providers in two_both.items()
           if len(providers)>=2 and name in same}
    if not set(twice)<=set(same)<=set(combined):
        raise AssertionError("Inconsistent external phenotype evidence hierarchy")
    return dict(same_published_corpus_both_traits=same,
                any_separate_corpora_consistent_both_traits=combined,
                two_distinct_named_corpora_both_traits_agree=twice),dict(
                n_conflicted_primary_colour_species=primary_conflicts,
                n_conflicted_form_species=form_conflicts,
                n_dual_corpus_primary_species=len(primary_dual),
                n_dual_corpus_form_species=len(form_dual))


def source_provenance_opportunity(donor_names:set,all_flora:dict,
                                  eligible_islands:list[str],
                                  target_pair:set,
                                  mainland:pd.DataFrame,
                                  assignment:pd.DataFrame):
    """Opportunity only, on entire preselected ORIGINAL target island population."""
    full_names=set().union(*all_flora.values())
    gift=set(mainland.work_species)
    valid=(set(donor_names)&gift)-full_names
    if valid&full_names:raise AssertionError("Original island accepted taxon leaked")
    source=mainland.loc[mainland.work_species.isin(valid)]
    by_region=source.groupby("entity_ID").work_species.agg(set).to_dict()
    near=assignment.loc[assignment.source_mode.eq("geo_k5")].copy()
    near["source_rank"]=pd.to_numeric(near.source_rank,errors="raise")
    near=near.loc[near.source_rank.between(1,5)]
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original geographical source ranking")
    by_island=near.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    n_islands_with_1=0
    n_islands_with_5=0
    n_islands_with_source=0
    n_islands_with_genus=0
    n_islands_missing_source=0
    target_total=0
    target_congeners=0
    source_donor_island_pairs=0
    for island in eligible_islands:
        orig=all_flora[island]
        target=orig&target_pair
        if len(target)<MIN_DIRECT:
            raise ValueError("Original pair-trait eligibility changed")
        target_total+=len(target)
        entities=by_island.get(island,())
        if not entities:
            n_islands_missing_source+=1
        else:
            n_islands_with_source+=1
        local=set().union(*(by_region.get(x,set()) for x in entities))
        local_count=len(local)
        source_donor_island_pairs+=local_count
        if local_count>=1:n_islands_with_1+=1
        if local_count>=5:n_islands_with_5+=1
        donor_genera={name.split()[0] for name in local}
        count=sum(x.split()[0] in donor_genera for x in target)
        if count:n_islands_with_genus+=1
        target_congeners+=count
    return dict(
        n_source_only_names_with_both_traits=len(valid),
        n_source_only_distinct_genera=len({x.split()[0] for x in valid}),
        n_GIFT_mainland_regions_with_any_source_only_joint_donor=len(by_region),
        n_target_islands_in_original_pairwise_population=len(eligible_islands),
        n_target_islands_with_geographic_source_assignment=n_islands_with_source,
        n_target_islands_missing_geographic_source_assignment=n_islands_missing_source,
        n_target_islands_with_at_least_one_external_joint_source_donor=n_islands_with_1,
        n_target_islands_with_at_least_five_external_joint_source_donors=n_islands_with_5,
        mean_source_only_joint_donor_names_per_target_island=(
            source_donor_island_pairs/len(eligible_islands)),
        n_direct_colour_form_target_plant_carrier_occurrences=target_total,
        n_target_plant_carriers_with_nearest_external_joint_congener=target_congeners,
        fraction_direct_joint_plant_carriers_with_source_genus_support=(
            target_congeners/target_total),
        n_target_islands_with_any_source_congeneric_external_joint_donor=n_islands_with_genus,
    )


def analyze(traits,flora,mainland,assignment,raw_corpora):
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original plant flora columns missing")
    if flora.duplicated(["island_id","accepted_species"]).any() or (
        not flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original species/status source invalid")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:
        raise ValueError("Original 4,505 botanical island floras missing")
    if not {"entity_ID","work_species"}<=set(mainland) or mainland.duplicated(
            ["entity_ID","work_species"]).any():
        raise ValueError("Original GIFT mainland native source corrupt")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assignment):
        raise ValueError("Original GIFT source assignment missing")
    names=set(traits["flower_primary_color"])&set(traits["floral_form"])
    target=sorted(i for i in full if len(full[i]&names)>=MIN_DIRECT)
    if len(target)!=2324:
        raise ValueError("Original colour/form 2324-island pair eligibility drift")
    provider={}
    rows=[]
    universe=set().union(*full.values())
    GIFT=set(mainland.work_species)
    for p,data in raw_corpora.items():
        states,details=nhm_external_categories(data)
        provider[p]=states
        shared=set(states["primary"])&set(states["form"])
        qualified=(shared&GIFT)-universe
        rows.append(dict(source_corpus=p,
                         n_original_text_descriptions=details["n_external_rows"],
                         n_original_with_primary_flower_colour=len(states["primary"]),
                         n_original_with_supported_flower_shape=len(states["form"]),
                         n_exact_same_source_species_primary_and_shape=len(shared),
                         n_exact_same_source_primary_and_shape_GIFT_mainland_source_ONLY=len(qualified),
                         n_source_only_actual_flower_form_genera=len({
                             x.split()[0] for x in qualified})))
    groups,audit=cross_corpus_joint(provider)
    counts=[]
    for tier,lab in groups.items():
        info=source_provenance_opportunity(
            set(lab),full,target,names,mainland,assignment)
        counts.append(dict(external_joint_evidence_tier=tier,
                           n_original_nhm_names_with_documented_joint_colour_form=len(lab),
                           **info))
    return dict(
        status="NHM_EXTERNAL_SOURCE_ONLY_SAME_SPECIES_JOINT_COLOR_FORM_SUPPORT",
        n_original_GIFT_botanical_island_floras=4505,
        n_original_flora_with_5plus_direct_colour_and_form=2324,
        n_original_direct_colour_form_species=len(names),
        n_original_independent_GIFT_native_mainland_source_regions=mainland.entity_ID.nunique(),
        n_original_separate_NHM_source_corpora=len(provider),
        original_each_NHM_text_source=sorted(rows,key=lambda x:x["source_corpus"]),
        source_only_external_colour_form_evidence_hierarchy=counts,
        cross_corpus_colour_form_quality_audit=audit,
        claims=dict(
            main_4505_botanical_island_universe_not_shrunk_to_external_morphology=True,
            original_direct_two_trait_mask_does_not_require_symmetry=True,
            no_external_trait_names_present_in_ANY_4505_island_botanical_list=True,
            both_traits_on_same_actual_accepted_source_species_no_imputation=True,
            separate_corpora_both_traits_more_stringent_than_across_provider_pair=True,
            gipt_native_mainland_name_presence_not_ancestral_source=True,
            external_botanical_text_not_independent_mainland_field_population=True,
            missing_text_flower_shape_not_monophorphism_or_absence=True,
            no_adaptive_flower_colour_form_pollen_fitness_or_pollinator_shift=True,
            no_private_plant_island_corpus_taxon_or_source_identifiers_exported=True,
        ))


def _sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):h.update(part)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for x in ("original-root","nhm-directory","output"):
        p.add_argument("--"+x,required=True,type=Path)
    args=p.parse_args()
    original={
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "gift":args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":args.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    sha={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "gift":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for k,path in original.items():
        if not path.is_file() or _sha(path)!=sha[k]:
            raise ValueError("Original GIFT source checksum changed "+k)
    frames={}
    for source,filename in PROVIDERS.items():
        path=args.nhm_directory/filename
        if not path.is_file() or nhm_git_blob_sha(path)!=BLOB_PINS[source]:
            raise ValueError("Original named NHM provider blob changed "+source)
        frames[source]=read_provider_file(path,source)
    traits=strict_traits(pd.read_csv(original["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    gift=pd.read_csv(original["gift"],dtype=str,usecols=[
        "entity_ID","work_species"])
    assign=pd.read_csv(original["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    output=analyze(traits,flora,gift,assign,frames)
    if output["n_original_GIFT_botanical_island_floras"]!=4505:
        raise AssertionError("Source original population changed")
    if output["n_original_flora_with_5plus_direct_colour_and_form"]!=2324:
        raise AssertionError("Full original paired endpoint changed")
    body=json.dumps(output,indent=2,sort_keys=True,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"work_species"',
                    '"entity_ID"','"source_signature"','"spatial_block"',
                    '"island_latitude"','"island_longitude"'):
        if private in body:
            raise ValueError("Original GIFT botanicals leaked into public evidence")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(status=output["status"],
                           n_paired_trait_eligible_original_islands=2324,
                           n_external_trait_evidence_tiers=len(counts) if False else len(
                               output["source_only_external_colour_form_evidence_hierarchy"]),
                           no_causal_flower_function_claim=True)))


if __name__=="__main__":main()
