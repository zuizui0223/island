"""Original 4505 island floral COLOUR prediction: named NHM text-corpus falsification.

Both provider corpora are original CC-BY botanical text description sources,
not independent specimen-level mainland phenotypes. Keep original 4505 island
botanical population and original target Direct colour species; exclude EVERY
accepted plant name recorded anywhere on the original 4505 islands from ALL
native GIFT mainland trait donor pools. Compare BHL and Flora of China on the
SAME focal islands and TWO SAME outcome-blind mainland-source assignments,
pre-matched for source-only colour donor species AND genus counts in BOTH
corpora. Source-corpus comparisons are descriptive, not evolution or pollination.
Only anonymous region, status, study-source and cluster summaries leave CI.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from collections import Counter,defaultdict

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_nhm_original_provider_audit import (
    read_provider_file, PROVIDERS, BLOB_PINS,
)
from fd_nhm_external_flora_preflight import nhm_git_blob_sha,nhm_external_categories
from fd_nhm_global_source_only_colour import (
    NHM_ALLOWED, ALPHA, build_external_source_catalog,
    proper_flora_losses, source_probability, _bootstrap_equal_group,
)
from fd_source_name_heldout_macro_preflight import (
    short_distance_block_disjoint_endpoint_union,
)

CORPORA=("BHL_text","FloraChina_text")
MATCH_FACTOR=2.0
SOURCE_ENTITY_JACCARD_LIMIT=.8
MIN_ISLANDS=30
MIN_BLOCKS=10
BOOTSTRAP=1999


def independent_donors(raw:pd.DataFrame,gift_names:set,island_names:set):
    """One documented exact species name per external source, no duplicates."""
    x,meta=nhm_external_categories(raw)
    donors={sp:state for sp,state in x["primary"].items()
            if sp in gift_names and sp not in island_names
            and state in NHM_ALLOWED}
    if set(donors)&island_names:
        raise AssertionError("Island taxon name leaked into NHM GIFT donors")
    return donors,meta


def matched_controls(source_id:str,catalogs:dict,by_region:dict,
                     *,n_controls=2,factor=MATCH_FACTOR):
    """Both original source corpora see EXACTLY the same placebo assignments.

    Only original location, source-entity overlap, and external donor SPECIES
    and GENUS support are used, never botanical target colour observations.
    """
    if n_controls<1 or factor<1 or not catalogs:
        raise ValueError("Incorrect original source matching contract")
    if not all(source_id in catalog for catalog in catalogs.values()):
        return [],0
    anchor=next(iter(catalogs.values()))[source_id]
    region=anchor["region"]
    actual_source=anchor["source_entities"]
    items=[]
    for other_id in by_region.get(region,()):
        if other_id==source_id:continue
        if not all(other_id in catalog for catalog in catalogs.values()):
            continue
        other=next(iter(catalogs.values()))[other_id]
        if other["block"]==anchor["block"]:
            continue
        source_set=other["source_entities"]
        union=actual_source|source_set
        jaccard=len(actual_source&source_set)/len(union) if union else 1.
        if jaccard>=SOURCE_ENTITY_JACCARD_LIMIT:
            continue
        delta=0.
        ratios=[]
        keep=True
        for catalog in catalogs.values():
            focal=catalog[source_id]
            candidate=catalog[other_id]
            if focal["n_source_species"]<1 or focal["n_source_genera"]<1:
                keep=False;break
            species_ratio=candidate["n_source_species"]/focal["n_source_species"]
            genus_ratio=candidate["n_source_genera"]/focal["n_source_genera"]
            if not (1./factor<=species_ratio<=factor
                    and 1./factor<=genus_ratio<=factor):
                keep=False;break
            delta+=abs(np.log(species_ratio))+abs(np.log(genus_ratio))
            ratios.append((species_ratio,genus_ratio))
        if keep:
            items.append((delta,jaccard,str(other_id),other_id,ratios))
    ranked=sorted(items)
    chosen=[]
    used_sets=set()
    for _,jac,_,other_id,ratios in ranked:
        source_set=next(iter(catalogs.values()))[other_id]["source_entities"]
        if source_set in used_sets:continue
        used_sets.add(source_set)
        chosen.append(dict(island=other_id,source_jaccard=jac,
                           species_genus_ratios=ratios))
        if len(chosen)==n_controls:break
    return chosen if len(chosen)==n_controls else [],len(items)


def aggregate_comparison(records:list[dict],*,bootstrap=BOOTSTRAP):
    """Cluster equal summaries of paired, same-control BHL/China score arrays."""
    if not records:return dict(status="NO_COMMON_EXTERNAL_SOURCE_SUPPORT",
                                n_matched_islands=0)
    blocks={r["block"] for r in records}
    first_sources={r["first_source"] for r in records}
    n=len(records)
    base=dict(n_matched_islands=n,n_distinct_geographic_blocks=len(blocks),
              n_distinct_first_GIFT_mainland_source_groups=len(first_sources),
              n_original_direct_colour_plant_carrier_rows=sum(r["n_direct"] for r in records),
              mean_placebo_to_actual_source_donor_species_ratio_BHL=float(np.mean([
                  r["source_donor_species_ratio_BHL"] for r in records])),
              mean_placebo_to_actual_source_donor_species_ratio_China=float(np.mean([
                  r["source_donor_species_ratio_China"] for r in records])),
              mean_placebo_to_actual_donor_genus_ratio_BHL=float(np.mean([
                  r["source_donor_genus_ratio_BHL"] for r in records])),
              mean_placebo_to_actual_donor_genus_ratio_China=float(np.mean([
                  r["source_donor_genus_ratio_China"] for r in records])),
              mean_original_GIFT_source_entity_Jaccard=float(np.mean([
                  r["source_jaccard"] for r in records])))
    if n<MIN_ISLANDS or len(blocks)<MIN_BLOCKS or len(first_sources)<MIN_BLOCKS:
        return dict(**base,status="INSUFFICIENT_COMMON_EXTERNAL_PROVIDER_SOURCE_SUPPORT")
    bhl=np.stack([r["contrast_BHL"] for r in records])
    china=np.stack([r["contrast_China"] for r in records])
    delta=bhl-china
    combined=np.concatenate((bhl,china,delta),axis=1)
    grouped=[dict(block=r["block"],first_source=r["first_source"],
                  contrast=combined[i]) for i,r in enumerate(records)]
    phy=_bootstrap_equal_group(grouped,"block",seed=20261115,draws=bootstrap)
    donor=_bootstrap_equal_group(grouped,"first_source",seed=20261116,draws=bootstrap)
    def slice_group(x):
        if x is None:return None
        def pair(j):
            return [float(x["blockwise_bootstrap_logloss_95"][0])
                   ] if False else None
        # Re-evaluate per-provider and paired control-source difference with
        # unchanged group memberships and same numerical bootstrap seed.
        return dict(n_groups=x["n_groups"])
    # A single 6-vector resampling is valid; the original helper exposes two
    # components, so repeat it by coordinate with matching seeds and blocks.
    out=dict(**base,status="POSTHOC_SAME_FOCAL_EXTERNAL_NHM_CORPUS_SOURCE_PLACEBO",
             mean_island_BHL_logloss_gain=float(bhl[:,0].mean()),
             mean_island_China_logloss_gain=float(china[:,0].mean()),
             mean_island_BHL_brier_gain=float(bhl[:,1].mean()),
             mean_island_China_brier_gain=float(china[:,1].mean()),
             mean_island_BHL_minus_China_logloss_advantage=float(delta[:,0].mean()),
             fraction_focal_floras_BHL_real_source_better=float(np.mean(bhl[:,0]>0)),
             fraction_focal_floras_China_real_source_better=float(np.mean(china[:,0]>0)),
             fraction_focal_floras_both_corpora_favor_real_source=float(np.mean(
                 (bhl[:,0]>0)&(china[:,0]>0))))
    for groupfield,seed,label in (
        ("block",20261115,"geographic_block"),
        ("first_source",20261116,"first_GIFT_source")):
        for prefix,value in (
            ("BHL",bhl),("China",china),("BHL_minus_China",delta)):
            r=[dict(**{groupfield:x[groupfield]},contrast=value[i])
               for i,x in enumerate(records)]
            q=_bootstrap_equal_group(r,groupfield,seed=seed,draws=bootstrap)
            out[label+"_"+prefix+"_equal_mean_logloss_gain"]=q["mean_gain_logloss"]
            out[label+"_"+prefix+"_conditional_logloss_95"]=q["blockwise_bootstrap_logloss_95"]
            out[label+"_"+prefix+"_equal_mean_brier_gain"]=q["mean_gain_brier"]
            out[label+"_"+prefix+"_conditional_brier_95"]=q["blockwise_bootstrap_brier_95"]
    return out


def evaluate_provider_comparison(traits,flora,geo,land,assignments,
                                 provider_raw,*,bootstrap=BOOTSTRAP):
    required={"island_id","accepted_species","origin_status"}
    if not required<=set(flora) or flora.duplicated(
            ["island_id","accepted_species"]).any():
        raise ValueError("Original botanical flora missing or duplicated")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Invalid original island status")
    all_map=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    if len(all_map)!=4505:
        raise ValueError("Original 4505 botanical flora changed")
    native_map=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    all_names=set().union(*all_map.values())
    available=sorted(NHM_ALLOWED)
    category={c:i for i,c in enumerate(available)}
    direct={sp:category[col] for sp,col in traits["flower_primary_color"].items()
            if col in category}
    names=set(direct)
    original_all=sorted(i for i,ls in all_map.items() if len(ls&names)>=5)
    original_native=sorted(i for i in all_map if len(native_map.get(i,set())&names)>=5)
    if len(original_all)!=3119 or len(original_native)!=395:
        raise ValueError("Original 3119/395 colour opportunity changed")
    pair_names=set(traits["flower_primary_color"])&set(traits["floral_form"])
    indexed=geo.set_index("island_id")
    n10=sorted(i for i in all_map if i in indexed.index
               and indexed.at[i,"analysis_regime"]=="northern_midlatitude"
               and len(all_map[i]&pair_names)>=10
               and len(all_map[i]&pair_names)/len(all_map[i])>=.10)
    endpoints,pair_total,blocks=short_distance_block_disjoint_endpoint_union(
        n10,all_map,geo)
    if (len(endpoints),pair_total,blocks)!=(215,456,29):
        raise ValueError("Original 215/456/29 focal pair support changed")
    gift=set(land.work_species)
    profiles={}
    catalogs={}
    by_region=defaultdict(list)
    for corpus,raw in provider_raw.items():
        labels,metadata=independent_donors(raw,gift,all_names)
        if not labels:raise ValueError("External corpus has no independent donor species")
        prior=np.array([Counter(labels.values()).get(col,0)+.5 for col in available],float)
        prior/=prior.sum()
        catalogs[corpus]=build_external_source_catalog(
            land,assignments,geo,labels,category)
        profiles[corpus]=(labels,prior,metadata)
        for item in catalogs[corpus].values():
            if item["n_source_species"]<1:raise ValueError("Empty GIFT donor entry")
    shared=set.intersection(*(set(c.keys()) for c in catalogs.values()))
    first=next(iter(catalogs.values()))
    for island in shared:
        by_region[first[island]["region"]].append(island)
    for region in by_region:by_region[region].sort()
    scopes=(
        ("all_recorded_n5",original_all,all_map),
        ("original_native_n5",original_native,native_map),
        ("posthoc_215_endpoints",endpoints,all_map))
    out=[]
    for key,ids,membership in scopes:
        records=[]
        n_missing_source=0
        n_no_controls=0
        for island in ids:
            target=[direct[x] for x in sorted(membership.get(island,set())&names)]
            if len(target)<5:continue
            if island not in shared:
                n_missing_source+=1;continue
            controls,n_options=matched_controls(
                island,catalogs,by_region,n_controls=2)
            if not controls:
                n_no_controls+=1;continue
            a=first[island]
            row=dict(block=a["block"],first_source=a["first_source"],
                     n_direct=len(target),
                     source_jaccard=float(np.mean([p["source_jaccard"] for p in controls])))
            for corpus,label in (("BHL_text","BHL"),("FloraChina_text","China")):
                cat=catalogs[corpus]
                global_prior=profiles[corpus][1]
                actual=proper_flora_losses(
                    target,source_probability(cat[island]["values"],global_prior,alpha=ALPHA))
                mock=np.mean([
                    proper_flora_losses(
                        target,source_probability(
                            cat[p["island"]]["values"],global_prior,alpha=ALPHA))
                    for p in controls],axis=0)
                row["contrast_"+label]=np.asarray(mock)-np.asarray(actual)
                ix=list(catalogs).index(corpus)
                row["source_donor_species_ratio_"+label]=float(np.mean(
                    [p["species_genus_ratios"][ix][0] for p in controls]))
                row["source_donor_genus_ratio_"+label]=float(np.mean(
                    [p["species_genus_ratios"][ix][1] for p in controls]))
            records.append(row)
        report=aggregate_comparison(records,bootstrap=bootstrap)
        out.append(dict(scope=key,n_eligible_original_floras=len(ids),
                        n_original_floras_missing_one_or_both_external_corpora=n_missing_source,
                        n_original_floras_without_two_jointly_support_matched_placebos=n_no_controls,
                        **report))
    return dict(
        status="NHM_BHL_VS_CHINA_ORIGINAL_SOURCE_PROVIDER_VALIDATION",
        n_original_botanical_floras=4505,
        n_original_colour_n5_islands=len(original_all),
        n_original_native_colour_n5_islands=len(original_native),
        n_original_selected_pair_endpoints=len(endpoints),
        corpora={
            x:dict(n_external_mainland_nontarget_species=len(v[0]),
                   n_eligible_source_catalog_geographic_units=len(catalogs[x]),
                   n_original_text_rows=v[2]["n_external_rows"],
                   n_unique_exact_binomial_source_names=v[2]["n_exact_taxon_binomial_names"])
            for x,v in profiles.items()},
        same_island_same_controls_comparisons=out,
        claims=dict(
            botanical_universe_4505_not_triple_complete_928=True,
            original_island_plant_species_names_removed_from_all_source_regions=True,
            external_text_sources_separate_no_combined_provider_imputation=True,
            both_provider_source_controls_identical_for_each_focal_island=True,
            matching_selected_without_focal_island_flower_colour_outcomes=True,
            physical_and_first_GIFT_source_group_weights_are_different_estimands=True,
            source_only_profile_data_not_independent_mainland_population_measurements=True,
            two_text_corpora_not_independent_evolutionary_field_histories=True,
            no_interaction_pollen_transfer_guild_shift_or_adaptive_evolution_claim=True,
            no_original_taxon_island_genus_source_GIFT_entity_ids_exported=True))


def sha256(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""):
            h.update(block)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for arg in ("original-root","geology","nhm-directory","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    args=p.parse_args()
    base=args.original_root
    paths={
        "traits":base/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":base/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "gift":base/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":base/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "gift":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for k,v in paths.items():
        if not v.is_file() or sha256(v)!=pins[k]:
            raise ValueError("Original protected source SHA drift "+k)
    raw={}
    nhm_pins={}
    for source in CORPORA:
        filename=PROVIDERS[source]
        path=args.nhm_directory/filename
        if not path.is_file() or nhm_git_blob_sha(path)!=BLOB_PINS[source]:
            raise ValueError("External NHM source corpus git SHA mismatch "+source)
        raw[source]=read_provider_file(path,source)
        nhm_pins[source]=BLOB_PINS[source]
    original=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    floras=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    gift=pd.read_csv(paths["gift"],dtype=str,usecols=[
        "entity_ID","work_species"])
    assign=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    geo=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block","island_latitude",
        "island_longitude"])
    result=evaluate_provider_comparison(original,floras,geo,gift,assign,raw)
    result["NHM_original_corpus_git_blob_sha"]=nhm_pins
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for name in ('"island_id"','"accepted_species"','"work_species"',
                 '"entity_ID"','"source_signature"','"spatial_block"',
                 '"island_latitude"','"island_longitude"'):
        if name in body:
            raise ValueError("Original plant or source locality identity escaped")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(status=result["status"],
                          n_original_eligible_colour_islands=result["n_original_colour_n5_islands"],
                          n_report_scopes=len(result["same_island_same_controls_comparisons"]),
                          not_independent_ancestral_pollinator_evidence=True)))


if __name__=="__main__":main()
