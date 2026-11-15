"""Original 4505 botanical GIFT candidate-source geography sensitivity.

Compare independently archived GIFT:
 - geo_k5: five nearest native GIFT mainland source regions by geographical km
 - geo50_climate10: ten climate-nearest source regions from geographically
   prefiltered nearest 50 source regions, all with actual geographic AND
   source-climate distances. The latter lacks source sets for 20 of the
   original 4505 accepted botanical island floras.

Both modes are audited by the SAME finite GIFT source-name
listing-breadth×actual-kilometre genus null, but source plant candidate
sets are DIFFERENT. This is not a causal effect of climate on establishment,
seed movement, adaptation or pollination. Pair by ORIGINAL island identity
for an explicit source-model sensitivity, and preserve all 4505 floras.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups,source_sectors
from fd_direct_source_breadth_matched import REGIONS
from fd_original_source_trait_free_genus_null import original_trait_island_eligibility
from fd_direct_source_distance_conditioned_genus_null import (
    closest_km_bin,exact_taxon_stratum_genus_reference,
    exact_same_names_breadth_baseline,group_equal
)
from fd_direct_source_breadth_matched import listing_breadth_bin

MODES=("geo_k5","geo50_climate10")
SCOPES=("strict_native_nonendemic_source_listed",
        "all_original_island_species_source_listed",
        "original_introduced_source_listed")
MIN_N=5
MIN_S=10
MIN_MATCH=5
MIN_ISLANDS=30
MIN_BLOCKS=10
BOOT=999
PAIRED=(
    "source_mode_climate_minus_geo5_original_genus_overlap_deficit",
    "source_mode_climate_minus_geo5_original_genus_concentration_residual"
)


def source_species_min_km(source_records:list[tuple[str,float]],
                          source_membership:dict[str,set[str]]):
    """Most proximal actual geographically measured native source REGION listing.

    Source records may be ranked by climate, not by km. Always choose
    true minimum source_distance_km among the SELECTED original entities.
    """
    names={}
    if not source_records:return names
    for entity,distance in source_records:
        if not isinstance(distance,(int,float,np.integer,np.floating)) or (
            not np.isfinite(distance) or distance<0):
            raise ValueError("GIFT original source regional distance missing")
        for sp in source_membership.get(entity,set()):
            if sp not in names or distance<names[sp]:
                names[sp]=float(distance)
    return names


def paired_mode_group_summary(records:list[dict],cluster:str,
                              *,n_boot=BOOT,seed=20261115):
    if not records:
        return dict(status="NO_BOTH_ORIGINAL_GIFT_SOURCE_MODES",n_islands=0,n_groups=0)
    df=pd.DataFrame(records).dropna(subset=[cluster,*PAIRED])
    n=len(df)
    G=df[cluster].nunique()
    result=dict(n_islands=n,n_groups=G)
    if n<MIN_ISLANDS or G<MIN_BLOCKS:
        result["status"]="INSUFFICIENT_SAME_ISLAND_SOURCE_PAIR_GROUPS"
        return result
    a=df.groupby(cluster,sort=True)[list(PAIRED)].mean().to_numpy(float)
    avg=a.mean(0)
    rng=np.random.default_rng(seed)
    boot=a[rng.integers(G,size=(n_boot,G))].mean(1)
    for j,name in enumerate(PAIRED):
        result[name+"_source_group_equal_mean"]=float(avg[j])
        result[name+"_conditional_group_95"]=[
            float(np.quantile(boot[:,j],.025)),
            float(np.quantile(boot[:,j],.975))]
        loo=[float(a[np.arange(G)!=i,j].mean()) for i in range(G)]
        result[name+"_leave_one_source_group_min"]=float(min(loo))
        result[name+"_leave_one_source_group_max"]=float(max(loo))
    result["status"]="ORIGINAL_SAME_ISLAND_GEO_VS_CLIMATE_SOURCE_DEFINITION_SENSITIVITY"
    result["mean_source_mode_union_candidate_name_jaccard"]=float(
        df.source_mode_original_gift_candidate_jaccard.mean())
    result["mean_original_geo_k5_mainland_candidate_species"]=float(
        df.geo5_source_candidate_species_count.mean())
    result["mean_original_climate10_mainland_candidate_species"]=float(
        df.climate10_source_candidate_species_count.mean())
    result["n_source_resamples_not_independent_colonizations"]=n_boot
    return result


def analyze(original,geo,assign,mainland,states,*,n_boot=BOOT):
    F=status_groups(original)
    if len(F)!=4505:raise ValueError("Original GIFT source botanical population changed")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or (
        geo.island_id.duplicated().any()):
        raise ValueError("Original independent GSHHG geography invalid")
    g=geo.set_index("island_id").reindex(sorted(F))
    if g.analysis_regime.notna().sum()!=4504:
        raise ValueError("Original one geography unresolved botanical flora changed")
    region={i:(str(g.at[i,"analysis_regime"])
               if pd.notna(g.at[i,"analysis_regime"])
               else "geography_unresolved")
            for i in F}
    reg_counts=Counter(region.values())
    if sum(reg_counts.values())!=4505 or set(reg_counts)!=set(REGIONS):
        raise ValueError("Original GIFT macro geography not complete")
    required={"island_id","source_mode","source_rank",
              "entity_ID","source_distance_km","climate_distance"}
    if not required<=set(assign):
        raise ValueError("Original GIFT geography/climate source assignment columns absent")
    A=assign.loc[assign.source_mode.isin(MODES)].copy()
    A.source_rank=pd.to_numeric(A.source_rank,errors="raise").astype(int)
    A.source_distance_km=pd.to_numeric(A.source_distance_km,errors="raise")
    A.climate_distance=pd.to_numeric(A.climate_distance,errors="coerce")
    if not np.isfinite(A.source_distance_km).all() or (A.source_distance_km<0).any():
        raise ValueError("Original source geographic distances invalid")
    if A.duplicated(["island_id","source_mode","source_rank"]).any():
        raise ValueError("Original candidate source mode ranked duplicates")
    if not A.loc[A.source_mode.eq("geo50_climate10"),"climate_distance"].notna().all():
        raise ValueError("Original climate-ranked GIFT sources have missing climate distances")
    if (not A.loc[A.source_mode.eq("geo_k5"),"source_rank"].between(1,5).all() or
        not A.loc[A.source_mode.eq("geo50_climate10"),"source_rank"].between(1,10).all()):
        raise ValueError("Original GIFT source rank 5/10 mode drift")
    expected={"geo_k5":41325,"geo50_climate10":82190}
    for mode,n in expected.items():
        if int(A.source_mode.eq(mode).sum())!=n:
            raise ValueError("Original geographic/climate source assignment opportunity changed")
    membership=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    breadth=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    source_modes=defaultdict(dict)
    for (id,mode),group in A.groupby(["island_id","source_mode"],sort=True):
        source_modes[id][mode]=list(zip(
            group.sort_values("source_rank").entity_ID.astype(str),
            group.sort_values("source_rank").source_distance_km.astype(float)))
    climate_original_covered=sum("geo50_climate10" in source_modes[id] for id in F)
    if climate_original_covered!=4485:
        raise ValueError("Original 4485/4505 climate candidate source flora coverage changed")
    sectors=source_sectors(states)
    if len(sectors)!=7:raise ValueError("Original source seven Direct sector mask definitions changed")
    all_panels=[]
    geo_support=[]
    for reg in REGIONS:
        ids=sorted(i for i in F if region[i]==reg)
        buffers=defaultdict(dict)
        mode_counts=Counter()
        direct_masks=defaultdict(set)
        for id in ids:
            item=F[id]
            potential_by_mode={}
            for mode in MODES:
                neighbor=source_modes[id].get(mode,[])
                if not neighbor:continue
                min_dist=source_species_min_km(neighbor,membership)
                all_names=set(min_dist)
                potential_by_mode[mode]=all_names
                outside=all_names-item["full"]
                # Original GEO mode defines pre-analytical eligible ISLAND
                # identity for each Direct trait, not the climate result.
                if mode=="geo_k5":
                    strict=all_names&item["native_nonendemic"]
                    for name,trait in sectors.items():
                        if original_trait_island_eligibility(
                            strict,outside,breadth,trait):
                            direct_masks[name].add(id)
                scope_names={
                    SCOPES[0]:all_names&item["native_nonendemic"],
                    SCOPES[1]:all_names&item["full"],
                    SCOPES[2]:all_names&item["introduced"],
                }
                mode_counts[mode+"_present_source_flora"]+=1
                for scope,ingroup in scope_names.items():
                    if len(ingroup)<MIN_N or len(outside)<MIN_S:continue
                    keys={sp:(listing_breadth_bin(int(breadth[sp])),
                              closest_km_bin(min_dist[sp]))
                          for sp in ingroup|outside}
                    test=exact_taxon_stratum_genus_reference(
                        ingroup,outside,keys)
                    if test is None:continue
                    q,names,_=test
                    if q["original_matched_effective_source_species"]<MIN_MATCH:
                        continue
                    q=exact_same_names_breadth_baseline(
                        q,*test[1:],breadth)
                    q.update(
                        physical_block=(str(g.at[id,"spatial_block"])
                            if pd.notna(g.at[id,"spatial_block"]) else None),
                        source_region_cluster=(
                            source_modes[id]["geo_k5"][0][0]
                            if source_modes[id].get("geo_k5") else None),
                        original_n_selected_source_regions=len(neighbor),
                        original_source_candidate_species_count=len(outside),
                        original_source_mode=mode)
                    buffers[(scope,mode)][id]=q
            if all(m in potential_by_mode for m in MODES):
                a_set=potential_by_mode["geo_k5"]
                b_set=potential_by_mode["geo50_climate10"]
                buffers[("private_source_jaccard","both")][id]=(
                    len(a_set&b_set)/len(a_set|b_set)
                    if a_set|b_set else 0.)
        for scope in SCOPES:
            scope_cohorts=(["all_original_source_taxa"]+
                ["same_geo5_islands_as_Direct_"+name for name in sectors]
                if scope==SCOPES[0] else
                ["all_original_source_taxa"])
            for label in scope_cohorts:
                qdict={}
                valid_names=set(ids) if label=="all_original_source_taxa" else (
                    direct_masks[label.replace("same_geo5_islands_as_Direct_","")])
                for mode in MODES:
                    data=[
                        v for key,v in buffers[(scope,mode)].items()
                        if key in valid_names]
                    row=dict(region=reg,scope=scope,cohort=label,source_mode=mode,
                             n_original_botanical_region_floras=len(ids),
                             n_original_current_mode_qualified_floras=len(data),
                             n_original_existing_source_mode_regions=len(
                                 set(key for key in ids if mode in source_modes[key])))
                    if data:
                        row["island_equal_genus_deficit_after_actual_source_distance"]=float(
                            np.mean([x["observed_minus_breadth_distance_expected_shared_genera"]
                                     for x in data]))
                    for kind,cname in (
                        ("physical_block","physical_10degree_block"),
                        ("source_region_cluster","first_original_GIFT_geo5_source")):
                        row[cname]=group_equal(
                            data,kind,boots=n_boot,
                            seed=int.from_bytes(hashlib.sha256(
                                (reg+scope+label+mode+kind).encode()).digest()[:8],"big"))
                    qdict[mode]=row
                    all_panels.append(row)
                common=(set(buffers[(scope,"geo_k5")])&
                        set(buffers[(scope,"geo50_climate10")])&
                        valid_names)
                paired=[]
                for id in sorted(common):
                    a=buffers[(scope,"geo_k5")][id]
                    b=buffers[(scope,"geo50_climate10")][id]
                    jaccard=buffers[("private_source_jaccard","both")].get(id,0.)
                    paired.append(dict(
                        physical_block=a["physical_block"],
                        source_region_cluster=a["source_region_cluster"],
                        source_mode_climate_minus_geo5_original_genus_overlap_deficit=(
                            b["observed_minus_breadth_distance_expected_shared_genera"]-
                            a["observed_minus_breadth_distance_expected_shared_genera"]),
                        source_mode_climate_minus_geo5_original_genus_concentration_residual=(
                            b["observed_minus_distance_expected_genus_concentration"]-
                            a["observed_minus_distance_expected_genus_concentration"]),
                        source_mode_original_gift_candidate_jaccard=jaccard,
                        geo5_source_candidate_species_count=a["original_source_candidate_species_count"],
                        climate10_source_candidate_species_count=b["original_source_candidate_species_count"],
                    ))
                p=dict(region=reg,scope=scope,cohort=label,
                       n_source_mode_geo5_qualified=qdict["geo_k5"]["n_original_current_mode_qualified_floras"],
                       n_source_mode_climate10_qualified=qdict["geo50_climate10"]["n_original_current_mode_qualified_floras"],
                       n_same_original_islands_qualified_under_both_source_modes=len(common),
                       original_source_modes_use_different_species_candidate_pools=True)
                for kind,key in (
                    ("physical_block","physical_10degree_block"),
                    ("source_region_cluster","first_original_GIFT_geo5_source")):
                    p[key]=paired_mode_group_summary(
                        paired,kind,n_boot=n_boot,
                        seed=int.from_bytes(hashlib.sha256(
                            (reg+scope+label+"PAIRED"+kind).encode()).digest()[:8],"big"))
                geo_support.append(p)
    return dict(
        status="ORIGINAL_4505_GEOGRAPHIC5_VS_CLIMATE10_GIFT_SOURCE_POOL_SENSITIVITY",
        n_original_botanical_GIFT_island_floras=4505,
        n_original_climate10_source_matched_botanical_island_floras=climate_original_covered,
        n_original_independent_gift_source_assignment_rows_by_mode={
            mode:expected[mode] for mode in MODES},
        n_original_Direct_trait_only_island_eligibility_masks=len(sectors),
        n_original_mainland_GIFT_source_entities=mainland.entity_ID.nunique(),
        original_complete_GSHHG_plant_region_population=[
            dict(region=reg,n_original_botanical_floras=int(reg_counts[reg]))
            for reg in REGIONS],
        original_alternate_mainland_source_genus_null_results=all_panels,
        original_same_island_geographic_vs_climate_source_mode_differences=geo_support,
        claims=dict(
            true_original_source_km_not_observed_seed_distance=True,
            source_pool_climate10_gift_regions_different_from_geographic5=True,
            original_climate_distances_recorded_for_selected_climate10_regions=True,
            no_false_climate_distance_for_geographic5_regions=True,
            each_source_mode_genus_null_conditions_on_real_GIFT_breadth_and_kilometers=True,
            paired_source_mode_comparison_uses_same_original_island_ID_but_different_candidate_taxa=True,
            same_original_flower_Direct_eligibility_mask_is_geo5_defined_not_result_selected=True,
            original_full_4505_botanical_population_retained=True,
            no_actual_arrival_success_plant_adaptation_or_pollinator_service_inferred=True,
            no_original_private_accepted_taxon_island_GIFT_entity_genus_coordinates_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for x in ("original-root","geology","output"):
        p.add_argument("--"+x,type=Path,required=True)
    args=p.parse_args()
    paths={
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "trait":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "mainland":args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":args.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "trait":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for key,name in paths.items():
        if sha(name)!=pins[key]:
            raise ValueError("Original GIFT botanical source SHA changed: "+key)
    original=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    trait=strict_traits(pd.read_csv(paths["trait"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    assign=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","entity_ID","source_mode","source_rank","source_distance_km","climate_distance"])
    mainland=pd.read_csv(paths["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    geo=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    if geo.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original 8264 GSHHG and 748 mainland source entity records drift")
    z=analyze(original,geo,assign,mainland,trait)
    assert z["n_original_botanical_GIFT_island_floras"]==4505
    blob=json.dumps(z,indent=2,sort_keys=True,allow_nan=False)
    for k in ('"island_id"','accepted_species"','entity_ID"','"work_species"',
              '"spatial_block"','"source_signature"','"island_latitude"','"island_longitude"'):
        if k in blob:raise ValueError("Original protected GIFT island/source taxonomy leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(blob+"\n")
    print(json.dumps(dict(status=z["status"],n_botanical_floras=4505,
       n_anonymous_source_mode_region_status_cohorts=len(z[
           "original_same_island_geographic_vs_climate_source_mode_differences"]),
       no_climate_caused_establishment_claim=True)))


if __name__=="__main__":
    main()
