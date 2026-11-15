"""Original 4505-island GIFT closest actual mainland-source kilometres conditioned genus null.

The ecological response is ACTUAL GIFT genus overlap and within-island genus
pair concentration on ALL accepted original source-listed plant names, no
flower trait label required. For each original species, source distance =
MIN independently archived source_distance_km among the island's five ranked
native mainland GIFT source REGION(S) genuinely listing this accepted taxon.
This is a geographic source proxy, NOT species occurrence edge, propagule
arrival or empirical seed dispersal kernel.

Exact sampling controls BOTH original accepted taxon SOURCE-region listing
breadth and candidate source distance (km) while preserving island-side
and source-only sample sizes and original source genus identity within
every stratum. Compare breadth-only and breadth×distance expectations on
the IDENTICAL original plant species cohort, not different selected names.

Native known, all-status and introduced source-listed floras are overlapping
documentary groups, not independent biological replications. Source-Direct
flower labels are ONLY used to name the exact previously phenotype-eligible
island masks, NEVER in genus response or conditional sampling null.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter, defaultdict
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import status_groups,source_sectors
from fd_direct_source_breadth_matched import BINS,REGIONS,listing_breadth_bin
from fd_direct_source_finite_sampling_null import (
    genus,prob_no_original_genus_in_native_or_candidate,
)
from fd_original_source_trait_free_genus_null import (
    exact_taxonomic_source_null,original_trait_island_eligibility
)

DISTANCE_BINS=("lt250km","250_750km","750_2000km","ge2000km")
SCOPES=("strict_native_nonendemic_source_listed",
        "all_original_island_species_source_listed",
        "original_introduced_source_listed")
MIN_NATIVE=5
MIN_CANDIDATE=10
MIN_MATCH=5
MIN_ISLANDS=30
MIN_BLOCKS=10
BOOT=999
METRICS=(
    "observed_minus_breadth_distance_expected_shared_genera",
    "observed_minus_breadth_only_expected_shared_genera_on_SAME_names",
    "additional_genus_sharing_expected_from_not_fixing_source_distance",
    "observed_minus_distance_expected_genus_concentration",
    "observed_minus_breadth_expected_genus_concentration_on_SAME_names"
)


def closest_km_bin(km:float)->str:
    if not isinstance(km,(int,float,np.integer,np.floating)) or (
            not np.isfinite(km) or km<0):
        raise ValueError("Original GIFT source distance must be finite nonnegative km")
    if km<250:return "lt250km"
    if km<750:return "250_750km"
    if km<2000:return "750_2000km"
    return "ge2000km"


def exact_taxon_stratum_genus_reference(focal:set[str],candidate:set[str],
                                       strata:dict[str,tuple[str,str]],
                                       *,min_effective:int=MIN_MATCH):
    """Exact independent finite name draw per fixed breadth-by-distance stratum.

    The returned two sets are private in-process botanical accepted names;
    callers must export ONLY the numerical anonymous result dictionary.
    """
    if not isinstance(min_effective,int) or min_effective<1:
        raise ValueError("Invalid source eligible effective species gate")
    if not focal or not candidate or focal&candidate:
        raise ValueError("Original focal/source accepted plants must be disjoint")
    if not focal|candidate<=set(strata):
        raise ValueError("Unverified real closest GIFT source distance for botanical name")
    groups=defaultdict(lambda:dict(focal=set(),candidate=set()))
    for sp in focal:groups[strata[sp]]["focal"].add(sp)
    for sp in candidate:groups[strata[sp]]["candidate"].add(sp)
    used=[]
    focal_used=set()
    candidate_used=set()
    for key,sets in groups.items():
        u,v=sets["focal"],sets["candidate"]
        if not u or not v:continue
        if key[0] not in BINS or key[1] not in DISTANCE_BINS:
            raise ValueError("Original geographic source bin unknown")
        focal_used.update(u)
        candidate_used.update(v)
        used.append((key,u,v))
    if not used:return None
    effective=sum(min(len(u),len(v)) for _,u,v in used)
    if effective<min_effective:return None
    actual_fgen={genus(s) for s in focal_used}
    actual_sgen={genus(s) for s in candidate_used}
    allg=actual_fgen|actual_sgen
    expected_absent_focal=defaultdict(lambda:1.)
    expected_absent_candidate=defaultdict(lambda:1.)
    expected_same_genus_within_pair=0.
    reference=[]
    n_focal=0
    for _,u,v in used:
        n,s=len(u),len(v)
        N=n+s
        n_focal+=n
        bygen=Counter(genus(x) for x in u|v)
        for g,m in bygen.items():
            nof,noc=prob_no_original_genus_in_native_or_candidate(
                int(N),int(n),int(m))
            expected_absent_focal[g]*=nof
            expected_absent_candidate[g]*=noc
        if n>1 and N>1:
            expected_same_genus_within_pair+=comb(n,2)*sum(
                m*(m-1) for m in bygen.values())/(N*(N-1))
        reference.append((n,N,bygen))
    between=0.
    for (n,N,counts),(m,M,other) in combinations(reference,2):
        between+=n*m*sum(v*other.get(g,0) for g,v in counts.items())/(N*M)
    expected_src_only=float(sum(expected_absent_focal[g] for g in allg))
    expected_focal_only=float(sum(expected_absent_candidate[g] for g in allg))
    expected_shared=float(len(allg)-expected_src_only-expected_focal_only)
    observed_shared=len(actual_fgen&actual_sgen)
    observed_src_only=len(actual_sgen-actual_fgen)
    observed_focal_only=len(actual_fgen-actual_sgen)
    if observed_shared+observed_src_only+observed_focal_only!=len(allg):
        raise AssertionError("Original real genus presence identity fails")
    if abs(expected_shared+expected_src_only+expected_focal_only-len(allg))>1e-8:
        raise AssertionError("Original exact finite source genus expected identity fails")
    counts=Counter(genus(x) for x in focal_used)
    observed_concentration=(
        sum(m*(m-1) for m in counts.values())/(n_focal*(n_focal-1))
        if n_focal>1 else None)
    expected_concentration=(
        (expected_same_genus_within_pair+between)/comb(n_focal,2)
        if n_focal>1 else None)
    numerical=dict(
        original_matched_effective_source_species=effective,
        n_original_joint_source_breadth_distance_strata=len(used),
        n_original_focal_source_names_in_joint_strata=len(focal_used),
        n_original_candidate_source_names_in_joint_strata=len(candidate_used),
        original_focal_source_taxon_retention_fraction=len(focal_used)/len(focal),
        original_candidate_source_taxon_retention_fraction=len(candidate_used)/len(candidate),
        original_actual_focal_shared_source_genera=observed_shared,
        original_actual_focal_only_source_genera=observed_focal_only,
        original_actual_candidate_only_source_genera=observed_src_only,
        expected_shared_genera_under_breadth_distance_exact_source_null=expected_shared,
        expected_focal_only_genera_under_breadth_distance_exact_source_null=expected_focal_only,
        expected_candidate_only_genera_under_breadth_distance_exact_source_null=expected_src_only,
        observed_minus_breadth_distance_expected_shared_genera=(
            observed_shared-expected_shared),
        original_observed_focal_genus_pair_concentration=float(observed_concentration),
        expected_focal_genus_pair_concentration_under_breadth_distance_exact_null=float(expected_concentration),
        observed_minus_distance_expected_genus_concentration=float(
            observed_concentration-expected_concentration),
    )
    # The exact same accepted species sets are returned for the caller to
    # recompute the previous original breadth-only null; never materialized
    # in the anonymous output artifact.
    return numerical,focal_used,candidate_used


def exact_same_names_breadth_baseline(
        strict_result:dict,used_focal:set[str],used_candidate:set[str],
        original_gift_breadth:dict[str,int]):
    benchmark=exact_taxonomic_source_null(
        used_focal,used_candidate,original_gift_breadth)
    if benchmark is None:
        raise AssertionError("Original distance-qualified source set lacks breadth reference")
    obs=int(strict_result["original_actual_focal_shared_source_genera"])
    if obs!=benchmark["observed_original_shared_genus_count"]:
        raise AssertionError("Compared different accepted flora names across source nulls")
    original_concentration=strict_result[
        "original_observed_focal_genus_pair_concentration"]
    if abs(original_concentration-
           benchmark["observed_original_island_genus_pair_concentration"])>1e-12:
        raise AssertionError("Original source richness differs across baseline and geography")
    strict_result.update(
        observed_minus_breadth_only_expected_shared_genera_on_SAME_names=float(
            obs-benchmark["exact_null_expected_original_shared_genus_count"]),
        original_breadth_only_expected_shared_genera_on_SAME_names=float(
            benchmark["exact_null_expected_original_shared_genus_count"]),
        additional_genus_sharing_expected_from_not_fixing_source_distance=float(
            benchmark["exact_null_expected_original_shared_genus_count"]-
            strict_result["expected_shared_genera_under_breadth_distance_exact_source_null"]),
        observed_minus_breadth_expected_genus_concentration_on_SAME_names=float(
            original_concentration-
            benchmark["exact_null_expected_original_island_genus_pair_concentration"]),
        original_breadth_only_expected_genus_concentration_on_SAME_names=float(
            benchmark["exact_null_expected_original_island_genus_pair_concentration"]),
        original_source_breadth_N_named_species_still_same_between_nulls=True,
    )
    return strict_result


def source_opportunity_distances(
        source_entity_distance:list[tuple[str,float]],
        original_gift_mainland_listings:dict[str,set[str]]):
    """Correct ORIGINAL source-listed species nearest candidate geography only."""
    near={}
    count_sources=Counter()
    prev=-float("inf")
    for rank,(entity,distance) in enumerate(source_entity_distance,start=1):
        if not np.isfinite(distance) or distance<0:
            raise ValueError("Unresolved native GIFT source distance")
        if distance+1e-7<prev:
            raise ValueError("Original ranked geo_k5 source distances not monotone")
        prev=distance
        for sp in original_gift_mainland_listings.get(entity,set()):
            # Original source entity membership, not inferred whole species range.
            if sp not in near:
                near[sp]=float(distance)
                count_sources[rank]+=1
    return near,count_sources


def group_equal(records:list[dict],group:str,*,boots=BOOT,seed=20261115):
    if not records:return dict(status="NO_ORIGINAL_BOTANICAL_DISTANCE_MATCHES",
                                n_original_island_floras=0,n_original_source_groups=0)
    df=pd.DataFrame(records).dropna(subset=[group,*METRICS])
    n=len(df);G=df[group].nunique()
    out=dict(n_original_island_floras=n,n_original_source_groups=G,
             n_original_islands_with_any_true_GIFT_source_distance_match=n)
    if n<MIN_ISLANDS or G<MIN_BLOCKS:
        return dict(**out,status="INSUFFICIENT_ORIGINAL_SOURCE_GEO_BLOCK_COVERAGE")
    scores=df.groupby(group,sort=True)[list(METRICS)].mean().to_numpy(float)
    observed=scores.mean(0)
    rng=np.random.default_rng(seed)
    boot=scores[rng.integers(G,size=(boots,G))].mean(1)
    for j,metric in enumerate(METRICS):
        out[metric+"_group_equal_mean"]=float(observed[j])
        out[metric+"_group_conditional_95"]=[
            float(np.quantile(boot[:,j],.025)),
            float(np.quantile(boot[:,j],.975))]
        loo=[float(scores[np.arange(G)!=i,j].mean()) for i in range(G)]
        out[metric+"_delete_one_group_min"]=float(min(loo))
        out[metric+"_delete_one_group_max"]=float(max(loo))
    out["status"]="EXACT_GIFT_DISTANCE_AND_LISTING_BREADTH_SOURCE_GENUS_NULL"
    out["mean_original_exact_genus_shared_observed"]=float(
        df.original_actual_focal_shared_source_genera.mean())
    out["mean_original_matched_taxon_source_distance_retention"]=float(
        df.original_focal_source_taxon_retention_fraction.mean())
    out["mean_original_matched_mainland_candidate_distance_retention"]=float(
        df.original_candidate_source_taxon_retention_fraction.mean())
    out["mean_original_source_distance_matched_breadth_effective_count"]=float(
        df.original_matched_effective_source_species.mean())
    out["n_conditional_original_cluster_resamples"]=boots
    out["cluster_resamples_are_not_evolutionary_independent_colonizations"]=True
    return out


def analyze(flora,geo,assign,mainland,states,*,boots=BOOT):
    source=status_groups(flora)
    if len(source)!=4505:raise ValueError("Original botanical 4505 island frame changed")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or (
        geo.island_id.duplicated().any()):
        raise ValueError("Original GSHHG 8264 geographic join invalid")
    g=geo.set_index("island_id").reindex(sorted(source))
    if int(g.analysis_regime.notna().sum())!=4504:
        raise ValueError("Original one geographic unresolved island changed")
    regions={i:(str(g.at[i,"analysis_regime"])
                if pd.notna(g.at[i,"analysis_regime"])
                else "geography_unresolved") for i in source}
    counts=Counter(regions.values())
    if set(counts)!=set(REGIONS) or sum(counts.values())!=4505:
        raise ValueError("Source original global region source support drift")
    needs={"island_id","source_mode","source_rank","entity_ID","source_distance_km"}
    if not needs<=set(assign):
        raise ValueError("Original independent mainland candidate source kilometres missing")
    ranks=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    ranks["source_rank"]=pd.to_numeric(ranks.source_rank,errors="raise").astype(int)
    ranks["source_distance_km"]=pd.to_numeric(ranks.source_distance_km,errors="raise")
    if (ranks.source_distance_km.isna().any() or
        not ranks.source_rank.between(1,5).all() or
        ranks.duplicated(["island_id","source_rank"]).any()):
        raise ValueError("Original GIFT five ranked actual source distances incomplete")
    neighbors=ranks.sort_values(["island_id","source_rank"]).groupby(
        "island_id").apply(
        lambda t:list(zip(t.entity_ID.astype(str),
                          t.source_distance_km.astype(float))),
        include_groups=False).to_dict()
    if not {"entity_ID","work_species"}<=set(mainland) or (
        mainland.duplicated(["entity_ID","work_species"]).any()):
        raise ValueError("Original source entity taxonomic GIFT names duplicated")
    byentity=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    gift_breadth=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    sector=source_sectors(states)
    if len(sector)!=7:raise ValueError("Original seven Direct source sector masks changed")
    panels=[]
    audit=[]
    for region in REGIONS:
        islands=sorted(i for i in source if regions[i]==region)
        cohorts=defaultdict(list)
        before=Counter()
        for island in islands:
            focal=source[island]
            five=neighbors.get(island,[])
            near,byrank=source_opportunity_distances(five,byentity)
            src_names=set(near)
            if not src_names:continue
            candidates=src_names-focal["full"]
            before["n_original_source_islands_with_five_nearest_regions"]+=int(len(five)==5)
            names_by_scope={
                "strict_native_nonendemic_source_listed":src_names & focal["native_nonendemic"],
                "all_original_island_species_source_listed":src_names & focal["full"],
                "original_introduced_source_listed":src_names & focal["introduced"]}
            tags=dict(
                physical_block=(str(g.at[island,"spatial_block"]) if
                                pd.notna(g.at[island,"spatial_block"]) else None),
                source_region_cluster=five[0][0] if five else None)
            for scope,actual_names in names_by_scope.items():
                if len(actual_names)<MIN_NATIVE or len(candidates)<MIN_CANDIDATE:
                    continue
                before[scope+"_n5_n10_unconditional_source_islands"]+=1
                strata={
                    sp:(listing_breadth_bin(int(gift_breadth[sp])),
                        closest_km_bin(near[sp]))
                    for sp in (actual_names|candidates)}
                null=exact_taxon_stratum_genus_reference(
                    actual_names,candidates,strata)
                if null is None:continue
                k,used_u,used_c=null
                if k["original_matched_effective_source_species"]<MIN_MATCH:
                    continue
                k=exact_same_names_breadth_baseline(
                    k,used_u,used_c,gift_breadth)
                k.update(tags)
                cohorts[(scope,"all_original_source_taxa")].append(k)
                # Real prior source-Direct masks on EXACT original island IDs.
                # Morphology only selects eligible islands; the GENUS outcome
                # above uses ALL botanical accepted source plant names.
                if scope==SCOPES[0]:
                    for category,trait in sector.items():
                        if original_trait_island_eligibility(
                            actual_names,candidates,gift_breadth,trait):
                            cohorts[(scope,"same_islands_as_Direct_"+category)].append(k)
        for scope in SCOPES:
            wanted=(
                ["all_original_source_taxa"]+
                ["same_islands_as_Direct_"+cat for cat in sector]
                if scope==SCOPES[0] else
                ["all_original_source_taxa"])
            for name in wanted:
                records=cohorts.get((scope,name),[])
                row=dict(region=region,scope=scope,cohort=name,
                         original_botanical_region_floras=len(islands),
                         n_original_island_floras_after_source_distance_and_listing_match=len(records),
                         n_islands_with_original_5source_n5n10_before_distance=before[
                             scope+"_n5_n10_unconditional_source_islands"])
                if records:
                    row.update(
                        island_equal_original_observed_genus_count=float(
                            np.mean([x["original_actual_focal_shared_source_genera"]
                                     for x in records])),
                        island_equal_original_distance_expected_genus_count=float(
                            np.mean([x["expected_shared_genera_under_breadth_distance_exact_source_null"]
                                     for x in records])),
                        island_equal_original_breadth_expected_genus_on_SAME_names=float(
                            np.mean([x["original_breadth_only_expected_shared_genera_on_SAME_names"]
                                     for x in records])))
                for cluster,label in (
                    ("physical_block","physical_10degree_block_equal"),
                    ("source_region_cluster","first_GIFT_mainland_source_equal")):
                    row[label]=group_equal(
                        records,cluster,boots=boots,seed=int.from_bytes(
                            hashlib.sha256(
                                (region+scope+name+label).encode()).digest()[:8],"big"))
                panels.append(row)
            audit.append(dict(
                region=region,scope=scope,
                n_original_region_botanical_islands=len(islands),
                n_original_source_listed_n5n10_before_distance_strata=before[
                    scope+"_n5_n10_unconditional_source_islands"],
                n_distance_and_breadth_matched_source_islands=len(
                    cohorts.get((scope,"all_original_source_taxa"),[]))))
    if sum(counts.values())!=4505:
        raise AssertionError("Global source flora population not preserved")
    return dict(
        status="ORIGINAL_4505_DIRECT_GIFT_KM_BREADTH_CONDITIONED_GENUS_SOURCE_NULL",
        n_original_botanical_GIFT_island_floras=4505,
        n_original_GIFT_mainland_native_source_regions=mainland.entity_ID.nunique(),
        n_independent_geo_k5_island_source_distance_records=len(ranks),
        n_real_original_Direct_trait_only_island_mask_types=len(sector),
        frozen_absolute_original_source_distance_km_bins=list(DISTANCE_BINS),
        original_GIFT_geography_flora_population=[
            dict(region=reg,n_original_botanical_floras=int(counts[reg]))
            for reg in REGIONS],
        original_geographical_source_distance_opportunity=audit,
        original_source_genus_null_distance_vs_breadth_same_taxa=panels,
        claims=dict(
            original_all_4505_botanical_floras_retained=True,
            source_distance_IS_real_GIFT_source_region_distance_not_dispersal_kernel=True,
            original_source_distance_only_assigned_if_source_entity_actually_lists_taxon=True,
            fixed_real_candidate_source_GIFT_listing_breadth_20_joint_bins=True,
            distance_and_breadth_models_evaluate_identical_accepted_plant_names=True,
            original_GIFT_source_rank_not_climate_or_source_ancestral_distance=True,
            no_floral_traits_used_in_primary_botanical_genus_response=True,
            strict_native_original_endemic_floristic_status_concordance_required=True,
            same_island_direct_trait_masks_overlap_and_are_not_new_replicates=True,
            original_source_taxon_records_not_actual_propagule_arrival_or_failure=True,
            group_bootstrap_not_independent_phylogenetic_histories=True,
            no_private_original_plant_island_genus_source_entity_coordinates_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for field in ("original-root","geology","output"):
        p.add_argument("--"+field,type=Path,required=True)
    a=p.parse_args()
    paths={
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "mainland":a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for key,path in paths.items():
        if sha(path)!=pins[key]:
            raise ValueError("Original protected source botanical file SHA drift: "+key)
    original=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    states=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    mainland=pd.read_csv(paths["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    assignment=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","entity_ID","source_mode","source_rank","source_distance_km"])
    geography=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    if geography.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original source geography 8264 and GIFT mainland 748 drift")
    z=analyze(original,geography,assignment,mainland,states)
    assert z["n_original_botanical_GIFT_island_floras"]==4505
    assert z["n_independent_geo_k5_island_source_distance_records"]==41325
    assert sum(x["n_original_botanical_floras"]
               for x in z["original_GIFT_geography_flora_population"])==4505
    output=json.dumps(z,indent=2,sort_keys=True,allow_nan=False)
    for forbidden in ('"island_id"','"accepted_species"','"entity_ID"',
                      '"work_species"','"spatial_block"','"source_signature"',
                      '"island_latitude"','"island_longitude"'):
        if forbidden in output:
            raise ValueError("Original private taxonomy/geographical source escaped")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(output+"\n")
    print(json.dumps(dict(
        status=z["status"],n_original_4505_floras=4505,
        n_source_geo_status_cohorts=len(z["original_source_genus_null_distance_vs_breadth_same_taxa"]),
        is_actual_seed_dispersal_kernel=False)))


if __name__=="__main__":main()
