"""Direct exact source-breadth-matched floristic category difference decomposition.

Full 4505 botanical island source frame, seven real Direct floral/reproductive
binary category sectors, strict source-listed native-nonendemic plants vs
candidate GIFT native mainland plants ABSENT from the whole original
island botanical checklist. Partition the SAME observed native-minus-source
frequency difference into four SIGNED contributions:
 (i)  shared-genus relative weight change,
 (ii) shared-genus between-species trait category frequency difference,
 (iii) native-only genus contribution,
 (iv)  negative source-only genus contribution.
All original species are real source-Direct accepted plants (no trait imputation).
GIFT source-listing breadth matched using original five fixed bins and
min(native,source) bin weights, with EXACT reconciliation to original
island native/mainland source-breadth-matched effect. A "source-exclusive"
genus is only absent from the observed candidate source or island checklist,
not truly absent from a seed rain. Nothing here is a causal colonization or
endemic island adaptive trait inference.
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
from fd_direct_three_stage_floristic_comparison import (
    status_groups,source_sectors,genus)
from fd_direct_source_breadth_matched import (
    BINS,REGIONS,listing_breadth_bin,breadth_matched_sector_difference
)

MIN_DIRECT_NATIVE=5
MIN_DIRECT_SOURCE=10
MIN_MATCHED_BREADTH=5
MIN_INDEPENDENT_ISLANDS=30
MIN_GEO_BLOCKS=10
SOURCE_GROUP_RESAMPLES=999

COMPONENTS=(
    "shared_genus_composition",
    "shared_genus_within_species_states",
    "native_exclusive_genus",
    "source_exclusive_genus_negative",
)
OUTCOMES=("original_breadth_matched_direct_difference",*COMPONENTS)


def exact_genus_bin_decomposition(native:list[str]|set[str],
                                  mainland:list[str]|set[str],
                                  trait:dict[str,int]):
    """One source breadth bin: all four signed terms add to actual category delta."""
    native=set(native);mainland=set(mainland)
    if not native or not mainland or native&mainland:
        raise ValueError("Both source lists must be present, distinct accepted names")
    if not native|mainland<=set(trait):
        raise ValueError("No invented original high/medium Direct species labels")
    a=defaultdict(list);b=defaultdict(list)
    for sp in native:a[genus(sp)].append(int(trait[sp]))
    for sp in mainland:b[genus(sp)].append(int(trait[sp]))
    ga=set(a);gb=set(b)
    total_a,total_b=len(native),len(mainland)
    out={x:0. for x in COMPONENTS}
    informative=0
    for g in sorted(ga&gb):
        wa,wb=len(a[g])/total_a,len(b[g])/total_b
        qa,qb=float(np.mean(a[g])),float(np.mean(b[g]))
        out["shared_genus_composition"]+=(wa-wb)*(qa+qb)/2
        out["shared_genus_within_species_states"]+=(wa+wb)*(qa-qb)/2
        if set(a[g]+b[g])=={0,1}:informative+=1
    for g in ga-gb:
        out["native_exclusive_genus"]+=len(a[g])/total_a*float(np.mean(a[g]))
    for g in gb-ga:
        out["source_exclusive_genus_negative"]-=len(b[g])/total_b*float(np.mean(b[g]))
    direct=float(np.mean([trait[x] for x in native])-
                 np.mean([trait[x] for x in mainland]))
    if abs(direct-sum(out.values()))>1e-12:
        raise AssertionError("Four signed genus components do not sum to original trait difference")
    return dict(
        original_direct_binary_sector_difference=direct,
        **out,
        n_both_status_shared_genera=len(ga&gb),
        n_native_unique_genera=len(ga-gb),
        n_source_unique_genera=len(gb-ga),
        n_distinct_shared_genera_with_actual_state_variation=informative,
        n_native_source_carriers_from_shared_genera=sum(len(a[g]) for g in ga&gb),
        n_mainland_source_carriers_from_shared_genera=sum(len(b[g]) for g in ga&gb),
        n_original_native_Direct_source_listed_species=total_a,
        n_original_mainland_candidate_source_only_Direct_species=total_b,
    )


def exact_breadth_genus_decomposition(native:set[str],mainland:set[str],
                                      trait:dict[str,int],
                                      GIFT_native_region_breadth:dict[str,int]):
    """Exact within-original-bin decomposition, all raw source bin counts.
    Native/listed-region species breadth ≥1 always source-verified. Only
    source bins with BOTH original groups exist in the matched contrast.
    """
    if not native or not mainland or native&mainland:
        raise ValueError("No original accepted native/source disjoint pair of groups")
    if not native|mainland<=set(trait):
        raise ValueError("Cannot impute missing Direct trait category")
    if not (native|mainland)<=set(GIFT_native_region_breadth):
        raise ValueError("Original source listing breadth unavailable for candidate accepted plant")
    pools={side:{p:set() for p in BINS} for side in ("native","candidate")}
    for sp in native:
        pools["native"][listing_breadth_bin(GIFT_native_region_breadth[sp])].add(sp)
    for sp in mainland:
        pools["candidate"][listing_breadth_bin(GIFT_native_region_breadth[sp])].add(sp)
    numerator=Counter()
    counts=Counter()
    matched=0
    used_bins=0
    shared_genus_names=set()
    informative_genus_names=set()
    n_shared_carriers_native=n_shared_carriers_source=0.
    genus_total_shared=genus_total_aonly=genus_total_bonly=0
    for p in BINS:
        a,b=pools["native"][p],pools["candidate"][p]
        if not a or not b:continue
        weight=min(len(a),len(b))
        if not weight:continue
        bin_result=exact_genus_bin_decomposition(a,b,trait)
        matched+=weight
        used_bins+=1
        counts["n_breadth_bins_with_zero_sharing_genera"]+=int(
            bin_result["n_both_status_shared_genera"]==0)
        for name in COMPONENTS:
            numerator[name]+=weight*bin_result[name]
        for g in set(genus(s) for s in a)&set(genus(s) for s in b):
            shared_genus_names.add(g)
            values={trait[s] for s in a|b if genus(s)==g}
            if values=={0,1}:informative_genus_names.add(g)
        n_shared_carriers_native+=weight*(
            bin_result["n_native_source_carriers_from_shared_genera"]/len(a))
        n_shared_carriers_source+=weight*(
            bin_result["n_mainland_source_carriers_from_shared_genera"]/len(b))
        genus_total_shared+=bin_result["n_both_status_shared_genera"]
        genus_total_aonly+=bin_result["n_native_unique_genera"]
        genus_total_bonly+=bin_result["n_source_unique_genera"]
    if not matched:return None
    result={name:float(numerator[name]/matched) for name in COMPONENTS}
    result["original_breadth_matched_direct_difference"]=float(sum(result.values()))
    old=breadth_matched_sector_difference(native,mainland,trait,
                                         GIFT_native_region_breadth)
    if (old is None or old["n_effectively_matched_source_listing_breadth_carriers"]!=matched
            or old["n_original_breadth_overlapping_bin_support"]!=used_bins
            or abs(result["original_breadth_matched_direct_difference"]-
                   old["source_breadth_matched_native_minus_candidate_sector_share"])>1e-12):
        raise AssertionError("New SAME-ISLAND exact decomposition changes original verified breadth matched source contrast")
    result.update(
        n_original_matched_effective_source_species=matched,
        n_independently_original_source_listing_breadth_bins_used=used_bins,
        n_bin_specific_shared_genus_units=genus_total_shared,
        n_bin_specific_native_exclusive_genus_units=genus_total_aonly,
        n_bin_specific_source_exclusive_genus_units=genus_total_bonly,
        n_distinct_shared_genera=len(shared_genus_names),
        n_distinct_SHARED_genera_with_REAL_category_variation=len(
            informative_genus_names),
        fraction_effective_native_carrier_weight_in_shared_genera=float(
            n_shared_carriers_native/matched),
        fraction_effective_source_candidate_weight_in_shared_genera=float(
            n_shared_carriers_source/matched),
        direct_unmatched_native_minus_source_difference=float(
            old["same_island_original_unmatched_native_minus_candidate_sector_share"]),
    )
    return result


def group_summary(records:list[dict],name:str,*,boot:int=SOURCE_GROUP_RESAMPLES,seed:int=20261115):
    """Physical-block or first-GIFT-source equal conditional original category diff."""
    if not records:
        return dict(status="NOT_IDENTIFIABLE_NO_DIRECT_SOURCE_COMPARISON",
                    n_original_comparable_island_floras=0,n_source_groups=0)
    df=pd.DataFrame(records)
    if not set(OUTCOMES)|{name}<=set(df):
        raise ValueError("Missing one signed true original source category component")
    x=df.dropna(subset=[name,*OUTCOMES])
    n=len(x)
    n_groups=x[name].nunique()
    out=dict(n_original_comparable_island_floras=n,n_source_groups=n_groups,
             island_equal_mean_direct_fraction=float(x[OUTCOMES[0]].mean()) if n else None)
    if n<MIN_INDEPENDENT_ISLANDS or n_groups<MIN_GEO_BLOCKS:
        return dict(**out,status="NOT_IDENTIFIABLE_FEW_INDEPENDENT_SOURCE_BLOCKS")
    a=x.groupby(name,sort=True)[list(OUTCOMES)].mean().to_numpy(float)
    result=a.mean(0)
    if abs(result[0]-result[1:].sum())>1e-10:
        raise AssertionError("Original physical block equal decomposition does not reconcile")
    rng=np.random.default_rng(seed)
    index=rng.integers(len(a),size=(boot,len(a)))
    draws=a[index].mean(1)
    if not np.allclose(draws[:,0],draws[:,1:].sum(1),atol=1e-10):
        raise AssertionError("Original conditional block reference loses signed summation")
    out.update(
        status="ORIGINAL_EXACT_SOURCE_BREADTH_BALANCED_GENUS_COMPONENTS",
        n_conditional_10degree_or_source_resamples=boot,
        mean_original_block_island_count=n/len(a),
        n_original_qualified_islands_with_SHARED_source_genus=int(
            (x.n_distinct_shared_genera>0).sum()),
        n_original_qualified_islands_with_ANY_REAL_variable_SHARED_genus=int(
            (x.n_distinct_SHARED_genera_with_REAL_category_variation>0).sum()),
        mean_effective_native_share_in_shared_genera=float(
            x.fraction_effective_native_carrier_weight_in_shared_genera.mean()),
        mean_effective_mainland_share_in_shared_genera=float(
            x.fraction_effective_source_candidate_weight_in_shared_genera.mean()),
        mean_n_qualifying_distinct_shared_genera=float(
            x.n_distinct_shared_genera.mean()),
        mean_n_qualifying_category_variable_genera=float(
            x.n_distinct_SHARED_genera_with_REAL_category_variation.mean()),
        fraction_original_islands_direct_difference_positive=float(
            (x[OUTCOMES[0]]>0).mean()),
    )
    for j,field in enumerate(OUTCOMES):
        out[field+"_physical_or_source_group_equal_mean"]=float(result[j])
        out[field+"_conditional_group_resampling_95"]=[
            float(np.quantile(draws[:,j],.025)),
            float(np.quantile(draws[:,j],.975))]
    out["four_components_add_exactly_to_direct_observed_delta"]=True
    out["signed_components_not_causal_percentages_of_evolution"]=True
    return out


def analyse(states:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            source_assignments:pd.DataFrame,mainland:pd.DataFrame,
            *,n_boot=SOURCE_GROUP_RESAMPLES):
    if not {"island_id","accepted_species","origin_status",
            "endemic_status","floristic_status"}<=set(flora):
        raise ValueError("Missing original separately audited nativeness/endemicity")
    full=status_groups(flora)
    if len(full)!=4505:raise ValueError("Original GIFT 4505 botanical island frame drift")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or (
        geo.island_id.duplicated().any()):
        raise ValueError("Original geography source invalid")
    g=geo.set_index("island_id").reindex(sorted(full))
    if g.analysis_regime.notna().sum()!=4504:
        raise ValueError("Original one unlocated GIFT botanical flora changed")
    region_of={i:(str(g.at[i,"analysis_regime"])
                  if pd.notna(g.at[i,"analysis_regime"])
                  else "geography_unresolved") for i in full}
    region_count=Counter(region_of.values())
    if set(region_count)!=set(REGIONS) or sum(region_count.values())!=4505:
        raise ValueError("The original entire GIFT source region frame changed")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(source_assignments):
        raise ValueError("Missing original five native GIFT mainland candidate rankings")
    ranks=source_assignments.loc[
        source_assignments.source_mode.eq("geo_k5")].copy()
    ranks.source_rank=pd.to_numeric(ranks.source_rank,errors="raise")
    ranks=ranks.loc[ranks.source_rank.between(1,5)]
    if ranks.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original source candidate rank")
    source_by_island=ranks.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if not {"entity_ID","work_species"}<=set(mainland) or mainland.duplicated(
            ["entity_ID","work_species"]).any():
        raise ValueError("Protected source GIFT mainland species list corrupted")
    source_names=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    listing=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    traits=source_sectors(states)
    summary=[]
    rawsupport=[]
    for region in REGIONS:
        ids=sorted(i for i in full if region_of[i]==region)
        for trait,traitlabel in traits.items():
            qualified=[]
            n_prequal=0
            original_native_direct=0
            original_source_candidate_direct=0
            for island in ids:
                q=full[island]
                five=source_by_island.get(island,())
                available=set().union(*(source_names.get(e,set()) for e in five))
                native=(q["native_nonendemic"]&available)&set(traitlabel)
                absent=(available-q["full"])&set(traitlabel)
                original_native_direct+=len(native)
                original_source_candidate_direct+=len(absent)
                if len(native)<MIN_DIRECT_NATIVE or len(absent)<MIN_DIRECT_SOURCE:
                    continue
                n_prequal+=1
                decomposition=exact_breadth_genus_decomposition(
                    native,absent,traitlabel,listing)
                if decomposition is None or decomposition[
                    "n_original_matched_effective_source_species"]<MIN_MATCHED_BREADTH:
                    continue
                region_block=g.at[island,"spatial_block"]
                qualified.append(dict(
                    **decomposition,
                    spatial_block=str(region_block) if pd.notna(region_block) else None,
                    first_source=five[0] if len(five) else None))
            rawsupport.append(dict(
                region=region,trait_sector=trait,
                n_original_botanical_region_island_floras=len(ids),
                n_original_source_prequalified_native_n5_candidate_n10_islands=n_prequal,
                n_original_effective_breadth_matched_islands=len(qualified),
                n_native_source_listed_Direct_island_x_species_occurrences=original_native_direct,
                n_candidate_mainland_source_only_Direct_island_x_species_occurrences=original_source_candidate_direct))
            panel=dict(region=region,trait_sector=trait,
                       n_original_qualified_islands=len(qualified),
                       original_matched_breadth_genus_composition_parts=True)
            for group,alias in (
                ("spatial_block","original_10degree_geographic_block_equal"),
                ("first_source","original_first_GIFT_mainland_source_equal")):
                panel[alias]=group_summary(
                    qualified,group,boot=n_boot,
                    seed=int.from_bytes(hashlib.sha256(
                        (region+trait+alias).encode()).digest()[:8],"big"))
            summary.append(panel)
    return dict(
        status="ORIGINAL_4505_SOURCE_BREADTH_MATCHED_EXACT_GENUS_FLORISTIC_PARTITION",
        n_original_GIFT_botanical_island_floras=4505,
        n_original_real_binary_source_trait_sectors=len(traits),
        n_original_typed_gift_native_mainland_region_entities=mainland.entity_ID.nunique(),
        n_original_geographical_population_source_regimes=len(REGIONS),
        original_4505_botanical_flora_source_region_census=[
            dict(region=r,n_original_botanical_floras=int(region_count[r]))
            for r in REGIONS],
        original_breadth_match_and_genus_support_by_trait_region=rawsupport,
        original_real_genus_partition_direct_flora_differences=summary,
        claims=dict(
            four_exact_signed_components_sum_to_observed_island_frequency_difference=True,
            actual_original_same_island_GIFT_listing_breadth_balance_preserved=True,
            all_original_focal_plant_names_excluded_from_candidate_mainland_only=True,
            actual_direct_source_category_states_no_unmeasured_imputation=True,
            within_shared_genus_different_species_not_local_population_adaptation=True,
            exclusive_genus_in_checklist_not_failed_true_seed_arrival=True,
            original_five_GIFT_mainland_region_source_not_ancestry=True,
            original_10degree_and_first_GIFT_source_groups_not_independent_evolution=True,
            no_genus_prediction_accuracy_fitted_as_primary_response=True,
            all_original_4505_botanical_floras_not_928_complete_cases_retained=True,
            no_historical_colonization_pollinator_loss_or_seed_fitness_claim=True,
            no_original_private_taxon_island_genus_source_entity_geo_ids_exported=True,
        ))


def sha(path):
    digest=hashlib.sha256()
    with path.open("rb") as f:
        for bit in iter(lambda:f.read(1<<20),b""):
            digest.update(bit)
    return digest.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for key in ("original-root","geology","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    source={
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "trait":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "gift":a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    expected={
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "trait":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "gift":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for key,name in source.items():
        if sha(name)!=expected[key]:
            raise ValueError("Protected GIFT botanical Direct source SHA drift "+key)
    state=strict_traits(pd.read_csv(source["trait"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(source["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status",
        "endemic_status","floristic_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assignment=pd.read_csv(source["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(source["gift"],dtype=str,usecols=[
        "entity_ID","work_species"])
    if geo.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original independently audited island and source GIFT frame changed")
    q=analyse(state,flora,geo,assignment,mainland)
    assert q["n_original_GIFT_botanical_island_floras"]==4505
    content=json.dumps(q,indent=2,sort_keys=True,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"work_species"',
                    '"entity_ID"','"spatial_block"','"source_signature"',
                    '"island_latitude"','"island_longitude"'):
        if private in content:
            raise ValueError("Protected accepted plant or GIFT source location leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(content+"\n")
    print(json.dumps(dict(
        status=q["status"],n_GIFT_original_botanical_islands=4505,
        n_direct_trait_region_panels=len(q[
            "original_real_genus_partition_direct_flora_differences"]),
        no_colonization_or_phenotype_adaptation_claim=True)))


if __name__=="__main__":
    main()
