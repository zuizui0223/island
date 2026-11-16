"""Original 4505 GIFT island WCVP v16 native-L3 breadth direct genus null.

No flower-category prediction. Reads independently published Kew WCVP
accepted SPECIES and TDWG-L3 native distributions from verified source ZIP,
calculates exact source sampling genus overlap under:
 A = fixed original GIFT mainland listing breadth + source-region distance,
 B = SAME EXACT accepted source names PLUS independent WCVP-native L3
     region-count category. All 100 possible source strata are fixed
before original geographic/category outcomes. Missing WCVP names/native
distribution never coded as 0 regions. Focal full flora names excluded from
source-only candidate list. Source-listed island group includes documented
strict native_nonendemic (primary) and all statuses (documentary sensitivity).
Potential source regions are NOT observed propagule arrivals or actual native
range-to-island spatial edge distances; counts of native L3 areas are not
natural distribution geographic area or measured dispersal.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import zipfile
from collections import Counter,defaultdict
from pathlib import Path
from math import comb
import numpy as np
import pandas as pd
from fd_island_functional_beta import strict_traits
from fd_direct_source_finite_sampling_null import genus,prob_no_original_genus_in_native_or_candidate
from fd_direct_source_breadth_matched import BINS,REGIONS,listing_breadth_bin
from fd_direct_source_distance_conditioned_genus_null import closest_km_bin,DISTANCE_BINS
from fd_direct_three_stage_floristic_comparison import status_groups,source_sectors
from fd_original_source_trait_free_genus_null import original_trait_island_eligibility

WCVP_BINS=("nativeL3_1","nativeL3_2_3","nativeL3_4_10",
           "nativeL3_11_30","nativeL3_31plus")
SOURCE_MODES=("strict_native_nonendemic_source_listed",
              "all_original_island_species_source_listed")
DIRECTIONAL_MASKS=("zygomorphic","tube_or_bell","autonomous_vs_absent",
                   "bilateral_tube_or_bell_joint")
MIN_FOCAL=5
MIN_SOURCE=10
MIN_MATCH=5
MIN_ISLANDS=30
MIN_BLOCKS=10
BOOT=999
METRICS=("observed_minus_range_and_source_expected_shared_genera",
         "observed_minus_source_km_breadth_only_expected_shared_genera_ON_SAME_SPECIES",
         "additional_expected_genus_sharing_after_WCVP_native_breadth_conditioning")


def sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()


def native_l3_bin(n:int)->str:
    if not isinstance(n,(int,np.integer)) or n<1:
        raise ValueError("WCVP native L3 region count unavailable")
    if n==1:return WCVP_BINS[0]
    if n<=3:return WCVP_BINS[1]
    if n<=10:return WCVP_BINS[2]
    if n<=30:return WCVP_BINS[3]
    return WCVP_BINS[4]


def accepted_wcvp_native_l3_range(zip_path:Path,targets:set[str]):
    """Exact-original-GIFT accepted name matching, zero introduced/doubtful/extinct.

    Raw accepted names, WCVP IDs and native L3 codes NEVER exported from
    original-protected Github Actions. The returned coverage is anonymous.
    """
    if not zipfile.is_zipfile(zip_path):
        raise ValueError("Official WCVP source is not a valid ZIP")
    with zipfile.ZipFile(zip_path) as z:
        paths=set(z.namelist())
        names_path=next((p for p in paths if p.endswith("wcvp_names.csv")),None)
        dist_path=next((p for p in paths if p.endswith("wcvp_distribution.csv")),None)
        if not names_path or not dist_path:
            raise ValueError("Official WCVP source ZIP lacks required native range tables")
        name_id={}
        with z.open(names_path) as f:
            for frame in pd.read_csv(
                    f,sep="|",dtype=str,usecols=[
                    "plant_name_id","taxon_rank","taxon_status","taxon_name"],
                    chunksize=250000,low_memory=False):
                accepted=frame.loc[
                    frame.taxon_rank.str.strip().str.casefold().eq("species")
                    & frame.taxon_status.str.strip().str.casefold().eq("accepted")
                    & frame.taxon_name.str.strip().isin(targets),
                    ["plant_name_id","taxon_name"]]
                for pid,name in accepted.itertuples(index=False,name=None):
                    name=name.strip()
                    if not pid or not name:continue
                    if pid in name_id and name_id[pid]!=name:
                        raise ValueError("Conflicting original WCVP accepted name ID")
                    name_id[pid]=name
        L3=defaultdict(set)
        with z.open(dist_path) as f:
            for chunk in pd.read_csv(
                    f,sep="|",dtype=str,usecols=[
                    "plant_name_id","area_code_l3","introduced",
                    "extinct","location_doubtful"],chunksize=400000,
                    low_memory=False):
                mask=chunk.plant_name_id.isin(name_id)
                if not mask.any():continue
                q=chunk.loc[mask].copy()
                for flag in ("introduced","extinct","location_doubtful"):
                    q=q.loc[pd.to_numeric(q[flag],errors="coerce").eq(0)]
                q=q.loc[q.area_code_l3.notna()
                        &q.area_code_l3.str.strip().ne("")]
                for pid,region in q[["plant_name_id","area_code_l3"]].itertuples(
                    index=False,name=None):
                    L3[name_id[pid]].add(region.strip())
    counts={k:len(v) for k,v in L3.items() if len(v)}
    names=set(name_id.values())
    return counts,dict(
        wcvp_official_archive_sha256=sha(zip_path),
        wcvp_accepted_exact_species_name_matches=len(names),
        wcvp_accepted_exact_species_with_positive_native_TDWG_L3=len(counts),
        wcvp_exact_accepted_names_without_positive_native_L3=len(names)-len(counts),
        n_input_original_GIFT_distinct_accepted_species=len(targets),
        native_distribution_filters=dict(introduced=0,extinct=0,location_doubtful=0),
        exact_accepted_only_no_synonym_guessing=True,
        WCVP_recording_native_region_count_NOT_true_range_area=True,
    )


def exact_genus_strata(focal:set[str],candidate:set[str],
                       strata:dict[str,tuple],*,
                       min_effective=MIN_MATCH):
    if not focal or not candidate or focal&candidate or not (
            focal|candidate)<=set(strata):
        raise ValueError("Original botanical names or WCVP source bins inconsistent")
    pools=defaultdict(lambda:[set(),set()])
    for s in focal:pools[strata[s]][0].add(s)
    for s in candidate:pools[strata[s]][1].add(s)
    # Validate the pre-frozen source categories even if effective sample
    # support later proves too sparse. Invalid bins must fail closed.
    for key in pools:
        if (len(key) not in (2,3) or key[0] not in BINS
            or key[1] not in DISTANCE_BINS
            or (len(key)==3 and key[2] not in WCVP_BINS)):
            raise ValueError("Undefined or posthoc WCVP/native GIFT source stratum")
    used=[(k,a,b) for k,(a,b) in pools.items() if a and b]
    if not used:return None
    effective=sum(min(len(a),len(b)) for _,a,b in used)
    if effective<min_effective:return None
    alln=set();allc=set()
    expected_none_n=defaultdict(lambda:1.)
    expected_none_c=defaultdict(lambda:1.)
    within_pairs=0.
    cross=[]
    for key,names,source in used:
        if len(key) not in (2,3) or key[0] not in BINS or (
                key[1] not in DISTANCE_BINS) or (
                len(key)==3 and key[2] not in WCVP_BINS):
            raise ValueError("Undefined or posthoc candidate native source range stratum")
        alln.update(names);allc.update(source)
        n=len(names);s=len(source);N=n+s
        genus_ct=Counter(genus(x) for x in names|source)
        for g,m in genus_ct.items():
            p_n,p_c=prob_no_original_genus_in_native_or_candidate(
                int(N),int(n),int(m))
            expected_none_n[g]*=p_n
            expected_none_c[g]*=p_c
        if n>1 and N>1:
            within_pairs+=comb(n,2)*sum(m*(m-1) for m in genus_ct.values())/(N*(N-1))
        cross.append((n,N,genus_ct))
    cross_pairs=0.
    for i,(n,N,counts) in enumerate(cross):
        for m,M,other in cross[i+1:]:
            cross_pairs+=(n*m)*sum(c*other.get(g,0) for g,c in counts.items())/(N*M)
    gnames=set(expected_none_n)
    actual_gen_n={genus(x) for x in alln}
    actual_gen_c={genus(x) for x in allc}
    actual_shared=len(actual_gen_n&actual_gen_c)
    expect_shared=len(gnames)-sum(expected_none_n.values())-sum(expected_none_c.values())
    if (actual_shared+len(actual_gen_n-actual_gen_c)+len(actual_gen_c-actual_gen_n)
            !=len(gnames)):
        raise AssertionError("Original actual source genus partition fails")
    if not -1e-7<=expect_shared<=len(gnames)+1e-7:
        raise AssertionError("Exact finite WCVP source genus expectation outside range")
    n=len(alln)
    actualcounts=Counter(genus(x) for x in alln)
    actual_concentration=sum(v*(v-1) for v in actualcounts.values())/(n*(n-1)) if n>1 else None
    exp_concentration=(within_pairs+cross_pairs)/comb(n,2) if n>1 else None
    return dict(
        original_actual_shared_botanical_genera=int(actual_shared),
        expected_shared_botanical_genera_exact_finite_source=float(expect_shared),
        observed_minus_finite_expected_shared_genera=float(actual_shared-expect_shared),
        observed_minus_expected_focal_genus_concentration=float(
            actual_concentration-exp_concentration),
        n_original_botanical_focal_names_in_eligible_source_strata=int(n),
        n_original_botanical_source_only_names_in_eligible_strata=int(len(allc)),
        n_effective_equal_status_source_strata=int(effective),
        n_original_source_native_breadth_km_WCVP_strata=int(len(used)),
    ),alln,allc


def geography_matched_one_island(own:set[str],source_only:set[str],
                                gift_listing:dict[str,int],
                                nearest_km:dict[str,float],
                                native_L3:dict[str,int]):
    known=set(native_L3)
    focal=own&known&set(nearest_km)
    candidate=source_only&known&set(nearest_km)
    if len(focal)<MIN_FOCAL or len(candidate)<MIN_SOURCE:return None
    allnames=focal|candidate
    triple={
        s:(listing_breadth_bin(gift_listing[s]),
           closest_km_bin(float(nearest_km[s])),
           native_l3_bin(int(native_L3[s])))
        for s in allnames}
    three=exact_genus_strata(focal,candidate,triple)
    if three is None:return None
    full,ret_f,ret_s=three
    two={s:triple[s][:2] for s in ret_f|ret_s}
    baseline=exact_genus_strata(ret_f,ret_s,two)
    if baseline is None:raise AssertionError("Triple-supported source pool fails two-stratum baseline")
    old,base_f,base_s=baseline
    if ret_f!=base_f or ret_s!=base_s or (
            old["original_actual_shared_botanical_genera"]!=full["original_actual_shared_botanical_genera"]):
        raise AssertionError("Not comparing identical original accepted source names")
    out=dict(
        observed_minus_range_and_source_expected_shared_genera=full[
            "observed_minus_finite_expected_shared_genera"],
        observed_minus_source_km_breadth_only_expected_shared_genera_ON_SAME_SPECIES=old[
            "observed_minus_finite_expected_shared_genera"],
        additional_expected_genus_sharing_after_WCVP_native_breadth_conditioning=(
            old["expected_shared_botanical_genera_exact_finite_source"]-
            full["expected_shared_botanical_genera_exact_finite_source"]),
        original_actual_shared_botanical_genera=full[
            "original_actual_shared_botanical_genera"],
        original_both_reference_expected_genus_sharing=full[
            "expected_shared_botanical_genera_exact_finite_source"],
        n_original_focal_with_WCVP_positive_native_range=len(focal),
        n_original_source_only_with_WCVP_positive_native_range=len(candidate),
        n_original_focal_matched_after_ALL_strata=full[
            "n_original_botanical_focal_names_in_eligible_source_strata"],
        n_original_source_matched_after_ALL_strata=full[
            "n_original_botanical_source_only_names_in_eligible_strata"],
        WCVP_source_focal_retention_after_triple_overlap=len(ret_f)/len(focal),
        WCVP_source_candidate_retention_after_triple_overlap=len(ret_s)/len(candidate),
        n_original_equal_breadth_km_WCVP_matched_effective=full[
            "n_effective_equal_status_source_strata"],
        n_original_triple_source_strata=full[
            "n_original_source_native_breadth_km_WCVP_strata"],
        observed_minus_expected_genus_concentration_WCVP_range_matched=full[
            "observed_minus_expected_focal_genus_concentration"],
        observed_minus_expected_genus_concentration_km_breadth_only_ON_SAME_SPECIES=old[
            "observed_minus_expected_focal_genus_concentration"],
    )
    return out


def grouped_summary(records,group:str,*,boots=BOOT,seed=20261116):
    if not records:
        return dict(status="NO_INDEPENDENT_WCVP_NATIVE_RANGES_FOR_SOURCE_COMPARISON",
                    n_original_matched_floras=0,n_original_groups=0)
    df=pd.DataFrame(records).dropna(subset=[group,*METRICS])
    n=len(df);G=df[group].nunique()
    summary=dict(n_original_matched_floras=n,n_original_groups=G)
    if n<MIN_ISLANDS or G<MIN_BLOCKS:
        return dict(**summary,status="INSUFFICIENT_NATIVE_RANGE_SOURCE_SUPPORT")
    v=df.groupby(group,sort=True)[list(METRICS)].mean().to_numpy(float)
    rng=np.random.default_rng(seed)
    draws=v[rng.integers(G,size=(boots,G))].mean(1)
    m=v.mean(0)
    summary["status"]="ORIGINAL_WCVP_NATIVE_RANGE_CONDITIONAL_SOURCE_GENUS_NULL"
    summary["n_conditional_original_geographic_bootstraps"]=boots
    for j,k in enumerate(METRICS):
        summary[k+"_group_equal"]=float(m[j])
        summary[k+"_conditional_group_95"]=[
            float(np.quantile(draws[:,j],.025)),
            float(np.quantile(draws[:,j],.975))]
        loo=[float(v[np.arange(G)!=i,j].mean()) for i in range(G)]
        summary[k+"_delete_one_group_min"]=float(min(loo))
        summary[k+"_delete_one_group_max"]=float(max(loo))
    return summary


def analyze(flora,geo,assignment,mainland,traits,native_L3,source_manifest,*,
            boots=BOOT):
    original=status_groups(flora)
    if len(original)!=4505:raise ValueError("Full 4505 GIFT original botanical flora changed")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or (
            geo.island_id.duplicated().any()):
        raise ValueError("GSHHG independent physical source crosswalk changed")
    gg=geo.set_index("island_id").reindex(sorted(original))
    if int(gg.analysis_regime.notna().sum())!=4504:
        raise ValueError("Original one GSHHG geography-unresolved botanical flora changed")
    region_of={i:str(gg.at[i,"analysis_regime"]) if pd.notna(
        gg.at[i,"analysis_regime"]) else "geography_unresolved" for i in original}
    if set(region_of.values())!=set(REGIONS) or len(region_of)!=4505:
        raise ValueError("Original region source accounting incomplete")
    if not {"island_id","source_mode","source_rank",
            "source_distance_km","entity_ID"}<=set(assignment):
        raise ValueError("Original native GIFT source geographic distances missing")
    a=assignment.loc[assignment.source_mode.eq("geo_k5")].copy()
    a.source_rank=pd.to_numeric(a.source_rank,errors="raise")
    a.source_distance_km=pd.to_numeric(a.source_distance_km,errors="raise")
    if a.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original geography ranked native source duplicate")
    a=a.loc[a.source_rank.between(1,5)]
    byisland=a.sort_values(["island_id","source_rank"]).groupby(
        "island_id").apply(lambda x:list(zip(x.entity_ID,x.source_distance_km)),
                          include_groups=False).to_dict()
    if not {"entity_ID","work_species"}<=set(mainland) or (
            mainland.duplicated(["entity_ID","work_species"]).any()):
        raise ValueError("Source GIFT native mainland taxa not distinct by region")
    regions=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    GIFTlisting=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    targets=set(GIFTlisting)|set(flora.accepted_species)
    root_status=dict(
        n_all_original_GIFT_candidate_mainland_species=len(GIFTlisting),
        n_original_GIFT_candidate_species_with_WCVP_positive_native_L3=len(
            set(GIFTlisting)&set(native_L3)),
        n_original_island_accepted_species=len(set(flora.accepted_species)),
        n_original_island_species_with_WCVP_positive_native_L3=len(
            set(flora.accepted_species)&set(native_L3)),
        source_taxonomic_names_set_used_for_WCVP_matching=len(targets),
    )
    labels=source_sectors(traits)
    cohort_codes={
       "all_GIFT_status_source_listed":None,
       **{"same_own_islands_as_original_Direct_"+k:v
          for k,v in labels.items() if k in DIRECTIONAL_MASKS}}
    result=[]
    for region in REGIONS:
        islands=[i for i in sorted(original) if region_of[i]==region]
        bags=defaultdict(list)
        cover=defaultdict(Counter)
        for i in islands:
            F=original[i]
            near=byisland.get(i,[])
            pool=set()
            closest={}
            for entity,km in near:
                if not np.isfinite(km) or km<0:
                    raise ValueError("Original GIFT native source geographic distance invalid")
                names=regions.get(entity,set())
                pool.update(names)
                for sp in names:
                    if sp not in closest or km<closest[sp]:closest[sp]=float(km)
            source_only=pool-F["full"]
            for status in SOURCE_MODES:
                own=(F["native_nonendemic"] if status.startswith("strict")
                     else F["full"])&pool
                matched_own=own&set(native_L3)
                matched_source=source_only&set(native_L3)
                count=cover[status]
                count["n_original_botanical_islands"]+=1
                count["n_candidate_only_source_accepted_island_taxon_records"]+=len(source_only)
                count["n_candidate_only_source_WCVP_native_L3_records"]+=len(matched_source)
                count["n_focal_accepted_source_listed_records"]+=len(own)
                count["n_focal_WCVP_native_L3_records"]+=len(matched_own)
                if len(own)<MIN_FOCAL or len(source_only)<MIN_SOURCE:continue
                count["n_original_basic_GIFT_eligible_source_floras"]+=1
                q=geography_matched_one_island(
                    own,source_only,GIFTlisting,closest,native_L3)
                if q is None:continue
                block=gg.at[i,"spatial_block"]
                record=dict(**q,
                    spatial_block=str(block) if pd.notna(block) else None,
                    first_original_GIFT_mainland_source=near[0][0] if near else None)
                bags[(status,"all_GIFT_status_source_listed")].append(record)
                for cohort,sector in cohort_codes.items():
                    if sector is None:continue
                    if status.startswith("strict") and original_trait_island_eligibility(
                            own,source_only,GIFTlisting,sector):
                        bags[(status,cohort)].append(record)
        for status in SOURCE_MODES:
            obs=cover[status]
            for cohort in cohort_codes:
                if cohort!="all_GIFT_status_source_listed" and not status.startswith("strict"):
                    continue
                rows=bags.get((status,cohort),[])
                v=dict(region=region,status_scope=status,
                       island_cohort=cohort,
                       n_original_botanical_floras_in_region=len(islands),
                       n_original_full_GIFT_source_n5_candidate_n10_islands=obs[
                           "n_original_basic_GIFT_eligible_source_floras"],
                       n_original_WCVP_native_L3_matched_and_100bin_qualified_islands=len(rows),
                       original_focal_GIFT_source_species_records=obs[
                           "n_focal_accepted_source_listed_records"],
                       original_focal_source_names_with_WCVP_native_L3=obs[
                           "n_focal_WCVP_native_L3_records"],
                       original_candidate_only_source_species_records=obs[
                           "n_candidate_only_source_accepted_island_taxon_records"],
                       original_candidate_only_source_names_with_WCVP_native_L3=obs[
                           "n_candidate_only_source_WCVP_native_L3_records"],
                       focal_WCVP_native_L3_positive_source_fraction=(
                           obs["n_focal_WCVP_native_L3_records"]/
                           obs["n_focal_accepted_source_listed_records"]
                           if obs["n_focal_accepted_source_listed_records"] else None),
                       candidate_WCVP_native_L3_positive_source_fraction=(
                           obs["n_candidate_only_source_WCVP_native_L3_records"]/
                           obs["n_candidate_only_source_accepted_island_taxon_records"]
                           if obs["n_candidate_only_source_accepted_island_taxon_records"] else None))
                if rows:
                    v.update(
                        mean_retained_WCVP_positive_native_focal_fraction=float(
                            np.mean([r["WCVP_source_focal_retention_after_triple_overlap"]
                                     for r in rows])),
                        mean_retained_WCVP_positive_native_candidate_fraction=float(
                            np.mean([r["WCVP_source_candidate_retention_after_triple_overlap"]
                                     for r in rows])),
                        n_original_groups_with_WCVP_taxon_source_Direct_consistency=len(rows))
                for group,label in (
                    ("spatial_block","physical_10degree_block"),
                    ("first_original_GIFT_mainland_source","first_GIFT_mainland_source_group")):
                    v[label]=grouped_summary(rows,group,boots=boots,
                         seed=int.from_bytes(hashlib.sha256((
                             region+status+cohort+label).encode()).digest()[:8],"big"))
                result.append(v)
    return dict(
        status="ORIGINAL_4505_INDEPENDENT_Kew_WCVP_NATIVE_RANGE_BREADTH_SOURCE_GENUS_NULL",
        n_original_botanical_GIFT_island_floras=4505,
        n_original_independent_GIFT_native_mainland_source_regions=len(regions),
        n_original_individually_Direct_sector_masks=len(cohort_codes)-1,
        original_WCVP_v16_source_manifest=source_manifest,
        original_source_name_coverage_WCVP_independent_native_L3=root_status,
        source_original_native_WCVP_range_genus_coverage_and_null=result,
        claims=dict(
            original_4505_population_preserved=True,
            all_original_plant_names_in_focal_flora_excluded_from_source_candidate=True,
            WCVP_native_distribution_metadata_not_actual_local_source_edge_distance=True,
            original_source_region_km_not_species_range_edge_km=True,
            original_GIFT_range_listing_WCVP_native_range_Direct_trait_statuses_separate=True,
            no_unmatched_WCVP_species_as_zero_native_range=True,
            same_actual_accepted_name_cohorts_used_for_two_vs_three_stratum_null=True,
            source_100stratum_support_and_original_candidate_attrition_reported=True,
            source_null_is_exact_finite_name_exchangeability_not_propagule_arrival=True,
            multi_region_climate_phylogenetic_and_dispersal_traits_NOT_controlled=True,
            original_native_flora_status_selected_and_WCVP_provider_nonindependence=True,
            no_flowering_adaptation_or_pollinator_delivery_inferred=True,
            no_private_GIFT_taxon_island_genus_WCVP_L3_or_location_exported=True))


def main():
    p=argparse.ArgumentParser()
    for k in ("original-root","geology","wcvp-zip","output"):
        p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    protected=dict(
        flora=a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        traits=a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        mainland=a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        assign=a.original_root/"fixed/source/result/island_source_assignments.csv.gz")
    pins=dict(
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        mainland="3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")
    for k,path in protected.items():
        if sha(path)!=pins[k]:
            raise ValueError("Original protected GIFT botanical direct source changed "+k)
    flora=pd.read_csv(protected["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status",
        "floristic_status"])
    traits=strict_traits(pd.read_csv(protected["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    if geo.island_id.nunique()!=8264:
        raise ValueError("Original 8264 independently mapped GSHHG islands changed")
    assignment=pd.read_csv(protected["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","source_distance_km","entity_ID"])
    mainland=pd.read_csv(protected["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    if mainland.entity_ID.nunique()!=748:
        raise ValueError("Original mainland GIFT candidate region source changed")
    target=set(mainland.work_species)|set(flora.accepted_species)
    native,manifest=accepted_wcvp_native_l3_range(a.wcvp_zip,target)
    result=analyze(flora,geo,assignment,mainland,traits,native,manifest)
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for secret in ('"island_id"','"accepted_species"','"work_species"',
                   '"entity_ID"','"spatial_block"','"native_l3_codes"',
                   '"area_code_l3"','"plant_name_id"','"island_latitude"',
                   '"island_longitude"'):
        if secret in body:raise ValueError("Protected original botanical names/area source leaked")
    if result["n_original_botanical_GIFT_island_floras"]!=4505:
        raise ValueError("Full original 4505 floral botanical population lost")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(
        status=result["status"],
        n_source_floras=4505,
        n_WCVP_native_region_accepted_species=manifest[
            "wcvp_accepted_exact_species_with_positive_native_TDWG_L3"],
        source_frozen_accepted_name_Direct_labels_preserved=True,
        historical_colonization_or_pollinator_evolution_claim=False)))


if __name__=="__main__":
    main()
