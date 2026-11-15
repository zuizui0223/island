"""Directional *which floral profile region is absent?* nulls on true Direct form×symmetry.

External original population: 4505 botanical island floras, true same-species
Direct floral form and symmetry (NO third-trait completion). The sector is
zygomorphic AND tubular/bell versus an actinomorphic/open-radial comparator.
Only real original observed paired profile categories are in the reference.

Null A (equal k) samples k unique ORIGINAL region reference joint profiles
uniformly: a negative-control visualization, not a genus/donor-balanced test.
Null B samples n actual Direct source species from the region's measured pool
without replacement, preserving category frequency. Null C samples EXACT
island genus carrier counts from within each in-region genus WITHOUT replacement
and swaps WHOLE observed paired profiles, preserving category pairing.
All three are in-sample regional descriptive source references, NOT
independently field-measured ancestors or causal selection.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from collections import Counter,defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from fd_island_functional_beta import strict_traits
from fd_directional_flower_isolation_gradients import weighted_cluster_fit,MIN_DIRECT

SECTORS={
    "bilateral_tube_or_bell":lambda x:x[1]=="zygomorphic" and x[0] in (
        "tubular","bell_campanulate"),
    "open_radial_actinomorphic":lambda x:x[1]=="actinomorphic" and x[0]=="open_radial",
}
REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "highlatitude_ge60_nested")
SCOPES=("all_recorded","original_native")
MIN_BLOCKS=10
MIN_ISLANDS=30


def no_category_draw_probability(N,K,n):
    """P(no success category on sample of n distinct names), exact finite pool."""
    if not all(isinstance(x,(int,np.integer)) for x in (N,K,n)):
        raise ValueError("Hypergeometric original counts must be integer")
    if N<0 or K<0 or n<0 or K>N or n>N:
        raise ValueError("Impossible original finite botanical source pool")
    if n==0 or K==0:return 1.
    if n>N-K:return 0.
    return float(hypergeom.pmf(0,N,K,n))


def null_expectation_one_island(known:set,joint:dict,
                                pool_names:set,sector_name:str):
    """Actual category breadth compared with fixed-k, prevalence, genus null."""
    if sector_name not in SECTORS or not known or not known<=pool_names:
        raise ValueError("Original real within-region category source missing")
    reference={sp:joint[sp] for sp in pool_names}
    categories=sorted(set(reference.values()))
    predicate=SECTORS[sector_name]
    specialist=[c for c in categories if predicate(c)]
    if len(specialist)==0:return None
    region_source=Counter(reference.values())
    obs=set(joint[s] for s in known)
    K=len(categories)
    R=len(specialist)
    k=len(obs)
    n=len(known)
    N=len(pool_names)
    observed_r=sum(c in obs for c in specialist)
    uniform_exp=k*R/K
    uniform_no=no_category_draw_probability(K,R,k)
    prev_exp=sum(1.-no_category_draw_probability(
        N,region_source[c],n) for c in specialist)
    prevalence_no=no_category_draw_probability(
        N,sum(region_source[c] for c in specialist),n)
    genus_names=defaultdict(list)
    for sp in pool_names:
        genus=sp.split()[0] if len(sp.split())>=2 else "fixed:"+sp
        genus_names[genus].append(sp)
    k_i=Counter((sp.split()[0] if len(sp.split())>=2 else "fixed:"+sp)
                for sp in known)
    # Exact finite genus profile permutations, NO Monte Carlo pseudo-replication.
    if any(k_i[g]>len(genus_names[g]) for g in k_i):
        raise AssertionError("Island genus has more original species than its region")
    def missing_probability(accept):
        ans=1.
        for g,n_take in k_i.items():
            members=genus_names[g]
            hits=sum(accept(joint[s]) for s in members)
            ans*=no_category_draw_probability(len(members),hits,n_take)
        return ans
    genus_exp=sum(1.-missing_probability(lambda c,target=c:c==target)
                  for c in specialist)
    genus_no=missing_probability(predicate)
    return dict(
        n_direct_species=n,n_observed_joint_profiles=k,
        n_region_reference_joint_profiles=K,
        n_region_reference_special_joint_profiles=R,
        n_region_reference_distinct_species=N,
        observed_special_region_occupied=observed_r,
        observed_sector_absence=int(observed_r==0),
        expected_occupied_special_uniform_given_k=uniform_exp,
        expected_absence_uniform_given_k=uniform_no,
        expected_occupied_special_category_prevalence_fixed_n=prev_exp,
        expected_absence_category_prevalence_fixed_n=prevalence_no,
        expected_occupied_special_regional_genus_swap=genus_exp,
        expected_absence_regional_genus_swap=genus_no,
        excess_special_missing_uniform_given_k=(uniform_exp-observed_r)/R,
        excess_special_missing_prevalence_fixed_n=(prev_exp-observed_r)/R,
        excess_special_missing_genus_fixed_n=(genus_exp-observed_r)/R,
        excess_complete_sector_absence_genus_fixed_n=(
            int(observed_r==0)-genus_no),
    )


def run(states,flora,geo,assign,*,n_boot=999):
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("No original source accepted species roster")
    if flora.duplicated(["island_id","accepted_species"]).any() or not (
            flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original botanically accepted island status invalid")
    required={"island_id","analysis_regime","spatial_block",
              "distance_to_continent_km","area_km2","island_latitude"}
    if not required<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Verified original independent island geography missing")
    joint_species=set(states["floral_form"])&set(states["floral_symmetry"])
    joint={s:(states["floral_form"][s],states["floral_symmetry"][s])
           for s in joint_species}
    by_island=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    if len(by_island)!=4505:
        raise ValueError("Missing original 4505 island floras")
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    islands=sorted(by_island)
    g=geo.set_index("island_id").reindex(islands)
    first=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    first.source_rank=pd.to_numeric(first.source_rank,errors="raise")
    if first.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT source ranked candidate repeats")
    first=first.loc[first.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    regions={}
    for region in REGIONS:
        if region=="highlatitude_ge60_nested":
            iids=[i for i in islands if pd.notna(g.at[i,"island_latitude"])
                  and abs(float(g.at[i,"island_latitude"]))>=60]
        else:
            iids=[i for i in islands if g.at[i,"analysis_regime"]==region]
        regions[region]=iids
    outputs=[]
    for scope in SCOPES:
        lookup=by_island if scope=="all_recorded" else native
        for region,ri in regions.items():
            source=[lookup.get(i,set()) for i in ri]
            source_names=set().union(*source)&joint_species
            if not source_names:continue
            for sector_name in SECTORS:
                rows=[]
                for island in ri:
                    orig=lookup.get(island,set())
                    known=orig&joint_species
                    if len(known)<MIN_DIRECT:continue
                    geo_i=g.loc[island]
                    try:
                        dist=float(geo_i.distance_to_continent_km)
                        area=float(geo_i.area_km2)
                        lat=float(geo_i.island_latitude)
                    except (TypeError,ValueError):continue
                    if not np.isfinite([dist,area,lat]).all() or (
                        dist<0 or area<=0 or abs(lat)>90 or pd.isna(geo_i.spatial_block)):
                        continue
                    out=null_expectation_one_island(
                        known,joint,source_names,sector_name)
                    if out is None:continue
                    rows.append(dict(
                        **out,
                        log2_isolation_plus1=math.log2(1+dist),
                        log1p_area=math.log1p(area),
                        direct_coverage=len(known)/len(orig),
                        log_direct_n=math.log(len(known)),
                        absolute_latitude=abs(lat),
                        spatial_block=str(geo_i.spatial_block),
                        source_cluster=first.get(island) or None))
                field_names={
                    "excess_special_missing_uniform_given_k":"fixed_k_uniform",
                    "excess_special_missing_prevalence_fixed_n":"frequency_fixed_n",
                    "excess_special_missing_genus_fixed_n":"genus_fixed_n",
                    "excess_complete_sector_absence_genus_fixed_n":"genus_absence_fixed_n"}
                report=dict(scope=scope,region=region,sector=sector_name,
                            n_original_source_floras_geographic_region=len(ri),
                            n_direct_eligible_islands_with_recorded_isolation=len(rows),
                            n_reference_real_direct_joint_species=len(source_names),
                            reference_region_is_in_sample_not_independent_ancestral_pool=True)
                if not rows:
                    report["status"]="NOT_IDENTIFIABLE_NO_DIRECT_SAME_SPECIES_SUPPORT"
                    outputs.append(report);continue
                frame=pd.DataFrame(rows)
                report.update(
                    n_actual_genus_available_special_sector_profiles=int(
                        frame.n_region_reference_special_joint_profiles.iloc[0]),
                    n_real_region_joint_profiles=int(
                        frame.n_region_reference_joint_profiles.iloc[0]),
                    mean_actual_special_sector_profile_count=float(
                        frame.observed_special_region_occupied.mean()),
                    fraction_original_islands_with_no_documented_special_profile=float(
                        frame.observed_sector_absence.mean()),
                    mean_predicted_no_special_profile_given_k=float(
                        frame.expected_absence_uniform_given_k.mean()),
                    mean_predicted_no_special_profile_given_n_and_genus=float(
                        frame.expected_absence_regional_genus_swap.mean()),
                    mean_unadjusted_missing_special_vs_genus_null=float(
                        frame.excess_special_missing_genus_fixed_n.mean()),
                )
                if len(frame)<MIN_ISLANDS or frame.spatial_block.nunique()<MIN_BLOCKS:
                    report["status"]="INSUFFICIENT_GEOGRAPHIC_SUPPORT"
                    outputs.append(report);continue
                for field,label in field_names.items():
                    model=frame.copy()
                    model["response"]=frame[field]
                    stats=weighted_cluster_fit(model,cluster="spatial_block",
                           n_boot=n_boot,
                           seed=int.from_bytes(hashlib.sha256(
                               (scope+region+sector_name+label).encode()).digest()[:8],"big"))
                    for k,v in stats.items():
                        report[label+"__physical_block__"+k]=v
                report["status"]=report.get("genus_fixed_n__physical_block__status")
                outputs.append(report)
    return dict(
        status="ORIGINAL_DIRECTIONAL_MISSING_FLOWER_PROFILE_SPACE_NULLS",
        n_original_botanical_island_floras=4505,
        n_original_global_Direct_form_and_symmetry_same_species=len(joint),
        n_original_islands_with_any_direct_source_form_symmetry=sum(
            len(v&joint_species)>0 for v in by_island.values()),
        sector_names=sorted(SECTORS),
        original_regional_missingness_results=outputs,
        claims=dict(
            actual_same_species_morphology_pairs_not_three_trait_completed=True,
            uniform_k_null_not_genus_aware_or_prevalence_aware=True,
            frequency_null_uses_region_unique_species_without_replacement=True,
            genus_null_uses_true_joint_profile_label_swap_without_replacement=True,
            all_nulls_in_sample_not_evolutionary_source_ancestry=True,
            no_category_absence_claim_beyond_documented_species=True,
            physical_blocks_not_independent_phylogenetic_history=True,
            no_original_species_island_GIFT_source_geographic_identity_export=True,
            no_causal_pollinator_loss_or_adaptive_floral_evolution_claim=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for t in ("original-root","geology","output"):
        p.add_argument("--"+t,type=Path,required=True)
    p.add_argument("--n-boot",type=int,default=999)
    a=p.parse_args()
    original={
       "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
       "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
       "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    known={
       "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
       "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
       "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for k,path in original.items():
        if sha(path)!=known[k]:raise ValueError("Original botanical source SHA drift: "+k)
    if a.n_boot<99:raise ValueError("Insufficient spatial-block sensitivity")
    states=strict_traits(pd.read_csv(original["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "distance_to_continent_km","area_km2","island_latitude"])
    assign=pd.read_csv(original["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if geo.island_id.nunique()!=8264:
        raise ValueError("Original independent physical island frame changed")
    z=run(states,flora,geo,assign,n_boot=a.n_boot)
    if z["n_original_global_Direct_form_and_symmetry_same_species"]!=4450:
        raise ValueError("Original two-trait Source Direct eligibility drift")
    body=json.dumps(z,indent=2,sort_keys=True,allow_nan=False)
    for term in ('"island_id"','"accepted_species"','"entity_ID"',
                 '"source_signature"','"spatial_block"',
                 '"island_latitude"','"island_longitude"'):
        if term in body:raise ValueError("Original private identity exported")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(status=z["status"],n_flora=4505,
              n_joint_morph_direct_species=4450,n_result_rows=len(
                z["original_regional_missingness_results"]),
              no_adaptive_selection_from_cross_section=True)))


if __name__=="__main__":main()
