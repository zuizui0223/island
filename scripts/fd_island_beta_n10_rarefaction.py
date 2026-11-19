"""Equal-10 directly observed flower carriers, same geographic block-disjoint dyads.

POSTHOC sensitivity of source-original northern-midlatitude 500–2000km
n>=10 and >=10% trait coverage cohort. The 4505 original botanical
floras remain the population frame. Source-region genus donor profiles
are fixed at the ORIGINAL ≥5-Direct regional reference across thinning
draws; no local donor species or phenotype states are inferred or exported.
"""
from __future__ import annotations
import argparse
import json
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from fd_island_functional_beta import (
    strict_traits,prepare_endpoint,compare_block,membership,
    geographical_distance)
from fd_island_beta_block_matching import match_blocks

REGION="northern_midlatitude"
MIN_CARRIERS=10
MIN_FRACTION=.10
N_THIN=29
N_RESELECT=199


def analyze(states:dict,flora:pd.DataFrame,geo:pd.DataFrame, *,
            n_thinnings:int=N_THIN,n_matches:int=N_RESELECT,
            min_pairs:int=100) -> dict:
    if n_thinnings<1 or n_matches<1:
        raise ValueError("Positive number of deterministic source sensitivities required")
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original flora fields missing")
    if flora.duplicated(["island_id","accepted_species"]).any() or (
        not flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original accepted plant status invalid")
    if not {"island_id","analysis_regime","spatial_block",
            "island_latitude","island_longitude"}<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Geography not independently joined")
    both=set(states["flower_primary_color"])&set(states["floral_form"])
    joint={sp:(states["flower_primary_color"][sp],states["floral_form"][sp])
           for sp in both}
    all_f=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    geo=geo.set_index("island_id")
    all_islands=sorted(all_f)
    qualified_n5=[
        i for i in all_islands
        if len(all_f[i]&both)>=5 and
        i in geo.index and geo.at[i,"analysis_regime"]==REGION]
    local_reference_names=set().union(*(all_f[i]&both for i in qualified_n5))
    local_source={sp:joint[sp] for sp in local_reference_names}
    # Freeze regional genus donors using all region-n5 source botany, *not*
    # only the n10 focal island carriers or thinned island sample.
    candidate=[
        i for i in qualified_n5
        if len(all_f[i]&both)>=MIN_CARRIERS and
        len(all_f[i]&both)/len(all_f[i])>=MIN_FRACTION]
    if len(candidate)<2:
        raise ValueError("Insufficient original 10-carrier verified focal support")
    full=[all_f[i] for i in candidate]
    direct=[all_f[i]&both for i in candidate]
    fA=membership(full,sorted(set().union(*full)))
    share=(fA@fA.T).toarray()
    g=geo.reindex(candidate)
    lat=pd.to_numeric(g.island_latitude,errors="coerce").to_numpy()
    lon=pd.to_numeric(g.island_longitude,errors="coerce").to_numpy()
    block_id,block_names=pd.factorize(g.spatial_block.to_numpy())
    distance=geographical_distance(lat,lon,0,len(candidate))
    I,J=np.triu_indices(len(candidate),1)
    valid=(
        (share[I,J]==0)&
        (block_id[I]>=0)&(block_id[J]>=0)&(block_id[I]!=block_id[J])&
        (distance[I,J]>=500)&(distance[I,J]<2000)&
        np.isfinite(distance[I,J])&
        np.isfinite(lat[I])&np.isfinite(lat[J])&
        np.isfinite(lon[I])&np.isfinite(lon[J])&
        (np.abs(lat[I])<=90)&(np.abs(lat[J])<=90)&
        (np.abs(lon[I])<=180)&(np.abs(lon[J])<=180))
    I,J=I[valid],J[valid]
    if len(I)<min_pairs:
        return dict(status="SUPPORT_INSUFFICIENT",n_full_botanical_islands=len(all_islands),
                    n_disjoint_pair_candidates=int(len(I)),no_biological_claim=True)
    b1,b2=block_id[I],block_id[J]
    pick=[
        match_blocks(b1,b2,seed=20261115+seed)
        for seed in range(n_matches)]
    n_picks={len(x) for x in pick}
    if len(n_picks)!=1 or min(n_picks)<1:
        raise AssertionError("Maximum block-cardinality not fixed")
    base=prepare_endpoint(direct,local_source)
    bmat=compare_block(base,0,len(candidate))
    full_residual=(bmat["squared_frequency_distance"]-
                   bmat["genus_null_squared_distance"])[I,J]
    full_sel=np.array([float(np.mean(full_residual[ix])) for ix in pick])
    original_pair_mean=float(np.mean(full_residual))
    original_matched_median=float(np.median(full_sel))
    # All thinning draws use exactly the SAME eligible island dyads, the SAME
    # externally unchanged region in-sample genus donors, and same matches.
    rng=np.random.default_rng(20261115)
    thin_pair=[]
    thin_match=[]
    n_fixed_reference=len(local_source)
    for draw in range(n_thinnings):
        drawn=[set(rng.choice(sorted(roster),size=MIN_CARRIERS,replace=False))
               for roster in direct]
        p=prepare_endpoint(drawn,local_source)
        if p["global_n"]!=n_fixed_reference:
            raise AssertionError("Regional reference collapsed during rarefaction")
        z=compare_block(p,0,len(candidate))
        residual=(z["squared_frequency_distance"]-
                  z["genus_null_squared_distance"])[I,J]
        thin_pair.append(float(np.mean(residual)))
        thin_match.append([float(np.mean(residual[ix])) for ix in pick])
    thin_pair=np.asarray(thin_pair)
    thin_match=np.asarray(thin_match)
    output=dict(
        status="ORIGINAL_N10_COVERAGE10_EQUAL10_SPATIAL_STRESS",
        n_full_botanical_islands=len(all_islands),
        n_original_region_n5_reference_islands=len(qualified_n5),
        n_original_region_genus_donor_species=n_fixed_reference,
        n_original_candidate_islands=len(candidate),
        n_original_fixed_full_flora_species_disjoint_pairs=len(I),
        n_distinct_candidate_spatial_blocks=len(set(np.r_[b1,b2])),
        n_maximally_block_disjoint_pairs=int(next(iter(n_picks))),
        original_full_carrier_pair_equal_mean_residual_D2=original_pair_mean,
        original_full_carrier_block_disjoint_median_residual_D2=original_matched_median,
        original_full_carrier_block_disjoint_positive_fraction=float(np.mean(full_sel>0)),
        n_deterministic_10species_thinnings=n_thinnings,
        n_outcome_blind_10degree_block_matchings=n_matches,
        equal10_pair_equal_median_residual_D2=float(np.median(thin_pair)),
        equal10_pair_equal_p05_D2=float(np.quantile(thin_pair,.05)),
        equal10_pair_equal_p95_D2=float(np.quantile(thin_pair,.95)),
        equal10_block_disjoint_median_residual_D2=float(np.median(thin_match)),
        equal10_block_disjoint_p05_D2=float(np.quantile(thin_match,.05)),
        equal10_block_disjoint_p95_D2=float(np.quantile(thin_match,.95)),
        equal10_block_disjoint_positive_fraction=float(np.mean(thin_match>0)),
        equal10_per_thinning_maximal_matching_median_min=float(
            np.min(np.median(thin_match,axis=1))),
        equal10_per_thinning_maximal_matching_median_max=float(
            np.max(np.median(thin_match,axis=1))),
        source_region=REGION,
        source_distance_bin_km=[500,2000],
        claims=dict(
            all_4505_original_flora_retained_in_population=True,
            selected_secondary_n10_and_fraction_gate_is_posthoc=True,
            exact_same_full_flora_species_disjoint_pair_graph=True,
            no_dyad_or_matching_reselection_as_independent_test=True,
            local_genus_reference_in_sample_fixed_across_rarefaction=True,
            no_adaptive_evolution_or_pollinator_shift_inference=True,
            no_taxon_island_provider_geography_identifiers_exported=True,
        ))
    return output


def main():
    a=argparse.ArgumentParser()
    for x in ("original-root","geology","output"):
        a.add_argument("--"+x,required=True,type=Path)
    args=a.parse_args()
    orig={
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"}
    sha={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"}
    for name,path in orig.items():
        h=hashlib.sha256()
        with path.open("rb") as fp:
            for chunk in iter(lambda:fp.read(1<<20),b""):h.update(chunk)
        if h.hexdigest()!=sha[name]:
            raise ValueError("Original protected source fingerprint changed")
    trait=pd.read_csv(orig["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(orig["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(args.geology,dtype=str)
    if flora.island_id.nunique()!=4505:raise ValueError("Source 4505 flora baseline changed")
    q=analyze(strict_traits(trait),flora,geo)
    if q["status"]!="ORIGINAL_N10_COVERAGE10_EQUAL10_SPATIAL_STRESS":
        raise ValueError("Original direct 10-species opportunity insufficient")
    if q["n_original_fixed_full_flora_species_disjoint_pairs"]!=456:
        raise ValueError("Original same-cohort 456 species-disjoint pair ledger changed")
    body=json.dumps(q,sort_keys=True,indent=2,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"spatial_block"',
                    '"island_latitude"','"island_longitude"','"source_signature"'):
        if private in body:
            raise ValueError("Private botanical source identity leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+'\n',encoding='utf-8')
    print(json.dumps(dict(status=q["status"],
                          n_eligible_original_plant_pairs=q["n_original_fixed_full_flora_species_disjoint_pairs"],
                          n_monte_carlo_resamples=q["n_deterministic_10species_thinnings"],
                          no_independent_evolutionary_p_value=True)))


if __name__=="__main__":
    main()
