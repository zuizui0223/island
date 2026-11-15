"""Original 4505-island directional near/far TRUE form×symmetry profile-set comparison.

The partial-subset null conditions k_far<k_near; the exact regional-genus
reference evaluates ALL outcome-blind eligible pairs, without selecting
on realized far k. The in-region genus reference holds actual near profiles
fixed and randomizes far genus profiles without replacement. Not a joint
historical source, adaptive loss or pollinator effectiveness model.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from collections import Counter,defaultdict
from functools import lru_cache
from math import comb
import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits,membership,geographical_distance
from fd_island_beta_block_matching import match_blocks,graph_support

REGIONS=("northern_midlatitude","tropical","southern_extratropical")
SCOPES=("all_recorded","original_native")
BANDS=(("all_distances",0.,float("inf")),
       ("less_than_2000km",0.,2000.),
       ("at_least_2000km",2000.,float("inf")))
MATCHES=199
ISOLATION_RATIO=2.
MIN_DIRECT=5
MIN_PAIRS,MIN_ISLANDS,MIN_BLOCKS=100,30,10


@lru_cache(maxsize=400000)
def finite_all_allowed_probability(N:int,allowed:int,n:int)->float:
    """Exact finite-population all-success probability."""
    if any(not isinstance(x,int) for x in (N,allowed,n)):
        raise ValueError("Original species source cardinalities must be integer")
    if not 0<=n<=N or not 0<=allowed<=N:
        raise ValueError("Invalid finite genus reference counts")
    if n>allowed:return 0.
    p=1.
    for j in range(n):p*=(allowed-j)/(N-j)
    return float(p)


def elementary_symmetric(weights:tuple[float,...],k:int)->float:
    """Exact weighted-k finite distinct-category elementary symmetric term."""
    if k<0:raise ValueError("Negative original number of profiles")
    if k>len(weights):return 0.
    a=np.zeros(k+1)
    a[0]=1.
    for i,w in enumerate(weights):
        if not np.isfinite(w) or w<0:
            raise ValueError("Negative/nonfinite recorded source category frequency")
        for j in range(min(k,i+1),0,-1):
            a[j]+=w*a[j-1]
    return float(a[k])


def uniform_fixed_k(K:int,k_near:int,k_far:int)->float:
    if not 0<=k_near<=K or not 0<=k_far<=K:
        raise ValueError("Original profile count beyond source reference")
    return float(comb(k_near,k_far)/comb(K,k_far)) if k_far<=k_near else 0.


def weighted_fixed_k(near_mask:int,k_far:int,weights:tuple[float,...])->float:
    if near_mask<0 or near_mask >= (1<<len(weights)) or k_far<0:
        raise ValueError("Impossible observed joint profile category mask")
    if near_mask.bit_count()<k_far:return 0.
    denominator=elementary_symmetric(weights,k_far)
    if denominator<=0:raise ValueError("Weighted original profile source null unsupported")
    subset=tuple(w for i,w in enumerate(weights) if near_mask&(1<<i))
    return float(elementary_symmetric(subset,k_far)/denominator)


def genus_subset_expectation(near_mask:int,far_genus:Counter,reference:dict,
                             cat_index:dict)->float:
    """Unconditional observed far subset probability; observed near set FIXED."""
    p=1.
    for genus,n in far_genus.items():
        counts=reference[genus]
        N=int(sum(counts.values()))
        allowable=int(sum(c for state,c in counts.items()
                          if near_mask&(1<<cat_index[state])))
        p*=finite_all_allowed_probability(N,allowable,int(n))
        if p==0.:break
    return float(p)


def outcome_blind_matching(rows:list[dict],draws:int=MATCHES)->dict:
    a=np.array([z["block_near"] for z in rows],int)
    b=np.array([z["block_far"] for z in rows],int)
    unique=set()
    selections=[]
    for seed in range(draws):
        chosen=match_blocks(a,b,seed=20261115+seed)
        selections.append(chosen)
        unique.update(map(int,chosen))
    return dict(selections=selections,unique=unique,graph=graph_support(a,b))


def report_matching(rows,matched)->dict:
    values=defaultdict(list)
    for chosen in matched["selections"]:
        if not len(chosen):continue
        subset=[rows[i] for i in chosen]
        low=[x for x in subset if x["kfar"]<x["knear"]]
        values["n_disjoint_block_pairs"].append(len(chosen))
        values["observed_all_pair_far_subset_probability"].append(float(np.mean([
            x["subset"] for x in subset])))
        values["genus_conditional_all_pair_expected_subset_probability"].append(float(
            np.mean([x["genus_null"] for x in subset])))
        values["observed_minus_genus_all_pairs"].append(float(np.mean([
            x["subset"]-x["genus_null"] for x in subset])))
        values["fraction_matched_pairs_kfar_lt_knear"].append(len(low)/len(subset))
        if low:
            values["observed_strict_subset_given_kfar_lt_knear"].append(float(
                np.mean([x["strict_subset"] for x in low])))
            values["exact_uniform_fixed_k_reference"].append(float(np.mean([
                x["uniform"] for x in low])))
            values["exact_weighted_fixed_k_reference"].append(float(np.mean([
                x["weighted"] for x in low])))
            values["observed_minus_uniform_fixed_k"].append(float(np.mean([
                x["strict_subset"]-x["uniform"] for x in low])))
            values["observed_minus_weighted_fixed_k"].append(float(np.mean([
                x["strict_subset"]-x["weighted"] for x in low])))
    results=dict(n_outcome_blind_block_pair_reselections=len(matched["selections"]),
                 **matched["graph"])
    for key,x in values.items():
        if key=="n_disjoint_block_pairs":
            results[key+"_minimum"]=int(min(x))
            results[key+"_maximum"]=int(max(x))
        else:
            results[key+"_median"]=float(np.median(x))
            results[key+"_p05"]=float(np.quantile(x,.05))
            results[key+"_p95"]=float(np.quantile(x,.95))
            results[key+"_positive_fraction"]=float(np.mean(np.array(x)>0))
    return results


def analyze(states:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            *,draws:int=MATCHES,min_pairs:int=MIN_PAIRS,
            min_islands:int=MIN_ISLANDS,min_blocks:int=MIN_BLOCKS)->dict:
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original accepted source botanical flora incomplete")
    if flora.duplicated(["island_id","accepted_species"]).any() or (
        not flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original GIFT plant statuses or unique names invalid")
    need={"island_id","analysis_regime","spatial_block",
          "distance_to_continent_km","island_latitude","island_longitude","area_km2"}
    if not need<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Original geographic isolation missing")
    both=set(states["floral_form"])&set(states["floral_symmetry"])
    joint={s:(states["floral_form"][s],states["floral_symmetry"][s])
           for s in both}
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("Original 4505 island population changed")
    geo=geo.set_index("island_id")
    coverage=[]
    results=[]
    for scope,source in (("all_recorded",full),("original_native",native)):
        for region in REGIONS:
            ids=[id for id in sorted(full)
                 if id in geo.index and geo.at[id,"analysis_regime"]==region]
            pool=set().union(*(source.get(i,set()) for i in ids))&both if ids else set()
            category_freq=Counter(joint[x] for x in pool)
            categories=sorted(category_freq,key=repr)
            if not categories:continue
            lookup={cat:i for i,cat in enumerate(categories)}
            weights=tuple(category_freq[cat]/len(pool) for cat in categories)
            genus_reference=defaultdict(Counter)
            for sp in pool:
                genus=sp.split()[0] if len(sp.split())>=2 else "fixed:"+sp
                genus_reference[genus][joint[sp]]+=1
            z=[]
            for id in ids:
                known=source.get(id,set())&both
                if len(known)<MIN_DIRECT:continue
                record=geo.loc[id]
                try:
                    d=float(record.distance_to_continent_km)
                    lat=float(record.island_latitude)
                    lon=float(record.island_longitude)
                    area=float(record.area_km2)
                except (TypeError,ValueError):continue
                if not np.isfinite([d,lat,lon,area]).all() or (
                    d<0 or area<=0 or abs(lat)>90 or abs(lon)>180 or
                    pd.isna(record.spatial_block)):continue
                m=0
                for sp in known:m|=1<<lookup[joint[sp]]
                gen=Counter(sp.split()[0] if len(sp.split())>=2
                            else "fixed:"+sp for sp in known)
                z.append(dict(id=id,full=full[id],mask=m,k=m.bit_count(),
                              n=len(known),genera=gen,d=d,lat=lat,lon=lon,
                              block=str(record.spatial_block)))
            coverage.append(dict(scope=scope,region=region,
                n_original_botanical_region_floras=len(ids),
                n_n5_geo_isolation_eligible_floras=len(z),
                n_original_region_Direct_species=len(pool),
                n_source_recorded_joint_profiles=len(categories)))
            if len(z)<2:continue
            n=len(z)
            m=membership([v["full"] for v in z],
                         sorted(set().union(*(v["full"] for v in z))))
            shared=(m@m.T).toarray()
            block,labels=pd.factorize([v["block"] for v in z])
            d=np.array([v["d"] for v in z])
            km=geographical_distance(
                np.array([v["lat"] for v in z]),
                np.array([v["lon"] for v in z]),0,n)
            i,j=np.triu_indices(n,1)
            far=np.where(d[i]>d[j],i,j)
            near=np.where(d[i]>d[j],j,i)
            eligible=(shared[i,j]==0)&(block[i]!=block[j])&(
                1+d[far]>=ISOLATION_RATIO*(1+d[near]))
            F=far[eligible]
            N=near[eligible]
            D=km[i[eligible],j[eligible]]
            for label,lower,upper in BANDS:
                take=np.flatnonzero((D>=lower)&(D<upper))
                far_i=F[take];near_i=N[take]
                n_pairs=len(take)
                active=len(set(np.r_[far_i,near_i]))
                nblocks=len(set(np.r_[block[far_i],block[near_i]]))
                base=dict(scope=scope,region=region,physical_pair_distance=label,
                    n_eligible_full_flora_species_disjoint_near_far_pairs=n_pairs,
                    n_unique_eligible_original_island_floras=active,
                    n_original_physical_10deg_blocks=nblocks,
                    n_original_reference_joint_categories=len(categories))
                if n_pairs<min_pairs or active<min_islands or nblocks<min_blocks:
                    results.append(dict(**base,status="INSUFFICIENT_ORIGINAL_PAIRS"));continue
                @lru_cache(maxsize=400000)
                def weighted(mask,k):return weighted_fixed_k(mask,k,weights)
                @lru_cache(maxsize=100000)
                def uniform(near_k,far_k):
                    return uniform_fixed_k(len(categories),near_k,far_k)
                rows=[]
                for fi,ni in zip(far_i,near_i):
                    f=z[fi];nearer=z[ni]
                    contained=int((f["mask"]&~nearer["mask"])==0)
                    strict=int(contained and f["mask"]!=nearer["mask"])
                    low=f["k"]<nearer["k"]
                    rows.append(dict(
                        block_near=int(block[ni]),block_far=int(block[fi]),
                        far=int(fi),near=int(ni),
                        kfar=f["k"],knear=nearer["k"],
                        subset=contained,strict_subset=strict,
                        uniform=uniform(nearer["k"],f["k"]) if low else 0.,
                        weighted=weighted(nearer["mask"],f["k"]) if low else 0.))
                small=[r for r in rows if r["kfar"]<r["knear"]]
                base.update(
                    fraction_far_fewer_profiles_than_near=len(small)/len(rows),
                    n_pairs_where_strict_subset_mathematically_possible=len(small),
                    actual_far_subset_near_fraction_all=float(np.mean([
                        x["subset"] for x in rows])),
                    actual_far_strict_subset_near_fraction_all=float(np.mean([
                        x["strict_subset"] for x in rows])),
                    actual_strict_subset_fraction_conditioned_on_lower_k=(
                        float(np.mean([x["strict_subset"] for x in small]))
                        if small else None),
                    exact_uniform_fixed_k_expected_subset_fraction=(
                        float(np.mean([x["uniform"] for x in small]))
                        if small else None),
                    exact_weighted_fixed_k_expected_subset_fraction=(
                        float(np.mean([x["weighted"] for x in small]))
                        if small else None))
                selection=outcome_blind_matching(rows,draws=draws)
                for ix in selection["unique"]:
                    r=rows[ix]
                    r["genus_null"]=genus_subset_expectation(
                        z[r["near"]]["mask"],z[r["far"]]["genera"],
                        genus_reference,lookup)
                base["outcome_blind_block_disjoint_sensitivity"]=report_matching(
                    rows,selection)
                results.append(dict(**base,status="ORIGINAL_EXPLORATORY_NEAR_FAR_SUBSET"))
    return dict(
        status="ORIGINAL_DIRECTIONAL_NEAR_FAR_CATEGORICAL_NESTEDNESS",
        n_original_botanical_island_floras=len(full),
        n_original_same_species_Direct_form_symmetry_taxa=len(joint),
        original_region_source_support=coverage,
        original_near_far_directional_tests=results,
        claims=dict(
            full_4505_source_floras_maintained=True,
            original_far_has_at_least_twice_1_plus_continent_distance=True,
            original_full_flora_zero_shared_accepted_species_required=True,
            original_uniform_k_null_only_conditional_on_far_lower_k=True,
            original_weighted_k_null_conditions_on_observed_profile_richness=True,
            genus_null_applied_to_all_outcome_blind_pairs_not_k_selected=True,
            genus_null_conditions_on_actual_near_profiles_and_far_genus_counts=True,
            same_region_original_source_genus_reference_is_in_sample=True,
            numerical_matched_block_draws_not_independent_replicates=True,
            no_independent_evolution_or_pollination_mechanism_claim=True,
            no_private_island_species_source_geo_identifiers_exported=True))


def sha(p):
    h=hashlib.sha256()
    with p.open("rb") as f:
        for x in iter(lambda:f.read(1<<20),b""):h.update(x)
    return h.hexdigest()


def main():
    parser=argparse.ArgumentParser()
    for x in ("original-root","geology","output"):
        parser.add_argument("--"+x,type=Path,required=True)
    a=parser.parse_args()
    original={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz"}
    checks={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a"}
    for name,p in original.items():
        if sha(p)!=checks[name]:
            raise ValueError("Protected original source archived fingerprint changed")
    states=strict_traits(pd.read_csv(original["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(original["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geography=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block",
        "distance_to_continent_km","island_latitude",
        "island_longitude","area_km2"])
    if geography.island_id.nunique()!=8264:
        raise ValueError("Original independently classified geography changed")
    report=analyze(states,flora,geography)
    if (report["n_original_botanical_island_floras"]!=4505 or
            report["n_original_same_species_Direct_form_symmetry_taxa"]!=4450):
        raise ValueError("Full original source Direct flora pair support changed")
    body=json.dumps(report,indent=2,sort_keys=True,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"entity_ID"',
                    '"spatial_block"','"source_signature"',
                    '"island_latitude"','"island_longitude"'):
        if private in body:raise ValueError("Private botanical original key leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(status=report["status"],n_original_floras=4505,
          n_original_result_strata=len(report["original_near_far_directional_tests"]),
          no_historical_trait_loss_or_pollinator_claim=True)))


if __name__=="__main__":main()
