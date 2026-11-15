"""Source-list breadth matched original 4505-island GENUS null, NO floral trait score.

This genuinely independent audit of the SOURCE trait-documentation gate asks
whether original island accepted-plant genus richness/overlap itself departs
from source-name exchangeability, when sampling all source-list botanical names
rather than just plants with Direct flower annotations. Original strict native
nonendemic and the broader all-status source-listed flora are different,
nonindependent source opportunities. Exactly the SAME island IDs previously
Direct-qualified for each of seven floral/reproductive phenotypes are audited
without using ANY of their phenotype category values in the genus outcome.

Original GIFT mainland candidate lists are NOT observed seed rain, individual
establishment trials or ancestors. Unknown botanical native status cannot be
treated as native in the all-status negative-control sample.
"""
from __future__ import annotations
import argparse,hashlib,json
from collections import Counter,defaultdict
from functools import lru_cache
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    status_groups,source_sectors)
from fd_direct_source_breadth_matched import BINS,REGIONS,listing_breadth_bin
from fd_direct_source_finite_sampling_null import (
    genus,prob_no_original_genus_in_native_or_candidate
)

BOOT=999
MIN_NATIVE=5
MIN_CANDIDATE=10
MIN_MATCH=5
MIN_ISLANDS=30
MIN_BLOCKS=10
SCOPES=("strict_native_nonendemic_source_listed",
        "all_original_island_species_source_listed",
        "original_introduced_source_listed")
PRIMARY_SCOPES=SCOPES[:2]
METRICS=(
    "shared_original_genus_count_minus_finite_source_expectation",
    "source_only_original_genus_count_minus_finite_source_expectation",
    "island_only_original_genus_count_minus_finite_source_expectation",
    "observed_minus_expected_original_island_genus_pair_concentration"
)


def exact_taxonomic_source_null(native:set[str],candidate:set[str],
                                listing:dict[str,int]):
    """Exact genotype-free genus richness and concentration expectation.

    Null draws the ORIGINAL number n native accepted named species in each
    original source listing breadth bin from ALL actual source taxa in the
    matched bin, without replacement. Does NOT require Direct trait records.
    """
    if not native or not candidate or native&candidate:
        raise ValueError("Original source and focal accepted names disjointness invalid")
    if not (native|candidate)<=set(listing):
        raise ValueError("Original GIFT native candidate listing breadth missing")
    bins={b:{"native":set(),"source":set()} for b in BINS}
    for sp in native:bins[listing_breadth_bin(listing[sp])]["native"].add(sp)
    for sp in candidate:bins[listing_breadth_bin(listing[sp])]["source"].add(sp)
    used=[]
    shared_genus=set()
    actual_native_genus=set()
    actual_source_genus=set()
    p_absent_from_native=defaultdict(lambda:1.)
    p_absent_from_source=defaultdict(lambda:1.)
    within_expected_pairs=0.
    source_by_bin_genus=[]
    n_drawn_total=0
    effective=0
    for breadth,parts in bins.items():
        nset,sset=parts["native"],parts["source"]
        if not nset or not sset:continue
        n,s=len(nset),len(sset)
        N=n+s
        effective+=min(n,s)
        n_drawn_total+=n
        counts=Counter(genus(x) for x in nset|sset)
        actual_native_genus.update(genus(x) for x in nset)
        actual_source_genus.update(genus(x) for x in sset)
        for g,m in counts.items():
            qnat,qsrc=prob_no_original_genus_in_native_or_candidate(
                int(N),int(n),int(m))
            p_absent_from_native[g]*=qnat
            p_absent_from_source[g]*=qsrc
        if N>1 and n>1:
            p_same_source=sum(m*(m-1) for m in counts.values())/(N*(N-1))
            within_expected_pairs+=comb(n,2)*p_same_source
        source_by_bin_genus.append((n,N,counts))
        used.append(breadth)
    if effective==0:return None
    expected_cross_pairs=0.
    for ix,(n,N,counts) in enumerate(source_by_bin_genus):
        for m,M,other_counts in source_by_bin_genus[ix+1:]:
            # Different GIFT source breadth bins are disjoint taxon pools;
            # drawing each sampled group independently yields this exact
            # expected same-genus cross-bin pair count.
            numerator=sum(count*other_counts.get(g,0)
                          for g,count in counts.items())
            expected_cross_pairs+=(n*m)*numerator/(N*M)
    G=set(p_absent_from_native)
    observed_shared=len(actual_native_genus&actual_source_genus)
    observed_native_only=len(actual_native_genus-actual_source_genus)
    observed_source_only=len(actual_source_genus-actual_native_genus)
    exp_source_only=float(sum(p_absent_from_native[g] for g in G))
    exp_native_only=float(sum(p_absent_from_source[g] for g in G))
    exp_shared=float(len(G)-exp_native_only-exp_source_only)
    if observed_shared+observed_source_only+observed_native_only!=len(G):
        raise AssertionError("Source-checklist actual genus set accounting not exact")
    if abs(exp_shared+exp_native_only+exp_source_only-len(G))>1e-8:
        raise AssertionError("Exact original source genus expectation not additive")
    n=n_drawn_total
    observed_native_genus_counts=Counter(
        genus(x) for breadth in used for x in bins[breadth]["native"])
    if sum(observed_native_genus_counts.values())!=n:
        raise AssertionError("Original source native richness changed")
    observed_pairs=sum(m*(m-1) for m in observed_native_genus_counts.values())/(
        n*(n-1)) if n>=2 else None
    expected_pairs=(within_expected_pairs+expected_cross_pairs)/comb(n,2) if n>=2 else None
    if expected_pairs is not None and not 0<=expected_pairs<=1+1e-9:
        raise AssertionError("Original null expected genus concentration impossible")
    return dict(
        n_original_effective_breadth_matched_species=effective,
        n_original_source_breadth_bins_with_accepted_names=len(used),
        n_original_island_accepted_names_in_matched_bins=n,
        n_original_genus_names_in_complete_source_union=len(G),
        observed_original_shared_genus_count=observed_shared,
        exact_null_expected_original_shared_genus_count=exp_shared,
        shared_original_genus_count_minus_finite_source_expectation=float(
            observed_shared-exp_shared),
        observed_original_source_only_genus_count=observed_source_only,
        exact_null_expected_original_source_only_genus_count=exp_source_only,
        source_only_original_genus_count_minus_finite_source_expectation=float(
            observed_source_only-exp_source_only),
        observed_original_island_only_genus_count=observed_native_only,
        exact_null_expected_original_island_only_genus_count=exp_native_only,
        island_only_original_genus_count_minus_finite_source_expectation=float(
            observed_native_only-exp_native_only),
        observed_original_island_genus_pair_concentration=float(
            observed_pairs) if observed_pairs is not None else None,
        exact_null_expected_original_island_genus_pair_concentration=float(
            expected_pairs) if expected_pairs is not None else None,
        observed_minus_expected_original_island_genus_pair_concentration=float(
            observed_pairs-expected_pairs) if n>=2 else None,
    )


def original_trait_island_eligibility(native_source:set[str],only_source:set[str],
                                      GIFT_listing:dict[str,int],
                                      trait_state:dict[str,int])->bool:
    """Exactly previous original 7-trait n5/n10/effective-breadth gate.

    The individual Direct category is NEVER used to calculate the
    unrestricted botanical genus null: only source trait DOCUMENTATION
    presence determines this same-island cohort mask.
    """
    D=set(trait_state)
    native=native_source&D
    external=only_source&D
    if len(native)<MIN_NATIVE or len(external)<MIN_CANDIDATE:return False
    n=Counter(listing_breadth_bin(GIFT_listing[x]) for x in native)
    m=Counter(listing_breadth_bin(GIFT_listing[x]) for x in external)
    return sum(min(n[k],m[k]) for k in BINS)>=MIN_MATCH


def source_group_summary(records,cluster,*,bootstrap=BOOT,seed=20261115):
    if not records:return dict(status="NO_ORIGINAL_BOTANICAL_SOURCE_COMPARISON",
                                n_islands=0,n_clusters=0)
    df=pd.DataFrame(records).dropna(subset=[cluster,*METRICS])
    n=len(df);g=df[cluster].nunique()
    result=dict(n_islands=n,n_clusters=g)
    if n<MIN_ISLANDS or g<MIN_BLOCKS:
        result["status"]="INSUFFICIENT_ORIGINAL_INDEPENDENT_GEOLOGICAL_BLOCK_SUPPORT"
        return result
    v=df.groupby(cluster,sort=True)[list(METRICS)].mean().to_numpy(float)
    avg=v.mean(0)
    rng=np.random.default_rng(seed)
    draws=v[rng.integers(g,size=(bootstrap,g))].mean(1)
    result["status"]="EXACT_SOURCE_LIST_BREADTH_ALL_BOTANICAL_GENUS_NULL"
    result["original_999_group_resamples_not_independent_plant_histories"]=True
    result["n_group_conditional_resamples"]=bootstrap
    for j,metric in enumerate(METRICS):
        result[metric+"_group_equal_mean"]=float(avg[j])
        result[metric+"_group_conditional_95"]=[
            float(np.quantile(draws[:,j],.025)),
            float(np.quantile(draws[:,j],.975))]
        loo=[v[np.arange(g)!=k,j].mean() for k in range(g)]
        result[metric+"_delete_one_group_min"]=float(min(loo))
        result[metric+"_delete_one_group_max"]=float(max(loo))
    result["mean_source_genus_overlap_actual"]=float(
        df.observed_original_shared_genus_count.mean())
    result["mean_source_genus_overlap_null"]=float(
        df.exact_null_expected_original_shared_genus_count.mean())
    result["mean_original_observed_island_genus_concentration"]=float(
        df.observed_original_island_genus_pair_concentration.mean())
    result["mean_original_expected_island_genus_concentration"]=float(
        df.exact_null_expected_original_island_genus_pair_concentration.mean())
    result["mean_original_accepted_source_species_per_island_in_matched_bins"]=float(
        df.n_original_island_accepted_names_in_matched_bins.mean())
    return result


def analyze(original:pd.DataFrame,geo:pd.DataFrame,assign:pd.DataFrame,
            mainland:pd.DataFrame,original_states:dict,
            *,nboot:int=BOOT):
    flora=status_groups(original)
    if len(flora)!=4505:raise ValueError("Original 4505 botanical island flora population drift")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Independently reconstructed original GSHHG source geographic frame invalid")
    g=geo.set_index("island_id").reindex(sorted(flora))
    if int(g.analysis_regime.notna().sum())!=4504:
        raise ValueError("Original one unresolved geography botanical island changed")
    regime={i:(str(g.at[i,"analysis_regime"]) if pd.notna(g.at[i,"analysis_regime"])
               else "geography_unresolved") for i in flora}
    counts=Counter(regime.values())
    if set(counts)!=set(REGIONS) or sum(counts.values())!=4505:
        raise ValueError("Original entire GIFT biological geography accounting incomplete")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assign):
        raise ValueError("Original candidate native mainland region rankings absent")
    A=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    A.source_rank=pd.to_numeric(A.source_rank,errors="raise")
    A=A.loc[A.source_rank.between(1,5)]
    if A.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original botanical source region rank")
    nearest=A.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if not {"entity_ID","work_species"}<=set(mainland) or mainland.duplicated(
            ["entity_ID","work_species"]).any():
        raise ValueError("Original source mainland native list contains repeated accepted taxa")
    byentity=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    source_breadth=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    sectors=source_sectors(original_states)
    if len(sectors)!=7:raise ValueError("Original source Direct trait sector contract changed")
    panels=[]
    opportunity=[]
    for region in REGIONS:
        ids=sorted(i for i in flora if regime[i]==region)
        rowbags=defaultdict(list)
        studied=Counter()
        for i in ids:
            F=flora[i]
            key=nearest.get(i,())
            potential=set().union(*(byentity.get(z,set()) for z in key))
            source_only=potential-F["full"]
            block=g.at[i,"spatial_block"]
            tags=dict(
                spatial_block=str(block) if pd.notna(block) else None,
                first_source=key[0] if key else None)
            studied["n_original_source_listed_island_floras_with_five_regions"]+=int(len(key)==5)
            for scope in SCOPES:
                if scope=="strict_native_nonendemic_source_listed":
                    on=F["native_nonendemic"]&potential
                elif scope=="all_original_island_species_source_listed":
                    on=F["full"]&potential
                else:
                    on=F["introduced"]&potential
                if len(on)<MIN_NATIVE or len(source_only)<MIN_CANDIDATE:
                    continue
                q=exact_taxonomic_source_null(on,source_only,source_breadth)
                if q is None or q["n_original_effective_breadth_matched_species"]<MIN_MATCH:
                    continue
                v=dict(**tags,**q)
                rowbags[(scope,"all_original_source_taxa")].append(v)
                if scope=="strict_native_nonendemic_source_listed":
                    for trait,states in sectors.items():
                        if original_trait_island_eligibility(
                            on,source_only,source_breadth,states):
                            rowbags[(scope,"same_islands_as_Direct_"+trait)].append(v)
        for scope in SCOPES:
            cohorts=(("all_original_source_taxa",) if scope==SCOPES[2]
                     else ("all_original_source_taxa",*(
                         "same_islands_as_Direct_"+z for z in sectors
                     )))
            for cohort in cohorts:
                data=rowbags.get((scope,cohort),[])
                out=dict(region=region,source_status_scope=scope,
                         island_cohort=cohort,
                         original_botanical_region_floras=len(ids),
                         n_original_source_listed_islands_with_valid_full_name_null=len(data))
                if data:
                    out.update(
                        island_equal_mean_actual_shared_genera=float(np.mean([
                            x["observed_original_shared_genus_count"] for x in data])),
                        island_equal_mean_expected_shared_genera=float(np.mean([
                            x["exact_null_expected_original_shared_genus_count"] for x in data])),
                        island_equal_mean_genus_overlap_minus_null=float(np.mean([
                            x["shared_original_genus_count_minus_finite_source_expectation"] for x in data])),
                        island_equal_mean_native_genus_concentration_minus_null=float(np.mean([
                            x["observed_minus_expected_original_island_genus_pair_concentration"] for x in data])),
                        mean_original_accepted_species_in_matched_bins=float(np.mean([
                            x["n_original_island_accepted_names_in_matched_bins"] for x in data])))
                for group,key in (("spatial_block","original_physical_10degree_block"),
                                  ("first_source","original_first_GIFT_mainland_source")):
                    out[key]=source_group_summary(
                        data,group,bootstrap=nboot,
                        seed=int.from_bytes(hashlib.sha256(
                            (region+scope+cohort+group).encode()).digest()[:8],"big"))
                panels.append(out)
            opportunity.append(dict(
                region=region,source_status_scope=scope,
                n_original_botanical_region_islands=len(ids),
                n_valid_name_only_source_floras=len(
                    rowbags.get((scope,"all_original_source_taxa"),[]))))
    for x in ("tropical","northern_midlatitude"):
        wanted=[
            (x,"strict_native_nonendemic_source_listed",
             "same_islands_as_Direct_zygomorphic"),
            (x,"strict_native_nonendemic_source_listed",
             "same_islands_as_Direct_tube_or_bell")]
        expected={("tropical","zygomorphic"):98,
                  ("tropical","tube_or_bell"):88,
                  ("northern_midlatitude","zygomorphic"):230,
                  ("northern_midlatitude","tube_or_bell"):224}
        for reg,scope,cohort in wanted:
            actual=sum(len(rowbags) for rowbags in []) if False else next(
                p["n_original_source_listed_islands_with_valid_full_name_null"]
                for p in panels if p["region"]==reg and
                p["source_status_scope"]==scope and p["island_cohort"]==cohort)
            if actual!=expected[(reg,cohort.replace("same_islands_as_Direct_",""))]:
                raise AssertionError("Original source trait-qualified SAME-island cohort mask changed")
    if sum(counts.values())!=4505:
        raise AssertionError("Global full flora taxonomic genus source source population changed")
    return dict(
        status="ORIGINAL_4505_ALL_BOTANICAL_NAME_FINITE_SOURCE_GENUS_NULL",
        n_original_source_island_botanical_floras=4505,
        n_original_source_native_GIFT_mainland_region_entities=mainland.entity_ID.nunique(),
        n_original_separate_existing_Direct_sector_island_masks=len(sectors),
        n_original_geographic_population_regimes=len(REGIONS),
        original_4505_botanical_source_region_population=[
            dict(region=r,n_original_botanical_island_floras=int(counts[r]))
            for r in REGIONS],
        original_status_source_name_null_opportunity=opportunity,
        original_botanical_all_source_genus_exact_finite_null_panels=panels,
        claims=dict(
            original_genus_null_calculated_without_original_flower_or_reproduction_category_values=True,
            same_islands_as_previous_source_Direct_matched_cohorts=True,
            source_rare_vs_common_original_GIFT_listing_breadth_counts_fixed=True,
            native_vs_allrecorded_original_status_are_nested_not_independent=True,
            source_original_island_accepted_name_presence_not_seeding_arrival=True,
            exact_matched_source_taxonomic_sampling_without_replacement=True,
            expected_original_genus_richness_and_pair_concentration_exact=True,
            full_4505_botanical_population_not_928_three_trait_complete_subsample=True,
            conditional_10degree_GIFT_source_group_resamples_not_independent_evolution=True,
            source_native_status_unresolved_taxa_not_assumed_indigenous=True,
            source_potential_genus_overlap_not_proof_of_failed_establishment=True,
            no_original_private_plant_genus_island_source_entity_or_coordinates_exported=True))


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
    source=dict(
        flora=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        traits=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        mainland=args.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        assign=args.original_root/"fixed/source/result/island_source_assignments.csv.gz")
    hashes=dict(
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        mainland="3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")
    for k,path in source.items():
        if sha(path)!=hashes[k]:
            raise ValueError("Original GIFT botanical source SHA drift: "+k)
    states=strict_traits(pd.read_csv(source["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    original=pd.read_csv(source["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    geo=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assign=pd.read_csv(source["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(source["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    if geo.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original GSHHG and native GIFT source entity universe changed")
    result=analyze(original,geo,assign,mainland,states)
    if result["n_original_source_island_botanical_floras"]!=4505:
        raise ValueError("Original global island botanical source universe changed")
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for priv in ('"island_id"','"accepted_species"','"entity_ID"',
                 '"work_species"','"spatial_block"','"source_signature"',
                 '"island_latitude"','"island_longitude"'):
        if priv in body:raise ValueError("Protected source species/island entity identity exported")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n")
    print(json.dumps(dict(status=result["status"],
                          n_original_botanical_floras=4505,
                          n_regional_status_source_genus_panels=len(result[
                              "original_botanical_all_source_genus_exact_finite_null_panels"]),
                          no_actual_seed_arrival_or_adaptive_evolution_inferred=True)))


if __name__=="__main__":main()
