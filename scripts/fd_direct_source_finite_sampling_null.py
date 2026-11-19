"""Original GIFT 4505 island native-source floral composition: exact finite source null.

This is NOT a genus prediction model. For each actual source-breadth-matched
island native_nonendemic versus absent-from-checklist GIFT mainland candidate
floral sector, test:
  A. source taxon-label exchangeability conditional on the actual number of
     native and candidate DIRECT species in each original GIFT listing
     breadth bin, including EXACT expected genus richness/presence;
  B. source taxon-label exchangeability additionally conditional on observed
     local-vs-candidate genus incidence in each matched breadth bin.
Both nulls condition on real original high/medium species-DIRECT category
documentation. Original distinct accepted plant identity is NOT an observed
seed arrival, independently observed establishment, historical ancestor,
adaptive island morphology or an effective pollinator interaction.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter,defaultdict
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_direct_three_stage_floristic_comparison import (
    genus as _parse_genus,source_sectors,status_groups
)
from fd_direct_source_breadth_matched import (
    BINS,REGIONS,listing_breadth_bin,breadth_matched_sector_difference
)

MIN_NATIVE=5
MIN_CANDIDATE=10
MIN_EFFECTIVE=5
MIN_ISLANDS=30
MIN_BLOCKS=10
BOOT=999
@lru_cache(maxsize=120000)
def genus(sp:str)->str:
    return _parse_genus(sp)


METRICS=(
    "original_direct_native_minus_candidate_category_difference",
    "expected_direct_difference_genus_fixed_original_source_null",
    "actual_minus_genus_fixed_original_source_null",
    "observed_minus_random_expected_native_only_genus_count",
    "observed_minus_random_expected_candidate_only_genus_count",
    "observed_minus_random_expected_shared_genus_count",
)


def hypergeom_count_variance(N:int,K:int,n:int)->float:
    """Var actual number of 1-states in n names drawn without replacement."""
    if not all(isinstance(z,(int,np.integer)) for z in (N,K,n)):
        raise ValueError("Finite original GIFT source counts must be integers")
    if N<0 or not 0<=K<=N or not 0<=n<=N:
        raise ValueError("Invalid original candidate source counts")
    if N<=1 or n==0 or n==N or K in (0,N):
        return 0.
    p=K/N
    return float(n*(N-n)/(N-1)*p*(1-p))


@lru_cache(maxsize=250000)
def prob_no_original_genus_in_native_or_candidate(
        N:int,n_native:int,m_genus:int)->tuple[float,float]:
    """P(native has NONE of genus), P(candidate has NONE); exact finite pool."""
    if not all(isinstance(z,(int,np.integer)) for z in (
            N,n_native,m_genus)) or not (0<=n_native<=N and
            0<=m_genus<=N):
        raise ValueError("Invalid original source genus-count opportunity")
    n_candidate=N-n_native
    if m_genus==0:return 1.,1.
    def source_side_none(side_count):
        if m_genus>side_count:return 0.
        result=1.
        for i in range(m_genus):
            result*=(side_count-i)/(N-i)
        return float(result)
    return (source_side_none(n_candidate),
            source_side_none(n_native))


def exact_original_breadth_null(native:set[str],candidate:set[str],
                                scored:dict[str,int],
                                listed_count:dict[str,int]):
    """All 4 EXACT genus/non-genus controls on SAME REAL source-balanced island."""
    if not native or not candidate or native&candidate:
        raise ValueError("No distinct original candidate/native names")
    if not (native|candidate)<=set(scored) or not (
            native|candidate)<=set(listed_count):
        raise ValueError("Original source Direct label or GIFT breadth absent")
    pools={b:{"native":set(),"candidate":set()} for b in BINS}
    for s in native:
        pools[listing_breadth_bin(listed_count[s])]["native"].add(s)
    for s in candidate:
        pools[listing_breadth_bin(listed_count[s])]["candidate"].add(s)
    bins=[]
    actual_genus_native=set()
    actual_genus_candidate=set()
    pool_total_genus=set()
    p_none_native=defaultdict(lambda:1.)
    p_none_candidate=defaultdict(lambda:1.)
    weighted_sum=0.
    expected_g=0.
    var_a=0.
    var_g=0.
    effective=0
    genus_variable_shared=set()
    for breadth,part in pools.items():
        a,b=part["native"],part["candidate"]
        if not a or not b:continue
        n,s=len(a),len(b)
        N=n+s
        names=a|b
        weight=min(n,s)
        actual=float(np.mean([scored[x] for x in a])-
                     np.mean([scored[x] for x in b]))
        K=sum(scored[x] for x in names)
        assert 0<=K<=N
        var_unconditional=hypergeom_count_variance(N,int(K),n)*(1/n+1/s)**2
        bygenus=defaultdict(lambda:dict(native=[],candidate=[]))
        for x in a:bygenus[genus(x)]["native"].append(scored[x])
        for x in b:bygenus[genus(x)]["candidate"].append(scored[x])
        expectation=0.
        variance=0.
        for g,v in bygenus.items():
            x,y=v["native"],v["candidate"]
            m=len(x)+len(y)
            k=sum(x)+sum(y)
            # Fixed source genus incidence, random accepted names and
            # their REAL original 0/1 category state within that genus.
            expectation+=(len(x)/n-len(y)/s)*(k/m)
            variance+=hypergeom_count_variance(
                m,int(k),len(x))*(1/n+1/s)**2
            qnat,qcand=prob_no_original_genus_in_native_or_candidate(
                N,n,m)
            p_none_native[g]*=qnat
            p_none_candidate[g]*=qcand
            pool_total_genus.add(g)
            if len(x) and len(y) and 0<k<m:
                genus_variable_shared.add(g)
        weighted_sum+=weight*actual
        expected_g+=weight*expectation
        var_a+=weight*weight*var_unconditional
        var_g+=weight*weight*variance
        effective+=weight
        actual_genus_native|={genus(x) for x in a}
        actual_genus_candidate|={genus(x) for x in b}
        bins.append(dict(
            bin=breadth,n_original_native_Direct=n,
            n_original_candidate_Direct=s,
            n_original_source_unique_Direct_species=N,
            n_documented_category_positive_source_species=int(K),
            observed_direct_difference=actual,
            genus_fixed_conditional_expected_difference=float(expectation),
            conditional_random_native_source_variance=float(var_unconditional),
            conditional_genus_preserving_source_variance=float(variance)))
    if not effective:return None
    matched=weighted_sum/effective
    expected_genus=expected_g/effective
    gen_total=len(pool_total_genus)
    E_source_only=float(sum(p_none_native[g] for g in pool_total_genus))
    E_native_only=float(sum(p_none_candidate[g] for g in pool_total_genus))
    E_shared=float(gen_total-E_source_only-E_native_only)
    observed_shared=len(actual_genus_native&actual_genus_candidate)
    observed_native_only=len(actual_genus_native-actual_genus_candidate)
    observed_candidate_only=len(actual_genus_candidate-actual_genus_native)
    expected_native_richness=gen_total-E_source_only
    expected_source_richness=gen_total-E_native_only
    if abs(E_shared+E_native_only+E_source_only-gen_total)>1e-9:
        raise AssertionError("Exact genus source availability richness null fails")
    if observed_shared+observed_native_only+observed_candidate_only!=gen_total:
        raise AssertionError("Original actual shared/exclusive genus states inconsistent")
    original=breadth_matched_sector_difference(
        native,candidate,scored,listed_count)
    if not original or abs(matched-original[
        "source_breadth_matched_native_minus_candidate_sector_share"])>1e-12:
        raise AssertionError("Original previously verified direct flora source difference changed")
    if not np.isfinite([var_a,var_g,expected_genus,matched]).all():
        raise AssertionError("Nonfinite original direct source category null")
    return dict(
        original_direct_native_minus_candidate_category_difference=float(matched),
        nullA_exact_exchangeable_source_trait_expected_difference=0.,
        nullA_exact_exchangeable_source_trait_sampling_sd=float(
            np.sqrt(var_a)/effective),
        expected_direct_difference_genus_fixed_original_source_null=float(
            expected_genus),
        nullB_exact_genus_fixed_original_source_trait_sampling_sd=float(
            np.sqrt(var_g)/effective),
        actual_minus_genus_fixed_original_source_null=float(
            matched-expected_genus),
        n_original_effectively_breadth_matched_Direct_species=int(effective),
        n_original_matched_breadth_bins=len(bins),
        n_original_Direct_source_names_in_matched_bins=int(
            sum(v["n_original_source_unique_Direct_species"] for v in bins)),
        n_original_actual_source_genera_with_variable_species_categories_and_both_memberships=len(
            genus_variable_shared),
        observed_native_genus_richness=len(actual_genus_native),
        nullA_expected_native_genus_richness=float(expected_native_richness),
        observed_candidate_genus_richness=len(actual_genus_candidate),
        nullA_expected_candidate_genus_richness=float(expected_source_richness),
        observed_native_exclusive_genus_count=observed_native_only,
        nullA_expected_native_exclusive_genus_count=E_native_only,
        observed_candidate_exclusive_genus_count=observed_candidate_only,
        nullA_expected_candidate_exclusive_genus_count=E_source_only,
        observed_shared_genus_count=observed_shared,
        nullA_expected_shared_genus_count=E_shared,
        observed_minus_random_expected_native_only_genus_count=float(
            observed_native_only-E_native_only),
        observed_minus_random_expected_candidate_only_genus_count=float(
            observed_candidate_only-E_source_only),
        observed_minus_random_expected_shared_genus_count=float(
            observed_shared-E_shared),
        full_exact_real_original_source_category_breadth_bins=bins,
    )


def aggregate_by_original_clusters(rows:list[dict],group:str,*,bootstrap=BOOT,
                                   seed=20261115):
    if not rows:
        return dict(status="NOT_IDENTIFIABLE_NO_ORIGINAL_DIRECT_MATCHES",
                    n_islands=0,n_groups=0)
    df=pd.DataFrame(rows).dropna(subset=[group,*METRICS])
    n=len(df)
    g=df[group].nunique()
    out=dict(n_islands=n,n_groups=g,
             island_equal_mean_true_direct_native_source_category_difference=(
                 float(df[METRICS[0]].mean()) if n else None))
    if n<MIN_ISLANDS or g<MIN_BLOCKS:
        return dict(**out,status="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_GROUPS")
    groups=df.groupby(group,sort=True)[list(METRICS)].mean().to_numpy(float)
    mean=groups.mean(0)
    rng=np.random.default_rng(seed)
    resampled=groups[rng.integers(len(groups),
        size=(bootstrap,len(groups)))].mean(1)
    loo=[groups[np.arange(len(groups))!=j].mean(0) for j in range(len(groups))]
    L=np.asarray(loo)
    out["status"]="SOURCE_BREADTH_FINITE_NULL_OBSERVATIONAL_CLUSTERED_DIFFERENCE"
    out["n_conditional_original_source_group_resamples"]=bootstrap
    for k,col in enumerate(METRICS):
        out[col+"_group_equal_mean"]=float(mean[k])
        out[col+"_original_group_resample_95"]=[
            float(np.quantile(resampled[:,k],.025)),
            float(np.quantile(resampled[:,k],.975))]
        out[col+"_delete_one_group_min"]=float(L[:,k].min())
        out[col+"_delete_one_group_max"]=float(L[:,k].max())
    out["mean_original_trait_annotation_coverage_native"]=float(
        df.original_trait_annotation_coverage_native.mean())
    out["mean_original_trait_annotation_coverage_mainland_candidate"]=float(
        df.original_trait_annotation_coverage_candidate.mean())
    out["n_islands_with_at_least_one_distinct_variable_shared_genus"]=int(
        (df.n_original_actual_source_genera_with_variable_species_categories_and_both_memberships>0).sum())
    out["mean_distinct_variable_shared_genera_per_island"]=float(
        df.n_original_actual_source_genera_with_variable_species_categories_and_both_memberships.mean())
    out["mean_original_nullA_exact_trait_sampling_sd_per_island"]=float(
        df.nullA_exact_exchangeable_source_trait_sampling_sd.mean())
    out["mean_original_nullB_exact_trait_sampling_sd_per_island"]=float(
        df.nullB_exact_genus_fixed_original_source_trait_sampling_sd.mean())
    out["resampled_source_blocks_not_independent_biological_colonization_trials"]=True
    return out


def run(states,flora,geo,assign,mainland,*,bootstrap=BOOT):
    if not {"island_id","accepted_species","origin_status",
            "endemic_status","floristic_status"}<=set(flora):
        raise ValueError("Incomplete original source accepted status labels")
    original=status_groups(flora)
    if len(original)!=4505:
        raise ValueError("GIFT original botanical 4505 population changed")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Independent GSHHG original 8264 physical island crosswalk absent")
    g=geo.set_index("island_id").reindex(sorted(original))
    if int(g.analysis_regime.notna().sum())!=4504:
        raise ValueError("One original geography-unresolved accepted island requirement changed")
    region={
        i:(str(g.at[i,"analysis_regime"]) if pd.notna(g.at[i,"analysis_regime"])
           else "geography_unresolved")
        for i in original}
    counts=Counter(region.values())
    if sum(counts.values())!=4505 or set(counts)!=set(REGIONS):
        raise ValueError("Full original four region + singleton botanical population drift")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assign):
        raise ValueError("Original 5-nearest native GIFT mainland sources missing")
    ranks=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    ranks.source_rank=pd.to_numeric(ranks.source_rank,errors="raise")
    ranks=ranks.loc[ranks.source_rank.between(1,5)]
    if ranks.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original native source entity rank")
    nearest=ranks.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if not {"entity_ID","work_species"}<=set(mainland) or mainland.duplicated(
            ["entity_ID","work_species"]).any():
        raise ValueError("Original GIFT region×source accepted species roster corrupt")
    bysource=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    breadth=mainland.groupby("work_species").entity_ID.nunique().to_dict()
    sectors=source_sectors(states)
    if len(sectors)!=7:raise ValueError("Original seven biological source sectors changed")
    panels=[]
    censuses=[]
    for r in REGIONS:
        ids=sorted(i for i in original if region[i]==r)
        for trait,phenotypes in sectors.items():
            direct=set(phenotypes)
            records=[]
            status=Counter()
            for i in ids:
                original_i=original[i]
                five=nearest.get(i,())
                source=set().union(*(bysource.get(x,set()) for x in five))
                native_source=source&original_i["native_nonendemic"]
                candidate_source=source-original_i["full"]
                native=native_source&direct
                candidate=candidate_source&direct
                status["n_source_islands_with_all_five_ranked_mainland_regions"]+=int(
                    len(five)==5)
                status["n_original_native_name_status_source_listed"]+=len(native_source)
                status["n_original_candidate_source_only_names"]+=len(candidate_source)
                status["n_original_native_source_Direct_carriers"]+=len(native)
                status["n_original_candidate_only_source_Direct_carriers"]+=len(candidate)
                if len(native)<MIN_NATIVE or len(candidate)<MIN_CANDIDATE:
                    continue
                status["n_original_islands_n5_native_n10_source_Direct"]+=1
                q=exact_original_breadth_null(native,candidate,phenotypes,breadth)
                if q is None or q["n_original_effectively_breadth_matched_Direct_species"]<MIN_EFFECTIVE:
                    continue
                loc=g.at[i,"spatial_block"]
                block=str(loc) if pd.notna(loc) else None
                rows=q.copy()
                rows.pop("full_exact_real_original_source_category_breadth_bins",None)
                rows["original_trait_annotation_coverage_native"]=len(native)/len(native_source)
                rows["original_trait_annotation_coverage_candidate"]=len(candidate)/len(candidate_source)
                rows["original_annotation_coverage_difference_absolute"]=abs(
                    rows["original_trait_annotation_coverage_native"]-
                    rows["original_trait_annotation_coverage_candidate"])
                rows["spatial_block"]=block
                rows["source_cluster"]=five[0] if len(five) else None
                records.append(rows)
            panel=dict(region=r,trait_sector=trait,
                original_region_population_floras=len(ids),
                n_original_source_native_n5_candidate_n10_islands=status[
                    "n_original_islands_n5_native_n10_source_Direct"],
                n_original_breadth_matched_source_floras=len(records),
                n_original_native_source_listed_all_traits_unknown_name_rows=status[
                    "n_original_native_name_status_source_listed"],
                n_original_candidate_source_only_all_traits_unknown_name_rows=status[
                    "n_original_candidate_source_only_names"],
                n_original_native_source_Direct_island_plant_rows=status[
                    "n_original_native_source_Direct_carriers"],
                n_original_candidate_source_Direct_island_plant_rows=status[
                    "n_original_candidate_only_source_Direct_carriers"])
            if records:
                df=pd.DataFrame(records)
                panel.update(
                    n_original_islands_with_Direct_annotation_coverage_gap_le5pp=int(
                        (df.original_annotation_coverage_difference_absolute<=.05).sum()),
                    n_original_islands_with_Direct_annotation_coverage_gap_le10pp=int(
                        (df.original_annotation_coverage_difference_absolute<=.10).sum()),
                    original_island_equal_mean_exact_native_minus_candidate_fraction=float(
                        df.original_direct_native_minus_candidate_category_difference.mean()),
                    original_island_equal_mean_nullB_genus_fixed_expectation=float(
                        df.expected_direct_difference_genus_fixed_original_source_null.mean()),
                    original_island_equal_mean_residual_after_genus_null=float(
                        df.actual_minus_genus_fixed_original_source_null.mean()),
                    original_island_equal_mean_observed_shared_genera=float(
                        df.observed_shared_genus_count.mean()),
                    original_island_equal_mean_expected_shared_genera_under_source_random=float(
                        df.nullA_expected_shared_genus_count.mean()),
                    original_island_equal_mean_observed_native_only_genera=float(
                        df.observed_native_exclusive_genus_count.mean()),
                    original_island_equal_mean_expected_native_only_genera_under_source_random=float(
                        df.nullA_expected_native_exclusive_genus_count.mean()),
                    original_island_equal_mean_observed_source_only_genera=float(
                        df.observed_candidate_exclusive_genus_count.mean()),
                    original_island_equal_mean_expected_source_only_genera_under_source_random=float(
                        df.nullA_expected_candidate_exclusive_genus_count.mean()),
                    mean_original_native_Direct_annotation_coverage=float(
                        df.original_trait_annotation_coverage_native.mean()),
                    mean_original_source_candidate_Direct_annotation_coverage=float(
                        df.original_trait_annotation_coverage_candidate.mean()),
                    mean_original_matched_breadth_effective_Direct_carriers=float(
                        df.n_original_effectively_breadth_matched_Direct_species.mean()),
                    mean_original_actual_state_variable_shared_genera=float(
                        df.n_original_actual_source_genera_with_variable_species_categories_and_both_memberships.mean()))
                for label,subset in (
                    ("all_matched_direct_source",records),
                    ("recording_coverage_gap_le5pp",[
                        x for x in records if
                        x["original_annotation_coverage_difference_absolute"]<=.05]),
                    ("recording_coverage_gap_le10pp",[
                        x for x in records if
                        x["original_annotation_coverage_difference_absolute"]<=.10])):
                    pane={}
                    for group,group_label in (
                        ("spatial_block","physical_10degree_block"),
                        ("source_cluster","first_GIFT_mainland_source")):
                        pane[group_label]=aggregate_by_original_clusters(
                            subset,group,bootstrap=bootstrap,
                            seed=int.from_bytes(hashlib.sha256(
                                (r+trait+label+group_label).encode()).digest()[:8],
                                "big"))
                    pane["n_original_recorded_island_floras"]=len(subset)
                    panel[label]=pane
            else:
                panel["status"]="NOT_IDENTIFIABLE_NO_ORIGINAL_MATCHED_SOURCE_TRAIT"
            panels.append(panel)
            censuses.append(dict(
                region=r,trait_sector=trait,
                n_original_botanical_region_floras=len(ids),
                n_original_source_listed_native_direct_source_carrier_rows=status[
                    "n_original_native_source_Direct_carriers"],
                n_original_candidate_mainland_direct_source_carrier_rows=status[
                    "n_original_candidate_only_source_Direct_carriers"]))
    if len(panels)!=35 or sum(counts.values())!=4505:
        raise AssertionError("Original 7 sector×5 geography full botanical source census lost")
    return dict(
        status="ORIGINAL_4505_DIRECT_BREADTH_FINITE_SOURCE_TAXONOMIC_NULLS",
        n_original_GIFT_botanical_island_floras=4505,
        n_real_original_direct_trait_sectors=7,
        n_original_native_GIFT_mainland_regions=mainland.entity_ID.nunique(),
        original_4505_botanical_GSHHG_region_population=[
            dict(region=r,n_original_botanical_islands=int(counts[r]))
            for r in REGIONS],
        original_full_source_region_trait_opportunity=censuses,
        original_finite_taxonomic_source_null_panels=panels,
        claims=dict(
            baseline_original_breadth_matched_island_difference_reconciled_exactly=True,
            source_NullA_exact_finite_species_draws_preserve_n_and_GIFT_breadth=True,
            source_NullB_preserves_same_real_direct_category_labels_and_genus_counts=True,
            source_nulls_condition_on_high_medium_Direct_annotation_not_unrecorded_trait_values=True,
            original_trait_recording_coverage_gaps_explicit_and_predefined_sensitivities=True,
            original_shared_exclusive_genus_null_exact_finite_probabilities=True,
            original_global_4505_botanical_population_retained=True,
            no_original_endemic_or_intro_inferred_from_native_status=True,
            source_candidates_not_actual_seed_arrival_or_biological_absence=True,
            no_ancestral_trait_transition_or_reproductive_plant_fitness=True,
            original_source_group_resampling_not_independent_plant_evolution=True,
            original_private_taxon_genus_island_GIFT_entity_or_location_not_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda:f.read(1<<20),b""):
            h.update(part)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for x in ("original-root","geology","output"):
        p.add_argument("--"+x,type=Path,required=True)
    a=p.parse_args()
    names={
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "mainland":a.original_root/"fixed/source/source/gift_native_mainland_flora.csv.gz",
        "assign":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "mainland":"3f106d1dfecfc4f79d8421999dcdb30e8d31992712f33e4aa133f3efdf8e2cef",
        "assign":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for key,path in names.items():
        if sha(path)!=pins[key]:
            raise ValueError("Original GIFT source SHA drift "+key)
    states=strict_traits(pd.read_csv(names["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(names["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status","endemic_status","floristic_status"])
    geo=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assign=pd.read_csv(names["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    mainland=pd.read_csv(names["mainland"],dtype=str,usecols=[
        "entity_ID","work_species"])
    if geo.island_id.nunique()!=8264 or mainland.entity_ID.nunique()!=748:
        raise ValueError("Original GIFT native mainland and independent source GSHHG changed")
    result=run(states,flora,geo,assign,mainland)
    if result["n_original_GIFT_botanical_island_floras"]!=4505 or len(
        result["original_finite_taxonomic_source_null_panels"])!=35:
        raise AssertionError("Original source global botanical/trait universe drift")
    body=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"work_species"',
                    '"entity_ID"','"spatial_block"','"source_signature"',
                    '"island_latitude"','"island_longitude"'):
        if private in body:
            raise ValueError("Private original GIFT taxon/plant geography source leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n")
    print(json.dumps(dict(
        status=result["status"],n_original_GIFT_botanical_island_floras=4505,
        n_exact_taxonomic_null_panels=len(result[
            "original_finite_taxonomic_source_null_panels"]),
        independent_seed_arrival_or_adaptive_evolution_claim=False)))


if __name__=="__main__":main()
