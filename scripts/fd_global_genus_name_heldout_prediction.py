"""Full-original-flora-name-heldout genus conditional prediction of floral category assemblages.

Original source: all 4505 GIFT botanical floras and 10 Direct traits + 3
true observed joint pairs. Each tested island is removed by EVERY accepted
plant species NAME appearing anywhere in its FULL original botanical
checklist, not just test-carrier species. Global donor = other original Direct
species; regional donor = different original Direct species recorded in other
source floras of same broad region. Other islands may share donors; this is
species-NAME heldout within the same trait ledger, NOT independent phylogeny,
mainland-specific phenotypes or independent island histories.

Known original genus identities and documented genus counts at test island
form a conditional oracle. Compare squared observed CATEGORY FREQUENCY
prediction loss to the exact SAME leave-out donor genus-blind source prior.
Unseen source genus returns prior. 5 pseudocount shrinkage fixed before result.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits,TRAITS,PAIRS
from fd_global_genus_assembly_geometry import source_endpoints

REGIONS=("northern_midlatitude","tropical","southern_extratropical")
SCOPES=("all_recorded","original_native")
MODE=("regional_source_only","global_Direct_ledger")
ALPHA=5.0
MIN_DIRECT=5
MIN_BLOCKS=10
MIN_ISLANDS=30
N_BOOT=999


def genus(name):
    terms=str(name).split()
    return terms[0] if len(terms)>=2 else "fixed_species:"+str(name)


def compile_donors(state:dict,source:set,category_to_idx:dict):
    """One source-Direct observation per distinct original accepted species name."""
    valid=set(state)&set(source)
    K=len(category_to_idx)
    total=np.zeros(K,dtype=int)
    by_genus=defaultdict(lambda:np.zeros(K,dtype=int))
    encoded={}
    for sp in valid:
        state_idx=category_to_idx[state[sp]]
        total[state_idx]+=1
        g=genus(sp)
        by_genus[g][state_idx]+=1
        encoded[sp]=(g,int(state_idx))
    if int(total.sum())!=len(valid):
        raise AssertionError("Original Direct source donor duplicated")
    return dict(names=valid,n=len(valid),total=total,genus=by_genus,
                code=encoded,n_categories=K,category_to_idx=category_to_idx,
                state_names=set(state))


def predict_original_flora(accepted_entire_full_flora:set,
                           directly_documented_carriers:set,
                           state:dict,donor:dict,
                           *,alpha:float=ALPHA):
    """Exclude every focal accepted source plant NAME from the donor universe."""
    if not directly_documented_carriers or not directly_documented_carriers<=donor["state_names"]:
        raise ValueError("True Direct states required, no imputed source traits")
    if not directly_documented_carriers<=accepted_entire_full_flora:
        raise ValueError("Original Direct target carrier not in original full flora")
    if alpha<=0:raise ValueError("Prior shrinkage must be positive")
    K=donor["n_categories"]
    cat_index=donor["category_to_idx"]
    observed=np.zeros(K)
    for s in directly_documented_carriers:
        observed[cat_index[state[s]]]+=1
    observed/=len(directly_documented_carriers)
    excluded=accepted_entire_full_flora&donor["names"]
    counts=donor["total"].astype(float).copy()
    excluded_genus=defaultdict(lambda:np.zeros(K))
    for sp in excluded:
        g,cat=donor["code"][sp]
        counts[cat]-=1
        excluded_genus[g][cat]+=1
    after=int(counts.sum())
    if after<1:return None
    if (counts<-1e-9).any() or after!=donor["n"]-len(excluded):
        raise AssertionError("Original focal botanical name still in donor training")
    prior=(counts+.5)/(counts.sum()+.5*K)
    target_genus=Counter(genus(x) for x in directly_documented_carriers)
    mean=np.zeros(K)
    n_congener=0
    n_congener_variable=0
    for g,n in target_genus.items():
        source_count=(donor["genus"].get(g,np.zeros(K)).astype(float)-
                      excluded_genus[g])
        if np.any(source_count<-.000000001):
            raise AssertionError("Source genus donor-name subtraction underflow")
        d=max(0.,source_count.sum())
        if d>=1:
            n_congener+=n
            if np.count_nonzero(source_count)>1:
                n_congener_variable+=n
        q=(source_count+alpha*prior)/(d+alpha)
        mean+=n*q
    mean/=len(directly_documented_carriers)
    if not np.allclose(mean.sum(),1,atol=1e-10) or (
            np.any(mean<-.000000001)):
        raise AssertionError("Invalid category frequency source model")
    blind=float(np.square(observed-prior).sum())
    predicted=float(np.square(observed-mean).sum())
    return dict(
        blind_category_frequency_D2=blind,
        name_heldout_genus_category_frequency_D2=predicted,
        name_heldout_genus_prediction_gain_D2=blind-predicted,
        n_direct_carriers=len(directly_documented_carriers),
        n_all_original_flora_names_excluded=len(excluded),
        n_unique_external_same_ledger_species_donors_after_name_holdout=after,
        share_carriers_with_independent_NAME_genus_donor=(
            n_congener/len(directly_documented_carriers)),
        share_carriers_with_at_least_two_distinct_genus_profiles=(
            n_congener_variable/len(directly_documented_carriers)),
        is_better_than_genus_blind_prior=bool(predicted<blind),
    )


def block_bootstrap_paired_prediction(frame:pd.DataFrame,*,group:str,
                                      n_boot:int=N_BOOT,seed:int=20261115):
    """Group-equal loss difference and paired-group resampling sensitivity."""
    cols=["name_heldout_genus_prediction_gain_D2",
          "blind_category_frequency_D2",
          "name_heldout_genus_category_frequency_D2"]
    required=set(cols)|{group}
    if not required<=set(frame):raise ValueError("Original group prediction incomplete")
    df=frame.dropna(subset=list(required))
    if len(df)<MIN_ISLANDS or df[group].nunique()<MIN_BLOCKS:
        return dict(status="INSUFFICIENT_ORIGINAL_INDEPENDENT_SOURCE_CLUSTERS",
                    n_scored_islands=len(df),n_groups=df[group].nunique())
    groups=df.groupby(group,sort=True)[cols].mean()
    values=groups.to_numpy(float)
    mean=values.mean(0)
    rng=np.random.default_rng(seed)
    index=rng.integers(len(values),size=(n_boot,len(values)))
    draws=values[index].mean(1)
    loo=np.array([values[np.arange(len(values))!=i].mean(0)
                  for i in range(len(values))])
    if len(values)<2:raise AssertionError("Cannot leave one source group out")
    return dict(
        status="ORIGINAL_GROUP_EQUAL_NAME_HELDOUT_GENUS_PREDICTION",
        n_scored_islands=len(df),n_groups=len(values),
        n_original_number_resampling_draws=n_boot,
        genus_vs_prior_group_equal_D2_gain=float(mean[0]),
        genus_vs_prior_group_equal_gain_95=[
            float(np.quantile(draws[:,0],.025)),
            float(np.quantile(draws[:,0],.975))],
        genus_vs_prior_group_equal_gain_leave_one_group_min=float(loo[:,0].min()),
        genus_vs_prior_group_equal_gain_leave_one_group_max=float(loo[:,0].max()),
        source_conditional_group_equal_relative_frequency_loss_gain=(
            float(mean[0]/mean[1]) if mean[1]>1e-12 else None),
        group_equal_genus_model_loss=float(mean[2]),
        group_equal_blind_model_loss=float(mean[1]),
        bootstrap_draws_not_independent_biological_replicates=True)


def analyse(states:dict,flora:pd.DataFrame,geo:pd.DataFrame,assign:pd.DataFrame,
            *,n_boot=N_BOOT):
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Full original flora species and status columns needed")
    if flora.duplicated(["island_id","accepted_species"]).any() or (
        not flora.origin_status.isin(("native","introduced","unresolved")).all()):
        raise ValueError("Original species member duplicates/status invalid")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo):
        raise ValueError("Original audited geography missing")
    if geo.island_id.duplicated().any():
        raise ValueError("Multiple GSHHG locations for same original island")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assign):
        raise ValueError("Original GIFT candidate source rank missing")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("Full 4505 original botanical island universe drift")
    geo=geo.set_index("island_id").reindex(sorted(full))
    source=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise")
    if source.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original mainland source rank")
    first=source.loc[source.source_rank.eq(1)].set_index(
        "island_id").entity_ID.to_dict()
    ep=source_endpoints(states)
    output=[]
    coverage={}
    opportunity=[]
    for scope,lookup in (("all_recorded",full),("original_native",native)):
        for endpoint,trait in ep.items():
            names=set(trait)
            coverage[scope+"__"+endpoint]=sum(
                len(lookup.get(i,set())&names)>=MIN_DIRECT for i in full)
            category={c:j for j,c in enumerate(sorted(set(trait.values()),key=repr))}
            global_donor=compile_donors(trait,names,category)
            for region in REGIONS:
                source_ids=[i for i in sorted(full)
                            if geo.at[i,"analysis_regime"]==region]
                pool=set().union(*(lookup.get(i,set()) for i in source_ids))&names
                regional_donor=compile_donors(trait,pool,category)
                for mode in MODE:
                    donor=(regional_donor if mode=="regional_source_only" else
                           global_donor)
                    scored=[]
                    for island in source_ids:
                        target=lookup.get(island,set())&names
                        if len(target)<MIN_DIRECT:continue
                        z=predict_original_flora(
                            full[island],target,trait,donor)
                        if z is None:continue
                        block=geo.at[island,"spatial_block"]
                        scored.append(dict(
                            **z,physical_block=str(block) if pd.notna(block) else None,
                            source_cluster=first.get(island) if first.get(island) else None))
                    opportunity.append(dict(scope=scope,trait_endpoint=endpoint,
                        region=region,reference_donor_mode=mode,
                        n_original_botanical_region_islands=len(source_ids),
                        n_n5_trait_opportunity_islands=sum(
                            len(lookup.get(i,set())&names)>=MIN_DIRECT for i in source_ids),
                        n_scored_source_original_islands=len(scored),
                        n_original_distinct_reference_Direct_donor_species=donor["n"],
                        n_original_reference_genera=len(donor["genus"])))
                    if not scored:continue
                    private=pd.DataFrame(scored)
                    gain=private.name_heldout_genus_prediction_gain_D2
                    mean_blind=float(private.blind_category_frequency_D2.mean())
                    row=dict(scope=scope,trait_endpoint=endpoint,region=region,
                        reference_donor_mode=mode,
                        n_original_source_floras_with_name_heldout_scores=len(private),
                        n_reference_direct_species_before_focal_name_removal=donor["n"],
                        n_original_islands_genus_prediction_better_than_prior=int(
                            (gain>0).sum()),
                        fraction_original_islands_genus_prediction_better_than_prior=float(
                            (gain>0).mean()),
                        island_equal_mean_name_heldout_genus_vs_prior_D2_gain=float(
                            gain.mean()),
                        island_equal_fraction_frequency_loss_reduction_NOT_causal=(
                            float(gain.mean()/mean_blind) if mean_blind>1e-12 else None),
                        island_equal_mean_blind_frequency_D2=mean_blind,
                        island_equal_mean_name_heldout_genus_frequency_D2=float(
                            private.name_heldout_genus_category_frequency_D2.mean()),
                        original_direct_carrier_weighted_fraction_with_NAME_disjoint_congener=float(
                            np.average(private.share_carriers_with_independent_NAME_genus_donor,
                                       weights=private.n_direct_carriers)),
                        original_direct_carrier_weighted_fraction_with_variable_donor_genus=float(
                            np.average(private.share_carriers_with_at_least_two_distinct_genus_profiles,
                                       weights=private.n_direct_carriers)),
                        n_original_direct_target_carrier_rows=int(private.n_direct_carriers.sum()),
                        mean_direct_carriers_per_qualified_island=float(
                            private.n_direct_carriers.mean()),
                        median_original_unique_train_species_after_name_exclusion=float(
                            private.n_unique_external_same_ledger_species_donors_after_name_holdout.median()),
                        all_test_original_species_names_excluded_from_both_source_modes=True)
                    for grouping,key in (("physical_block","physical_block_equal"),
                                         ("source_cluster","first_GIFT_source_group_equal")):
                        stats=block_bootstrap_paired_prediction(
                            private,group=grouping,n_boot=n_boot,
                            seed=int.from_bytes(hashlib.sha256(
                                (scope+endpoint+region+mode+key).encode()).digest()[:8],"big"))
                        row[key]=stats
                    output.append(row)
    if (coverage["all_recorded__flower_primary_color__floral_form"]!=2324 or
            coverage["original_native__flower_primary_color__floral_form"]!=327):
        raise ValueError("Original true Direct paired trait opportunity changed")
    # Cross-region equal-source comparisons have to use their ORIGINAL region
    # denominators rather than pretending the regions are independent taxa.
    summaries=[]
    for scope in SCOPES:
        for endpoint in ep:
            for mode in MODE:
                x=[r for r in output if r["scope"]==scope and
                   r["trait_endpoint"]==endpoint and
                   r["reference_donor_mode"]==mode]
                if not x:continue
                n=sum(r["n_original_source_floras_with_name_heldout_scores"] for r in x)
                if n==0:continue
                summaries.append(dict(scope=scope,trait_endpoint=endpoint,
                   reference_donor_mode=mode,
                   n_total_nonoverlapping_three_region_scored_island_floras=n,
                   n_original_source_regional_groups=len(x),
                   island_equal_mean_genus_minus_blind_loss_gain=float(sum(
                       r["island_equal_mean_name_heldout_genus_vs_prior_D2_gain"]*
                       r["n_original_source_floras_with_name_heldout_scores"]
                       for r in x)/n),
                   globally_pooled_island_weighted_not_phylogenetically_independent=True))
    return dict(
       status="ORIGINAL_GLOBAL_13_ENDPOINT_FULL_FLORA_NAME_HELDOUT_GENUS_PREDICTION",
       n_original_source_botanical_island_floras=4505,
       n_original_direct_trait_endpoints=13,
       source_full_4505_island_endpoint_eligibility=coverage,
       source_original_region_name_disjoint_donor_support=opportunity,
       original_name_heldout_genus_vs_prior_island_scores=output,
       pooled_nonindependent_region_weighted_summary=summaries,
       claims=dict(
           full_all_status_focal_original_flora_species_names_excluded_from_training=True,
           global_ledger_donors_still_not_independent_mainland_plant_measurement=True,
           independently_recorded_genus_identity_known_oracle_not_colonization=True,
           twenty_six_endpoint_status_analyses_with_distinct_island_support=True,
           original_paired_traits_require_true_same_species_states=True,
           no_928_triple_complete_case_population_shrink=True,
           negative_genus_prediction_gains_retained_not_filtered=True,
           conditional_frequency_loss_gain_not_causal_adaptation_percentage=True,
           original_region_source_donors_have_original_source_trait_provider_overlap=True,
           original_spatial_and_GIFT_source_cluster_draws_not_independent_evolution=True,
           no_private_original_accepted_plant_island_genus_GIFT_source_coordinates_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):h.update(chunk)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for key in ("original-root","geology","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    paths={
        "traits":a.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":a.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        "assignment":a.original_root/"fixed/source/result/island_source_assignments.csv.gz"}
    pins={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        "assignment":"8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44"}
    for key,path in paths.items():
        if sha(path)!=pins[key]:
            raise ValueError("Original protected GIFT botanical source changed")
    state=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    floras=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geology=pd.read_csv(a.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assignment=pd.read_csv(paths["assignment"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if geology.island_id.nunique()!=8264:
        raise ValueError("Original GSHHG source frame changed")
    z=analyse(state,floras,geology,assignment)
    body=json.dumps(z,indent=2,sort_keys=True,allow_nan=False)
    for bad in ('"island_id"','"accepted_species"','"spatial_block"',
                '"entity_ID"','"source_signature"',
                '"island_latitude"','"island_longitude"'):
        if bad in body:raise ValueError("Original botanical source identity leaked")
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(status=z["status"],n_botanical_floras=4505,
        n_endpoint_source_strata=len(z["original_name_heldout_genus_vs_prior_island_scores"]),
        actual_independent_mainland_study=False)))


if __name__=="__main__":
    main()
