"""SOURCE-ORIGINAL global 4505-island across-region within-genus Direct-trait transfer.

A region-local donor plant name must appear ONLY in its target GSHHG region.
An external donor name must appear ONLY in other GSHHG regions, NEVER in
the entire target region. Shared-region plant names are removed from BOTH
donor sources, regardless of individual island. All target island accepted
plant names are additionally excluded from BOTH source donor corpora.
For each known target genus sample the SAME number (max 20) of two
distinct exclusively-region/direct donor species using a deterministic
outcome-blind name hash. Same Jeffreys external source prior, same genus
and same species count. Scores on real original Direct category states,
NOT observed local trait evolution or independent continental phenotypes.

The full 4505 source botanical population and 13 separate Direct endpoints
remain; no old 928 triple-trait complete-case gate, 17D KDE, pollinator
classification, ancestral source, or unverified source geography.
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
from fd_global_genus_assembly_geometry import source_endpoints

REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude")
SCOPES=("all_recorded","original_native")
MIN_DIRECT=5
MIN_MATCHED_PER_ISLAND=5
MIN_ISLANDS,MIN_BLOCKS=30,10
MAX_DONOR_PER_GENUS=20
N_CLUSTER_BOOT=999
ALPHA=5.0


def genus(name:str)->str:
    bits=str(name).split()
    return bits[0] if len(bits)>=2 else "original_single_species:"+str(name)


def donor_source_region_ownership(source_by_region:dict,region:str):
    """No shared accepted species name in either donor source pool."""
    if region not in source_by_region:
        raise ValueError("Target GSHHG region missing")
    other=set().union(*(x for r,x in source_by_region.items() if r!=region))
    local=set(source_by_region[region])
    exclusively_local=local-other
    exclusively_external=set().union(*(
        x for r,x in source_by_region.items()
        if r not in (region,"geography_unresolved")))-local
    if exclusively_external&local or exclusively_local&other or (
            exclusively_external&exclusively_local):
        raise AssertionError("Named original plant occurs in both donor geographies")
    return exclusively_local,exclusively_external


def hashed_donors_by_genus(donor_names:set,*,seed:int):
    """Outcome-BLIND deterministic species-name re-ordering, not trait selection."""
    bins=defaultdict(list)
    for name in donor_names:
        bins[genus(name)].append(name)
    for g,names in bins.items():
        names.sort(key=lambda sp:hashlib.sha256(
            f"fd:20261115:{seed}:{sp}".encode()).hexdigest())
    return bins


def reference_profile_prior(reference:set,state:dict,indices:dict):
    """Same external-only independent-NAME category prior for both sources."""
    if not reference or not reference<=set(state):
        raise ValueError("Strict external-region named trait donors unavailable")
    x=np.full(len(indices),.5,float)
    for sp in reference:x[indices[state[sp]]]+=1.
    x/=x.sum()
    return x


def balanced_probabilities(g:str,original_entire_focal_flora:set,
                           local_group:dict,external_group:dict,
                           state:dict,indices:dict,prior:np.ndarray,
                           *,seed:int=20261115):
    """The same matching n for the same documented genus, no focal name leakage."""
    a=[sp for sp in local_group.get(g,()) if sp not in original_entire_focal_flora]
    b=[sp for sp in external_group.get(g,()) if sp not in original_entire_focal_flora]
    m=min(len(a),len(b),MAX_DONOR_PER_GENUS)
    if m==0:return prior,prior,dict(
        match_n=0,local_names=0,external_names=0,
        local_multicategory=False,external_multicategory=False)
    la=a[:m];ex=b[:m]
    if (set(la)|set(ex))&original_entire_focal_flora or set(la)&set(ex):
        raise AssertionError("Original island target trait name leaked")
    K=len(prior)
    lx=np.zeros(K);ey=np.zeros(K)
    for sp in la:lx[indices[state[sp]]]+=1.
    for sp in ex:ey[indices[state[sp]]]+=1.
    lp=(lx+ALPHA*prior)/(m+ALPHA)
    ep=(ey+ALPHA*prior)/(m+ALPHA)
    if not np.isclose(lp.sum(),1) or not np.isclose(ep.sum(),1):
        raise AssertionError("Original genus source probability simplex invalid")
    return lp,ep,dict(
        match_n=int(m),local_names=int(len(a)),external_names=int(len(b)),
        local_multicategory=bool(np.count_nonzero(lx)>=2),
        external_multicategory=bool(np.count_nonzero(ey)>=2))


def island_heldout_comparison(whole_flora:set,target_roster:set,
                              state:dict,indices:dict,prior:np.ndarray,
                              local_group:dict,external_group:dict):
    """Genuine species Direct profile scores with same island and category labels."""
    if not target_roster or not target_roster<=whole_flora or not (
            target_roster<=set(state)):
        raise ValueError("Original target Direct plant record not in source flora")
    if len(prior)!=len(indices):raise ValueError("External source prior mismatch")
    K=len(indices)
    observed=np.zeros(K)
    genus_targets=defaultdict(list)
    for sp in target_roster:
        j=indices[state[sp]]
        observed[j]+=1.
        genus_targets[genus(sp)].append(j)
    observed/=len(target_roster)
    local=np.zeros(K);external=np.zeros(K)
    comparable_n=0
    comparable_ll_local=comparable_ll_external=0.
    comparable_brier_local=comparable_brier_external=0.
    variable_local=variable_external=0
    nmatched_total=0
    focal_n=len(target_roster)
    for g,codes in genus_targets.items():
        l,e,support=balanced_probabilities(
            g,whole_flora,local_group,external_group,state,indices,prior)
        n=len(codes)
        local+=n*l/focal_n
        external+=n*e/focal_n
        if support["match_n"]:
            comparable_n+=n
            nmatched_total+=support["match_n"]
            if support["local_multicategory"]:variable_local+=n
            if support["external_multicategory"]:variable_external+=n
            for j in codes:
                comparable_ll_local-=np.log(l[j])
                comparable_ll_external-=np.log(e[j])
                comparable_brier_local+=1+np.dot(l,l)-2*l[j]
                comparable_brier_external+=1+np.dot(e,e)-2*e[j]
    if not np.allclose(local.sum(),1) or not np.allclose(external.sum(),1):
        raise AssertionError("Original matched outside-region genus compositional sum failed")
    result=dict(
        n_original_real_direct_carriers=focal_n,
        n_matched_genus_direct_carriers=comparable_n,
        matched_carrier_fraction=comparable_n/focal_n,
        matched_genus_donor_names_sampled_per_side=int(nmatched_total),
        fraction_matched_carriers_in_variable_local_genus=(
            variable_local/comparable_n if comparable_n else 0.),
        fraction_matched_carriers_in_variable_external_genus=(
            variable_external/comparable_n if comparable_n else 0.),
        locally_exclusive_minus_externally_exclusive_frequency_D2_gain=float(
            np.square(observed-external).sum() -
            np.square(observed-local).sum()),
        local_exclusive_frequency_D2=float(np.square(observed-local).sum()),
        external_exclusive_frequency_D2=float(np.square(observed-external).sum()),
    )
    if comparable_n:
        result.update(
            local_minus_external_proper_logloss_gain_per_comparable_plant=float(
                (comparable_ll_external-comparable_ll_local)/comparable_n),
            local_minus_external_proper_Brier_gain_per_comparable_plant=float(
                (comparable_brier_external-comparable_brier_local)/comparable_n))
    else:
        result.update(local_minus_external_proper_logloss_gain_per_comparable_plant=None,
                      local_minus_external_proper_Brier_gain_per_comparable_plant=None)
    return result


def group_equal_pair_stats(frame:pd.DataFrame,group:str,*,seed:int,nboot:int=N_CLUSTER_BOOT):
    cols=("locally_exclusive_minus_externally_exclusive_frequency_D2_gain",
          "local_minus_external_proper_logloss_gain_per_comparable_plant",
          "local_minus_external_proper_Brier_gain_per_comparable_plant")
    q=frame.dropna(subset=[group,*cols])
    n=len(q)
    groups=q[group].nunique()
    if n<MIN_ISLANDS or groups<MIN_BLOCKS:
        return dict(status="NOT_IDENTIFIABLE_INDEPENDENT_ORIGINAL_GROUP_SUPPORT",
                    n_scored_islands=n,n_clusters=groups)
    G=q.groupby(group,sort=True)[list(cols)].mean().to_numpy()
    rng=np.random.default_rng(seed)
    draws=G[rng.integers(groups,size=(nboot,groups))].mean(1)
    loo=[G[np.arange(groups)!=k].mean(0) for k in range(groups)]
    a=np.array(loo)
    result=dict(
        status="ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON",
        n_scored_islands=n,n_clusters=groups,
        n_outcome_blind_group_resamples=nboot)
    for i,name in enumerate(cols):
        result[name+"_group_equal_mean"]=float(G[:,i].mean())
        result[name+"_group_resample_95"]=[
            float(np.quantile(draws[:,i],.025)),
            float(np.quantile(draws[:,i],.975))]
        result[name+"_leave_one_group_min"]=float(a[:,i].min())
        result[name+"_leave_one_group_max"]=float(a[:,i].max())
    return result


def analyse(traits:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            assign:pd.DataFrame,*,nboot:int=N_CLUSTER_BOOT):
    if not {"island_id","accepted_species","origin_status"}<=set(flora) or (
        flora.duplicated(["island_id","accepted_species"]).any()):
        raise ValueError("Invalid original accepted botanical source species roster")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Unrecognized original botanical status")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or (
            geo.island_id.duplicated().any()):
        raise ValueError("Independent original GSHHG classification missing")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assign):
        raise ValueError("Original geographic nearest GIFT donor source rank missing")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:raise ValueError("Original 4505 botanical floras changed")
    geography=geo.set_index("island_id").reindex(sorted(full))
    actual_region={
        i:(str(geography.at[i,"analysis_regime"])
           if pd.notna(geography.at[i,"analysis_regime"])
           and str(geography.at[i,"analysis_regime"]).strip()
           else "geography_unresolved")
        for i in full}
    region_counts=Counter(actual_region.values())
    if (sum(region_counts.values())!=4505 or
            sum(region_counts[r] for r in REGIONS)!=4504):
        raise ValueError("All original 4505 regional source islands not accounted for")
    others=sorted(set(region_counts)-set(REGIONS))
    if others!=["geography_unresolved"]:
        raise ValueError("Unexpected original geographic source regime")
    near=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    near.source_rank=pd.to_numeric(near.source_rank,errors="raise")
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original source rank duplicate")
    first=near.loc[near.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    ep=source_endpoints(traits)
    results=[]
    census=[]
    coverage={}
    for scope,mapping in (("all_recorded",full),("original_native",native)):
        by_region={
            r:set().union(*(mapping.get(i,set())
                            for i in full if actual_region[i]==r))
            for r in (*REGIONS,"geography_unresolved")}
        for endpoint,category_state in ep.items():
            category_state_names=set(category_state)
            allowed={v:i for i,v in enumerate(sorted(
                set(category_state.values()),key=repr))}
            coverage[scope+"__"+endpoint]=sum(
                len(mapping.get(i,set())&category_state_names)>=MIN_DIRECT
                for i in full)
            for region in REGIONS:
                target_ids=[i for i in sorted(full) if actual_region[i]==region]
                in_only,out_only=donor_source_region_ownership(by_region,region)
                local_candidates=in_only&category_state_names
                external_candidates=out_only&category_state_names
                original_n5=sum(
                    len(mapping.get(i,set())&category_state_names)>=MIN_DIRECT
                    for i in target_ids)
                meta=dict(
                    scope=scope,trait_endpoint=endpoint,region=region,
                    n_original_region_botanical_floras=len(target_ids),
                    n_original_n5_direct_trait_islands=original_n5,
                    n_regional_exclusive_direct_donor_species=len(local_candidates),
                    n_other_regions_exclusive_direct_donor_species=len(external_candidates),
                    n_regional_exclusive_donor_genera=len({
                        genus(s) for s in local_candidates}),
                    n_other_regions_exclusive_donor_genera=len({
                        genus(s) for s in external_candidates}))
                if not external_candidates or not local_candidates:
                    census.append(dict(**meta,
                        status="NOT_IDENTIFIABLE_EXCLUSIVE_REGION_DIRECT_DONORS"))
                    continue
                prior=reference_profile_prior(external_candidates,category_state,allowed)
                local_groups=hashed_donors_by_genus(
                    local_candidates,seed=20261115)
                outside_groups=hashed_donors_by_genus(
                    external_candidates,seed=20261115)
                paired_genus=set(local_groups)&set(outside_groups)
                meta["n_distinct_original_genera_with_both_name_disjoint_sources"]=len(
                    paired_genus)
                scored=[]
                known_comparable=0
                eligible_with_5=0
                for island in target_ids:
                    target=mapping.get(island,set())&category_state_names
                    if len(target)<MIN_DIRECT:continue
                    n_comparable=sum(genus(s) in paired_genus for s in target)
                    if n_comparable>=5:eligible_with_5+=1
                    known_comparable+=n_comparable
                    z=island_heldout_comparison(
                        full[island],target,category_state,allowed,prior,
                        local_groups,outside_groups)
                    block=geography.at[island,"spatial_block"]
                    z["spatial_block"]=(str(block) if pd.notna(block) else None)
                    z["source_cluster"]=first.get(island)
                    scored.append(z)
                meta.update(
                    n_target_carrier_records_in_nonoverlapping_geographic_region=sum(
                        r["n_original_real_direct_carriers"] for r in scored),
                    n_target_carriers_with_potential_matched_congeners=known_comparable,
                    n_eligible_islands_with_at_least_5_matched_carriers=eligible_with_5,
                    n_n5_islands_considered=len(scored),
                    n_regions_in_reference=3)
                census.append(meta)
                if not scored:continue
                full_df=pd.DataFrame(scored)
                data=full_df.loc[
                    full_df.n_matched_genus_direct_carriers>=MIN_MATCHED_PER_ISLAND].copy()
                if data.empty:
                    results.append(dict(**meta,
                        status="NOT_IDENTIFIABLE_NO_FOCAL_FLORA_WITH_5_MATCHED_CONGENERS"))
                    continue
                q=dict(**meta,
                       n_same_focal_islands_matched_carrier_n5=len(data),
                       n_original_matched_focal_Direct_carrier_records=int(
                           data.n_matched_genus_direct_carriers.sum()),
                       mean_original_matched_carrier_fraction=float(
                           data.matched_carrier_fraction.mean()),
                       mean_outcome_blind_matched_training_donor_names_per_side=float(
                           data.matched_genus_donor_names_sampled_per_side.mean()),
                       mean_matched_direct_carrier_fraction_in_variable_LOCAL_genus=float(
                           np.average(
                               data.fraction_matched_carriers_in_variable_local_genus,
                               weights=data.n_matched_genus_direct_carriers)),
                       mean_matched_direct_carrier_fraction_in_variable_EXTERNAL_genus=float(
                           np.average(
                               data.fraction_matched_carriers_in_variable_external_genus,
                               weights=data.n_matched_genus_direct_carriers)))
                for col in (
                    "locally_exclusive_minus_externally_exclusive_frequency_D2_gain",
                    "local_minus_external_proper_logloss_gain_per_comparable_plant",
                    "local_minus_external_proper_Brier_gain_per_comparable_plant"):
                    q[col+"_island_equal_mean"]=float(data[col].mean())
                    q[col+"_positive_island_fraction"]=float(np.mean(data[col]>0))
                for key,out_name in (
                    ("spatial_block","physical_10degree_block"),
                    ("source_cluster","first_GIFT_mainland_source")):
                    q[out_name]=group_equal_pair_stats(
                        data,key,nboot=nboot,seed=int.from_bytes(
                            hashlib.sha256((scope+endpoint+region+key).encode()
                                           ).digest()[:8],"big"))
                q["status"]=q["physical_10degree_block"]["status"]
                results.append(q)
    if (coverage["all_recorded__flower_primary_color__floral_form"]!=2324 or
            coverage["original_native__flower_primary_color__floral_form"]!=327):
        raise ValueError("Original 2324/327 Direct real-pair trait support changed")
    if sum(region_counts.values())!=4505:
        raise AssertionError("Original botanical source population not complete")
    return dict(
        status="ORIGINAL_4505_REGION_EXCLUSIVE_EQUAL_GENUS_DONOR_TRANSFER",
        n_original_botanical_island_floras=4505,
        n_original_accepted_source_trait_endpoints=13,
        n_frozen_source_regional_main_comparisons=len(REGIONS),
        n_frozen_source_analysis_scopes=len(SCOPES),
        original_complete_botanical_region_universe=[
            dict(region=r,n_original_botanical_floras=int(region_counts[r]),
                 region_geo_evaluable=r in REGIONS)
            for r in (*REGIONS,"geography_unresolved")],
        source_original_n5_endpoint_coverage=coverage,
        original_exclusive_training_taxonomic_source_support=census,
        original_equal_genus_donor_transfer_results=results,
        claims=dict(
            entire_4505_original_botanical_population_preserved=True,
            target_original_all_accepted_plant_names_removed_from_both_donor_sets=True,
            names_occurring_in_target_and_other_regions_removed_from_both_source_sets=True,
            external_donor_names_never_recorded_in_original_target_geographical_region=True,
            local_and_external_same_original_genus_and_exact_equal_donor_species_count=True,
            donor_species_selected_by_original_name_hash_not_by_trait_state=True,
            same_exclusive_external_species_category_prior_for_both_models=True,
            unmatched_target_genera_neutral_same_prior_for_both_models=True,
            strict_native_source_status_is_nested_nonindependent_scope=True,
            different_regions_and_traits_not_independent_evolutionary_replications=True,
            no_independent_mainland_or_phylogenetic_trait_source_or_adaptation_test=True,
            no_real_pollinator_service_guild_replacement_or_fitness_claim=True,
            no_original_private_taxon_island_genus_entity_location_data_exported=True))


def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as fp:
        for c in iter(lambda:fp.read(1<<20),b""):h.update(c)
    return h.hexdigest()


def main():
    p=argparse.ArgumentParser()
    for s in ("original-root","geology","output"):
        p.add_argument("--"+s,type=Path,required=True)
    args=p.parse_args()
    paths=dict(
        traits=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        flora=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        assign=args.original_root/"fixed/source/result/island_source_assignments.csv.gz")
    pins=dict(
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")
    for k,v in paths.items():
        if sha(v)!=pins[k]:raise ValueError("Protected original Direct flora source SHA drift: "+k)
    trait=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    region=pd.read_csv(args.geology,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assignment=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if region.island_id.nunique()!=8264:
        raise ValueError("Original 8264 GSHHG source geography changed")
    results=analyse(trait,flora,region,assignment)
    if results["n_original_botanical_island_floras"]!=4505:
        raise ValueError("Original 4505 botanical macro population changed")
    output=json.dumps(results,sort_keys=True,indent=2,allow_nan=False)
    for x in ('"island_id"','"accepted_species"','"entity_ID"',
              '"source_signature"','"spatial_block"','"island_latitude"',
              '"island_longitude"'):
        if x in output:raise ValueError("Original private source flora identity leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(output+"\n")
    print(json.dumps(dict(
        status=results["status"],n_original_botanical_floras=4505,
        n_anonymous_trait_support=len(results["original_exclusive_training_taxonomic_source_support"]),
        n_supported_original_region_transfer_rows=len(results["original_equal_genus_donor_transfer_results"]),
        no_mainland_independent_anthropogenic_or_evolutionary_claim=True)))


if __name__=="__main__":main()
