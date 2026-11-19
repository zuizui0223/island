"""GLOBAL 4,505 original GIFT floral taxa: exclude local neighbor data from TRAINING.

The existing region-exclusive genus prediction removes the focal island's
species names, but can still train on DIFFERENT accepted species observed
on other islands in the SAME 10-degree physical block or first potential
GIFT mainland source entity group. This is a stricter negative-control:
  - block_out: exclude ANY accepted species name recorded in ANY original
    GIFT island flora of the focal island's physical 10-degree block.
  - source_out: exclude ANY name in ANY original GIFT island flora whose
    original geo_k5 first candidate mainland source entity matches focal.
  - block_and_source_out: exclude the UNION of both sets.

Exclusion happens *before* learning both local/external genus profiles AND
the common externally region-exclusive blind category prior. Species
occurring in both distant and excluded islands are wholly excluded.
The original source-region-exclusive name split is computed from all
original floras, so region eligibility does not use target phenotype.
For each focal island, original baseline and every stricter model are
compared on the EXACT SAME accepted Direct target species/plant genera
still supported by the stricter model, NOT different trait/plant cohorts.
Original 4,505 botanically accepted island floras, 13 Direct endpoints,
all provenance and independent geographic regimes remain visible as
ANONYMOUS support censuses; no 928/927 complete trait cohort.

Actual focal genus composition is supplied to both models: CONDITIONAL
taxonomic prediction, not prediction of arriving taxa. Other source plant
names are species-wide Direct plant traits in GIFT island floras, not
independently phenotyped ancestral mainland populations. Source group
clustering and resampling are exploratory, not independent adaptation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from collections import defaultdict,Counter
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_global_genus_assembly_geometry import source_endpoints
from fd_interregional_genus_name_exclusive_transfer import (
    donor_source_region_ownership,hashed_donors_by_genus,
    reference_profile_prior,island_heldout_comparison,
    group_equal_pair_stats,sha,REGIONS,SCOPES,MIN_DIRECT,
    MIN_MATCHED_PER_ISLAND,MIN_ISLANDS,MIN_BLOCKS,
)

MODES=("block_out","first_source_out","block_and_first_source_out")
GROUPS=(("spatial_block","physical_10degree_block"),
        ("source_cluster","first_candidate_source_group"))
MAIN_COLS=(
    "strict_local_vs_external_Brier",
    "baseline_on_identical_carriers_local_vs_external_Brier",
    "strict_minus_baseline_local_vs_external_Brier",
    "strict_local_vs_external_logloss",
    "baseline_on_identical_carriers_local_vs_external_logloss",
    "strict_minus_baseline_local_vs_external_logloss",
    "strict_other_region_genus_vs_blind_Brier",
)
NB_BOOT=499


def _source_group(v):
    if not isinstance(v,str) or not v.strip() or v.strip().lower() in (
            "nan","none","null","unknown","na"):
        return None
    return v.strip()


def original_full_flora_group_blacklists(
    full:dict[str,set[str]],
    block_by:dict[str,str|None],
    first_by:dict[str,str|None],
)->dict:
    """Original geographic TRAINING ban set includes whole accepted floras.

    Any listed species in any co-block or same-first-source GIFT island is
    forbidden even if some different island also contains that species.
    """
    if len(full)!=4505 and len(full)<2:
        raise ValueError("Original accepted GIFT island source must be nonempty")
    if set(block_by)!=set(full) or set(first_by)!=set(full):
        raise ValueError("Original geographic/source assignments not aligned to full flora")
    block=defaultdict(set)
    source=defaultdict(set)
    nblk=Counter()
    nsrc=Counter()
    for island,original_names in full.items():
        if not original_names:
            raise ValueError("Original island flora may not be empty")
        b=_source_group(block_by[island])
        g=_source_group(first_by[island])
        if b is not None:
            block[b].update(original_names)
            nblk[b]+=1
        if g is not None:
            source[g].update(original_names)
            nsrc[g]+=1
    return dict(block=dict(block),first_source=dict(source),
        n_islands_per_block=dict(nblk),n_islands_per_source=dict(nsrc))


def source_training_ban(
    focal:str,whole_focal_flora:set[str],
    block_id:str|None,first_id:str|None,
    blacklists:dict,mode:str,
)->set[str]|None:
    """Never fallback to unblocked training on missing cluster identities."""
    if mode not in MODES:
        raise ValueError("Unrecognized original source-training exclusion")
    b=_source_group(block_id);g=_source_group(first_id)
    if mode in ("block_out","block_and_first_source_out") and (
            b is None or b not in blacklists["block"]):
        return None
    if mode in ("first_source_out","block_and_first_source_out") and (
            g is None or g not in blacklists["first_source"]):
        return None
    forbidden=set(whole_focal_flora)
    if mode in ("block_out","block_and_first_source_out"):
        forbidden.update(blacklists["block"][b])
    if mode in ("first_source_out","block_and_first_source_out"):
        forbidden.update(blacklists["first_source"][g])
    if not whole_focal_flora<=forbidden:
        raise AssertionError("All target names need to be excluded")
    return forbidden


def make_donor_models(
    local_source:set[str],external_source:set[str],
    forbidden:set[str],
    state:dict,indices:dict,
)->dict|None:
    """Re-fit all source-conditioned quantities AFTER training exclusions."""
    if local_source&external_source:
        raise ValueError("Regional exclusive donor names were not disjoint")
    if not (local_source|external_source)<=set(state):
        raise ValueError("Original candidate plant names have no Direct trait")
    local=local_source-forbidden
    external=external_source-forbidden
    if not local or not external:
        return None
    prior=reference_profile_prior(external,state,indices)
    lb=hashed_donors_by_genus(local,seed=20261115)
    eb=hashed_donors_by_genus(external,seed=20261115)
    genera=set(lb)&set(eb)
    return dict(
        prior=prior,local_groups=lb,external_groups=eb,
        matched_genera=genera,
        local_names=len(local),external_names=len(external),
        n_matched_genera=len(genera),
        n_local_excluded=len(local_source)-len(local),
        n_external_excluded=len(external_source)-len(external),
    )


def score_paired_identical_original_species(
    original_whole_focal_flora:set[str],original_direct_target:set[str],
    state:dict,indices:dict,
    baseline:dict|None,strict:dict|None,
    *,min_carriers:int=MIN_MATCHED_PER_ISLAND,
)->dict|None:
    """Old and strict models scored on EXACT SAME accepted target species."""
    if baseline is None or strict is None:
        return None
    actual=set(original_direct_target)
    if not actual or not actual<=original_whole_focal_flora:
        raise ValueError("Target Direct species not on ORIGINAL focal flora")
    shared=baseline["matched_genera"]&strict["matched_genera"]
    target={sp for sp in actual if sp.split()[0] in shared}
    if len(target)<min_carriers:
        return None
    scores={}
    for name,model in (("before",baseline),("strict",strict)):
        z=island_heldout_comparison(
            original_whole_focal_flora,target,state,indices,
            model["prior"],model["local_groups"],model["external_groups"])
        if z["n_matched_genus_direct_carriers"]!=len(target):
            raise AssertionError("Unmatched focal genus leaked into comparison")
        scores[name]=z
    a=scores["before"];b=scores["strict"]
    def dd(key):
        if a.get(key) is None or b.get(key) is None:
            raise ValueError("Matched source proper scoring unavailable")
        return float(b[key]-a[key])
    return dict(
        common_target_direct_carrier_count=len(target),
        n_different_target_genera=len({s.split()[0] for s in target}),
        original_n_target_direct_carriers=len(actual),
        strict_local_vs_external_Brier=float(
            b["local_minus_external_proper_Brier_gain_per_comparable_plant"]),
        baseline_on_identical_carriers_local_vs_external_Brier=float(
            a["local_minus_external_proper_Brier_gain_per_comparable_plant"]),
        strict_minus_baseline_local_vs_external_Brier=dd(
            "local_minus_external_proper_Brier_gain_per_comparable_plant"),
        strict_local_vs_external_logloss=float(
            b["local_minus_external_proper_logloss_gain_per_comparable_plant"]),
        baseline_on_identical_carriers_local_vs_external_logloss=float(
            a["local_minus_external_proper_logloss_gain_per_comparable_plant"]),
        strict_minus_baseline_local_vs_external_logloss=dd(
            "local_minus_external_proper_logloss_gain_per_comparable_plant"),
        strict_other_region_genus_vs_blind_Brier=float(
            b["external_genus_vs_blind_proper_Brier_gain_per_comparable_plant"]),
        baseline_on_identical_carriers_other_region_genus_vs_blind_Brier=float(
            a["external_genus_vs_blind_proper_Brier_gain_per_comparable_plant"]),
        strict_source_matched_genus_names=(
            len(strict["matched_genera"])),
        strict_local_source_names=int(strict["local_names"]),
        strict_external_source_names=int(strict["external_names"]),
        n_original_local_donor_species_removed=int(strict["n_local_excluded"]),
        n_original_external_donor_species_removed=int(strict["n_external_excluded"]),
    )


def original_group_equal_stress(frame:pd.DataFrame,cluster:str, *,
                               seed:int,nboot:int=NB_BOOT)->dict:
    """Comparable donor-replacement effect at physical/source group units."""
    x=frame.dropna(subset=[cluster,*MAIN_COLS]).copy()
    n=len(x);nclusters=x[cluster].nunique()
    if n<MIN_ISLANDS or nclusters<MIN_BLOCKS:
        return dict(status="NOT_IDENTIFIABLE_ORIGINAL_GROUP_SUPPORT",
                    n_scored_island_floras=int(n),
                    n_distinct_original_source_clusters=int(nclusters),
                    group_resampling_not_performed=True)
    grouped=x.groupby(cluster,sort=True)[list(MAIN_COLS)].mean().to_numpy()
    rng=np.random.default_rng(seed)
    sampled=grouped[rng.integers(nclusters,size=(nboot,nclusters))].mean(axis=1)
    row=dict(
        status="ORIGINAL_STRICT_TRAINING_EXCLUSION_PAIRED_FOCAL_NAME_HOLDOUT",
        n_scored_island_floras=int(n),
        n_distinct_original_source_clusters=int(nclusters),
        n_outcome_blind_original_group_resamples=nboot,
        n_distinct_heldout_floras_source_reused_within_cluster_not_independent=True)
    for i,c in enumerate(MAIN_COLS):
        row[c+"_group_equal_mean"]=float(grouped[:,i].mean())
        row[c+"_group_resample_95"]=[
            float(np.quantile(sampled[:,i],.025)),
            float(np.quantile(sampled[:,i],.975))]
        q=np.array([grouped[np.arange(nclusters)!=j,i].mean()
                    for j in range(nclusters)])
        row[c+"_leave_one_group_min"]=float(q.min())
        row[c+"_leave_one_group_max"]=float(q.max())
    return row


def analyse(original_direct:dict,flora:pd.DataFrame,geo:pd.DataFrame,
            assign:pd.DataFrame, *,
            modes=MODES,nboot:int=NB_BOOT,
            endpoint_filter:set[str]|None=None)->dict:
    if nboot<1:raise ValueError("Positive source sampling resamples needed")
    if not {"island_id","accepted_species","origin_status"}<=set(flora) or (
            flora.duplicated(["island_id","accepted_species"]).any()):
        raise ValueError("Original accepted floras duplicated/incomplete")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original accepted botanical status drift")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Original independent physical regional geography missing")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(assign):
        raise ValueError("Original GIFT candidate mainland first-source identity missing")
    if any(mode not in MODES for mode in modes) or len(set(modes))!=len(modes):
        raise ValueError("Unknown or duplicated frozen training exclusion modes")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:
        raise ValueError("Original entire GIFT 4505 botanical flora population drifted")
    geo=geo.set_index("island_id").reindex(sorted(full))
    near=assign.loc[assign.source_mode.eq("geo_k5")].copy()
    near.source_rank=pd.to_numeric(near.source_rank,errors="raise")
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original source first ranks duplicate")
    original_first=near.loc[near.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    block={i:_source_group(geo.at[i,"spatial_block"]) for i in full}
    first={i:_source_group(original_first.get(i)) for i in full}
    region={i:(_source_group(geo.at[i,"analysis_regime"]) or "geography_unresolved")
            for i in full}
    regional_counts=Counter(region.values())
    if sum(regional_counts.values())!=4505:
        raise AssertionError("The global GIFT flora total changed")
    if set(regional_counts)-set(REGIONS)-{"geography_unresolved"}:
        raise ValueError("Original supported region unexpected")
    blacklist=original_full_flora_group_blacklists(full,block,first)
    endpoint_source=source_endpoints(original_direct)
    if len(endpoint_source)!=13:
        raise ValueError("Original 10 single + 3 REAL same species pair endpoints changed")
    filter_set=set(endpoint_source) if endpoint_filter is None else set(endpoint_filter)
    if not filter_set<=set(endpoint_source):
        raise ValueError("Original source requested trait endpoint nonexistent")
    if not filter_set:
        raise ValueError("At least one Direct endpoint required")
    result=[]
    denominator=[]
    for scope,roster in (("all_recorded",full),("original_native",native)):
        byregion={r:set().union(*(roster.get(i,set()) for i in full
                if region[i]==r))
             for r in (*REGIONS,"geography_unresolved")}
        for endpoint,state in endpoint_source.items():
            if endpoint not in filter_set:continue
            indices={v:j for j,v in enumerate(sorted(set(state.values()),key=repr))}
            name_set=set(state)
            for reg in REGIONS:
                local,external=donor_source_region_ownership(byregion,reg)
                local &= name_set
                external &= name_set
                focal_ids=[i for i in sorted(full) if region[i]==reg]
                eligible=[i for i in focal_ids if
                          len(roster.get(i,set())&name_set)>=MIN_DIRECT]
                base_model=make_donor_models(local,external,set(),state,indices)
                for mode in modes:
                    rows=[]
                    missing_group=0
                    no_donors=0
                    not_comparable=0
                    cache={}
                    for island in eligible:
                        forbidden=source_training_ban(
                            island,full[island],block[island],first[island],
                            blacklist,mode)
                        if forbidden is None:
                            missing_group+=1
                            continue
                        if base_model is None:
                            no_donors+=1;continue
                        cachekey=(
                            mode,
                            block[island] if mode in ("block_out","block_and_first_source_out") else None,
                            first[island] if mode in ("first_source_out","block_and_first_source_out") else None,
                        )
                        if cachekey not in cache:
                            cache[cachekey]=make_donor_models(
                                local,external,forbidden,state,indices)
                        strict=cache[cachekey]
                        if strict is None:
                            no_donors+=1;continue
                        observed=score_paired_identical_original_species(
                            full[island],roster.get(island,set())&name_set,
                            state,indices,base_model,strict)
                        if observed is None:
                            not_comparable+=1;continue
                        observed["spatial_block"]=block[island]
                        observed["source_cluster"]=first[island]
                        rows.append(observed)
                    entry=dict(
                        original_scope=scope,
                        original_trait_endpoint=endpoint,
                        original_source_region=reg,
                        source_training_exclusion_mode=mode,
                        n_original_botanical_floras_in_region=len(focal_ids),
                        n_original_true_direct_n5_eligible_focal_islands=len(eligible),
                        n_missing_required_original_physical_or_first_source_group=int(missing_group),
                        n_original_islands_lacking_both_strict_training_sources=int(no_donors),
                        n_original_islands_lacking_5_common_direct_comparable_species=int(not_comparable),
                        n_original_original_GIFT_islands_with_strict_common_focal_scores=len(rows),
                        source_original_common_carrier_n5_vs_original_denominator_is_conditional=True,
                        group_by_training_location_excluded_before_source_prior_and_genus_refit=True,
                        original_target_actual_genus_composition_KNOWN=True,
                    )
                    if rows:
                        df=pd.DataFrame(rows)
                        entry.update(
                            original_mean_direct_species_compared_per_focal=float(
                                df.common_target_direct_carrier_count.mean()),
                            original_n_distinct_strict_local_donor_species_average=float(
                                df.strict_local_source_names.mean()),
                            original_n_distinct_strict_external_donor_species_average=float(
                                df.strict_external_source_names.mean()),
                            n_original_compared_real_Direct_species_records=int(
                                df.common_target_direct_carrier_count.sum()),
                            island_equal_strict_local_minus_external_Brier=float(
                                df.strict_local_vs_external_Brier.mean()),
                            island_equal_same_focal_baseline_local_minus_external_Brier=float(
                                df.baseline_on_identical_carriers_local_vs_external_Brier.mean()),
                            island_equal_strict_minus_baseline_Brier=float(
                                df.strict_minus_baseline_local_vs_external_Brier.mean()),
                            island_equal_strict_local_minus_external_logloss=float(
                                df.strict_local_vs_external_logloss.mean()),
                            island_equal_same_focal_baseline_local_minus_external_logloss=float(
                                df.baseline_on_identical_carriers_local_vs_external_logloss.mean()),
                            original_source_local_removed_species_mean=float(
                                df.n_original_local_donor_species_removed.mean()),
                            original_source_external_removed_species_mean=float(
                                df.n_original_external_donor_species_removed.mean()),
                        )
                        for group,key in GROUPS:
                            entry[key]=original_group_equal_stress(
                                df,group,nboot=nboot,seed=int.from_bytes(hashlib.sha256(
                                    (scope+endpoint+reg+mode+key).encode()
                                ).digest()[:8],"big"))
                    else:
                        entry["status"]="NOT_IDENTIFIABLE_NO_SAME_COHORT_MATCHED_DIRECT_FOCALS"
                    if entry["n_original_original_GIFT_islands_with_strict_common_focal_scores"]+(
                            missing_group+no_donors+not_comparable)!=len(eligible):
                        raise AssertionError("Training-exclusion island eligibility waterfall not exhaustive")
                    denominator.append(entry)
                    if rows:
                        if (entry["physical_10degree_block"]["status"]==
                             "ORIGINAL_STRICT_TRAINING_EXCLUSION_PAIRED_FOCAL_NAME_HOLDOUT" and
                            entry["first_candidate_source_group"]["status"]==
                             "ORIGINAL_STRICT_TRAINING_EXCLUSION_PAIRED_FOCAL_NAME_HOLDOUT"):
                            entry["status"]="ORIGINAL_STRICT_TRAINING_EXCLUSION_BOTH_GROUPS_QUALIFIED"
                        else:
                            entry["status"]="NOT_IDENTIFIABLE_INDEPENDENT_ORIGINAL_CLUSTER_SUPPORT"
                    result.append(entry)
    return dict(
        status="ORIGINAL_GIFT_4505_GLOBAL_GEO_SOURCE_TRAINING_EXCLUDED_REGION_TRAIT_TRANSFER",
        n_original_botanical_GIFT_island_floras=4505,
        n_original_floral_direct_endpoint_families=13,
        n_actual_trait_endpoint_families_tested=len(filter_set),
        trait_endpoints_tested=sorted(filter_set),
        n_true_original_regional_botanical_flora_islands=[
            dict(region=r,n_original_floras=int(regional_counts.get(r,0)))
            for r in (*REGIONS,"geography_unresolved")],
        n_original_floras_with_valid_physical_10degree_block=int(
            sum(v is not None for v in block.values())),
        n_original_floras_with_geo5_first_GIFT_candidate_mainland_source=int(
            sum(v is not None for v in first.values())),
        n_distinct_original_10degree_blocks=len(blacklist["block"]),
        n_distinct_original_first_candidate_source_entity_groups=len(
            blacklist["first_source"]),
        original_source_training_exclusion_modes=list(modes),
        original_source_training_eligibility_and_scoring=denominator,
        original_source_training_scored_results=result,
        scientific_limitations=dict(
            focal_accepted_full_flora_names_excluded_from_all_original_source_donor_profiles=True,
            SAME_physical_block_OTHER_island_species_names_removed_from_local_and_external_training=True,
            SAME_original_first_mainland_source_OTHER_island_species_names_removed_from_local_and_external_training=True,
            all_original_source_animal_or_pollen_not_studied=True,
            original_training_species_borrowed_from_GIFT_island_floras_NOT_independently_phenotyped_MAINLAND=True,
            original_plant_genus_composition_given_as_ORACLE_not_predicted=True,
            rescore_baseline_on_IDENTICAL_strict_source_supported_accepted_species=True,
            no_observed_absence_from_source_training_unavailability=True,
            original_4505_botanical_population_retained_despite_endpoint_qualification=True,
            one_analytic_region_endpoint_count_not_independent_evolutionary_replication=True,
            no_adaptation_phylogenetic_convergence_or_guild_shifting_causality_inferred=True,
            source_anonymous_aggregates_only_no_original_island_species_genus_source_region_ID_exported=True,
        ),
    )


def main():
    p=argparse.ArgumentParser()
    for arg in ("original-root","geography","output"):
        p.add_argument("--"+arg,type=Path,required=True)
    p.add_argument("--nboot",type=int,default=NB_BOOT)
    args=p.parse_args()
    paths=dict(
        traits=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        flora=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        assign=args.original_root/"fixed/source/result/island_source_assignments.csv.gz")
    expected=dict(
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44")
    for name,path in paths.items():
        if sha(path)!=expected[name]:
            raise ValueError("Original GIFT source authenticity changed: "+name)
    traits=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(args.geography,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    assign=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    if geo.island_id.nunique()!=8264:
        raise ValueError("Independent GSHHG geography full source coverage drift")
    result=analyse(traits,flora,geo,assign,nboot=args.nboot)
    if result["n_original_botanical_GIFT_island_floras"]!=4505:
        raise AssertionError("Original GIFT 4505 parent dropped")
    body=json.dumps(result,sort_keys=True,indent=2,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"entity_ID"',
                    '"source_region_id"','"spatial_block"','"source_cluster"',
                    '"source_genus"','"island_latitude"','"island_longitude"'):
        if private in body:
            raise ValueError("Original source accepted flora identity leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    rows=result["original_source_training_scored_results"]
    print("ANONYMOUS_GIFT_TRAINING_HOLDOUT_SUPPORT="+json.dumps(dict(
        status=result["status"],full_global_floras=4505,
        original_trait_endpoints=13,
        original_physical_blocks=result["n_distinct_original_10degree_blocks"],
        original_first_source_groups=result[
            "n_distinct_original_first_candidate_source_entity_groups"],
        n_original_floral_region_scope_train_holdout_panels=len(rows),
        n_both_source_qualifying_panels=sum(x.get("status")==
            "ORIGINAL_STRICT_TRAINING_EXCLUSION_BOTH_GROUPS_QUALIFIED"
            for x in rows),original_source_treatment_NOT_original_mainland_phenotypes=True
    ),sort_keys=True))
    for row in rows:
        if row["original_scope"]!="all_recorded" or row[
            "original_trait_endpoint"] not in (
                "floral_form","flower_primary_color__floral_symmetry",
                "flower_primary_color__floral_form","tube_depth_class",
                "floral_symmetry"):
            continue
        if row["original_source_region"] not in (
            "northern_midlatitude","tropical","southern_extratropical"):
            continue
        out={k:row.get(k) for k in (
            "original_scope","original_trait_endpoint","original_source_region",
            "source_training_exclusion_mode",
            "n_original_true_direct_n5_eligible_focal_islands",
            "n_original_original_GIFT_islands_with_strict_common_focal_scores",
            "n_missing_required_original_physical_or_first_source_group",
            "n_original_islands_lacking_both_strict_training_sources",
            "n_original_islands_lacking_5_common_direct_comparable_species",
            "island_equal_strict_local_minus_external_Brier",
            "island_equal_same_focal_baseline_local_minus_external_Brier",
            "status","physical_10degree_block","first_candidate_source_group")}
        print("ANONYMOUS_GIFT_TRAINING_HOLDOUT_FOCUS="+json.dumps(
            out,sort_keys=True))
if __name__=="__main__":
    main()
