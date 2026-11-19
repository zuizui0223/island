"""Original GIFT 4,505 island floral profile prediction WITHOUT focal GENUS oracle.

Input: each original accepted GIFT island flora and 13 species-Direct
endpoints, independent GSHHG macroregions and original first GIFT
candidate mainland-source entity. Donor training bans entire accepted
species rosters from focal 10-degree block and/or first-source group.

Original other-region genus-specific P(trait | genus) trained on exact
different accepted species plus a SHARED external blind category prior.
Target-region genus frequency weights learned ONLY from distant GIFT
island floras' SOURCE-ELIGIBLE distinct accepted species with a Direct
state. No focal genus identity, abundance or trait labels are used
to choose model weights. Predict categorical distribution for ALL
observed Direct target species on the focal island (not only source
genus-supported plants) and compare three different predictions:
  external source BLIND category prior,
  source-region blind phenotype frequency baseline,
  TARGET-region TRAINING genus-mixture × EXTERNAL-region genus traits.
An additional oracle-known-focal-genus conditional comparator quantifies
the original model-information advantage, never hidden in positive claims.

This is NO prospective colonization SDM: broad target macroregion is known,
genus weights represent regional reported Direct species names in other
GIFT ISLAND floras, and original species-wide Direct literature phenotypes
are shared across source islands. No animal/pollen/adaptive claim.
"""
from __future__ import annotations
from collections import Counter,defaultdict
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pandas as pd

from fd_island_functional_beta import strict_traits
from fd_global_genus_assembly_geometry import source_endpoints
from fd_interregional_genus_name_exclusive_transfer import (
    donor_source_region_ownership,hashed_donors_by_genus,
    reference_profile_prior,genus,sha,REGIONS,SCOPES,
    MIN_DIRECT,MIN_ISLANDS,MIN_BLOCKS,ALPHA,MAX_DONOR_PER_GENUS
)
from fd_global_geo_source_training_holdout import (
    original_full_flora_group_blacklists,source_training_ban,
    _source_group,MODES
)

NBOOT=499
COLUMNS=(
    "genus_mix_vs_external_blind_Brier",
    "genus_mix_vs_external_blind_logloss",
    "genus_mix_vs_local_blind_Brier",
    "genus_mix_vs_local_blind_logloss",
    "oracle_known_genus_vs_genus_mix_Brier",
    "oracle_known_genus_vs_genus_mix_logloss",
    "oracle_known_genus_vs_external_blind_Brier",
    "oracle_known_genus_vs_external_blind_logloss",
)
QUALIFIED="ORIGINAL_GIFT_SOURCE_POOL_GENUS_ORACLE_FREE_GROUP_HELDOUT_PROFILE"
MIN_LOCAL_SOURCE_SPECIES=20
MIN_EXTERNAL_SOURCE_SPECIES=20
MIN_SHARED_TAXONOMIC_GENERA=2
MIN_LOCAL_NAMES_WITH_EXTERNAL_CONGENER_FRACTION=.10


def _loss_to_plant_states(prob:np.ndarray,category_codes:list[int])->tuple[float,float]:
    p=np.asarray(prob,dtype=float)
    if (not len(category_codes) or not np.isfinite(p).all()
        or (p<=0).any() or not np.isclose(p.sum(),1)):
        raise ValueError("Original source probability distribution/target state incomplete")
    codes=np.asarray(category_codes,dtype=int)
    if (codes.min()<0 or codes.max()>=len(p)):
        raise ValueError("Source actual Direct target category outside training ontology")
    brier=float(1.+np.dot(p,p)-2.*p[codes].mean())
    logloss=float(-np.log(p[codes]).mean())
    return brier,logloss


def refit_no_oracle_regional_genus_source(
    local_region_only:set[str],external_region_only:set[str],
    forbidden:set[str],actual_states:dict,category_indices:dict,
)->dict|None:
    """Region-level expected GENUS mixture never sees focal island genus roster.

    Weights are number of DISTANT target-region GIFT source species with
    directly recorded source trait category by genus, not occupancy or
    unobserved GIFT island genus prediction. This is SOURCE-POOL taxonomic
    expected category mixture; do not label it empirical island colonization.
    """
    if local_region_only&external_region_only:
        raise ValueError("Original region-exclusive local/other species labels overlap")
    if not local_region_only|external_region_only<=set(actual_states):
        raise ValueError("Original donor species lacks High/Medium Direct trait state")
    if not forbidden:
        raise ValueError("Original focal/co-block/source names must be excluded")
    local=local_region_only-forbidden
    external=external_region_only-forbidden
    if not local or not external:
        return None
    if (local|external)&forbidden or local&external:
        raise AssertionError("Donor contains excluded focal/source plant")
    p_external=reference_profile_prior(external,actual_states,category_indices)
    p_local=reference_profile_prior(local,actual_states,category_indices)
    local_counts=Counter(genus(sp) for sp in local)
    outside=hashed_donors_by_genus(external,seed=20261119)
    local_training=hashed_donors_by_genus(local,seed=20261119)
    common_genera=set(local_counts)&set(outside)
    K=len(p_external)
    mix=np.zeros(K,dtype=float)
    genus_profiles={}
    for g,count in local_counts.items():
        eligible=outside.get(g,())
        n=min(len(eligible),len(local_training[g]),MAX_DONOR_PER_GENUS)
        if not n:
            profile=p_external
        else:
            obs=np.zeros(K,dtype=float)
            for name in eligible[:n]:
                obs[category_indices[actual_states[name]]]+=1.
            profile=(obs+ALPHA*p_external)/(n+ALPHA)
        genus_profiles[g]=profile
        mix+=count/len(local)*profile
    if (not np.isclose(mix.sum(),1.) or (mix<=0).any()
        or not np.isclose(p_local.sum(),1.)):
        raise ValueError("Source-inferred floral category probability not normalized")
    comparable_species=sum(local_counts[g] for g in common_genera)
    return dict(
        source_genus_weighted_external_flower_profile=mix,
        external_genus_blind_prior=p_external,
        regional_local_genus_blind_prior=p_local,
        external_conditional_genus_profiles=genus_profiles,
        n_distinct_local_direct_source_species=len(local),
        n_distinct_external_direct_source_species=len(external),
        n_distinct_local_genus_profiles=len(local_counts),
        n_distinct_genus_profiles_with_external_direct_congeners=len(common_genera),
        n_local_direct_source_species_in_genera_with_outside_trait_donors=comparable_species,
        local_source_species_fraction_with_external_same_genus_donor=comparable_species/len(local),
        n_disjoint_source_excluded_local_direct_species=len(local_region_only)-len(local),
        n_disjoint_source_excluded_external_direct_species=len(external_region_only)-len(external),
        original_focal_island_genus_identity_NOT_USED_IN_MODELLING=True,
        source_local_genus_weights_not_focal_island_occupancy=True,
    )


def score_original_island_without_genus_input(
    original_target_plant_roster:set[str],whole_target_flora:set[str],
    real_direct_states:dict,indices:dict,source:dict,
    *,min_direct:int=MIN_DIRECT,
)->dict|None:
    """All focal Direct plants, including genera without external donors."""
    target=set(original_target_plant_roster)&set(real_direct_states)
    if not target<=whole_target_flora:
        raise ValueError("Focal original Direct target not on accepted flora")
    if len(target)<min_direct:
        return None
    codes=[indices[real_direct_states[sp]] for sp in sorted(target)]
    profile=source["source_genus_weighted_external_flower_profile"]
    blind=source["external_genus_blind_prior"]
    local=source["regional_local_genus_blind_prior"]
    profile_l=_loss_to_plant_states(profile,codes)
    blind_l=_loss_to_plant_states(blind,codes)
    local_l=_loss_to_plant_states(local,codes)
    # Oracle comparator explicitly KNOWS every focal genus and is NOT used
    # in the no-oracle source mixture or pooled category prior.
    genus_map=source["external_conditional_genus_profiles"]
    oracle_brier=oracle_logloss=0.
    for species in sorted(target):
        p=genus_map.get(genus(species),blind)
        b,l=_loss_to_plant_states(p,[indices[real_direct_states[species]]])
        oracle_brier+=b
        oracle_logloss+=l
    oracle=(oracle_brier/len(target),oracle_logloss/len(target))
    return dict(
        n_original_direct_plant_carriers=len(target),
        n_original_focal_direct_carriers_in_SOURCE_known_genera=sum(
            genus(sp) in genus_map for sp in target),
        genus_mix_vs_external_blind_Brier=blind_l[0]-profile_l[0],
        genus_mix_vs_external_blind_logloss=blind_l[1]-profile_l[1],
        genus_mix_vs_local_blind_Brier=local_l[0]-profile_l[0],
        genus_mix_vs_local_blind_logloss=local_l[1]-profile_l[1],
        oracle_known_genus_vs_genus_mix_Brier=profile_l[0]-oracle[0],
        oracle_known_genus_vs_genus_mix_logloss=profile_l[1]-oracle[1],
        oracle_known_genus_vs_external_blind_Brier=blind_l[0]-oracle[0],
        oracle_known_genus_vs_external_blind_logloss=blind_l[1]-oracle[1],
        mean_source_regional_genus_coverage=source[
            "local_source_species_fraction_with_external_same_genus_donor"],
    )


def group_equal_original_predictive_gains(frame:pd.DataFrame,group:str,
                                          *,seed:int,nboot:int=NBOOT)->dict:
    n=len(frame)
    if n and (not set(COLUMNS)<=set(frame) or group not in frame):
        raise ValueError("Source group loss columns missing")
    x=frame.dropna(subset=[group,*COLUMNS])
    n=x.shape[0]
    groups=x[group].nunique()
    if n<MIN_ISLANDS or groups<MIN_BLOCKS:
        return dict(status="NOT_IDENTIFIABLE_ORIGINAL_SOURCE_BLOCK_SUPPORT",
                    n_source_scored_floras=n,
                    n_source_independent_groups=groups)
    G=x.groupby(group,sort=True)[list(COLUMNS)].mean().to_numpy()
    if not np.isfinite(G).all():
        raise ValueError("Source original group trait scores not finite")
    rng=np.random.default_rng(seed)
    drawn=G[rng.integers(0,groups,size=(nboot,groups))].mean(axis=1)
    result=dict(
        status=QUALIFIED,
        n_source_scored_floras=int(n),
        n_source_independent_groups=int(groups),
        original_group_bootstrap_draws=nboot,
        original_group_bootstrap_exploratory_not_multiplicity_adjusted=True)
    for j,name in enumerate(COLUMNS):
        z=np.array([G[np.arange(groups)!=i,j].mean()
                    for i in range(groups)],dtype=float)
        result[name+"_group_equal_mean"]=float(G[:,j].mean())
        result[name+"_group_resample_95"]=[
            float(np.quantile(drawn[:,j],.025)),
            float(np.quantile(drawn[:,j],.975))]
        result[name+"_leave_one_source_group_min"]=float(z.min())
        result[name+"_leave_one_source_group_max"]=float(z.max())
    return result


def analyse(original_direct:dict,flora:pd.DataFrame,
            geography:pd.DataFrame,source_assignments:pd.DataFrame,
            *,modes=MODES,nboot:int=NBOOT,
            endpoint_filter:set[str]|None=None)->dict:
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original GIFT flora source schema incomplete")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicated original accepted GIFT flora plant")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Unexpected original GIFT accepted plant status")
    if not {"island_id","analysis_regime","spatial_block"}<=set(geography):
        raise ValueError("Original GSHHG physical-block macroregion missing")
    if geography.island_id.duplicated().any():
        raise ValueError("Original GSHHG island physical block duplicated")
    if not {"island_id","source_mode","source_rank","entity_ID"}<=set(source_assignments):
        raise ValueError("Original first candidate mainland source assignment absent")
    full=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    if len(full)!=4505:
        raise ValueError("Original 4505 GIFT botanical island flora parent changed")
    geo=geography.set_index("island_id").reindex(sorted(full))
    source=source_assignments[source_assignments.source_mode.eq("geo_k5")].copy()
    source.source_rank=pd.to_numeric(source.source_rank,errors="raise")
    if source.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original GIFT source rank duplicate")
    first=source.loc[source.source_rank.eq(1)].set_index("island_id").entity_ID.to_dict()
    block={i:_source_group(geo.at[i,"spatial_block"]) for i in full}
    first_by={i:_source_group(first.get(i)) for i in full}
    regions={i:_source_group(geo.at[i,"analysis_regime"]) or "geography_unresolved"
             for i in full}
    actual_counts=Counter(regions.values())
    if (sum(actual_counts.values())!=4505 or
        sum(actual_counts.get(reg,0) for reg in REGIONS)!=4504
        or actual_counts.get("geography_unresolved",0)!=1):
        raise ValueError("Original 4505 GIFT geographic study-region denominator drift")
    blacklists=original_full_flora_group_blacklists(full,block,first_by)
    endpoints=source_endpoints(original_direct)
    if len(endpoints)!=13:
        raise ValueError("Original true ten single and three paired Direct endpoints changed")
    filter_set=set(endpoints) if endpoint_filter is None else set(endpoint_filter)
    if not filter_set or not filter_set<=set(endpoints):
        raise ValueError("Unrecognized original source endpoint filter")
    if len(set(modes))!=len(modes) or any(m not in MODES for m in modes):
        raise ValueError("Source region group exclusion mode invalid")

    all_rows=[]
    for scope,mapping in (("all_recorded",full),("original_native",native)):
        sources_by_region={reg:set().union(*(mapping.get(i,set()) for i in full
            if regions[i]==reg))
            for reg in (*REGIONS,"geography_unresolved")}
        for endpoint,states in endpoints.items():
            if endpoint not in filter_set:continue
            indices={v:j for j,v in enumerate(sorted(set(states.values()),key=repr))}
            known=set(states)
            for region in REGIONS:
                local_src,external_src=donor_source_region_ownership(sources_by_region,region)
                local_src &= known
                external_src &= known
                targets=[i for i in sorted(full) if regions[i]==region]
                elig=[i for i in targets if len(mapping.get(i,set())&known)>=MIN_DIRECT]
                for mode in modes:
                    missing=missing_model=low_coverage=0
                    scored=[]
                    trained_cache={}
                    for island in elig:
                        forbidden=source_training_ban(
                            island,full[island],block[island],first_by[island],
                            blacklists,mode)
                        if forbidden is None:
                            missing+=1
                            continue
                        cache_id=(mode,
                            block[island] if mode in ("block_out","block_and_first_source_out") else None,
                            first_by[island] if mode in ("first_source_out","block_and_first_source_out") else None,
                        )
                        if cache_id not in trained_cache:
                            trained_cache[cache_id]=refit_no_oracle_regional_genus_source(
                                local_src,external_src,forbidden,states,indices)
                        model=trained_cache[cache_id]
                        if model is None:
                            missing_model+=1
                            continue
                        if (model["n_distinct_local_direct_source_species"]<MIN_LOCAL_SOURCE_SPECIES
                            or model["n_distinct_external_direct_source_species"]<MIN_EXTERNAL_SOURCE_SPECIES
                            or model["n_distinct_genus_profiles_with_external_direct_congeners"]<MIN_SHARED_TAXONOMIC_GENERA
                            or model["local_source_species_fraction_with_external_same_genus_donor"]<
                                MIN_LOCAL_NAMES_WITH_EXTERNAL_CONGENER_FRACTION):
                            low_coverage+=1
                            continue
                        q=score_original_island_without_genus_input(
                            mapping.get(island,set())&known,full[island],
                            states,indices,model)
                        if q is None:
                            raise AssertionError("Source original prefiltered Direct focal n>=5 vanished")
                        q["spatial_block"]=block[island]
                        q["source_cluster"]=first_by[island]
                        scored.append(q)
                    if len(scored)+missing+missing_model+low_coverage!=len(elig):
                        raise AssertionError("Original 4505 source no-oracle eligibility waterfall missing floras")
                    row=dict(
                        scope=scope,trait_endpoint=endpoint,region=region,
                        training_source_group_deletion_mode=mode,
                        original_full_parent_GIFT_floras=4505,
                        n_original_botanical_islands_in_region=len(targets),
                        n_original_islands_with_direct_trait_n5=len(elig),
                        n_missing_source_block_or_source_group=missing,
                        n_no_original_local_or_external_name_disjoint_donors=missing_model,
                        n_original_eligible_floras_insufficient_source_pool_genus_overlap=low_coverage,
                        n_original_island_floras_scored_without_actual_focal_genus_roster=len(scored),
                        original_source_pool_genus_support_rule=dict(
                            min_distinct_target_region_direct_species=MIN_LOCAL_SOURCE_SPECIES,
                            min_distinct_external_region_direct_species=MIN_EXTERNAL_SOURCE_SPECIES,
                            min_shared_botanical_genera=MIN_SHARED_TAXONOMIC_GENERA,
                            min_local_direct_species_fraction_with_external_congener=
                                MIN_LOCAL_NAMES_WITH_EXTERNAL_CONGENER_FRACTION),
                        target_original_genus_identity_NOT_USED_TO_SELECT_PREDICTIVE_MIXTURE=True,
                    )
                    if scored:
                        df=pd.DataFrame(scored)
                        row.update(
                            n_scored_actual_source_Direct_plant_records=int(
                                df.n_original_direct_plant_carriers.sum()),
                            n_mean_original_plant_carriers_per_island=float(
                                df.n_original_direct_plant_carriers.mean()),
                            mean_fraction_of_focal_direct_species_with_known_regional_external_genus=float(
                                df.n_original_focal_direct_carriers_in_SOURCE_known_genera.sum()/
                                df.n_original_direct_plant_carriers.sum()),
                            mean_original_SOURCE_pool_genus_coverage=float(
                                df.mean_source_regional_genus_coverage.mean()),
                            island_equal_genus_mix_vs_external_blind_Brier=float(
                                df.genus_mix_vs_external_blind_Brier.mean()),
                            island_equal_genus_mix_vs_external_blind_logloss=float(
                                df.genus_mix_vs_external_blind_logloss.mean()),
                            island_equal_genus_mix_vs_local_blind_Brier=float(
                                df.genus_mix_vs_local_blind_Brier.mean()),
                            island_equal_genus_mix_vs_local_blind_logloss=float(
                                df.genus_mix_vs_local_blind_logloss.mean()),
                            island_equal_oracle_known_genus_vs_genus_mix_Brier=float(
                                df.oracle_known_genus_vs_genus_mix_Brier.mean()),
                            island_equal_oracle_known_genus_vs_genus_mix_logloss=float(
                                df.oracle_known_genus_vs_genus_mix_logloss.mean()),
                        )
                        row["physical_10degree_block"]=group_equal_original_predictive_gains(
                            df,"spatial_block",seed=int.from_bytes(
                                hashlib.sha256((scope+endpoint+region+mode+"block").encode()).digest()[:8],"big"),
                            nboot=nboot)
                        row["first_candidate_mainland_source_group"]=group_equal_original_predictive_gains(
                            df,"source_cluster",seed=int.from_bytes(
                                hashlib.sha256((scope+endpoint+region+mode+"source").encode()).digest()[:8],"big"),
                            nboot=nboot)
                        row["status"]=(
                            "ORIGINAL_NO_ORACLE_SOURCE_GENUS_MIXTURE_BOTH_CLUSTER_SUPPORT"
                            if row["physical_10degree_block"]["status"]==QUALIFIED
                            and row["first_candidate_mainland_source_group"]["status"]==QUALIFIED
                            else "NOT_IDENTIFIABLE_SOURCE_CLUSTER_SUPPORT"
                        )
                    else:
                        row["status"]="NOT_IDENTIFIABLE_ORIGINAL_DIRECT_SOURCE_POOL_GENUS_OVERLAP"
                    all_rows.append(row)
    if len(filter_set)==13 and len(all_rows)!=312:
        raise AssertionError("Full 4,505 flora 13 endpoint 2scope 4region 3mode census incomplete")
    return dict(
        status="ORIGINAL_GIFT_4505_NO_FOCAL_GENUS_ORACLE_GEOGRAPHIC_HOLDOUT_FLORAL_PROFILES",
        n_original_botanical_island_floras=4505,
        n_original_direct_endpoint_families=13,
        n_actual_original_direct_endpoints_tested=len(filter_set),
        original_direct_endpoints_tested=sorted(filter_set),
        region_original_floristic_denominators=[
            dict(region=r,n_original_accepted_botanical_floras=actual_counts.get(r,0))
            for r in (*REGIONS,"geography_unresolved")],
        source_training_modes=list(modes),
        n_original_geographical_10degree_blocks=len(blacklists["block"]),
        n_first_GIFT_candidate_mainland_source_groups=len(blacklists["first_source"]),
        region_trait_status_source_group_results=all_rows,
        conceptual_limits=dict(
            original_all_4505_botanical_parent_retained=True,
            original_observed_focal_genus_composition_NOT_used_in_no_oracle_predictor=True,
            original_observed_focal_genus_identity_used_ONLY_in_explicit_oracle_comparator=True,
            original_training_flora_local_genus_weights_learned_from_source_DIRECT_SPECIES_COUNTS_not_actual_island_colonization=True,
            original_candidate_first_mainland_source_group_used_ONLY_FOR_TRAINING_BLACKLIST_and_RESAMPLING=True,
            all_local_and_external_training_accepted_species_names_in_focal_sameblock_firstsource_groups_removed=True,
            actual_local_flora_genus_mixture_not_estimated_from_target_island_trait_distribution=True,
            source_plant_direct_states_from_other_GIFT_ISLANDS_not_independent_MAINLAND_ancestors=True,
            genus_blind_local_source_trait_frequency_is_positive_control_not_phylogenetic_source=True,
            original_13_direct_trait_categories_and_native_statuses_not_independent_replicates=True,
            no_effective_pollinator_guild_transition_pollen_fitness_or_island_local_adaptation_claim=True,
            missing_source_genus_overlap_does_not_imply_no_biological_effect=True,
            only_source_anonymous_results_no_original_island_species_genus_coordinates_exported=True,
        ),
    )


def main():
    ap=argparse.ArgumentParser()
    for arg in ("original-root","geography","output"):
        ap.add_argument("--"+arg,type=Path,required=True)
    ap.add_argument("--nboot",type=int,default=NBOOT)
    args=ap.parse_args()
    paths=dict(
        traits=args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        flora=args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
        assign=args.original_root/"fixed/source/result/island_source_assignments.csv.gz",
    )
    pins=dict(
        traits="9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        flora="6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
        assign="8a9d79c698afa20e611663adb884c8caf648e1a3836853fad94e6a0cd6facb44",
    )
    for k,p in paths.items():
        if sha(p)!=pins[k]:
            raise ValueError("Original global GIFT full botanic trait/source identity drift: "+k)
    traits=strict_traits(pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"]))
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    source=pd.read_csv(paths["assign"],dtype=str,usecols=[
        "island_id","source_mode","source_rank","entity_ID"])
    geography=pd.read_csv(args.geography,dtype=str,usecols=[
        "island_id","analysis_regime","spatial_block"])
    if geography.island_id.nunique()!=8264:
        raise ValueError("Original independent source GSHHG region geography altered")
    result=analyse(traits,flora,geography,source,nboot=args.nboot)
    if result["n_original_botanical_island_floras"]!=4505:
        raise AssertionError("Parent 4,505 GIFT botanical island flora denominator drift")
    data=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"entity_ID"',
                    '"spatial_block"','"source_cluster"','"source_genus"',
                    '"island_latitude"','"island_longitude"'):
        if private in data:
            raise ValueError("Original protected taxon/island/source identity leaked")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(data+"\n",encoding="utf-8")
    print("ANONYMOUS_ORACLE_FREE_4505_SUMMARY="+json.dumps(dict(
        status=result["status"],original_full_island_floras=4505,
        n_direct_endpoints=result["n_actual_original_direct_endpoints_tested"],
        n_original_region_trait_source_modes=len(
            result["region_trait_status_source_group_results"]),
        n_supported_both_independent_geographical_group_panels=sum(
            r.get("status")=="ORIGINAL_NO_ORACLE_SOURCE_GENUS_MIXTURE_BOTH_CLUSTER_SUPPORT"
            for r in result["region_trait_status_source_group_results"]),
        no_focal_genus_composition_oracle=True),sort_keys=True))
    for q in result["region_trait_status_source_group_results"]:
        if (q["scope"]!="all_recorded" or
            q["region"] not in (
            "northern_midlatitude","tropical","southern_extratropical") or
            q["trait_endpoint"] not in (
            "floral_form","floral_symmetry","flower_primary_color",
            "flower_primary_color__floral_form",
            "flower_primary_color__floral_symmetry","floral_form__floral_symmetry")):
            continue
        print("ANONYMOUS_ORACLE_FREE_FLOWER="+json.dumps(
            {k:q.get(k) for k in (
            "scope","trait_endpoint","region",
            "training_source_group_deletion_mode","status",
            "n_original_islands_with_direct_trait_n5",
            "n_original_island_floras_scored_without_actual_focal_genus_roster",
            "island_equal_genus_mix_vs_external_blind_Brier",
            "island_equal_genus_mix_vs_external_blind_logloss",
            "island_equal_genus_mix_vs_local_blind_Brier",
            "island_equal_genus_mix_vs_local_blind_logloss",
            "island_equal_oracle_known_genus_vs_genus_mix_Brier",
            "physical_10degree_block","first_candidate_mainland_source_group")},
            sort_keys=True))
if __name__=="__main__":
    main()
