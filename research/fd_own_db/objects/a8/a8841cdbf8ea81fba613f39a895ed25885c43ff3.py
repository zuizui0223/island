"""GIFT original 4505 botanical floras: test no-oracle genus source transport.

Reads ONLY anonymized original protected original genus-mix source output.
The focal island's true genus roster is excluded from no-oracle model
construction. Compare external genus mix against two sensible baselines:
external source-blind categorical trait prior, and strong target-region
botanical trait frequencies learned from OTHER heldout-group islands.
Both proper Brier and log losses must agree in BOTH original spatial and
first-candidate source group bootstraps.

A separate oracle-known-focal-genus comparator quantifies the degree of
taxonomic foreknowledge (not a fixed positivity expectation). Negative
source-transfer gains and source-insufficient panels must remain visible.
"""
from __future__ import annotations
from itertools import product
import math

REGIONS=("northern_midlatitude","tropical",
         "southern_extratropical","northern_high_latitude")
SCOPES=("all_recorded","original_native")
MODES=("block_out","first_source_out","block_and_first_source_out")
GROUPS=("physical_10degree_block","first_candidate_mainland_source_group")
QUALIFIED="ORIGINAL_NO_ORACLE_SOURCE_GENUS_MIXTURE_BOTH_CLUSTER_SUPPORT"
GOOD_GROUP="ORIGINAL_GIFT_SOURCE_POOL_GENUS_ORACLE_FREE_GROUP_HELDOUT_PROFILE"
METRICS=(
    ("external_blind","genus_mix_vs_external_blind"),
    ("local_blind","genus_mix_vs_local_blind"),
    ("oracle_gap","oracle_known_genus_vs_genus_mix"),
)
TRAIT_NAMES=(
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","inflorescence_display","tube_depth_class",
    "self_incompatibility","mating_system","autonomous_selfing_capacity",
    "cleistogamy","flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
)


def _number(v):
    return type(v) in (float,int) and math.isfinite(v)


def _score(source:dict|None,key:str):
    if (not isinstance(source,dict) or source.get("status")!=GOOD_GROUP
        or source.get("n_source_independent_groups",0)<10):
        return None
    middle=source.get(key+"_group_equal_mean")
    bounds=source.get(key+"_group_resample_95")
    if (not _number(middle) or not isinstance(bounds,list) or len(bounds)!=2
        or not all(_number(z) for z in bounds) or bounds[0]>bounds[1]):
        return None
    return dict(value=float(middle),lower=float(bounds[0]),
        upper=float(bounds[1]),source_groups=int(source["n_source_independent_groups"]),
        source_scored_floras=int(source["n_source_scored_floras"]))


def _direction(row:dict,prefix:str):
    collected={}
    for method in ("Brier","logloss"):
        for group in GROUPS:
            name=prefix+"_"+method
            collected[(method,group)]=_score(row.get(group),name)
    if row.get("status")!=QUALIFIED or not all(collected.values()):
        status="not_identifiable_source"
    elif all(v["lower"]>0 for v in collected.values()):
        status="positive"
    elif all(v["upper"]<0 for v in collected.values()):
        status="negative"
    else:
        status="ambiguous_or_loss_disagreement"
    return dict(
        status=status,
        physical_block_Brier=collected[("Brier",GROUPS[0])],
        physical_block_logloss=collected[("logloss",GROUPS[0])],
        first_source_Brier=collected[("Brier",GROUPS[1])],
        first_source_logloss=collected[("logloss",GROUPS[1])],
    )


def audit_no_oracle_source_flora(x:dict)->dict:
    if (x.get("status")!=
        "ORIGINAL_GIFT_4505_NO_FOCAL_GENUS_ORACLE_GEOGRAPHIC_HOLDOUT_FLORAL_PROFILES"
        or x.get("n_original_botanical_island_floras")!=4505
        or x.get("n_original_direct_endpoint_families")!=13
        or x.get("n_actual_original_direct_endpoints_tested")!=13
        or x.get("source_training_modes")!=list(MODES)):
        raise ValueError("Original 4,505 whole-flora/no-oracle 13-endpoint provenance failed")
    limits=x.get("conceptual_limits",{})
    for key in (
        "original_all_4505_botanical_parent_retained",
        "original_observed_focal_genus_composition_NOT_used_in_no_oracle_predictor",
        "original_observed_focal_genus_identity_used_ONLY_in_explicit_oracle_comparator",
        "original_training_flora_local_genus_weights_learned_from_source_DIRECT_SPECIES_COUNTS_not_actual_island_colonization",
        "all_local_and_external_training_accepted_species_names_in_focal_sameblock_firstsource_groups_removed",
        "source_plant_direct_states_from_other_GIFT_ISLANDS_not_independent_MAINLAND_ancestors",
        "original_13_direct_trait_categories_and_native_statuses_not_independent_replicates",
    ):
        if limits.get(key) is not True:
            raise ValueError("Original GIFT source/true observational claim ceiling drift: "+key)
    rows=x.get("region_trait_status_source_group_results")
    if not isinstance(rows,list) or len(rows)!=312:
        raise ValueError("Original 13x4x2x3 support tensor incomplete")
    keys=set()
    obj=[]
    for r in rows:
        key=(r.get("scope"),r.get("trait_endpoint"),r.get("region"),
             r.get("training_source_group_deletion_mode"))
        if key in keys:
            raise ValueError("Repeated original source region/status/trait/mode")
        keys.add(key)
        if r.get("original_full_parent_GIFT_floras")!=4505:
            raise ValueError("Original 4505 parent not retained in every source result")
        original_n=r.get("n_original_islands_with_direct_trait_n5")
        scored=r.get("n_original_island_floras_scored_without_actual_focal_genus_roster")
        x1=r.get("n_missing_source_block_or_source_group")
        x2=r.get("n_no_original_local_or_external_name_disjoint_donors")
        x3=r.get("n_original_eligible_floras_insufficient_source_pool_genus_overlap")
        if not all(type(z)==int and z>=0 for z in (original_n,scored,x1,x2,x3)):
            raise ValueError("Original source support waterfall is malformed")
        if scored+x1+x2+x3!=original_n:
            raise ValueError("Original no-oracle study support attrition does not reconcile")
        summary={}
        for name,prefix in METRICS:
            summary[name]=_direction(r,prefix)
        obj.append(dict(
            scope=key[0],trait_endpoint=key[1],region=key[2],
            source_mode=key[3],
            n_original_eligible_direct_trait_island_floras=original_n,
            n_matched_full_target_genus_oracle_free_island_floras=scored,
            original_source_group_status=r.get("status"),
            gains=summary,
        ))
    expected=set(product(SCOPES,TRAIT_NAMES,REGIONS,MODES))
    if keys!=expected:
        raise ValueError("Original 13 endpoint source keys incomplete or unknown")
    aggreg=[]
    for scope,mode,reg in product(SCOPES,MODES,REGIONS):
        matched=[r for r in obj if r["scope"]==scope and
                 r["source_mode"]==mode and r["region"]==reg]
        row=dict(scope=scope,source_mode=mode,region=reg,
            n_original_direct_endpoints=13,
            n_qualifying_both_source_group_panels=sum(
                r["original_source_group_status"]==QUALIFIED for r in matched))
        for name,_ in METRICS:
            row[name+"_positive_on_both_losses_both_groupings"]=sum(
                r["gains"][name]["status"]=="positive" for r in matched)
            row[name+"_negative_on_both_losses_both_groupings"]=sum(
                r["gains"][name]["status"]=="negative" for r in matched)
        aggreg.append(row)
    by={(r["scope"],r["trait_endpoint"],r["region"],r["source_mode"]):r
        for r in obj}
    triple={}
    allmode={}
    for scope in SCOPES:
        for name,_ in METRICS:
            supported={}
            for mode in MODES:
                positive=set()
                for trait in TRAIT_NAMES:
                    c=[by[(scope,trait,reg,mode)]["gains"][name]["status"]
                       for reg in REGIONS[:3]]
                    if all(s=="positive" for s in c):
                        positive.add(trait)
                supported[mode]=sorted(positive)
            triple[scope+"__"+name]=supported
            intersection=set(TRAIT_NAMES).intersection(
                *(set(supported[m]) for m in MODES))
            allmode[scope+"__"+name]=sorted(intersection)
    no_source_conclusion=dict(
        original_4505_global_botanical_population_preserved=True,
        original_target_genus_roster_never_available_to_no_oracle_predictor=True,
        regional_source_genus_frequency_weights_NOT_target_genus_occupancy=True,
        cross_region_source_species_are_other_GIFT_ISLAND_records_not_measured_MAINLAND=True,
        native_labelled_scope_nested_not_separate_ecological_replication=True,
        direct_traits_and_regional_groupings_correlated_not_13_independent_evolutionary_tests=True,
        no_original_mainland_ancestral_selection_or_effective_pollinator_guild_shifts_identified=True,
        source_no_oracle_plant_name_exclusion_depends_on_observed_focal_flora_for_holdout_validation_only=True,
        bootstrap_resampled_original_source_not_multiplicity_adjusted_pvalue=True,
    )
    return dict(
        status="ORIGINAL_GIFT_4505_ORACLE_FREE_REGIONAL_GENUS_MIXTURE_VS_TWO_BLIND_BASELINES",
        n_original_source_botanical_floras=4505,
        n_original_Direct_endpoints=13,
        source_original_13_endpoint_region_scope_training_mode_results=obj,
        source_region_scope_mode_oracle_free_prediction_counts=aggreg,
        three_region_supported_direct_endpoints_by_holdout=triple,
        all_three_geography_holdout_intersection=allmode,
        scientific_ceilings=no_source_conclusion,
    )
