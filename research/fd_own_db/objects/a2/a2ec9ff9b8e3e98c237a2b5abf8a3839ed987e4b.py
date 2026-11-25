"""Audit SHA-verified original GIFT 4505 regional genus trait transportability.

Inputs are ONLY anonymous historical original-source CI aggregates, no
protected island/plant/genus names. Full original 4505 accepted island
floras, all 13 Direct endpoints, all 4 observed broad geography regions,
and both recorded-provenance scopes are retained. This audit does NOT
re-fit an original individual-level model, and paired resampling intervals
are descriptive, not independently replicated evolutionary experiments.

Regional LOCAL/E XTERNAL donor names are mutually region-EXCLUSIVE,
different original accepted plant species but SAME botanical genus, with
outcome-blind equal numbers (up to 20 per genus). These donor names come
from GIFT region-specific ISLAND source floras, not verified ancestral
mainland plant phenotypes. The focal island's genus composition is GIVEN:
this is a taxonomic conditional prediction, not a colonization SDM.
"""
from __future__ import annotations

import math

REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude")
SCOPES=("all_recorded","original_native")
TRAITS=(
    "flower_primary_color","floral_form","floral_symmetry",
    "flower_size_class","inflorescence_display","tube_depth_class",
    "self_incompatibility","mating_system",
    "autonomous_selfing_capacity","cleistogamy",
    "flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
)
METRIC={
    "Brier":"local_minus_external_proper_Brier_gain_per_comparable_plant",
    "logloss":"local_minus_external_proper_logloss_gain_per_comparable_plant",
    "external_genus_Brier":"external_genus_vs_blind_proper_Brier_gain_per_comparable_plant",
}
N_FOCAL=4505
VALID="ORIGINAL_CONDITIONAL_REGION_EXCLUSIVE_DONOR_COMPARISON"


def _finite(x)->bool:
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def _summary(r:dict,metric:str,group:str)->dict|None:
    node=r.get(group)
    if not isinstance(node,dict):
        return None
    name=METRIC[metric]
    center=node.get(name+"_group_equal_mean")
    bounds=node.get(name+"_group_resample_95")
    clusters=node.get("n_clusters")
    if (not _finite(center) or not isinstance(bounds,list)
            or len(bounds)!=2 or not all(_finite(q) for q in bounds)
            or not _finite(clusters) or clusters<10 or
            bounds[0]>bounds[1]):
        return None
    return dict(estimate=float(center),
        lower=float(bounds[0]),upper=float(bounds[1]),
        n_clusters=int(clusters))


def directions(row:dict)->dict:
    """Robust sign if two proper losses AND block/source CIs agree.

    An exploratory CI excluding zero is not a preregistered inferential p.
    Original per-region trait endpoints and island floras are dependent.
    """
    output={}
    for metric in ("Brier","logloss","external_genus_Brier"):
        block=_summary(row,metric,"physical_10degree_block")
        source=_summary(row,metric,"first_GIFT_mainland_source")
        sign="not_qualified"
        if block and source and row.get("status")==VALID:
            intervals=(block,source)
            if all(x["lower"]>0 for x in intervals):sign="local_better" if metric!="external_genus_Brier" else "cross_region_genus_better"
            elif all(x["upper"]<0 for x in intervals):sign="external_better" if metric!="external_genus_Brier" else "cross_region_genus_worse"
            else:sign="source_or_block_ambiguous"
        output[metric]=dict(
            direction=sign,
            physical_block_equal=block,
            first_mainland_source_group_equal=source)
    b=output["Brier"]["direction"];l=output["logloss"]["direction"]
    output["two_proper_loss_direction"]=(
        b if b==l and b in ("local_better","external_better")
        else "not_robust_two_metrics")
    return output


def _key(row:dict)->tuple:
    return row["scope"],row["trait_endpoint"],row["region"]


def validate_source(a:dict)->dict:
    if (a.get("status")!="ORIGINAL_4505_REGION_EXCLUSIVE_EQUAL_GENUS_DONOR_TRANSFER"
        or a.get("n_original_botanical_island_floras")!=N_FOCAL
        or a.get("n_original_accepted_source_trait_endpoints")!=len(TRAITS)):
        raise ValueError("Original full GIFT population/Direct endpoint source drift")
    coverage=a.get("original_complete_botanical_region_universe")
    if not isinstance(coverage,list) or len(coverage)!=5:
        raise ValueError("Original four regions plus unresolved source denominator missing")
    region_counts={x["region"]:x["n_original_botanical_floras"] for x in coverage}
    if (set(region_counts)!=set(REGIONS)|{"geography_unresolved"}
        or sum(region_counts.values())!=N_FOCAL
        or region_counts["geography_unresolved"]!=1):
        raise ValueError("Original four-region GIFT geographic population changed")
    expected=set((scope,trait,region)
                 for scope in SCOPES for trait in TRAITS for region in REGIONS)
    supported=a.get("original_exclusive_training_taxonomic_source_support")
    if not isinstance(supported,list) or len(supported)!=len(expected):
        raise ValueError("Original source donor availability 104 panels changed")
    seen_supported={_key(r) for r in supported}
    if seen_supported!=expected:
        raise ValueError("Original source regional trait eligibility duplicated/incomplete")
    rows=a.get("original_equal_genus_donor_transfer_results")
    if not isinstance(rows,list) or len(rows)!=100:
        raise ValueError("Original source scored endpoint count changed")
    seen=set()
    count={}
    for row in rows:
        k=_key(row)
        if k not in expected or k in seen:
            raise ValueError("Original source phenotypic panel duplicated or outside 4-region cohort")
        seen.add(k)
        status=row.get("status")
        count[status]=count.get(status,0)+1
        if status not in (VALID,
            "NOT_IDENTIFIABLE_INDEPENDENT_ORIGINAL_GROUP_SUPPORT",
            "NOT_IDENTIFIABLE_NO_FOCAL_FLORA_WITH_5_MATCHED_CONGENERS"):
            raise ValueError("Unreported source train/test status")
        match=row.get("n_same_focal_islands_matched_carrier_n5")
        if match is not None and _finite(match):
            if match>N_FOCAL or match<0:
                raise ValueError("Impossible source original focal sample")
        if status==VALID:
            if not _finite(match) or match<30:
                raise ValueError("Supported source original region has fewer than 30 islands")
            if not isinstance(row.get("physical_10degree_block"),dict):
                raise ValueError("Qualified physical block support must be reported")
    if count.get(VALID)!=63 or len(expected-seen)!=4:
        raise ValueError("Original source independent-donor panel support changed")
    original_claims=a.get("claims",{})
    for key in (
        "entire_4505_original_botanical_population_preserved",
        "target_original_all_accepted_plant_names_removed_from_both_donor_sets",
        "names_occurring_in_target_and_other_regions_removed_from_both_source_sets",
        "local_and_external_same_original_genus_and_exact_equal_donor_species_count",
        "donor_species_selected_by_original_name_hash_not_by_trait_state",
        "strict_native_source_status_is_nested_nonindependent_scope",
        "no_independent_mainland_or_phylogenetic_trait_source_or_adaptation_test",
        "no_real_pollinator_service_guild_replacement_or_fitness_claim",
    ):
        if original_claims.get(key) is not True:
            raise ValueError("Original source region-exclusivity or noncausal contract missing: "+key)
    return dict(original_region_counts=region_counts,
                source_qualified_count_by_status=count,
                n_total_target_original_scored_cells=len(rows))


def verify_two_original_runs(first:dict,second:dict)->dict:
    """Check historical raw taxonomic results replicated, no pseudoindependence.

    The older source artifact omitted the later external-genus-vs-blind score;
    check the ORIGINAL COMMON 100 cells, not artificially compare new fields.
    """
    validate_source(first)
    validate_source(second)
    left={_key(x):x for x in first["original_equal_genus_donor_transfer_results"]}
    right={_key(x):x for x in second["original_equal_genus_donor_transfer_results"]}
    if left.keys()!=right.keys():
        raise ValueError("Original 4505 source-run eligibility does not replicate")
    n_comparisons=0
    worst=0.
    for key in left:
        a=left[key];b=right[key]
        for name in (
            "status","n_same_focal_islands_matched_carrier_n5",
            "n_n5_islands_considered","n_original_n5_direct_trait_islands"):
            x=a.get(name);y=b.get(name)
            if x!=y and not (isinstance(x,float) and isinstance(y,float)
                and math.isnan(x) and math.isnan(y)):
                raise ValueError("Original source cohort drift across independent CI replay")
        for name in (
            "local_minus_external_proper_Brier_gain_per_comparable_plant_island_equal_mean",
            "local_minus_external_proper_logloss_gain_per_comparable_plant_island_equal_mean",
            "locally_exclusive_minus_externally_exclusive_frequency_D2_gain_island_equal_mean",
        ):
            x=a.get(name);y=b.get(name)
            if _finite(x) and _finite(y):
                worst=max(worst,abs(x-y));n_comparisons+=1
                if abs(x-y)>1e-10:raise ValueError("Source score drift above numeric tolerance")
        for group in ("physical_10degree_block","first_GIFT_mainland_source"):
            ga=a.get(group);gb=b.get(group)
            if not isinstance(ga,dict) or not isinstance(gb,dict):continue
            for name in (
                "local_minus_external_proper_Brier_gain_per_comparable_plant_group_equal_mean",
                "local_minus_external_proper_logloss_gain_per_comparable_plant_group_equal_mean",
            ):
                x=ga.get(name);y=gb.get(name)
                if _finite(x) and _finite(y):
                    worst=max(worst,abs(x-y));n_comparisons+=1
                    if abs(x-y)>1e-10:raise ValueError("Original spatial/source grouping drift")
    if n_comparisons<100:
        raise ValueError("Too few independently replayed source metrics")
    return dict(n_original_matched_cells=100,
                n_replayed_common_metric_checks=n_comparisons,
                max_absolute_source_difference=float(worst),
                deterministic_same_input_two_CI_runs_not_independent_field_replication=True)


def audit_actual_candidate_mainland_source_novelty(a:dict)->dict:
    """STRICTLY separate GIFT mainland candidate lists from region-exclusive
    OTHER-ISLAND donor traits.

    An originally source-listed name is not independent phenotyping of a
    mainland ancestor: both source panels inherit the same species-wide
    original Direct trait ledger. A completely new name would still require
    independent actual mainland phenotype measurement.
    """
    if (a.get("status")!="ORIGINAL_GIFT_4505_MAINLAND_SOURCE_NOVEL_SPECIES_HOLDOUT_VS_INTERREGIONAL_GENUS_ORACLE"
        or a.get("parent_original_global_GIFT_botanical_island_floras")!=4505):
        raise ValueError("Original actual mainland-vs-GIFT source null not original full population")
    if (a.get("n_original_mainland_GIFT_native_source_region_lists")!=748
        or a.get("n_original_global_all_recorded_direct_colour_form_test_islands")!=2324
        or a.get("n_original_direct_colour_form_source_species_before_strict_holdout")!=3606):
        raise ValueError("GIFT original source candidate mainland coverage changed")
    if (a.get("n_original_source_colour_form_species_outside_ALL_4505_GIFT_island_floras")!=0
        or a.get("n_original_candidate_mainland_donors_after_2324_trait_eligible_island_full_name_holdout")!=9):
        raise ValueError("Original source independent mainland names claim changed")
    claims=a.get("original_source_provenance_limit",{})
    for key in (
        "global_4505_all_flora_name_holdout_leaves_ZERO_independent_source_direct_pair_donors",
        "cohort_2324_name_holdout_leaves_NINE_source_direct_pair_donors_but_novelty_not_ancestry",
        "original_single_species_wide_Direct_trait_ledger_shared_by_both_sources",
        "other_region_conditional_genus_predictor_KNOWS_focal_island_genus_composition",
        "independent_candidate_mainland_GIFT_lists_not_ancestral_population_phenotypes",
    ):
        if claims.get(key) is not True:
            raise ValueError("Independent mainland and island source training distinction absent: "+key)
    cohorts=a.get("original_actual_candidate_mainland_GIFT_source_holdout_cohorts")
    if not isinstance(cohorts,list) or not cohorts:
        raise ValueError("Original different-cohort source denominator not available")
    all_cohort=[r for r in cohorts if r.get("source_cohort")=="global_all_n5"]
    if len(all_cohort)!=1 or (
        all_cohort[0].get("n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed")!=9
        or all_cohort[0].get("n_original_source_donor_genera_with_two_plus_joint_category_states")!=0):
        raise ValueError("Original 2324 GIFT mainland candidate donor source status drift")
    return dict(
        parent_global_GIFT_island_botanical_floras=4505,
        n_candidate_mainland_GIFT_native_source_regions=748,
        n_original_direct_colour_form_candidate_mainland_species_names=3606,
        n_original_direct_colour_form_mainland_species_names_outside_ALL_4505_island_floras=0,
        n_candidate_source_species_after_2324_full_flora_names_heldout=9,
        n_candidate_source_genera_with_true_pair_variation_after_2324_holdout=0,
        result_independent_mainland_only_cannot_test_prospective_4505_global_flower_colour_form_prediction=True,
        region_exclusive_other_GIFT_island_species_source_transfer_is_NOT_independent_mainland_source=True,
        no_source_lineage_phylogenetic_or_guild_adaptation_inference=True,
        source_candidate_cohorts_overlap_and_cannot_directly_subtract_effects=True,
    )


def analyze_original_global_regionality(source:dict,older:dict,
                                        mainland_candidate:dict|None=None)->dict:
    provenance=validate_source(source)
    replay=verify_two_original_runs(source,older)
    mainland_independence=(audit_actual_candidate_mainland_source_novelty(
        mainland_candidate) if mainland_candidate is not None else None)
    rows=source["original_equal_genus_donor_transfer_results"]
    aggregated=[]
    qualified={}
    for row in rows:
        scope,trait,region=_key(row)
        z=directions(row)
        evidence=dict(scope=scope,trait_endpoint=trait,region=region,
            original_source_panel_status=row["status"],
            n_focal_islands_matched_original_distinct_genus_donors=(
                int(row["n_same_focal_islands_matched_carrier_n5"])
                if _finite(row.get("n_same_focal_islands_matched_carrier_n5"))
                else 0),
            n_original_source_region_floras=int(row["n_original_region_botanical_floras"]),
            source_matched_carrier_fraction=(
                float(row["mean_original_matched_carrier_fraction"])
                if _finite(row.get("mean_original_matched_carrier_fraction")) else None),
            confidence=directions(row))
        aggregated.append(evidence)
        qualified[(scope,trait,region)]=evidence
    totals=[]
    for scope in SCOPES:
        for region in REGIONS:
            x=[a for a in aggregated if a["scope"]==scope and a["region"]==region]
            b=[a for a in x if a["original_source_panel_status"]==VALID]
            totals.append(dict(scope=scope,region=region,
                n_all_13_original_direct_endpoints=len(x),
                n_region_endpoints_qualified_after_original_source_eligibility=len(b),
                n_robust_local_better_both_proper_losses=sum(
                    a["confidence"]["two_proper_loss_direction"]=="local_better" for a in b),
                n_robust_cross_region_better_both_proper_losses=sum(
                    a["confidence"]["two_proper_loss_direction"]=="external_better" for a in b),
                n_cross_region_genus_better_than_external_genus_blind_brier=sum(
                    a["confidence"]["external_genus_Brier"]["direction"]=="cross_region_genus_better"
                    for a in b),
                note_interdependent_endpoints_not_independent_replications=True))
    region_reversals=[]
    for scope in SCOPES:
        for trait in TRAITS:
            n=qualified.get((scope,trait,"northern_midlatitude"))
            t=qualified.get((scope,trait,"tropical"))
            if not n or not t:continue
            one=n["confidence"]["two_proper_loss_direction"]
            two=t["confidence"]["two_proper_loss_direction"]
            if one=="local_better" and two=="external_better":
                region_reversals.append(dict(
                    scope=scope,trait_endpoint=trait,
                    north_n_focal_islands=n["n_focal_islands_matched_original_distinct_genus_donors"],
                    tropical_n_focal_islands=t["n_focal_islands_matched_original_distinct_genus_donors"],
                    northern_block_Brier=n["confidence"]["Brier"]["physical_block_equal"],
                    northern_first_source_Brier=n["confidence"]["Brier"]["first_mainland_source_group_equal"],
                    tropical_block_Brier=t["confidence"]["Brier"]["physical_block_equal"],
                    tropical_first_source_Brier=t["confidence"]["Brier"]["first_mainland_source_group_equal"],
                    brier_and_logloss_both_source_and_block_sign_reversal=True))
    return dict(
        status="ORIGINAL_GIFT_4505_CONDITIONAL_INTERREGIONAL_GENUS_TRAIT_TRANSPORT_AUDIT",
        original_n_global_botanical_floras=N_FOCAL,
        n_original_highmedium_direct_endpoints=len(TRAITS),
        source_original_four_region_and_missing_denominator=provenance,
        independent_execution_same_original_data_replay=replay,
        independent_mainland_candidate_species_provenance_gate=mainland_independence,
        all_original_qualifying_trait_endpoint_scope_region_results=aggregated,
        source_robust_13_endpoint_panel_counts=totals,
        original_north_to_tropical_double_metric_sign_reversal=region_reversals,
        scientific_ceilings=dict(
            original_regional_local_donors_from_ISLAND_regional_floristic_lists_not_independent_mainland_ancestral_observations=True,
            original_focal_botanical_genus_composition_is_ORACLE_in_prediction=True,
            true_original_focal_accepted_plant_species_names_excluded_globally=True,
            two_CI_replays_are_not_two_independent_ecological_datasets=True,
            paired_group_bootstrap_intervals_exploratory_not_multiplicity_adjusted=True,
            source_regions_and_endpoint_traits_not_independent_evolutionary_replicates=True,
            observed_region_specific_conditional_prediction_not_adaptive_island_trait_evolution=True,
            original_Pollinator_guild_shift_selection_fitness_not_observed=True,
            main_global_4505_botanical_population_not_restricted_to_qualifying_subsets=True,
        ),
    )
