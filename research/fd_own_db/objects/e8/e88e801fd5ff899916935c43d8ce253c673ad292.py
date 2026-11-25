"""Strict donor-origin reality check for GIFT 4,505 genus transfer claims.

Compare separate, signed ORIGINAL anonymous source artifacts:
(A) regional GENUS prediction trained on name-exclusive species of
    OTHER GIFT ISLAND floras (13 Direct endpoints, known focal genus);
(B) actual GIFT MAINLAND candidate source-entity species lists for
    real same-species Direct flower colour×form, with all focal
    accepted-plant names excluded from source regions.

These have different traits/test cohorts/source rosters, so NO direct
difference or pooled p-value is computed. This is an identification
and source-data novelty gate, not ancestral colonization/selection.
"""
from __future__ import annotations
import math

COHORTS=(
    "global_all_n5","global_original_native_n5",
    "northmid_n10_coverage10",
    "northmid_500_2000km_exact456_endpoint_union",
)
SOURCE_STATUS="ORIGINAL_STRICT_GIFT_SOURCE_NAME_HOLDOUT_PREFLIGHT"
TRANSFER_STATUS="ORIGINAL_4505_GIFT_THREE_REGION_TWO_LOSS_TWO_GROUP_GENUS_TRANSFER_REPLICATION"


def _cohort(q:dict)->dict:
    req={
       "cohort","status","n_original_focal_islands",
       "n_unique_source_direct_pair_donor_species_before_holdout",
       "n_unique_source_direct_pair_donor_species_after_exact_name_holdout",
       "n_donor_genera_with_multiple_joint_profiles",
       "n_original_focal_full_flora_names",
       "ratio_test_carriers_with_any_global_same_genus_donor",
       "ratio_test_carriers_with_same_genus_local_five_source_donor",
       "source_local_minus_global_genus_improvement_island_brier",
       "source_local_minus_global_genus_improvement_island_log_loss",
    }
    if not req<=set(q):
        raise ValueError("Missing original candidate-mainland source holdout fields")
    n=q["n_original_focal_islands"]
    before=q["n_unique_source_direct_pair_donor_species_before_holdout"]
    after=q["n_unique_source_direct_pair_donor_species_after_exact_name_holdout"]
    variable=q["n_donor_genera_with_multiple_joint_profiles"]
    if (not all(isinstance(a,int) and a>=0 for a in (n,before,after,variable))
        or after>before):
        raise ValueError("Original native source opportunities invalid")
    for c in ("ratio_test_carriers_with_any_global_same_genus_donor",
              "ratio_test_carriers_with_same_genus_local_five_source_donor"):
        v=q[c]
        if not isinstance(v,(int,float)) or not (0<=v<=1):
            raise ValueError("Original candidate source genus coverage invalid")
    if q["status"] not in (
          "POSTHOC_OUT_OF_TARGET_NAME_SOURCE_PREDICTION",
          "INSUFFICIENT_INDEPENDENT_GEOGRAPHY",
          "NOT_IDENTIFIABLE_NO_DIRECT_PAIR_GIFT_DONORS_AFTER_HOLDOUT"):
        raise ValueError("Original candidate source status unqualified")
    for k in ("source_local_minus_global_genus_improvement_island_brier",
              "source_local_minus_global_genus_improvement_island_log_loss"):
        if not isinstance(q[k],(int,float)) or not math.isfinite(q[k]):
            raise ValueError("Original proper loss finite holdout score missing")
    return dict(
       source_cohort=q["cohort"],
       original_test_botanical_islands=n,
       original_focal_full_flora_accepted_species_names_excluded=int(
          q["n_original_focal_full_flora_names"]),
       n_candidate_mainland_colour_form_direct_species_before_name_holdout=before,
       n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed=after,
       n_original_source_donor_genera_with_two_plus_joint_category_states=variable,
       share_direct_focal_carriers_with_any_donor_genus=float(
          q["ratio_test_carriers_with_any_global_same_genus_donor"]),
       share_direct_focal_carriers_with_local_five_mainland_donor_genus=float(
          q["ratio_test_carriers_with_same_genus_local_five_source_donor"]),
       local_mainland_genus_minus_global_source_genus_Brier_gain_descriptive=float(
          q["source_local_minus_global_genus_improvement_island_brier"]),
       local_mainland_genus_minus_global_source_genus_logloss_gain_descriptive=float(
          q["source_local_minus_global_genus_improvement_island_log_loss"]),
       source_model_original_status=q["status"],
       no_direct_cross_cohort_p_value_or_confidence_interval=True)


def audit_candidate_mainland_reality(mainland:dict,transfer:dict)->dict:
    if (mainland.get("status")!=SOURCE_STATUS or
        mainland.get("n_original_botanical_islands")!=4505 or
        mainland.get("n_original_distinct_source_regions")!=748 or
        mainland.get("n_original_mainland_native_source_rows")!=1666186):
        raise ValueError("Not the SHA-verified full GIFT 4,505 original mainland candidate source")
    claims=mainland.get("claims",{})
    if not all(claims.get(k) is True for k in (
        "full_4505_botanical_universe_retained",
        "direct_species_flower_colour_and_form_joint_not_imputed",
        "same_name_removed_from_ALL_mainland_regions_before_training",
        "genus_known_is_a_posthoc_oracle_not_colonization_prediction",
        "no_double_counting_of_local_donors_overlapping_global_genus",
    )):
        raise ValueError("Original potential mainland-source holdout contract broken")
    if (transfer.get("status")!=TRANSFER_STATUS
        or transfer.get("n_original_botanical_global_island_floras")!=4505
        or transfer.get("n_original_direct_endpoints")!=13
        or transfer.get("n_all_recorded_endpoints_replicating_across_all_three_geographies")!=8
        or not transfer.get("strict_approach",{}).get(
            "source_training_other_regions_island_floras_not_proven_ancestral_mainland_species")):
        raise ValueError("Original regional GIFT-island genus transfer source identity drift")
    raw=mainland.get("cohorts")
    if not isinstance(raw,list) or len(raw)!=4 or {q["cohort"] for q in raw}!=set(COHORTS):
        raise ValueError("Original four geography/native candidate mainland cohorts changed")
    rows=[_cohort(q) for q in raw]
    by={x["source_cohort"]:x for x in rows}
    full=by["global_all_n5"]
    if (full["original_test_botanical_islands"]!=2324 or
        full["n_candidate_mainland_colour_form_direct_species_before_name_holdout"]!=3606 or
        full["n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed"]!=9
        or full["n_original_source_donor_genera_with_two_plus_joint_category_states"]!=0):
        raise ValueError("Original GIFT full-parent mainland source-name holdout ceiling changed")
    if mainland.get("n_source_direct_pair_species_after_ALL_4505_island_name_exclusion")!=0:
        raise ValueError("All 4505 accepted-flora names leave ZERO source-only Direct colour/form mainland donors")
    strict_native=by["global_original_native_n5"]
    if (strict_native["original_test_botanical_islands"]!=327 or
        strict_native["n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed"]!=679):
        raise ValueError("Known-native source data support changed")
    return dict(
        status="ORIGINAL_GIFT_4505_MAINLAND_SOURCE_NOVEL_SPECIES_HOLDOUT_VS_INTERREGIONAL_GENUS_ORACLE",
        parent_original_global_GIFT_botanical_island_floras=4505,
        n_original_mainland_GIFT_native_source_region_lists=748,
        n_original_mainland_native_region_species_occurrence_rows=1666186,
        n_original_direct_colour_form_source_species_before_strict_holdout=3606,
        n_original_source_colour_form_species_outside_ALL_4505_GIFT_island_floras=0,
        n_original_candidate_mainland_donors_after_2324_trait_eligible_island_full_name_holdout=9,
        n_original_global_all_recorded_direct_colour_form_test_islands=2324,
        original_actual_candidate_mainland_GIFT_source_holdout_cohorts=rows,
        conditional_other_region_GIFT_island_genus_predictive_endpoints_across_three_regions=int(
            transfer["n_all_recorded_endpoints_replicating_across_all_three_geographies"]),
        conditional_other_region_GIFT_island_genus_predictive_endpoints_native_only_three_region=int(
            transfer["n_native_only_endpoints_replicating_across_all_three_geographies"]),
        original_source_provenance_limit=dict(
            original_GIFT_flora_species_not_used_as_mainland_training_same_name=True,
            complete_target_full_flora_name_exclusion_hardens_against_memorization=True,
            all_floras_global_parent_4505_not_just_928_complete_three_trait=True,
            independent_candidate_mainland_GIFT_lists_not_ancestral_population_phenotypes=True,
            original_single_species_wide_Direct_trait_ledger_shared_by_both_sources=True,
            other_region_conditional_genus_predictor_KNOWS_focal_island_genus_composition=True,
            actual_source_mainland_holdout_compares_only_colour_and_form_pair_not_all_13_traits=True,
            distinct_cohort_and_endpoint_scores_cannot_be_directly_subtracted=True,
            global_4505_all_flora_name_holdout_leaves_ZERO_independent_source_direct_pair_donors=True,
            cohort_2324_name_holdout_leaves_NINE_source_direct_pair_donors_but_novelty_not_ancestry=True,
            native_only_and_short_distance_source_results_are_selective_subsets_no_global_generalization=True,
            no_original_source_ancestral_seed_dispersal_or_pollinator_causal_interpretation=True,
            no_private_flora_island_source_region_genus_species_or_location_id_exported=True,
        ))
