"""Mainland candidate source vs island-to-island genus oracle must not be conflated."""
import copy
import pytest
from fd.mainland_source_novel_donor_gate import (
    audit_candidate_mainland_reality,COHORTS,
)


def _fixtures():
    testsize={
        "global_all_n5":(2324,9,0),
        "global_original_native_n5":(327,679,68),
        "northmid_n10_coverage10":(1055,1561,177),
        "northmid_500_2000km_exact456_endpoint_union":(215,2899,330),
    }
    rows=[]
    for k in COHORTS:
        n,remain,variable=testsize[k]
        rows.append(dict(
           cohort=k,status="POSTHOC_OUT_OF_TARGET_NAME_SOURCE_PREDICTION",
           n_original_focal_islands=n,
           n_original_focal_full_flora_names=n*3,
           n_unique_source_direct_pair_donor_species_before_holdout=3606,
           n_unique_source_direct_pair_donor_species_after_exact_name_holdout=remain,
           n_donor_genera_with_multiple_joint_profiles=variable,
           ratio_test_carriers_with_any_global_same_genus_donor=0.30 if remain>9 else 0.0038,
           ratio_test_carriers_with_same_genus_local_five_source_donor=0.09 if remain>9 else 0.0004,
           source_local_minus_global_genus_improvement_island_brier=0.001 if remain>9 else 0.0,
           source_local_minus_global_genus_improvement_island_log_loss=0.003 if remain>9 else 0.0))
    mainland=dict(
        status="ORIGINAL_STRICT_GIFT_SOURCE_NAME_HOLDOUT_PREFLIGHT",
        n_original_botanical_islands=4505,
        n_original_distinct_source_regions=748,
        n_original_mainland_native_source_rows=1666186,
        n_source_direct_pair_species_after_ALL_4505_island_name_exclusion=0,
        claims={
          "full_4505_botanical_universe_retained":True,
          "direct_species_flower_colour_and_form_joint_not_imputed":True,
          "same_name_removed_from_ALL_mainland_regions_before_training":True,
          "genus_known_is_a_posthoc_oracle_not_colonization_prediction":True,
          "no_double_counting_of_local_donors_overlapping_global_genus":True,
        },cohorts=rows)
    transfer=dict(
       status="ORIGINAL_4505_GIFT_THREE_REGION_TWO_LOSS_TWO_GROUP_GENUS_TRANSFER_REPLICATION",
       n_original_botanical_global_island_floras=4505,
       n_original_direct_endpoints=13,
       n_all_recorded_endpoints_replicating_across_all_three_geographies=8,
       n_native_only_endpoints_replicating_across_all_three_geographies=0,
       strict_approach={
           "source_training_other_regions_island_floras_not_proven_ancestral_mainland_species":True})
    return mainland,transfer


def test_original_global_full_flora_name_exclusion_leaves_no_variable_mainland_genus():
    m,t=_fixtures()
    r=audit_candidate_mainland_reality(m,t)
    assert r["parent_original_global_GIFT_botanical_island_floras"]==4505
    assert r["n_original_source_colour_form_species_outside_ALL_4505_GIFT_island_floras"]==0
    assert r["n_original_candidate_mainland_donors_after_2324_trait_eligible_island_full_name_holdout"]==9
    assert r["conditional_other_region_GIFT_island_genus_predictive_endpoints_across_three_regions"]==8
    rows={q["source_cohort"]:q for q in r[
        "original_actual_candidate_mainland_GIFT_source_holdout_cohorts"]}
    assert rows["global_all_n5"][
       "n_original_source_donor_genera_with_two_plus_joint_category_states"]==0
    assert rows["global_all_n5"][
       "local_mainland_genus_minus_global_source_genus_Brier_gain_descriptive"]==0
    assert rows["global_original_native_n5"][
       "n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed"]==679
    assert r["original_source_provenance_limit"][
        "distinct_cohort_and_endpoint_scores_cannot_be_directly_subtracted"]


def test_unpinned_ancestral_source_and_same_species_training_cannot_pass():
    m,t=_fixtures()
    bad=copy.deepcopy(m)
    bad["claims"]["same_name_removed_from_ALL_mainland_regions_before_training"]=False
    with pytest.raises(ValueError,match="holdout contract"):
        audit_candidate_mainland_reality(bad,t)
    bad=copy.deepcopy(m)
    bad["n_source_direct_pair_species_after_ALL_4505_island_name_exclusion"]=1
    with pytest.raises(ValueError,match="ZERO"):
        audit_candidate_mainland_reality(bad,t)
    bad=copy.deepcopy(t)
    bad["strict_approach"][
        "source_training_other_regions_island_floras_not_proven_ancestral_mainland_species"]=False
    with pytest.raises(ValueError,match="source identity"):
        audit_candidate_mainland_reality(m,bad)


def test_mismatched_population_or_invalid_source_denominator_blocked():
    m,t=_fixtures()
    bad=copy.deepcopy(m)
    bad["n_original_botanical_islands"]=928
    with pytest.raises(ValueError,match="4,505"):
        audit_candidate_mainland_reality(bad,t)
    bad=copy.deepcopy(m)
    bad["cohorts"][0]["n_unique_source_direct_pair_donor_species_after_exact_name_holdout"]=3900
    with pytest.raises(ValueError,match="opportunities invalid"):
        audit_candidate_mainland_reality(bad,t)
