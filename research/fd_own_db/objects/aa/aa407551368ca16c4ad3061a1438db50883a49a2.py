"""Synthetic checks for full-parent interregional genus transfer audit."""
import copy
from itertools import product
import pytest

from fd.interregional_transfer_robustness import (
    REGIONS,SCOPES,TRAITS,VALID,directions,
    analyze_original_global_regionality,validate_source,
    verify_two_original_runs,
)


def _stats(direction=1):
    return {
        "n_clusters":18,
        "n_scored_islands":30,
        "local_minus_external_proper_Brier_gain_per_comparable_plant_group_equal_mean":
            0.05*direction,
        "local_minus_external_proper_Brier_gain_per_comparable_plant_group_resample_95":
            [0.03,0.07] if direction>0 else [-0.07,-0.03],
        "local_minus_external_proper_logloss_gain_per_comparable_plant_group_equal_mean":
            0.12*direction,
        "local_minus_external_proper_logloss_gain_per_comparable_plant_group_resample_95":
            [0.09,0.16] if direction>0 else [-0.16,-0.09],
        "external_genus_vs_blind_proper_Brier_gain_per_comparable_plant_group_equal_mean":
            0.02,
        "external_genus_vs_blind_proper_Brier_gain_per_comparable_plant_group_resample_95":
            [0.01,0.03],
    }


def original_synthetic_global_fixture():
    rows=[]
    source_support=[]
    for scope,trait,region in product(SCOPES,TRAITS,REGIONS):
        source_support.append(dict(scope=scope,trait_endpoint=trait,region=region))
        if len(rows)>=100:continue
        qualified=len(rows)<63
        sign=-1 if region=="tropical" else 1
        rows.append(dict(scope=scope,trait_endpoint=trait,region=region,
            status=VALID if qualified else "NOT_IDENTIFIABLE_INDEPENDENT_ORIGINAL_GROUP_SUPPORT",
            n_original_region_botanical_floras=2205 if region=="northern_midlatitude" else 1558,
            n_same_focal_islands_matched_carrier_n5=35,
            n_n5_islands_considered=45,
            n_original_n5_direct_trait_islands=80,
            n_original_matched_focal_Direct_carrier_records=120,
            mean_original_matched_carrier_fraction=0.40,
            physical_10degree_block=_stats(sign) if qualified else None,
            first_GIFT_mainland_source=_stats(sign) if qualified else None,
            locally_exclusive_minus_externally_exclusive_frequency_D2_gain_island_equal_mean=0.004,
            local_minus_external_proper_Brier_gain_per_comparable_plant_island_equal_mean=0.05*sign,
            local_minus_external_proper_logloss_gain_per_comparable_plant_island_equal_mean=0.12*sign,
        ))
    claims={k:True for k in (
        "entire_4505_original_botanical_population_preserved",
        "target_original_all_accepted_plant_names_removed_from_both_donor_sets",
        "names_occurring_in_target_and_other_regions_removed_from_both_source_sets",
        "local_and_external_same_original_genus_and_exact_equal_donor_species_count",
        "donor_species_selected_by_original_name_hash_not_by_trait_state",
        "strict_native_source_status_is_nested_nonindependent_scope",
        "no_independent_mainland_or_phylogenetic_trait_source_or_adaptation_test",
        "no_real_pollinator_service_guild_replacement_or_fitness_claim",
    )}
    return dict(status="ORIGINAL_4505_REGION_EXCLUSIVE_EQUAL_GENUS_DONOR_TRANSFER",
        n_original_botanical_island_floras=4505,
        n_original_accepted_source_trait_endpoints=13,
        original_complete_botanical_region_universe=[
            dict(region="northern_midlatitude",n_original_botanical_floras=2205),
            dict(region="tropical",n_original_botanical_floras=1558),
            dict(region="southern_extratropical",n_original_botanical_floras=317),
            dict(region="northern_high_latitude",n_original_botanical_floras=424),
            dict(region="geography_unresolved",n_original_botanical_floras=1)],
        original_exclusive_training_taxonomic_source_support=source_support,
        original_equal_genus_donor_transfer_results=rows,claims=claims)


def test_source_provenance_and_nonindependent_complete_parent():
    q=original_synthetic_global_fixture()
    check=validate_source(q)
    assert sum(check["original_region_counts"].values())==4505
    assert check["source_qualified_count_by_status"][VALID]==63
    with pytest.raises(ValueError,match="source drift"):
        validate_source(dict(q,n_original_botanical_island_floras=928))
    b=copy.deepcopy(q)
    b["claims"]["target_original_all_accepted_plant_names_removed_from_both_donor_sets"]=False
    with pytest.raises(ValueError,match="noncausal contract"):
        validate_source(b)


def test_direction_must_agree_across_proper_losses_and_both_independent_groupings():
    row=dict(status=VALID,
        physical_10degree_block=_stats(-1),
        first_GIFT_mainland_source=_stats(-1))
    q=directions(row)
    assert q["two_proper_loss_direction"]=="external_better"
    assert q["external_genus_Brier"]["direction"]=="cross_region_genus_better"
    bad=copy.deepcopy(row)
    bad["first_GIFT_mainland_source"]["local_minus_external_proper_logloss_gain_per_comparable_plant_group_resample_95"]=[-0.2,+0.01]
    assert directions(bad)["two_proper_loss_direction"]=="not_robust_two_metrics"
    bad=copy.deepcopy(row)
    bad["first_GIFT_mainland_source"]["n_clusters"]=7
    assert directions(bad)["two_proper_loss_direction"]=="not_robust_two_metrics"


def test_compare_two_replays_same_original_data_not_two_independent_field_studies():
    a=original_synthetic_global_fixture()
    b=copy.deepcopy(a)
    assert verify_two_original_runs(a,b)["n_original_matched_cells"]==100
    b["original_equal_genus_donor_transfer_results"][0][
        "local_minus_external_proper_Brier_gain_per_comparable_plant_island_equal_mean"]+=1e-5
    with pytest.raises(ValueError,match="drift"):
        verify_two_original_runs(a,b)


def test_full_regional_temperature_or_butterfly_conclusion_not_inferred_from_genus_oracle():
    a=original_synthetic_global_fixture()
    q=analyze_original_global_regionality(a,copy.deepcopy(a))
    assert q["original_n_global_botanical_floras"]==4505
    assert q["n_original_highmedium_direct_endpoints"]==13
    assert q["scientific_ceilings"][
        "original_regional_local_donors_from_ISLAND_regional_floristic_lists_not_independent_mainland_ancestral_observations"]
    assert q["scientific_ceilings"]["original_focal_botanical_genus_composition_is_ORACLE_in_prediction"]
    assert q["independent_execution_same_original_data_replay"][
        "deterministic_same_input_two_CI_runs_not_independent_field_replication"]
    assert q["original_north_to_tropical_double_metric_sign_reversal"]


def original_synthetic_mainland_candidate():
    return dict(
        status="ORIGINAL_GIFT_4505_MAINLAND_SOURCE_NOVEL_SPECIES_HOLDOUT_VS_INTERREGIONAL_GENUS_ORACLE",
        parent_original_global_GIFT_botanical_island_floras=4505,
        n_original_mainland_GIFT_native_source_region_lists=748,
        n_original_global_all_recorded_direct_colour_form_test_islands=2324,
        n_original_direct_colour_form_source_species_before_strict_holdout=3606,
        n_original_source_colour_form_species_outside_ALL_4505_GIFT_island_floras=0,
        n_original_candidate_mainland_donors_after_2324_trait_eligible_island_full_name_holdout=9,
        original_source_provenance_limit=dict(
            global_4505_all_flora_name_holdout_leaves_ZERO_independent_source_direct_pair_donors=True,
            cohort_2324_name_holdout_leaves_NINE_source_direct_pair_donors_but_novelty_not_ancestry=True,
            original_single_species_wide_Direct_trait_ledger_shared_by_both_sources=True,
            other_region_conditional_genus_predictor_KNOWS_focal_island_genus_composition=True,
            independent_candidate_mainland_GIFT_lists_not_ancestral_population_phenotypes=True,
        ),
        original_actual_candidate_mainland_GIFT_source_holdout_cohorts=[
            dict(source_cohort="global_all_n5",
                n_candidate_mainland_direct_species_after_ALL_cohort_full_flora_names_removed=9,
                n_original_source_donor_genera_with_two_plus_joint_category_states=0),
        ],
    )


def test_mainland_candidate_flora_is_a_different_reference_from_other_island_regional_genus_training():
    from fd.interregional_transfer_robustness import (
        audit_actual_candidate_mainland_source_novelty)
    a=original_synthetic_mainland_candidate()
    q=audit_actual_candidate_mainland_source_novelty(a)
    assert q["parent_global_GIFT_island_botanical_floras"]==4505
    assert q["n_original_direct_colour_form_mainland_species_names_outside_ALL_4505_island_floras"]==0
    assert q["n_candidate_source_species_after_2324_full_flora_names_heldout"]==9
    assert q["n_candidate_source_genera_with_true_pair_variation_after_2324_holdout"]==0
    b=original_synthetic_global_fixture()
    r=analyze_original_global_regionality(b,copy.deepcopy(b),a)
    assert r["independent_mainland_candidate_species_provenance_gate"][
        "region_exclusive_other_GIFT_island_species_source_transfer_is_NOT_independent_mainland_source"]
    a["n_original_source_colour_form_species_outside_ALL_4505_GIFT_island_floras"]=1
    with pytest.raises(ValueError,match="independent mainland names"):
        audit_actual_candidate_mainland_source_novelty(a)
