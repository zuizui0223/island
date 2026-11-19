"""Original GIFT strict training-exclusion cross-regional signal controls (synthetic)."""
from itertools import product
from copy import deepcopy
import pytest

from fd_strict_external_genus_blind_dual_score_audit import (
    MODES,REGIONS,SCOPES,MAIN_REGIONS,
    QUALIFIED,GROUP_QUALIFIED,
    audit_genus_blind_dual_scores,
)


ENDPOINTS=(
    "floral_form","floral_symmetry","flower_primary_color",
    "flower_size_class","inflorescence_display","tube_depth_class",
    "self_incompatibility","mating_system","autonomous_selfing_capacity",
    "cleistogamy","flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
)


def _group(brier=.04,logloss=.08,*,n=18):
    return dict(
        status=GROUP_QUALIFIED,
        n_scored_island_floras=55,
        n_distinct_original_source_clusters=n,
        strict_other_region_genus_vs_blind_Brier_group_equal_mean=brier,
        strict_other_region_genus_vs_blind_Brier_group_resample_95=(
            [.02,.06] if brier>0 else [-.06,-.02]),
        strict_other_region_genus_vs_blind_logloss_group_equal_mean=logloss,
        strict_other_region_genus_vs_blind_logloss_group_resample_95=(
            [.04,.12] if logloss>0 else [-.12,-.04]),
        baseline_on_identical_carriers_other_region_genus_vs_blind_Brier_group_equal_mean=.03,
        baseline_on_identical_carriers_other_region_genus_vs_blind_Brier_group_resample_95=[.01,.05],
        baseline_on_identical_carriers_other_region_genus_vs_blind_logloss_group_equal_mean=.06,
        baseline_on_identical_carriers_other_region_genus_vs_blind_logloss_group_resample_95=[.03,.09],
    )


def _fixture():
    names=list(product(MODES,SCOPES,REGIONS,ENDPOINTS))
    assert len(names)==312
    rows=[]
    for i,(mode,scope,region,trait) in enumerate(names):
        qualified=i<178
        row=dict(
            original_scope=scope,original_trait_endpoint=trait,
            original_source_region=region,
            source_training_exclusion_mode=mode,
            status=QUALIFIED if qualified else "NOT_IDENTIFIABLE_INDEPENDENT_ORIGINAL_CLUSTER_SUPPORT",
            n_original_original_GIFT_islands_with_strict_common_focal_scores=55 if qualified else 7,
            n_original_true_direct_n5_eligible_focal_islands=80,
        )
        if qualified:
            row["physical_10degree_block"]=_group()
            row["first_candidate_source_group"]=_group()
        rows.append(row)
    output=dict(
        status="ORIGINAL_GIFT_4505_GLOBAL_GEO_SOURCE_TRAINING_EXCLUDED_REGION_TRAIT_TRANSFER",
        n_original_botanical_GIFT_island_floras=4505,
        n_original_floral_direct_endpoint_families=13,
        n_actual_trait_endpoint_families_tested=13,
        original_source_training_exclusion_modes=list(MODES),
        original_source_training_eligibility_and_scoring=rows,
        original_source_training_scored_results=[r for r in rows if
            r["status"]==QUALIFIED],
        scientific_limitations={x:True for x in (
            "focal_accepted_full_flora_names_excluded_from_all_original_source_donor_profiles",
            "SAME_physical_block_OTHER_island_species_names_removed_from_local_and_external_training",
            "SAME_original_first_mainland_source_OTHER_island_species_names_removed_from_local_and_external_training",
            "rescore_baseline_on_IDENTICAL_strict_source_supported_accepted_species",
            "original_training_species_borrowed_from_GIFT_island_floras_NOT_independently_phenotyped_MAINLAND",
            "original_plant_genus_composition_given_as_ORACLE_not_predicted")},
    )
    return output


def test_4505_full_source_masks_and_strict_provenance():
    q=_fixture()
    z=audit_genus_blind_dual_scores(q)
    assert z["original_global_botanical_source"]["original_full_parent_island_floras"]==4505
    assert z["original_global_botanical_source"]["original_region_status_trait_mode_eligibility_panels"]==312
    assert z["original_global_botanical_source"]["n_original_source_qualified_region_status_modes"]==178
    assert z["scientific_claim_ceiling"]["actual_focal_island_genus_roster_PROVIDED_as_oracle"]
    bad=deepcopy(q)
    bad["n_original_botanical_GIFT_island_floras"]=928
    with pytest.raises(ValueError,match="Original 4,505"):
        audit_genus_blind_dual_scores(bad)
    bad=deepcopy(q)
    bad["scientific_limitations"]["rescore_baseline_on_IDENTICAL_strict_source_supported_accepted_species"]=False
    with pytest.raises(ValueError,match="source rigorous training exclusion"):
        audit_genus_blind_dual_scores(bad)


def test_both_proper_losses_and_both_group_bootstraps_required():
    q=_fixture()
    x=audit_genus_blind_dual_scores(q)
    row=x["full_region_status_trait_mode_scored_results"][0]
    assert row["strict_other_region_genus_vs_blind_dual_loss_direction"]=="external_genus_better_than_blind"
    a=q["original_source_training_scored_results"][0]
    a["first_candidate_source_group"]["strict_other_region_genus_vs_blind_logloss_group_resample_95"]=[-.02,.10]
    q2=audit_genus_blind_dual_scores(q)
    assert q2["full_region_status_trait_mode_scored_results"][0][
        "strict_other_region_genus_vs_blind_dual_loss_direction"
    ]=="not_robust_across_two_proper_losses_and_two_geographic_cluster_methods"
    a["first_candidate_source_group"]["n_distinct_original_source_clusters"]=6
    q3=audit_genus_blind_dual_scores(q)
    assert q3["full_region_status_trait_mode_scored_results"][0][
        "strict_other_region_genus_vs_blind_dual_loss_direction"
    ]=="not_robust_across_two_proper_losses_and_two_geographic_cluster_methods"


def test_repeated_13_endpoints_are_not_independent_evolution_experiments():
    q=_fixture()
    out=audit_genus_blind_dual_scores(q)
    counts={(x["mode"],x["scope"],x["region"]):x for x in
            out["source_geo_region_status_mode_dual_proper_score_counts"]}
    north=counts[("block_out","all_recorded","northern_midlatitude")]
    assert north["n_original_source_qualified_direct_endpoints"]==13
    assert north["n_both_proper_scores_and_both_groupings_positive"]==13
    triple=out["original_three_region_strict_double_proper_score_transfers"]
    assert any(x["original_direct_endpoint"]=="floral_form" and
               x["mode"]=="block_out" and
               x["scope"]=="all_recorded" for x in triple)
    assert out["scientific_claim_ceiling"][
        "positive_resamples_are_exploratory_not_multiplicity_corrected"]
    assert out["scientific_claim_ceiling"][
        "original_4505_botanical_parent_not_restricted_to_matched_source_subset"]


def test_external_genus_below_blind_is_retained_as_negative_result():
    q=_fixture()
    a=q["original_source_training_scored_results"][0]
    a["physical_10degree_block"]=_group(-.04,-.08)
    a["first_candidate_source_group"]=_group(-.04,-.08)
    x=audit_genus_blind_dual_scores(q)
    assert x["full_region_status_trait_mode_scored_results"][0][
        "strict_other_region_genus_vs_blind_dual_loss_direction"
    ]=="external_genus_worse_than_blind"
    bad=deepcopy(q)
    bad["original_source_training_eligibility_and_scoring"].pop()
    with pytest.raises(ValueError,match="312-panel"):
        audit_genus_blind_dual_scores(bad)
