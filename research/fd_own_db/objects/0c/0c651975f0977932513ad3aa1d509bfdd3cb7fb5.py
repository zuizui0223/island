"""Synthetic original source oracle-free global no-genera-placebo support."""
from itertools import product
from copy import deepcopy
import pytest

from fd.oracle_free_assembly_evidence import (
    REGIONS,SCOPES,MODES,TRAIT_NAMES,QUALIFIED,GOOD_GROUP,
    audit_no_oracle_source_flora,
)


def _group(sign=1):
    out=dict(status=GOOD_GROUP,
        n_source_independent_groups=15,
        n_source_scored_floras=35,
        original_group_bootstrap_draws=499)
    for prefix in ("genus_mix_vs_external_blind",
                   "genus_mix_vs_local_blind",
                   "oracle_known_genus_vs_genus_mix"):
        for metric in ("Brier","logloss"):
            mid=(.03 if metric=="Brier" else .08)*sign
            out[prefix+"_"+metric+"_group_equal_mean"]=mid
            out[prefix+"_"+metric+"_group_resample_95"]=(
                [.02,.04] if sign>0 and metric=="Brier" else
                [.06,.11] if sign>0 else
                [-.04,-.02] if metric=="Brier" else [-.11,-.06])
    return out


def _fixture():
    rows=[]
    for scope,trait,reg,mode in product(SCOPES,TRAIT_NAMES,REGIONS,MODES):
        rows.append(dict(
            scope=scope,trait_endpoint=trait,region=reg,
            training_source_group_deletion_mode=mode,
            original_full_parent_GIFT_floras=4505,
            n_original_islands_with_direct_trait_n5=50,
            n_original_island_floras_scored_without_actual_focal_genus_roster=35,
            n_missing_source_block_or_source_group=1,
            n_no_original_local_or_external_name_disjoint_donors=2,
            n_original_eligible_floras_insufficient_source_pool_genus_overlap=12,
            status=QUALIFIED,
            physical_10degree_block=_group(),
            first_candidate_mainland_source_group=_group(),
        ))
    return dict(
        status="ORIGINAL_GIFT_4505_NO_FOCAL_GENUS_ORACLE_GEOGRAPHIC_HOLDOUT_FLORAL_PROFILES",
        n_original_botanical_island_floras=4505,
        n_original_direct_endpoint_families=13,
        n_actual_original_direct_endpoints_tested=13,
        source_training_modes=list(MODES),
        region_trait_status_source_group_results=rows,
        conceptual_limits={k:True for k in (
            "original_all_4505_botanical_parent_retained",
            "original_observed_focal_genus_composition_NOT_used_in_no_oracle_predictor",
            "original_observed_focal_genus_identity_used_ONLY_in_explicit_oracle_comparator",
            "original_training_flora_local_genus_weights_learned_from_source_DIRECT_SPECIES_COUNTS_not_actual_island_colonization",
            "all_local_and_external_training_accepted_species_names_in_focal_sameblock_firstsource_groups_removed",
            "source_plant_direct_states_from_other_GIFT_ISLANDS_not_independent_MAINLAND_ancestors",
            "original_13_direct_trait_categories_and_native_statuses_not_independent_replicates",
        )},
    )


def test_source_4505_parent_and_all_13_true_direct_endpoints():
    x=_fixture()
    q=audit_no_oracle_source_flora(x)
    assert q["status"]=="ORIGINAL_GIFT_4505_ORACLE_FREE_REGIONAL_GENUS_MIXTURE_VS_TWO_BLIND_BASELINES"
    assert q["n_original_source_botanical_floras"]==4505
    assert len(q["source_original_13_endpoint_region_scope_training_mode_results"])==312
    assert len(q["source_region_scope_mode_oracle_free_prediction_counts"])==24
    assert len(q["all_three_geography_holdout_intersection"]["all_recorded__external_blind"])==13
    assert q["scientific_ceilings"]["original_target_genus_roster_never_available_to_no_oracle_predictor"]
    assert q["scientific_ceilings"]["bootstrap_resampled_original_source_not_multiplicity_adjusted_pvalue"]


def test_both_groupings_and_both_proper_losses_are_required():
    x=_fixture()
    r=x["region_trait_status_source_group_results"][0]
    # first source logloss ambiguous -> no positive effect, even if Brier positive
    r["first_candidate_mainland_source_group"][
        "genus_mix_vs_external_blind_logloss_group_resample_95"]=[-.01,.11]
    q=audit_no_oracle_source_flora(x)
    rec=q["source_original_13_endpoint_region_scope_training_mode_results"][0]
    assert rec["gains"]["external_blind"]["status"]=="ambiguous_or_loss_disagreement"
    assert rec["gains"]["local_blind"]["status"]=="positive"
    x2=_fixture()
    row=x2["region_trait_status_source_group_results"][0]
    row["physical_10degree_block"]=_group(-1)
    row["first_candidate_mainland_source_group"]=_group(-1)
    r2=audit_no_oracle_source_flora(x2)
    assert r2["source_original_13_endpoint_region_scope_training_mode_results"][0][
        "gains"]["external_blind"]["status"]=="negative"


def test_original_population_and_source_group_attrition_do_not_silently_change():
    x=_fixture()
    x["n_original_botanical_island_floras"]=928
    with pytest.raises(ValueError,match="4,505"):
        audit_no_oracle_source_flora(x)
    y=_fixture()
    y["region_trait_status_source_group_results"][0][
        "n_original_eligible_floras_insufficient_source_pool_genus_overlap"]=13
    with pytest.raises(ValueError,match="attrition"):
        audit_no_oracle_source_flora(y)
    z=_fixture()
    z["region_trait_status_source_group_results"].pop()
    with pytest.raises(ValueError,match="incomplete"):
        audit_no_oracle_source_flora(z)
    bad=_fixture()
    bad["conceptual_limits"][
        "original_observed_focal_genus_composition_NOT_used_in_no_oracle_predictor"]=False
    with pytest.raises(ValueError,match="claim ceiling"):
        audit_no_oracle_source_flora(bad)
