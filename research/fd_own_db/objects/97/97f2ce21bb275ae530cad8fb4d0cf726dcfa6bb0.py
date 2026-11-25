"""GIFT 4505 geographic genus source provenance and negative scores, synthetic."""
from copy import deepcopy
from itertools import product
import pytest

from fd.environment_genus_source_audit import (
    audit_geo_source_genus_prediction,MODELS,PAIRED,REGIONS,SCOPES,QUALIFIED,
)


def fixture():
    region_n=dict(zip(REGIONS,[2205,1558,317,424]))
    cells=[]
    for scope,region in product(SCOPES,REGIONS):
        count=region_n[region]
        sample=60 if scope=="all_recorded" else 40
        group={"status":QUALIFIED,
             "n_original_source_scored_floras":sample,"n_source_groups":12}
        for m in MODELS:
            for metric in ("Brier","logloss"):
                val=(-.02 if m=="geo_nearest" else .03)
                group[m+"_minus_region_equal_"+metric+"_group_equal_mean"]=val
                group[m+"_minus_region_equal_"+metric+"_group_resample_95"]=(
                    [-.04,-.01] if val<0 else [.01,.05])
        for label,ref,alternative in PAIRED:
            for loss in ("Brier","logloss"):
                score=lambda model: -.02 if model=="geo_nearest" else .03
                delta=score(alternative)-score(ref)
                group[label+"_"+loss+"_group_equal_mean"]=delta
                group[label+"_"+loss+"_group_resample_95"]=(
                    [delta-.01,delta+.01])
        cells.append(dict(
            scope=scope,region=region,status=
                "ORIGINAL_GIFT_GENUS_OCCUPANCY_COMPOSITION_ENVIRONMENTAL_TEST_BOTH_SOURCE_GROUPS",
            n_original_region_botanical_floras=count,
            n_original_source_islands_with_full_geo_and_genus_baseline=sample+10,
            n_original_floras_scored_for_genus_community_environmental_prediction=sample,
            n_original_focal_islands_no_complete_geographic_or_full_genus_record=count-sample-10,
            n_unscored_donor_coverage_by_reason=dict(
                n_focal_islands_lacking_10_donor_source_groups=10),
            physical_10degree_block=deepcopy(group),
            first_candidate_source_entity_group=deepcopy(group),
        ))
    controls={s:True for s in (
      "original_4505_island_botanical_flora_population_retained",
      "original_focal_true_genus_community_not_used_to_predict",
      "original_true_accepted_species_genus_used_ONLY_FOR_SCORING",
      "original_training_GIFT_islands_in_focal_10degree_block_EXCLUDED",
      "original_training_GIFT_islands_with_same_first_GIFT_candidate_mainland_source_EXCLUDED",
      "GIFT_source_donor_island_taxonomy_NOT_independent_mainland_ancestor",
      "source_species_normalized_by_island_and_UNSEEN_GENUS_held_in_common_unknown_category",
      "source_taxon_genus_nonlisting_is_NOT_confirmed_ecological_absence",
      "no_true_island_genus_colonization_probability_OR_independent_trait_adaptation",
      "no_pollinator_guild_selection_or_seedset_causality_claim",
      "source_group_bootstrap_exploratory_not_multiplicity_adjusted",
    )}
    return dict(
      status="ORIGINAL_GIFT_4505_NO_FOCAL_GENUS_SOURCE_ISLAND_ENVIRONMENTAL_GENUS_ASSEMBLY",
      original_GIFT_parent_island_floras=4505,
      source_genus_models=list(MODELS),
      controls=controls,
      original_no_oracle_original_island_genus_composition_source_predictions=cells,
    )


def test_four_regions_two_nested_statuses_and_no_genus_oracle():
    x=audit_geo_source_genus_prediction(fixture())
    assert x["n_original_global_botanical_floras"]==4505
    assert len(x["original_full_geographic_composition_cells"])==8
    assert len(x["original_macroregion_proper_loss_source_count_summary"])==8
    assert x["limitations"]["original_focal_botanical_genus_input_never_used"]
    assert x["limitations"]["original_genus_community_source_not_yet_flora_trait_prediction"]
    assert x["limitations"]["original_source_similarity_increment_TESTED_WITH_SAME_TARGET_AND_DIRECT_PAIRED_GROUP_RESAMPLING"]
    assert len(x["original_paired_model_direct_source_increment_summary"])==8


def test_proper_brier_logloss_and_BOTH_geographical_groupings_needed():
    z=fixture()
    q=audit_geo_source_genus_prediction(z)
    first=q["original_full_geographic_composition_cells"][0]
    assert first["original_source_models"]["geo_nearest"]["status"]=="predictor_worse_than_region_equal"
    assert first["original_source_models"]["environment_nearest"]["status"]=="predictor_better_than_region_equal"
    bad=deepcopy(z)
    bad["original_no_oracle_original_island_genus_composition_source_predictions"][0][
        "first_candidate_source_entity_group"]["environment_nearest_minus_region_equal_logloss_group_resample_95"]=[-.01,.04]
    out=audit_geo_source_genus_prediction(bad)
    assert out["original_full_geographic_composition_cells"][0][
        "original_source_models"]["environment_nearest"]["status"]=="ambiguous_or_disagreeing_losses_or_source_groups"


def test_no_unsupported_source_panels_become_positive():
    x=fixture()
    x["original_GIFT_parent_island_floras"]=928
    with pytest.raises(ValueError,match="4505"):
        audit_geo_source_genus_prediction(x)
    q=fixture()
    q["original_no_oracle_original_island_genus_composition_source_predictions"][0][
        "n_original_floras_scored_for_genus_community_environmental_prediction"]+=2
    with pytest.raises(ValueError,match="waterfall"):
        audit_geo_source_genus_prediction(q)
    q=fixture()
    q["controls"]["source_taxon_genus_nonlisting_is_NOT_confirmed_ecological_absence"]=False
    with pytest.raises(ValueError,match="ceiling missing"):
        audit_geo_source_genus_prediction(q)


def test_direct_PAIRED_original_source_increment_both_losses_and_source_groups_and_algebra():
    q=audit_geo_source_genus_prediction(fixture())
    x=q["original_full_geographic_composition_cells"][0]
    src=x["original_source_paired_predictor_increments"]
    assert src["source_similarity_increment_over_environment"]["status"]=="paired_ambiguous_or_loss_disagreement"
    assert src["environment_increment_over_geographic_nearness"]["status"]=="paired_predictor_better"
    assert src["environment_increment_over_geographic_nearness"]["proper_paired_difference_source_groups"][
        "Brier__physical_10degree_block"]["value"]==pytest.approx(.05)
    fail=fixture()
    z=fail["original_no_oracle_original_island_genus_composition_source_predictions"][0][
        "physical_10degree_block"]
    z["source_similarity_increment_over_environment_Brier_group_equal_mean"]+=.01
    with pytest.raises(ValueError,match="paired source environmental improvement failed"):
        audit_geo_source_genus_prediction(fail)
