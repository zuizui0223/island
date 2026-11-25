"""Synthetic-only tests for 78 original GIFT residual source partitions."""
import copy
import hashlib
import json
import zipfile

import pytest

from scripts.verify_gift_flower_residual_split import (
    STATUS,ENDPOINTS,REGIONS,SCOPES,read_original,verify
)


def synthetic_input():
    source=[]
    before=[]
    for reg in REGIONS:
        for scope in SCOPES:
            for ep in sorted(ENDPOINTS):
                geos={}
                old_geos={}
                for label in ("physical_block","first_source_group"):
                    part={"supported_for_inference":True}
                    prev={"supported_for_inference":True}
                    for baseline in ("regional_blind","external_blind","regional_genus_prior"):
                        for loss in ("Brier","logloss"):
                            metric=f"environment_source_nearest_gain_vs_{baseline}_{loss}"
                            part[metric]={"group_equal_mean":.03,"bootstrap_95":[.01,.05]}
                            prev[metric]={"group_equal_mean":.03,"bootstrap_95":[.01,.05]}
                    for stem,squared,proj in [
                        ("taxonomic_genus_forecast_residual",.10,.05),
                        ("direct_trait_documentation_residual",.08,.06),
                        ("within_genus_documented_trait_residual",.12,.09)
                    ]:
                        part[stem+"_squared_norm"]={"group_equal_mean":squared,
                                                   "bootstrap_95":[squared/2,squared*2]}
                        part[stem+"_signed_projection_on_total_error"]={
                            "group_equal_mean":proj,"bootstrap_95":[proj/2,proj*2]}
                    part["observed_minus_predicted_direct_trait_squared_norm"]={
                        "group_equal_mean":.20,"bootstrap_95":[.15,.25]}
                    part["observed_direct_trait_plugin_Q"]={
                        "group_equal_mean":.45,"bootstrap_95":[.4,.5]}
                    geos[label]=part;old_geos[label]=prev
                common={
                    "scope":scope,"region":reg,"trait_endpoint":ep,
                    "n_original_target_islands_scored":50,
                    "original_parent_botanical_floras":4505,
                    "original_full_accepted_name_group_ban_on_both_models":True,
                    "original_focal_genus_oracle_NOT_used_in_predictor":True,
                    "support_both_independent_source_groups":True,
                }
                source.append({**common,**geos,
                    "source_residual_projection_is_not_a_causal_selection_partition":True,
                    "whole_flora_and_direct_carrier_genus_used_only_after_prediction":True})
                before.append({**common,**old_geos})
    root={
        "original_parent_GIFT_botanical_floras":4505,
        "frozen_pilot_endpoints":sorted(ENDPOINTS),
    }
    return (
        {**root,"status":STATUS,
         "direct_residual_vector_identity_verified_per_original_target_island":True,
         "residual_projection_signed_not_causal_explained_variance":True,
         "observed_true_focal_genera_used_only_in_post_prediction_audit":True,
         "regional_source_heldout_results":source},
        {**root,"regional_source_heldout_results":before},
    )


def test_source_original_complete_78_and_all_group_component_projections():
    x,old=synthetic_input()
    verified=verify(x,old)
    assert verified["n_anonymous_source_panels"]==78
    assert verified["old_original_proper_category_scores_unchanged"]
    assert verified["counts"][
        "all_recorded|tropical|physical_block|largest_component|within_genus_documented_trait_residual"]==13
    assert all(a["physical_block_squared_error"]==pytest.approx(.2)
               for a in verified["source_anonymous_panels"])


def test_nonadditive_component_squares_must_not_be_wrongly_interpreted_as_variance_parts():
    x,old=synthetic_input()
    a=x["regional_source_heldout_results"][0]["physical_block"]
    total=a["observed_minus_predicted_direct_trait_squared_norm"]["group_equal_mean"]
    sum_norm=sum(a[f"{k}_squared_norm"]["group_equal_mean"] for k in (
        "taxonomic_genus_forecast_residual",
        "direct_trait_documentation_residual",
        "within_genus_documented_trait_residual"))
    assert sum_norm == pytest.approx(.30)
    assert total == pytest.approx(.20)
    assert verify(x,old)["n_anonymous_source_panels"]==78


def test_falsified_total_signed_projection_and_impossible_values_raise():
    x,old=synthetic_input()
    x["regional_source_heldout_results"][0]["physical_block"][
        "taxonomic_genus_forecast_residual_signed_projection_on_total_error"
    ]["group_equal_mean"]=-2
    with pytest.raises(ValueError,match="do not add"):
        verify(x,old)
    x,old=synthetic_input()
    x["regional_source_heldout_results"][0]["physical_block"][
        "observed_direct_trait_plugin_Q"]["group_equal_mean"]=1.2
    with pytest.raises(ValueError,match="outside valid interval"):
        verify(x,old)


def test_source_nooracle_model_or_target_cohort_change_rejected():
    x,old=synthetic_input()
    old["regional_source_heldout_results"][0][
        "n_original_target_islands_scored"]=51
    with pytest.raises(ValueError,match="cohort"):
        verify(x,old)
    x,old=synthetic_input()
    old["regional_source_heldout_results"][0]["first_source_group"][
        "environment_source_nearest_gain_vs_regional_blind_Brier"][
            "group_equal_mean"]=100.
    with pytest.raises(ValueError,match="proper scores changed"):
        verify(x,old)


def test_missing_source_panel_and_protected_botanical_identity_fail_closed(tmp_path):
    x,old=synthetic_input()
    x["regional_source_heldout_results"].pop()
    with pytest.raises(ValueError,match="13×3×2"):
        verify(x,old)
    x,old=synthetic_input()
    z=tmp_path/"source.zip"
    with zipfile.ZipFile(z,"w") as archive:
        archive.writestr("anonymous.json",json.dumps(x))
    digest=hashlib.sha256(z.read_bytes()).hexdigest()
    assert read_original(z,sha256=digest)["status"]==STATUS
    with pytest.raises(ValueError,match="digest"):
        read_original(z,sha256="0"*64)
    x["accidental"]={"accepted_species":"Private genus and species must never export"}
    raw=tmp_path/"bad.json"
    raw.write_text(json.dumps(x))
    with pytest.raises(ValueError,match="Protected"):
        read_original(raw)
