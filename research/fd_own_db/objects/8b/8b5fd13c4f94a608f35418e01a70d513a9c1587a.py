"""Synthetic CI contracts for GIFT original target accepted-species/genus sensitivity."""
import copy
import hashlib
import json
import zipfile

import pytest

from scripts.verify_gift_target_taxon_source import (
    ENDPOINTS,REGIONS,SCOPES,read_original,audit,signed_pair
)


def synthetic_original():
    allrows=[]
    for scope in SCOPES:
        for region in REGIONS:
            for trait in sorted(ENDPOINTS):
                score={"supported_for_inference":True,"n_original_group_units":37,
                       "resampling_unit":"distinct_accepted_species",
                       "nboot":99}
                genus_score={"supported_for_inference":True,
                             "n_original_group_units":12,
                             "resampling_unit":"distinct_accepted_genus",
                             "nboot":99}
                for model in (score,genus_score):
                    for prior in ("regional_blind","external_blind",
                                  "regional_genus_prior"):
                        for loss in ("Brier","logloss"):
                            model[f"environment_source_nearest_gain_vs_{prior}_{loss}"]={
                                "group_equal_mean":.02,
                                "bootstrap_95":[.01,.03]}
                original_scores={}
                for prior in ("regional_blind","external_blind",
                              "regional_genus_prior"):
                    for loss in ("Brier","logloss"):
                        original_scores[f"environment_source_nearest_gain_vs_{prior}_{loss}"]={
                            "group_equal_mean":.02,"bootstrap_95":[.01,.03]}
                allrows.append({
                    "scope":scope,"region":region,"trait_endpoint":trait,
                    "n_original_target_islands_scored":10,
                    "original_parent_botanical_floras":4505,
                    "original_full_accepted_name_group_ban_on_both_models":True,
                    "original_focal_genus_oracle_NOT_used_in_predictor":True,
                    "support_both_independent_source_groups":True,
                    "target_species_genus_reweighting":{
                        "n_distinct_target_species":37,
                        "n_distinct_target_genera":12,
                        "n_species_island_occurrences":100,
                        "mean_island_occurrences_per_distinct_species":100/37,
                        "target_accepted_names_absent_from_export":True,
                        "species_equal":score,"genus_equal":genus_score},
                    "species_equal_positive_local_external_pooled":{
                        "regional_blind":True,"external_blind":True,
                        "regional_genus_prior":True},
                    "genus_equal_positive_local_external_pooled":{
                        "regional_blind":True,"external_blind":True,
                        "regional_genus_prior":True},
                    "physical_block":copy.deepcopy(original_scores),
                    "first_source_group":copy.deepcopy(original_scores),
                })
    return {
        "status":"ORIGINAL_GIFT_STAGE2_TARGET_ACCEPTED_SPECIES_AND_GENUS_CLUSTER_SENSITIVITY",
        "original_parent_GIFT_botanical_floras":4505,
        "original_distinct_Direct_endpoints_total":13,
        "frozen_pilot_endpoints":sorted(ENDPOINTS),
        "unique_target_accepted_species_repetitions_not_independent":True,
        "species_equal_and_genus_equal_are_new_weighting_estimands":True,
        "source_island_based_original_proper_scores_unchanged":True,
        "regional_source_heldout_results":allrows,
    }


def test_original_gift_78_complete_and_both_species_genus_gates_supported():
    d=synthetic_original()
    result=audit(d)
    assert result["original_parent"]==4505
    assert result["original_panels"]==78
    assert len(result["all_anonymous_panels"])==78
    assert result["source_verdict_counts"][
        "all_recorded|tropical|species_equal|regional_blind|positive"]==13
    assert result["source_verdict_counts"][
        "original_native|southern_extratropical|genus_equal|external_blind|positive"]==13
    assert signed_pair(d["regional_source_heldout_results"][0],"species_equal","regional_blind")=="positive"


def test_repeated_target_specimen_count_has_to_match_unique_accepted_names():
    d=synthetic_original()
    d["regional_source_heldout_results"][0]["target_species_genus_reweighting"][
        "mean_island_occurrences_per_distinct_species"]=2
    with pytest.raises(ValueError,match="duplication ratio"):
        audit(d)
    d=synthetic_original()
    d["regional_source_heldout_results"][0]["target_species_genus_reweighting"][
        "n_species_island_occurrences"]=20
    with pytest.raises(ValueError,match="incidence"):
        audit(d)


def test_all_proper_species_vs_genus_reweighting_flags_have_to_be_correct():
    d=synthetic_original()
    r=d["regional_source_heldout_results"][0]
    a=r["target_species_genus_reweighting"]["species_equal"][
        "environment_source_nearest_gain_vs_regional_blind_Brier"]
    a["bootstrap_95"]=[-.01,.02]
    assert signed_pair(r,"species_equal","regional_blind")=="ambiguous"
    with pytest.raises(ValueError,match="positive gate"):
        audit(d)
    r["species_equal_positive_local_external_pooled"]["regional_blind"]=False
    assert audit(d)["original_panels"]==78


def test_prior_78_same_original_target_islands_and_group_scores_required():
    d=synthetic_original()
    old=copy.deepcopy(d)
    new=audit(d,old)
    assert new["original_island_equal_loss_forecasts_unchanged"]
    old["regional_source_heldout_results"][0]["n_original_target_islands_scored"]+=1
    with pytest.raises(ValueError,match="cohort denominator drift"):
        audit(d,old)
    old=copy.deepcopy(d)
    old["regional_source_heldout_results"][0]["physical_block"][
        "environment_source_nearest_gain_vs_regional_blind_Brier"]["group_equal_mean"]=.04
    with pytest.raises(ValueError,match="proper scores changed"):
        audit(d,old)


def test_missing_species_source_and_ontological_parent_fail_closed():
    d=synthetic_original()
    d["regional_source_heldout_results"].pop()
    with pytest.raises(ValueError,match="13×3×2"):
        audit(d)
    d=synthetic_original()
    d["original_parent_GIFT_botanical_floras"]=928
    with pytest.raises(ValueError,match="4505"):
        audit(d)
    d=synthetic_original()
    d["frozen_pilot_endpoints"][0]="made_up"
    with pytest.raises(ValueError,match="all13"):
        audit(d)


def test_privacy_deep_json_and_sha256_archive_guard(tmp_path):
    d=synthetic_original()
    f=tmp_path/"original.zip"
    with zipfile.ZipFile(f,"w") as z:
        z.writestr("source_anonymous.json",json.dumps(d))
    hashcode=hashlib.sha256(f.read_bytes()).hexdigest()
    assert read_original(f,sha256=hashcode)["original_parent_GIFT_botanical_floras"]==4505
    with pytest.raises(ValueError,match="SHA256"):
        read_original(f,sha256="0"*64)
    d["bad_nested"]={"accepted_species":"restricted_original_species"}
    bad=tmp_path/"bad.json"
    bad.write_text(json.dumps(d),encoding="utf-8")
    with pytest.raises(ValueError,match="Protected"):
        read_original(bad)



def test_unidentifiable_zero_direct_cohort_does_not_need_nonexistent_repeat_counts():
    # Original source emits a true UNSUPPORTED flora as an empty taxon
    # accumulator, with no occurrence count field at all.
    d=synthetic_original()
    row=d["regional_source_heldout_results"][0]
    row["n_original_target_islands_scored"]=0
    row["support_both_independent_source_groups"]=False
    taxon=row["target_species_genus_reweighting"]
    taxon.clear()
    taxon.update({
        "n_distinct_target_species":0,
        "n_distinct_target_genera":0,
        "target_accepted_names_absent_from_export":True,
        "species_equal":{"supported_for_inference":False},
        "genus_equal":{"supported_for_inference":False},
    })
    for u in ("species_equal","genus_equal"):
        row[f"{u}_positive_local_external_pooled"]={
            x:False for x in ("regional_blind","external_blind","regional_genus_prior")
        }
    a=audit(d)
    assert a["original_panels"]==78
    eligible=next(z for z in a["all_anonymous_panels"]
                  if (z["scope"],z["region"],z["trait_endpoint"]) ==
                  (row["scope"],row["region"],row["trait_endpoint"]))
    assert eligible["n_unique_direct_target_species"]==0
    assert eligible["n_direct_species_island_occurrences"]==0
    assert eligible["species_equal_gain_over_regional_blind"]=="unsupported"
    taxon["n_species_island_occurrences"]=7
    with pytest.raises(ValueError,match="invented occurrences"):
        audit(d)
