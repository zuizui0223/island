"""Synthetic tests for independent 78-panel GIFT proper-score vs Rao-Q gates."""
import copy
import hashlib
import json
import zipfile

import pytest

from scripts.verify_gift_category_vs_raoq_crossgate import (
    REGIONS, SCOPES, REFERENCES, GROUPS,
    read_source, evaluate, four_way_sign, rao_sign,
)


ENDPOINTS = [
    "flower_primary_color", "floral_form", "floral_symmetry",
    "flower_size_class", "tube_depth_class", "inflorescence_display",
    "self_incompatibility", "autonomous_selfing_capacity",
    "mating_system", "cleistogamy",
    "flower_primary_color__floral_form",
    "flower_primary_color__floral_symmetry",
    "floral_form__floral_symmetry",
]


def make_anonymous_sources():
    first,second=[],[]
    for region in REGIONS:
        for scope in SCOPES:
            for trait in ENDPOINTS:
                def group(with_q=False):
                    o={"supported_for_inference":True}
                    for baseline in REFERENCES:
                        for loss in ("Brier","logloss"):
                            o[f"environment_source_nearest_gain_vs_{baseline}_{loss}"]={
                                "bootstrap_95":[.01,.03],
                                "group_equal_mean":.02,
                            }
                    if with_q:
                        for estimator in ("plugin","unbiased"):
                            o[f"environment_source_nearest_gain_vs_regional_blind_categorical_simpson_{estimator}_abs_error"]={
                                "bootstrap_95":[.005,.025],
                                "group_equal_mean":.015,
                            }
                    return o
                common={
                    "region":region, "scope":scope,"trait_endpoint":trait,
                    "original_parent_botanical_floras":4505,
                    "original_full_accepted_name_group_ban_on_both_models":True,
                    "n_original_target_islands_scored":45,
                    "n_stage1_original_region_geography_eligible":100,
                    "n_original_region_flora":200,
                    "support_both_independent_source_groups":True,
                }
                first.append({**common,
                    "physical_block":group(), "first_source_group":group(),
                    "positive_vs_local_blind_both_scores_and_groups":True,
                    "positive_vs_external_blind_both_scores_and_groups":True,
                })
                second.append({**common,
                    "physical_block":group(with_q=True),
                    "first_source_group":group(with_q=True),
                    "categorical_simpson_plugin_abs_error_improves_both_groups":True,
                    "categorical_simpson_unbiased_abs_error_improves_both_groups":True,
                })
    base={"original_parent_GIFT_botanical_floras":4505,
          "frozen_pilot_endpoints":ENDPOINTS}
    return ({**base,"regional_source_heldout_results":first},
            {**base,
             "categorical_simpson_secondary_not_multivariate_hypervolume":True,
             "regional_source_heldout_results":second})


def test_complete_original_gift_two_source_artifacts_and_full_triple_gate():
    first,second=make_anonymous_sources()
    result=evaluate(first,second)
    assert result["n_original_source_panels"]==78
    assert len(result["panels"])==78
    assert result["summary_counts"]["all_recorded|joint_support|True"]==39
    assert result["summary_counts"]["original_native|triple_and_diversity|True"]==39
    assert result["panels"][0]["category_and_unbiased_Q_both_positive"]


def test_category_sign_and_unbiased_q_sign_can_disagree():
    first,second=make_anonymous_sources()
    c,r=first["regional_source_heldout_results"][0],second["regional_source_heldout_results"][0]
    for group in GROUPS:
        group_values=r[group]["environment_source_nearest_gain_vs_regional_blind_categorical_simpson_unbiased_abs_error"]
        group_values["bootstrap_95"]=[-.02,-.01]
    r["categorical_simpson_unbiased_abs_error_improves_both_groups"]=False
    assert rao_sign(r,"unbiased")=="negative"
    assert four_way_sign(c,"environment_source_nearest_gain_vs_regional_blind")=="positive"
    result=evaluate(first,second)
    row=next(x for x in result["panels"] if x["region"]==r["region"]
             and x["scope"]==r["scope"] and x["trait_endpoint"]==r["trait_endpoint"])
    assert row["category_vs_local"]=="positive"
    assert row["nominal_RaoQ_unbiased_vs_local"]=="negative"
    assert not row["category_and_unbiased_Q_both_positive"]


def test_cohort_drift_and_proper_score_drift_are_rejected():
    a,b=make_anonymous_sources()
    b["regional_source_heldout_results"][0]["n_original_target_islands_scored"]+=1
    with pytest.raises(ValueError,match="cohort"):
        evaluate(a,b)
    a,b=make_anonymous_sources()
    b["regional_source_heldout_results"][0]["physical_block"][
        "environment_source_nearest_gain_vs_regional_blind_Brier"]["group_equal_mean"]=.04
    with pytest.raises(ValueError,match="proper-score"):
        evaluate(a,b)


def test_13_traits_source_parent_and_duplicate_panel_guard():
    a,b=make_anonymous_sources()
    a["original_parent_GIFT_botanical_floras"]=928
    with pytest.raises(ValueError,match="4505"):
        evaluate(a,b)
    a,b=make_anonymous_sources()
    a["regional_source_heldout_results"][0]["trait_endpoint"]=ENDPOINTS[1]
    with pytest.raises(ValueError,match="Duplicate"):
        evaluate(a,b)
    a,b=make_anonymous_sources()
    b["frozen_pilot_endpoints"].pop()
    with pytest.raises(ValueError,match="13"):
        evaluate(a,b)


def test_source_zip_digest_and_private_original_identity_check(tmp_path):
    a,b=make_anonymous_sources()
    p=tmp_path/"anonymous_original.zip"
    with zipfile.ZipFile(p,"w") as z:
        z.writestr("all13.json",json.dumps(a))
    digest=hashlib.sha256(p.read_bytes()).hexdigest()
    assert read_source(p,digest)["original_parent_GIFT_botanical_floras"]==4505
    with pytest.raises(ValueError,match="digest"):
        read_source(p,"0"*64)
    with zipfile.ZipFile(p,"w") as z:
        z.writestr("all13.json",json.dumps(a))
        z.writestr("copy.json",json.dumps(a))
    with pytest.raises(ValueError,match="exactly one"):
        read_source(p)
    a["metadata"]={"accepted_species":"Never export source protected trait name"}
    raw = tmp_path/"prohibited_raw.json"
    raw.write_text(json.dumps(a))
    with pytest.raises(ValueError,match="source identity"):
        read_source(raw)
