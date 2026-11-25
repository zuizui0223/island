"""Synthetic original 4505-GIFT equal-source-N audit contracts."""
import copy
import hashlib
import json
import zipfile

import pytest

from scripts.verify_gift_equalN_source_reference import (
    DIRECT_ENDPOINTS, REGIONS, SCOPES, GROUPS, SEEDS, PARENT,
    audit, read_source, score_gate,
)


def synthetic_original():
    rows=[]
    for scope in SCOPES:
        for region in REGIONS:
            for trait in sorted(DIRECT_ENDPOINTS):
                group={"supported_for_inference":True}
                for seed in SEEDS:
                    for loss in ("Brier","logloss"):
                        forward=f"balanced_external_seed{seed}_gain_vs_balanced_local_seed{seed}_{loss}"
                        reverse=f"balanced_local_seed{seed}_gain_vs_balanced_external_seed{seed}_{loss}"
                        group[forward]={"bootstrap_95":[.01,.03],"group_equal_mean":.02}
                        group[reverse]={"bootstrap_95":[-.03,-.01],"group_equal_mean":-.02}
                rows.append({
                    "region":region,"scope":scope,"trait_endpoint":trait,
                    "original_parent_botanical_floras":PARENT,
                    "n_original_target_islands_scored":45,
                    "original_full_accepted_name_group_ban_on_both_models":True,
                    "original_focal_genus_oracle_NOT_used_in_predictor":True,
                    "support_both_independent_source_groups":True,
                    "balanced_direct_source_names_per_side_mean":200.,
                    "external_local_original_source_n_ratio_mean":1.2,
                    "balanced_external_better_all_three_source_samples_both_losses_groups":True,
                    "balanced_local_better_all_three_source_samples_both_losses_groups":False,
                    "physical_block":copy.deepcopy(group),
                    "first_source_group":copy.deepcopy(group),
                })
    return {
        "status":"ORIGINAL_GIFT_STAGE2_BALANCED_SOURCE_LOCAL_EXTERNAL_BLIND_SENSITIVITY",
        "original_parent_GIFT_botanical_floras":PARENT,
        "original_distinct_Direct_endpoints_total":13,
        "matched_source_direct_name_sample_seeds":[0,1,2],
        "matched_source_direct_name_per_side_cap":2000,
        "original_focal_trait_cohorts_and_nonoracle_predictions_unchanged":True,
        "direct_local_vs_external_blind_losses_bootstrapped_on_identical_source_carriers":True,
        "frozen_pilot_endpoints":sorted(DIRECT_ENDPOINTS),
        "regional_source_heldout_results":rows,
    }


def test_all_original_78_panels_three_hash_seed_double_loss_double_group():
    r=audit(synthetic_original())
    assert r["original_parent_GIFT_islands"]==4505
    assert r["n_original_region_status_panels"]==78
    assert r["n_source_hash_samples"]==3
    assert len(r["panels"])==78
    assert r["source_verdict_counts"]["all_recorded|tropical|external_advantage|positive"]==13
    assert r["source_verdict_counts"]["original_native|southern_extratropical|local_advantage|negative"]==13


def test_one_counter_seed_in_one_group_prevents_overclaim():
    x=synthetic_original()
    row=x["regional_source_heldout_results"][0]
    key="balanced_external_seed2_gain_vs_balanced_local_seed2_logloss"
    row["first_source_group"][key]["bootstrap_95"]=[-.001,.04]
    row["balanced_external_better_all_three_source_samples_both_losses_groups"]=False
    row["balanced_local_better_all_three_source_samples_both_losses_groups"]=False
    assert score_gate(row,beneficiary="external")=="ambiguous"
    r=audit(x)
    assert r["panels"][0]["n_original_scored_islands"]==45
    row["balanced_external_better_all_three_source_samples_both_losses_groups"]=True
    with pytest.raises(ValueError,match="Stored equal-N"):
        audit(x)


def test_original_unsupported_cells_keep_denominator_and_no_sign():
    x=synthetic_original()
    row=x["regional_source_heldout_results"][0]
    row["support_both_independent_source_groups"]=False
    row["balanced_external_better_all_three_source_samples_both_losses_groups"]=False
    row["balanced_local_better_all_three_source_samples_both_losses_groups"]=False
    row.pop("balanced_direct_source_names_per_side_mean")
    row.pop("external_local_original_source_n_ratio_mean")
    row["physical_block"]={}
    row["first_source_group"]={}
    assert score_gate(row,beneficiary="local")=="unsupported"
    r=audit(x)
    assert r["n_original_region_status_panels"]==78


def test_cannot_replace_original_direct_13_endpoint_family_or_botanical_parent():
    x=synthetic_original()
    x["frozen_pilot_endpoints"][0]="fake_trait"
    with pytest.raises(ValueError,match="13×3×2"):
        audit(x)
    x=synthetic_original()
    x["original_parent_GIFT_botanical_floras"]=928
    with pytest.raises(ValueError,match="identity"):
        audit(x)
    x=synthetic_original()
    x["regional_source_heldout_results"][1]["trait_endpoint"]=x["regional_source_heldout_results"][0]["trait_endpoint"]
    with pytest.raises(ValueError,match="duplicated"):
        audit(x)


def test_equalN_anonymous_source_zip_digest_and_no_botanical_id_leakage(tmp_path):
    x=synthetic_original()
    p=tmp_path/"source_original.zip"
    with zipfile.ZipFile(p,"w") as z:
        z.writestr("anonymous_original.json",json.dumps(x))
    h=hashlib.sha256(p.read_bytes()).hexdigest()
    assert read_source(p,h)["original_parent_GIFT_botanical_floras"]==4505
    with pytest.raises(ValueError,match="digest"):
        read_source(p,"0"*64)
    x["metadata"]={"accepted_species":"forbidden"}
    raw=tmp_path/"original.json"
    raw.write_text(json.dumps(x))
    with pytest.raises(ValueError,match="identifying"):
        read_source(raw)


def test_malformed_source_geography_intervals_fail_closed():
    x=synthetic_original()
    row=x["regional_source_heldout_results"][0]
    key="balanced_external_seed0_gain_vs_balanced_local_seed0_Brier"
    row["physical_block"][key]["bootstrap_95"]=[0,float("nan")]
    with pytest.raises(ValueError,match="malformed"):
        audit(x)
