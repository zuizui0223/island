"""Synthetic tests for independent anonymous original GIFT stage-2 gate."""
import json
import zipfile

import pytest

from scripts.verify_stage2_source_anonymous import load_anon, verify, direction

REGIONS = ("northern_midlatitude", "tropical", "southern_extratropical")
SCOPES = ("all_recorded", "original_native")
EPS = ("floral_form", "flower_primary_color__floral_form")


def panel(region, scope, endpoint, *, positive=True):
    def entries():
        return {
            f"{m}_gain_vs_{base}_{loss}": {
                "bootstrap_95": ([0.01, 0.03] if positive else [-0.03, -0.01]),
                "group_equal_mean": (0.02 if positive else -0.02)
            }
            for m in ("environment_source_nearest",)
            for base in ("regional_blind", "external_blind", "regional_genus_prior")
            for loss in ("Brier", "logloss")
        }
    group = lambda: {"supported_for_inference": True, "n_groups": 20, **entries()}
    return {
        "region": region, "scope": scope, "trait_endpoint": endpoint,
        "original_parent_botanical_floras": 4505,
        "original_full_accepted_name_group_ban_on_both_models": True,
        "n_original_target_islands_scored": 100,
        "positive_vs_local_blind_both_scores_and_groups": positive,
        "positive_vs_external_blind_both_scores_and_groups": positive,
        "physical_block": group(), "first_source_group": group()
    }


def complete():
    return {"original_parent_GIFT_botanical_floras": 4505,
            "status": "SYNTHETIC_ONLY",
            "frozen_pilot_endpoints": list(EPS),
            "regional_source_heldout_results": [
                panel(reg, scope, endpoint)
                for reg in REGIONS for scope in SCOPES for endpoint in EPS]}


def test_independent_dual_proper_loss_dual_cluster_signs():
    data = complete()
    result = verify(data)
    assert result["n_panels"] == 12
    assert result["three_region_supported_endpoints"]["all_recorded:local_blind"] == sorted(EPS)
    assert result["status_counts"]["all_recorded|local_blind|positive"] == 6
    assert result["n_endpoints"] == 2
    bad = data["regional_source_heldout_results"][0]
    bad["physical_block"]["environment_source_nearest_gain_vs_regional_blind_Brier"]["bootstrap_95"] = [-0.01, 0.02]
    bad["positive_vs_local_blind_both_scores_and_groups"] = False
    result = verify(data)
    assert result["status_counts"]["all_recorded|local_blind|ambiguous"] == 1
    bad["positive_vs_local_blind_both_scores_and_groups"] = True
    with pytest.raises(ValueError, match="Stored source positive gate"):
        verify(data)


def test_negative_and_insufficient_signs_remain_visible():
    p = panel("tropical", "original_native", "floral_form", positive=False)
    assert direction(p, model="environment_source_nearest", baseline="regional_blind") == "negative"
    p["first_source_group"]["supported_for_inference"] = False
    assert direction(p, model="environment_source_nearest", baseline="regional_blind") == "insufficient"
    p["first_source_group"]["supported_for_inference"] = True
    del p["physical_block"]["environment_source_nearest_gain_vs_regional_blind_Brier"]
    with pytest.raises(ValueError, match="Missing or invalid"):
        direction(p, model="environment_source_nearest", baseline="regional_blind")


def test_denominator_and_endpoint_panels_cannot_disappear():
    x = complete()
    x["regional_source_heldout_results"].pop()
    with pytest.raises(ValueError, match="Missing original"):
        verify(x)
    x = complete()
    x["original_parent_GIFT_botanical_floras"] = 928
    with pytest.raises(ValueError, match="4505"):
        verify(x)
    x = complete()
    x["regional_source_heldout_results"][2]["region"] = "tropical"
    with pytest.raises(ValueError, match="Repeated"):
        verify(x)


def test_reads_exact_one_public_anonymous_source_zip_without_sensitive_ids(tmp_path):
    original = complete()
    z = tmp_path/"original-anonymous.zip"
    with zipfile.ZipFile(z, "w") as archive:
        archive.writestr("anonymous.json", json.dumps(original))
    assert load_anon(z)["original_parent_GIFT_botanical_floras"] == 4505
    with zipfile.ZipFile(z, "w") as archive:
        archive.writestr("anonymous.json", json.dumps(original))
        archive.writestr("duplicate.json", json.dumps(original))
    with pytest.raises(ValueError, match="exactly one"):
        load_anon(z)
    plain = tmp_path/"raw.json"
    original["private_leakage"] = {"island_id": "hidden"}
    plain.write_text(json.dumps(original))
    with pytest.raises(ValueError, match="Private original"):
        load_anon(plain)
