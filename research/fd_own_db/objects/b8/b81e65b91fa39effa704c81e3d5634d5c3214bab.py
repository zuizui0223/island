"""Synthetic contracts for independently auditing original GIFT joint-source panels."""
import json
import zipfile

import pytest

from scripts.verify_gift_architecture_assurance_source import (
    ARCH, ASSURANCE, REGIONS, SCOPES,
    read_anonymous, audit_original_joint, sign,
)


def fake_panels(*, conditional=False):
    rows = []
    for region in sorted(REGIONS):
        for scope in sorted(SCOPES):
            for a in ARCH:
                for b in ASSURANCE:
                    frame = {
                        "supported_for_inference": True,
                        "predicted_joint_gain_vs_own_independence_Brier":
                            {"bootstrap_95": [.01, .03]},
                        "predicted_joint_gain_vs_own_independence_logloss":
                            {"bootstrap_95": [.03, .07]},
                        "regional_blind_joint_gain_vs_own_independence_Brier":
                            {"bootstrap_95": [-.05, -.01]},
                        "regional_blind_joint_gain_vs_own_independence_logloss":
                            {"bootstrap_95": [-.08, -.02]},
                    }
                    if conditional:
                        frame.update({
                            "predicted_joint_gain_vs_genus_conditional_independence_Brier":
                                {"bootstrap_95": [-.01, .03]},
                            "predicted_joint_gain_vs_genus_conditional_independence_logloss":
                                {"bootstrap_95": [.01, .04]},
                        })
                    row = {
                        "region": region, "scope": scope, "trait_endpoint": a+"__"+b,
                        "original_parent_botanical_floras":4505,
                        "original_full_accepted_name_group_ban_on_both_models":True,
                        "no_pollinator_guild_or_within_species_island_adaptation_measured":True,
                        "support_both_independent_source_groups":True,
                        "n_original_target_islands_scored":64,
                        "external_genus_joint_dependency_over_own_marginals_supported":True,
                        "regional_blind_joint_dependency_over_own_marginals_supported":False,
                        "physical_block":dict(frame),
                        "first_source_group":dict(frame),
                    }
                    if conditional:
                        row["external_genus_joint_dependency_over_genus_conditional_marginals_supported"]=False
                    rows.append(row)
    data={
        "original_parent_GIFT_botanical_floras":4505,
        "frozen_pilot_endpoints":[a+"__"+b for a in ARCH for b in ASSURANCE],
        "regional_source_heldout_results":rows,
        "status":"SYNTHETIC_ORIGINAL_JOINT_GATE_TEST_ONLY",
    }
    if conditional:
        data["genus_conditional_independence_retains_between_genus_compositional_covariance"]=True
    return data


def test_all_36_true_joint_panels_preserved_and_all_gates_verified():
    r=audit_original_joint(fake_panels())
    assert r["n_source_panels"]==36
    assert len(r["panels"])==36
    assert r["source_verdict_counts"]["all_recorded|tropical|joint_vs_global_independence|positive"]==6
    assert r["source_verdict_counts"]["original_native|northern_midlatitude|local_prior_joint_vs_independence|negative"]==6


def test_original_genus_conditional_null_does_not_inherit_global_independence_sign():
    r=audit_original_joint(fake_panels(conditional=True))
    assert r["stronger_genus_conditional_null_included"]
    assert r["source_verdict_counts"]["all_recorded|tropical|joint_vs_within_genus_independence|ambiguous"]==6
    assert r["source_verdict_counts"]["all_recorded|tropical|joint_vs_global_independence|positive"]==6


def test_failed_original_group_or_missing_interval_cannot_claim_positive():
    d=fake_panels()
    first=d["regional_source_heldout_results"][0]
    first["support_both_independent_source_groups"]=False
    first["external_genus_joint_dependency_over_own_marginals_supported"]=False
    first["regional_blind_joint_dependency_over_own_marginals_supported"]=False
    assert sign(first,"predicted_joint_gain_vs_own_independence")=="unsupported"
    assert audit_original_joint(d)["n_source_panels"]==36
    first["support_both_independent_source_groups"]=True
    first["physical_block"]["predicted_joint_gain_vs_own_independence_Brier"]["bootstrap_95"]=[-.01,.02]
    assert sign(first,"predicted_joint_gain_vs_own_independence")=="ambiguous"
    # A deliberately contradictory stored-positive flag must fail the audit.
    first["external_genus_joint_dependency_over_own_marginals_supported"]=True
    with pytest.raises(ValueError,match="stored positive flag"):
        audit_original_joint(d)
    first["external_genus_joint_dependency_over_own_marginals_supported"]=False
    del first["physical_block"]["predicted_joint_gain_vs_own_independence_Brier"]
    with pytest.raises(ValueError,match="interval missing"):
        audit_original_joint(d)


def test_original_botanical_parent_and_identity_privacy_must_be_explicit(tmp_path):
    d=fake_panels()
    d["regional_source_heldout_results"].pop()
    with pytest.raises(ValueError,match="Missing original"):
        audit_original_joint(d)
    d=fake_panels()
    d["original_parent_GIFT_botanical_floras"]=928
    with pytest.raises(ValueError,match="4505"):
        audit_original_joint(d)
    d=fake_panels()
    a=tmp_path/"anonymous.zip"
    with zipfile.ZipFile(a,"w") as z:z.writestr("source.json",json.dumps(d))
    assert read_anonymous(a)["original_parent_GIFT_botanical_floras"]==4505
    d["meta"]={"accepted_species":"Protected name should not appear"}
    with zipfile.ZipFile(a,"w") as z:z.writestr("source.json",json.dumps(d))
    with pytest.raises(ValueError,match="identity key"):
        read_anonymous(a)
