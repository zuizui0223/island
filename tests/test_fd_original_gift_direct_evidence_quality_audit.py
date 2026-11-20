"""Synthetic-only original GIFT species-Direct high/medium evidence sensitivity."""
import numpy as np
import pandas as pd
import pytest

from fd_original_gift_direct_evidence_quality_audit import (
    TRAITS, PAIRS, ENDPOINTS, direct_grade_states,
    variation_panel, bh_qvalues, REGIONS, SCOPES
)


def synthetic_ledger(n_genera=15):
    rows=[]
    # This constructs 15 genera with distinct accepted plant species,
    # each with 5 high and 5 medium species. We do NOT invent two real
    # plant trait observations for the same accepted species.
    for g in range(n_genera):
        for k in range(10):
            name=f"Genus{g} species{k}"
            for trait in TRAITS:
                rows.append({
                    "accepted_species":name,"trait_name":trait,
                    "evidence_scope":"species_direct",
                    "evidence_quality":"high" if k<5 else "medium",
                    "normalized_value":("classA" if k<5 else "classB"),
                })
    return pd.DataFrame(rows)


def test_true_species_direct_quality_and_genuine_same_species_pairs_only():
    data=direct_grade_states(synthetic_ledger(15))
    assert len(data)==len(ENDPOINTS)==13
    assert set(data)==set(TRAITS)|{a+"__"+b for a,b in PAIRS}
    key="flower_primary_color__floral_symmetry"
    assert len(data[key]["state"])==150
    assert data[key]["state"]["Genus0 species0"]==("classA","classA")
    assert data[key]["quality"]["Genus0 species0"]=="high"
    assert data[key]["quality"]["Genus0 species8"]=="medium"
    ledger=synthetic_ledger(15)
    ledger.loc[
        (ledger.accepted_species=="Genus0 species0")
        & (ledger.trait_name=="floral_symmetry"),
        "evidence_quality"
    ]="medium"
    changed=direct_grade_states(ledger)
    assert changed[key]["quality"]["Genus0 species0"]=="medium"
    assert changed[key]["state"]["Genus0 species0"]==("classA","classA")


def test_within_genus_category_quality_differences_permute_grade_labels():
    data=direct_grade_states(synthetic_ledger())
    source=set(data["floral_form"]["state"])
    r=variation_panel(
        source,data["floral_form"],universe_count=len(source),
        seed=20261121,nperm=99)
    assert r["n_source_distinct_genera_with_Direct"]==15
    assert r["n_within_genus_quality_matched_genera"]==15
    assert r["n_true_species_Direct_coverage"]==150
    assert r["n_direct_quality_high"]==75
    assert r["n_direct_quality_medium"]==75
    assert r["pooled_high_medium_category_total_variation"]==pytest.approx(1)
    assert r["within_genus_equal_mean_high_medium_total_variation"]==pytest.approx(1)
    assert r["within_genus_permutation_qualified"]
    assert r["within_genus_quality_label_permutation_upper_tail_p_exploratory"]<=.02
    assert r["within_genus_observed_minus_permutation_mean_TV"]>0
    assert r["biological_within_genus_variation_vs_annotation_not_identifiable"]


def test_identical_high_medium_distributions_have_no_quality_gradient():
    ledger=synthetic_ledger()
    ledger["normalized_value"]="classA"
    data=direct_grade_states(ledger)
    r=variation_panel(set(data["floral_form"]["state"]),data["floral_form"],
                      universe_count=150,seed=17,nperm=99)
    assert r["within_genus_permutation_qualified"]
    assert r["within_genus_equal_mean_high_medium_total_variation"]==pytest.approx(0)
    assert r["within_genus_quality_label_permutation_upper_tail_p_exploratory"]==1


def test_a_few_within_genus_matched_groups_must_report_unsupported():
    data=direct_grade_states(synthetic_ledger(4))
    r=variation_panel(set(data["floral_form"]["state"]),data["floral_form"],
                      universe_count=40,seed=19,nperm=99)
    assert r["n_within_genus_quality_matched_genera"]==4
    assert not r["within_genus_permutation_qualified"]
    assert r["within_genus_quality_label_permutation_upper_tail_p_exploratory"] is None


def test_unmatched_and_uncertified_botanical_name_not_made_fake_genus():
    data=direct_grade_states(synthetic_ledger(12))
    ep=data["floral_form"]
    ep["state"]["Sp."]="classA"
    ep["quality"]["Sp."]="high"
    source=set(ep["state"])
    r=variation_panel(source,ep,universe_count=len(source),seed=11,nperm=99)
    assert r["n_direct_names_without_certified_genus"]==1
    assert r["n_within_genus_quality_matched_genera"]==12
    assert r["n_true_species_Direct_coverage"]==len(source)


def test_source_original_direct_quality_duplicate_name_trait_fails():
    ledger=synthetic_ledger(2)
    ledger=pd.concat([ledger,ledger.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="duplicated"):
        direct_grade_states(ledger)


def test_quality_graded_unresolved_and_genus_incomplete_not_invented():
    data=direct_grade_states(synthetic_ledger(12))
    ep=data["floral_form"]
    assert len(ep["state"])==120
    empty=variation_panel(set(),ep,universe_count=0,seed=12,nperm=99)
    assert empty["n_true_species_Direct_coverage"]==0
    assert not empty["within_genus_permutation_qualified"]
    assert empty["pooled_high_medium_category_total_variation"] is None
    with pytest.raises(ValueError,match="cohort"):
        variation_panel(set(ep["state"]),ep,universe_count=1,seed=1,nperm=99)


def test_bh_quality_p_values_report_nonindependent_nested_family():
    x=[
        {"within_genus_permutation_qualified":True,
         "within_genus_quality_label_permutation_upper_tail_p_exploratory":.01},
        {"within_genus_permutation_qualified":True,
         "within_genus_quality_label_permutation_upper_tail_p_exploratory":.03},
        {"within_genus_permutation_qualified":True,
         "within_genus_quality_label_permutation_upper_tail_p_exploratory":.20},
        {"within_genus_permutation_qualified":False}
    ]
    bh_qvalues(x)
    assert [round(v["within_genus_quality_label_permutation_BH_q_exploratory"],3)
            for v in x[:3]]==[.04,.06,.20]
    assert "within_genus_quality_label_permutation_BH_q_exploratory" not in x[3]
