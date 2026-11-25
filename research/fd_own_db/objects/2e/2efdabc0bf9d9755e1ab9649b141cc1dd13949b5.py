"""Synthetic contracts for GIFT original quality vs floral categories source auditor."""
import copy
import hashlib
import json
import zipfile

import pytest

from scripts.verify_gift_original_direct_quality_source import (
    ENDPOINTS,REGIONS,SCOPES,SOURCE_STATUS,read_anon,audit
)


def source_record():
    panels=[]
    for region in REGIONS:
        for scope in SCOPES:
            for ep in sorted(ENDPOINTS):
                panels.append({
                    "region":region,"scope":scope,"trait_endpoint":ep,
                    "original_GIFT_botanical_parent_island_floras":4505,
                    "regional_original_parent_botanical_floras":200,
                    "source_region_original_distinct_botanical_species":100,
                    "n_original_distinct_accepted_species":100,
                    "n_true_species_Direct_coverage":80,
                    "original_distinct_species_Direct_coverage_fraction":.8,
                    "n_direct_quality_high":40,
                    "n_direct_quality_medium":40,
                    "high_grade_fraction_among_Direct":.5,
                    "n_source_distinct_genera_with_Direct":20,
                    "n_direct_names_without_certified_genus":0,
                    "n_within_genus_quality_matched_genera":12,
                    "n_different_category_states":3,
                    "pooled_high_medium_category_total_variation":.6,
                    "within_genus_equal_mean_high_medium_total_variation":.4,
                    "within_genus_permutation_qualified":True,
                    "within_genus_null_mean_total_variation":.2,
                    "within_genus_null_95th_total_variation":.3,
                    "within_genus_observed_minus_permutation_mean_TV":.2,
                    "within_genus_quality_label_permutation_upper_tail_p_exploratory":.01,
                    "within_genus_quality_label_permutation_BH_q_exploratory":.01,
                    "high_vs_medium_are_different_plant_species_not_repeated_measurements":True,
                    "biological_within_genus_variation_vs_annotation_not_identifiable":True
                })
    return {
        "status":SOURCE_STATUS,
        "original_parent_GIFT_botanical_island_floras":4505,
        "original_parent_geography_islands":8264,
        "n_genuine_Direct_endpoints":13,
        "n_region_scope_endpoint_panels":78,
        "quality_ontology":"species_direct_high_vs_medium_once_per_species",
        "joint_endpoint_high_grade_requires_both_same_species_traits_high":True,
        "original_species_name_and_genus_IDs_never_exported":True,
        "no_independent_within_species_or_local_island_phenotype_replicates":True,
        "BH_q_exploratory_across_all_qualified_78_nested_tests":True,
        "n_permutation_qualified_panels":78,
        "panels":panels
    }


def test_complete_78_panels_nested_noncausal_quality_test_and_bh_all():
    x=source_record()
    r=audit(x)
    assert r["GIFT_botanical_parent"]==4505
    assert r["geography_parent"]==8264
    assert r["n_full_region_status_panels"]==78
    assert r["n_qualified_within_genus_quality_tests"]==78
    assert len(r["anonymous_panels"])==78
    assert r["scope_region_source_support_counts"][
        "all_recorded|tropical|exploratory_BH_q_below_0_05"]==13
    assert "not repeated measures" in r["nonclaims"] or "DIFFERENT species" in r["nonclaims"]


def test_bh_q_recomputed_over_all_qualified_trait_family():
    x=source_record()
    row=x["panels"][0]
    row["within_genus_quality_label_permutation_BH_q_exploratory"]=.0001
    with pytest.raises(ValueError,match="BH q mismatch"):
        audit(x)
    x=source_record()
    row=x["panels"][0]
    row["within_genus_permutation_qualified"]=False
    for k in ("within_genus_null_mean_total_variation",
              "within_genus_null_95th_total_variation",
              "within_genus_observed_minus_permutation_mean_TV",
              "within_genus_quality_label_permutation_upper_tail_p_exploratory",
              "within_genus_quality_label_permutation_BH_q_exploratory"):
        row[k]=None
    x["n_permutation_qualified_panels"]=77
    assert audit(x)["n_qualified_within_genus_quality_tests"]==77


def test_direct_quality_counts_and_true_species_coverages_not_invented():
    x=source_record()
    x["panels"][0]["n_direct_quality_medium"]=39
    with pytest.raises(ValueError,match="denominators"):
        audit(x)
    x=source_record()
    x["panels"][0]["n_direct_quality_high"]=60
    x["panels"][0]["n_direct_quality_medium"]=20
    with pytest.raises(ValueError,match="proportions"):
        audit(x)
    x=source_record()
    x["panels"][0]["n_source_distinct_genera_with_Direct"]=4
    with pytest.raises(ValueError,match="denominators"):
        audit(x)


def test_missing_source_panels_and_duplicate_study_keys_rejected():
    x=source_record()
    x["panels"].pop()
    with pytest.raises(ValueError,match="13×3×2"):
        audit(x)
    x=source_record()
    x["panels"][1]["trait_endpoint"]=x["panels"][0]["trait_endpoint"]
    with pytest.raises(ValueError,match="Duplicate"):
        audit(x)
    x=source_record()
    x["original_parent_GIFT_botanical_island_floras"]=928
    with pytest.raises(ValueError,match="provenance"):
        audit(x)


def test_hash_pinned_original_source_archive_and_private_name_keys_rejected(tmp_path):
    x=source_record()
    f=tmp_path/"source.zip"
    with zipfile.ZipFile(f,"w") as archive:
        archive.writestr("original_anonymous.json",json.dumps(x))
    digest=hashlib.sha256(f.read_bytes()).hexdigest()
    assert read_anon(f,digest)["original_parent_GIFT_botanical_island_floras"]==4505
    with pytest.raises(ValueError,match="SHA256"):
        read_anon(f,"0"*64)
    x["accidental_original"]={"accepted_species":"Protected original name"}
    raw=tmp_path/"accidental.json"
    raw.write_text(json.dumps(x))
    with pytest.raises(ValueError,match="Protected"):
        read_anon(raw)


def test_nonqualified_genus_comparison_not_conflated_with_biological_absence():
    x=source_record()
    r=x["panels"][0]
    r.update({
        "n_source_distinct_genera_with_Direct":1,
        "n_within_genus_quality_matched_genera":0,
        "within_genus_equal_mean_high_medium_total_variation":None,
        "within_genus_permutation_qualified":False,
        "within_genus_null_mean_total_variation":None,
        "within_genus_null_95th_total_variation":None,
        "within_genus_observed_minus_permutation_mean_TV":None,
        "within_genus_quality_label_permutation_upper_tail_p_exploratory":None,
        "within_genus_quality_label_permutation_BH_q_exploratory":None,
    })
    x["n_permutation_qualified_panels"]=77
    z=audit(x)
    assert z["n_qualified_within_genus_quality_tests"]==77
    assert z["anonymous_panels"][0]["within_genus_permutation_qualified"] is False
