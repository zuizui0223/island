"""Synthetic strictly documented floristic native/endemic/introduced and source contrasts."""
import pandas as pd
import pytest
from fd_direct_three_stage_floristic_comparison import (
    status_groups,source_sectors,original_shared_genus,grouped_interval
)


def _statuses():
    return pd.DataFrame([
       dict(island_id="A",accepted_species="Gen1 e",
            origin_status="native",endemic_status="endemic",floristic_status="endemic"),
       dict(island_id="A",accepted_species="Gen1 n",
            origin_status="native",endemic_status="nonendemic",floristic_status="native_nonendemic"),
       dict(island_id="A",accepted_species="Gen2 u",
            origin_status="native",endemic_status="unresolved",floristic_status="native_endemism_unresolved"),
       dict(island_id="A",accepted_species="Gen3 i",
            origin_status="introduced",endemic_status="nonendemic",floristic_status="introduced"),
       dict(island_id="A",accepted_species="Gen4 z",
            origin_status="unresolved",endemic_status="unresolved",floristic_status="unresolved"),
       dict(island_id="B",accepted_species="Gen1 e",
            origin_status="native",endemic_status="endemic",floristic_status="endemic")
    ])


def test_endemic_is_independent_from_native_and_intro_with_true_statuses():
    q=status_groups(_statuses())
    assert len(q)==2
    assert q["A"]["endemic"]=={"Gen1 e"}
    assert q["A"]["native_nonendemic"]=={"Gen1 n"}
    assert q["A"]["introduced"]=={"Gen3 i"}
    assert q["A"]["native_endemism_unresolved"]=={"Gen2 u"}
    assert q["A"]["unresolved_origin"]=={"Gen4 z"}
    assert len(q["A"]["full"])==5
    assert q["B"]["endemic"]=={"Gen1 e"}


def test_endemicity_contradictions_and_duplicate_accepted_rows_fail_closed():
    g=_statuses()
    with pytest.raises(ValueError):
        status_groups(pd.concat([g,g.iloc[:1]],ignore_index=True))
    bad=g.copy()
    bad.loc[bad.accepted_species.eq("Gen1 e"),"origin_status"]="introduced"
    with pytest.raises(ValueError):status_groups(bad)
    bad=g.copy()
    bad.loc[bad.accepted_species.eq("Gen1 n"),"endemic_status"]="endemic"
    with pytest.raises(ValueError):status_groups(bad)


def test_each_real_flower_trait_stays_separate_with_strict_binary_states():
    state={
      "floral_symmetry":{"G e":"zygomorphic","G n":"actinomorphic"},
      "floral_form":{"G e":"tubular","G n":"open_radial"},
      "tube_depth_class":{"G e":"deep","G n":"shallow"},
      "self_incompatibility":{"G e":"SC","G n":"SI","G x":"mixed_or_variable"},
      "autonomous_selfing_capacity":{
         "G e":"absent","G n":"autonomous","G x":"delayed"}
    }
    s=source_sectors(state)
    assert len(s)==7
    assert s["zygomorphic"]=={"G e":1,"G n":0}
    assert s["tube_or_bell"]=={"G e":1,"G n":0}
    assert s["open_radial"]=={"G e":0,"G n":1}
    assert s["deep_tube"]=={"G e":1,"G n":0}
    assert s["SC_vs_SI"]=={"G e":1,"G n":0}
    assert s["autonomous_vs_absent"]=={"G e":0,"G n":1}
    assert s["bilateral_tube_or_bell_joint"]=={"G e":1,"G n":0}
    assert "G x" not in s["SC_vs_SI"]
    assert "G x" not in s["autonomous_vs_absent"]


def test_same_genus_comparison_counts_actually_variable_direct_categories():
    E={"GenA endemic1","GenB endemic2","GenC endemic3"}
    N={"GenA native1","GenB native2","GenD native4"}
    traits={
      "GenA endemic1":1,"GenA native1":0,
      "GenB endemic2":1,"GenB native2":1,
      "GenC endemic3":1,"GenD native4":0}
    out=original_shared_genus(E,N,traits)
    assert out["n_matched_same_island_genus_units"]==2
    assert out["n_matched_same_island_genus_with_actual_category_variation"]==1
    assert out["mean_genus_equal_endemic_minus_native_nonendemic_share"]==pytest.approx(.5)
    q=original_shared_genus(
        {"GenB endemic2"},{"GenB native2"},traits)
    assert q["n_matched_same_island_genus_with_actual_category_variation"]==0
    assert q["mean_genus_equal_endemic_minus_native_nonendemic_share"]==0
    r=original_shared_genus({"GenC endemic3"},{"GenD native4"},traits)
    assert r["n_matched_same_island_genus_units"]==0
    assert r["mean_genus_equal_endemic_minus_native_nonendemic_share"] is None


def test_not_well_supported_does_not_invent_zero_effect():
    q=grouped_interval(
        [dict(block="B1",difference=.2),dict(block="B1",difference=.3)],
        "block",boots=99)
    assert q["status"]=="NOT_IDENTIFIABLE_ORIGINAL_SOURCE_CLUSTER_SUPPORT"
    assert q["original_island_equal_difference"]==pytest.approx(.25)
    assert grouped_interval([],"block",boots=99)["status"].startswith("NOT_IDENTIFIABLE")


def test_block_equal_direct_difference_not_plant_duplicated_replicates():
    vals=[
        dict(block="B"+str(i//2),difference=.11+(-1)**i*.02)
        for i in range(80)]
    q=grouped_interval(vals,"block",boots=99)
    assert q["status"]=="DESCRIPTIVE_ORIGINAL_BLOCK_EQUAL_FLORISTIC_SECTOR_CONTRAST"
    assert q["n_comparable_floras"]==80
    assert q["n_source_clusters"]==40
    assert q["original_cluster_equal_difference"]==pytest.approx(.11)
    assert q["original_cluster_conditional_resample_95"]==pytest.approx([.11,.11])
