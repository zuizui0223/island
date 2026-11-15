"""Pure synthetic exact DIRECT genus decomposition of actual native/source trait differences.

No biological GIFT botanical IDs; no target phenotype from predicted genus.
"""
import numpy as np
import pandas as pd
import pytest
from fd_direct_source_genus_accounting import (
    COMPONENTS,OUTCOMES,exact_genus_bin_decomposition,
    exact_breadth_genus_decomposition,group_summary
)


def _synthetic():
    # Same original source-listing breadth bin, overlapping GenA, and
    # one source-only + one native-only genus. Actual original category state,
    # including negative signed contribution, must be preserved.
    native={"GenA native1","GenA native2","GenB native3"}
    source={"GenA source1","GenC source2"}
    phenotype={
        "GenA native1":1,"GenA native2":1,"GenB native3":0,
        "GenA source1":0,"GenC source2":1}
    listing={s:1 for s in phenotype}
    return native,source,phenotype,listing


def test_source_native_exact_four_signed_components():
    native,source,phenotype,_=_synthetic()
    q=exact_genus_bin_decomposition(native,source,phenotype)
    assert q["original_direct_binary_sector_difference"]==pytest.approx(1/6)
    assert q["shared_genus_composition"]==pytest.approx(1/12)
    assert q["shared_genus_within_species_states"]==pytest.approx(7/12)
    assert q["native_exclusive_genus"]==0
    assert q["source_exclusive_genus_negative"]==pytest.approx(-.5)
    assert q["n_both_status_shared_genera"]==1
    assert q["n_native_unique_genera"]==1
    assert q["n_source_unique_genera"]==1
    assert q["n_distinct_shared_genera_with_actual_state_variation"]==1
    assert sum(q[k] for k in COMPONENTS)==pytest.approx(1/6)


def test_no_shared_genus_reproduces_signed_exclusive_contributions():
    na={"A one","A two"}
    so={"B one","C one"}
    state={"A one":1,"A two":0,"B one":1,"C one":0}
    q=exact_genus_bin_decomposition(na,so,state)
    assert q["n_both_status_shared_genera"]==0
    assert q["shared_genus_composition"]==0
    assert q["shared_genus_within_species_states"]==0
    assert q["native_exclusive_genus"]==pytest.approx(.5)
    assert q["source_exclusive_genus_negative"]==pytest.approx(-.5)
    assert q["original_direct_binary_sector_difference"]==0


def test_fixed_source_listing_bins_reconcile_existing_verified_breadth_difference():
    native,source,state,listing=_synthetic()
    # Second real GIFT source breadth stratum has opposite direction.
    native|={"GenD native4"}
    source|={"GenE source3","GenF source4"}
    state.update({"GenD native4":0,"GenE source3":0,"GenF source4":1})
    listing.update({"GenD native4":12,"GenE source3":12,"GenF source4":12})
    q=exact_breadth_genus_decomposition(native,source,state,listing)
    assert q["n_independently_original_source_listing_breadth_bins_used"]==2
    assert q["n_original_matched_effective_source_species"]==3
    assert sum(q[k] for k in COMPONENTS)==pytest.approx(
        q["original_breadth_matched_direct_difference"])
    assert q["original_breadth_matched_direct_difference"]==pytest.approx(1/9)
    assert q["n_distinct_shared_genera"]==1
    assert q["n_distinct_SHARED_genera_with_REAL_category_variation"]==1
    assert 0<=q["fraction_effective_native_carrier_weight_in_shared_genera"]<=1
    assert 0<=q["fraction_effective_source_candidate_weight_in_shared_genera"]<=1


def test_no_matched_breadth_is_unidentified_not_biological_zero():
    native={"A first"};other={"B first"}
    states={"A first":1,"B first":0}
    listing={"A first":1,"B first":31}
    assert exact_breadth_genus_decomposition(native,other,states,listing) is None


def test_multiple_status_bin_consistency_without_negative_component_censorship():
    # A genus dominated by a different category state can create opposing
    # shared-genus weight and within-genus composition contributions.
    native={"A n1","A n2","B n3"}
    mainland={"A s1","A s2","C s3"}
    states={"A n1":0,"A n2":0,"B n3":1,"A s1":1,"A s2":1,"C s3":0}
    q=exact_genus_bin_decomposition(native,mainland,states)
    assert q["shared_genus_within_species_states"]<0
    assert sum(q[k] for k in COMPONENTS)==pytest.approx(
        q["original_direct_binary_sector_difference"])
    assert q["original_direct_binary_sector_difference"]<0


def test_group_equal_all_components_and_bootstrap_reconcile():
    raw=[]
    native,source,state,listing=_synthetic()
    q=exact_breadth_genus_decomposition(native,source,state,listing)
    for i in range(60):
        raw.append(dict(**q,spatial_block=f"B{i//2}",first_source=f"S{i//3}"))
    z=group_summary(raw,"spatial_block",boot=199,seed=3)
    assert z["status"]=="ORIGINAL_EXACT_SOURCE_BREADTH_BALANCED_GENUS_COMPONENTS"
    assert z["n_original_comparable_island_floras"]==60
    assert z["n_source_groups"]==30
    assert z["four_components_add_exactly_to_direct_observed_delta"]
    assert z["original_breadth_matched_direct_difference_physical_or_source_group_equal_mean"]==pytest.approx(1/6)
    assert z["shared_genus_within_species_states_conditional_group_resampling_95"]==pytest.approx([7/12,7/12])
    assert z["source_exclusive_genus_negative_physical_or_source_group_equal_mean"]==pytest.approx(-.5)
    no=group_summary(raw[:12],"spatial_block",boot=99)
    assert no["status"]=="NOT_IDENTIFIABLE_FEW_INDEPENDENT_SOURCE_BLOCKS"


def test_no_duplicate_original_target_species_and_no_imputation():
    native,source,trait,_=_synthetic()
    with pytest.raises(ValueError):
        exact_genus_bin_decomposition({"GenA native1"},{"GenA native1"},trait)
    with pytest.raises(ValueError):
        exact_genus_bin_decomposition({"GenA new"},{"GenA source1"},trait)
    with pytest.raises(ValueError):
        exact_breadth_genus_decomposition(native,source,trait,
             {"GenA native1":1})
