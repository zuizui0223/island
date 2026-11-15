"""Exact finite GIFT source draw null: synthetic random sampling and genus preservation.

No protected botanical record names in the test; actual original SOURCE
data always only via SHA-verified CI.
"""
import itertools
from collections import Counter
import math

import numpy as np
import pandas as pd
import pytest

from fd_direct_source_finite_sampling_null import (
    hypergeom_count_variance,prob_no_original_genus_in_native_or_candidate,
    exact_original_breadth_null,aggregate_by_original_clusters,
    METRICS)


def test_closed_form_exact_hypergeometric_native_category_variance():
    assert hypergeom_count_variance(4,2,2)==pytest.approx(1/3)
    assert hypergeom_count_variance(10,0,4)==0
    assert hypergeom_count_variance(10,10,4)==0
    assert hypergeom_count_variance(10,4,0)==0
    assert hypergeom_count_variance(10,4,10)==0
    with pytest.raises(ValueError):
        hypergeom_count_variance(10,11,3)
    with pytest.raises(ValueError):
        hypergeom_count_variance(10,5,11)


def test_exact_single_genus_probability_no_selected_native_or_source():
    p=prob_no_original_genus_in_native_or_candidate(4,2,2)
    assert p==(pytest.approx(1/6),pytest.approx(1/6))
    p=prob_no_original_genus_in_native_or_candidate(4,2,1)
    assert p==(pytest.approx(.5),pytest.approx(.5))
    assert prob_no_original_genus_in_native_or_candidate(4,2,4)==(0.,0.)
    with pytest.raises(ValueError):
        prob_no_original_genus_in_native_or_candidate(4,2,5)


def test_actual_original_island_source_diff_and_genus_null_exact_reconciliation():
    # Source A is shared and variable; B is native exclusive, C source exclusive.
    native={"A n","B n"}
    candidate={"A s","C s"}
    states={"A n":1,"A s":0,"B n":0,"C s":1}
    breadth={x:1 for x in states}
    q=exact_original_breadth_null(native,candidate,states,breadth)
    assert q["original_direct_native_minus_candidate_category_difference"]==pytest.approx(0.)
    assert q["nullA_exact_exchangeable_source_trait_expected_difference"]==0
    assert q["nullA_exact_exchangeable_source_trait_sampling_sd"]==pytest.approx(math.sqrt(1/3))
    assert q["expected_direct_difference_genus_fixed_original_source_null"]==pytest.approx(-.5)
    assert q["nullB_exact_genus_fixed_original_source_trait_sampling_sd"]==pytest.approx(.5)
    assert q["actual_minus_genus_fixed_original_source_null"]==pytest.approx(.5)
    assert q["observed_shared_genus_count"]==1
    assert q["observed_native_exclusive_genus_count"]==1
    assert q["observed_candidate_exclusive_genus_count"]==1
    assert q["nullA_expected_shared_genus_count"]==pytest.approx(2/3)
    assert q["nullA_expected_native_exclusive_genus_count"]==pytest.approx(7/6)
    assert q["nullA_expected_candidate_exclusive_genus_count"]==pytest.approx(7/6)
    assert q["n_original_actual_source_genera_with_variable_species_categories_and_both_memberships"]==1
    assert q["n_original_effectively_breadth_matched_Direct_species"]==2


def test_exact_enumeration_of_source_partition_null_and_genus_richness():
    """Enumerate every equiprobable status assignment; compare analytic moments."""
    species=("A n","A s","B n","C s")
    real_native={"A n","B n"}
    real_source={"A s","C s"}
    states={"A n":1,"A s":0,"B n":0,"C s":1}
    q=exact_original_breadth_null(real_native,real_source,states,
                                  {s:1 for s in species})
    possibilities=[]
    for chosen in itertools.combinations(species,2):
        n=set(chosen);s=set(species)-n
        # This null intentionally does not restrict assignments to originally
        # source-native labels; the original status label is the permuted event.
        d=np.mean([states[x] for x in n])-np.mean([states[x] for x in s])
        ng={x.split()[0] for x in n}
        sg={x.split()[0] for x in s}
        possibilities.append((d,len(ng&sg),len(ng-sg),len(sg-ng)))
    a=np.asarray(possibilities,float)
    assert a[:,0].mean()==pytest.approx(0)
    assert a[:,0].var()==pytest.approx(
        q["nullA_exact_exchangeable_source_trait_sampling_sd"]**2)
    assert a[:,1].mean()==pytest.approx(q["nullA_expected_shared_genus_count"])
    assert a[:,2].mean()==pytest.approx(
        q["nullA_expected_native_exclusive_genus_count"])
    assert a[:,3].mean()==pytest.approx(
        q["nullA_expected_candidate_exclusive_genus_count"])


def test_two_source_listing_breadth_bins_keep_earlier_original_direct_response():
    native={"A native1","B native2","D native3"}
    source={"A source1","C source2","E source3","F source4"}
    states={
        "A native1":1,"B native2":0,"D native3":1,
        "A source1":0,"C source2":1,"E source3":0,"F source4":0}
    breadth={
        "A native1":1,"B native2":1,"D native3":12,
        "A source1":1,"C source2":1,"E source3":12,"F source4":12}
    q=exact_original_breadth_null(native,source,states,breadth)
    assert q["n_original_matched_breadth_bins"]==2
    assert q["n_original_effectively_breadth_matched_Direct_species"]==3
    # First bin: .5-.5=0. Second bin: 1 - 0=1. min-count weights 2 and 1
    assert q["original_direct_native_minus_candidate_category_difference"]==pytest.approx(1/3)
    assert q["original_direct_native_minus_candidate_category_difference"]==pytest.approx(
        q["expected_direct_difference_genus_fixed_original_source_null"]+
        q["actual_minus_genus_fixed_original_source_null"])


def test_genus_invariance_is_null_no_information_not_absent_adaptation():
    nat={"A first","B first"}
    can={"A second","C second"}
    states={s:1 for s in nat|can}
    q=exact_original_breadth_null(nat,can,states,{x:1 for x in states})
    assert q["nullA_exact_exchangeable_source_trait_sampling_sd"]==0
    assert q["nullB_exact_genus_fixed_original_source_trait_sampling_sd"]==0
    assert q["actual_minus_genus_fixed_original_source_null"]==pytest.approx(0)
    assert q["n_original_actual_source_genera_with_variable_species_categories_and_both_memberships"]==0


def test_source_errors_and_missing_trait_labels_never_turned_into_zeros():
    states={"A n":1,"A s":0}
    with pytest.raises(ValueError):
        exact_original_breadth_null({"A n"},{"A n"},states,{"A n":1})
    with pytest.raises(ValueError):
        exact_original_breadth_null({"A n"},{"A s"},states,{"A n":1})
    with pytest.raises(ValueError):
        exact_original_breadth_null({"A n"},{"A unknown"},states,
                                   {"A n":1,"A unknown":1})


def test_group_mean_and_repeated_block_resamples_reconcile():
    nat={"A n","B n"};can={"A s","C s"}
    state={"A n":1,"B n":0,"A s":0,"C s":1}
    v=exact_original_breadth_null(nat,can,state,{x:1 for x in state})
    v.pop("full_exact_real_original_source_category_breadth_bins")
    rows=[]
    for k in range(60):
        rows.append(dict(**v,spatial_block="Block"+str(k//2),
                         original_trait_annotation_coverage_native=.7,
                         original_trait_annotation_coverage_candidate=.6))
    r=aggregate_by_original_clusters(rows,"spatial_block",bootstrap=199)
    assert r["status"]=="SOURCE_BREADTH_FINITE_NULL_OBSERVATIONAL_CLUSTERED_DIFFERENCE"
    assert r["n_islands"]==60
    assert r["n_groups"]==30
    assert r["actual_minus_genus_fixed_original_source_null_group_equal_mean"]==pytest.approx(.5)
    assert r["actual_minus_genus_fixed_original_source_null_original_group_resample_95"]==pytest.approx([.5,.5])
    assert r["original_direct_native_minus_candidate_category_difference_group_equal_mean"]==0
    assert r["mean_original_trait_annotation_coverage_native"]==pytest.approx(.7)
    assert r["mean_original_trait_annotation_coverage_mainland_candidate"]==pytest.approx(.6)
    assert aggregate_by_original_clusters(rows[:12],"spatial_block")[
        "status"]=="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_GROUPS"
