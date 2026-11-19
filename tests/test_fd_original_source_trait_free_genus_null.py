"""Independent synthetic original botanical genus availability source null, no floral traits.

Check exact genus overlap, taxonomic clustering, real source listing bins
and SAME original direct-qualified island masks.
"""
from __future__ import annotations
from collections import Counter
import itertools

import numpy as np
import pandas as pd
import pytest

from fd_original_source_trait_free_genus_null import (
    exact_taxonomic_source_null,original_trait_island_eligibility,
    source_group_summary
)


def test_original_taxonomic_genus_overlap_and_simpson_exact_enumeration():
    native={"A one","B one"}
    candidate={"A two","C one"}
    breadth={sp:1 for sp in native|candidate}
    q=exact_taxonomic_source_null(native,candidate,breadth)
    assert q["n_original_effective_breadth_matched_species"]==2
    assert q["n_original_island_accepted_names_in_matched_bins"]==2
    assert q["observed_original_shared_genus_count"]==1
    assert q["exact_null_expected_original_shared_genus_count"]==pytest.approx(2/3)
    assert q["shared_original_genus_count_minus_finite_source_expectation"]==pytest.approx(1/3)
    assert q["observed_original_island_genus_pair_concentration"]==0.
    assert q["exact_null_expected_original_island_genus_pair_concentration"]==pytest.approx(1/6)
    assert q["observed_minus_expected_original_island_genus_pair_concentration"]==pytest.approx(-1/6)
    u=tuple(native|candidate)
    drawn=[]
    for subset in itertools.combinations(u,2):
        f=set(subset)
        s=set(u)-f
        g={x.split()[0] for x in f}
        h={x.split()[0] for x in s}
        p=1. if len(g)==1 else 0.
        drawn.append((len(g&h),p))
    a=np.asarray(drawn)
    assert a[:,0].mean()==pytest.approx(q["exact_null_expected_original_shared_genus_count"])
    assert a[:,1].mean()==pytest.approx(q["exact_null_expected_original_island_genus_pair_concentration"])


def test_original_breadth_independent_bins_cross_genus_pair_probability():
    nat={"A first","A second"}
    src={"B first","C second"}
    bins={"A first":1,"B first":1,
          "A second":12,"C second":12}
    q=exact_taxonomic_source_null(nat,src,bins)
    assert q["n_original_source_breadth_bins_with_accepted_names"]==2
    assert q["n_original_effective_breadth_matched_species"]==2
    assert q["observed_original_island_genus_pair_concentration"]==1
    assert q["exact_null_expected_original_island_genus_pair_concentration"]==pytest.approx(.25)
    assert q["observed_minus_expected_original_island_genus_pair_concentration"]==pytest.approx(.75)
    assert q["observed_original_shared_genus_count"]==0
    assert q["exact_null_expected_original_shared_genus_count"]==pytest.approx(.5)


def test_original_source_all_unique_genera_has_no_concentration():
    n={"A first","B first"}
    s={"C first","D first"}
    q=exact_taxonomic_source_null(n,s,{k:1 for k in n|s})
    assert q["observed_original_shared_genus_count"]==0
    assert q["exact_null_expected_original_shared_genus_count"]==0
    assert q["observed_original_island_genus_pair_concentration"]==0
    assert q["exact_null_expected_original_island_genus_pair_concentration"]==0
    assert q["shared_original_genus_count_minus_finite_source_expectation"]==0


def test_trait_status_only_selects_original_island_mask_not_null_outcome():
    n={f"A species{k}" for k in range(5)}
    m={f"B species{k}" for k in range(10)}
    breadth={sp:1 for sp in n|m}
    real={sp:int(k%2) for k,sp in enumerate(sorted(n|m))}
    assert original_trait_island_eligibility(n,m,breadth,real)
    assert not original_trait_island_eligibility(
        n,m,breadth,{k:v for k,v in real.items() if k not in list(n)[:1]})
    other={sp:int(k%3==1) for k,sp in enumerate(sorted(n|m))}
    assert original_trait_island_eligibility(n,m,breadth,other)
    q=exact_taxonomic_source_null(n,m,breadth)
    assert q["n_original_island_accepted_names_in_matched_bins"]==5


def test_original_island_name_overlap_and_missing_source_breadth_do_not_count_zero():
    a={"A one"};b={"A two"}
    with pytest.raises(ValueError):
        exact_taxonomic_source_null(a,a,{"A one":1})
    with pytest.raises(ValueError):
        exact_taxonomic_source_null(a,b,{"A one":1})
    assert exact_taxonomic_source_null(a,b,{"A one":1,"A two":12}) is None


def test_original_group_equal_real_botanical_genus_and_source_statistics():
    a={"A one","B one"}
    b={"A two","C one"}
    q=exact_taxonomic_source_null(a,b,{sp:1 for sp in a|b})
    rows=[dict(**q,spatial_block=f"B{i//2}") for i in range(60)]
    r=source_group_summary(rows,"spatial_block",bootstrap=99)
    assert r["status"]=="EXACT_SOURCE_LIST_BREADTH_ALL_BOTANICAL_GENUS_NULL"
    assert r["n_islands"]==60 and r["n_clusters"]==30
    assert r["shared_original_genus_count_minus_finite_source_expectation_group_equal_mean"]==pytest.approx(1/3)
    assert r["shared_original_genus_count_minus_finite_source_expectation_group_conditional_95"]==pytest.approx([1/3,1/3])
    assert r["observed_minus_expected_original_island_genus_pair_concentration_group_equal_mean"]==pytest.approx(-1/6)
    assert source_group_summary(rows[:12],"spatial_block")["status"]==(
        "INSUFFICIENT_ORIGINAL_INDEPENDENT_GEOLOGICAL_BLOCK_SUPPORT")
