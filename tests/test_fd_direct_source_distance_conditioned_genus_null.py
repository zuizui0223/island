"""Synthetic exact closest-GIFT-source kilometre and source-list-breadth nulls."""
from __future__ import annotations
from itertools import combinations
import numpy as np
import pandas as pd
import pytest

from fd_direct_source_distance_conditioned_genus_null import (
    closest_km_bin,source_opportunity_distances,
    exact_taxon_stratum_genus_reference,exact_same_names_breadth_baseline,
    group_equal,DISTANCE_BINS,
)


def test_original_source_real_geographic_distance_frozen_band_boundaries():
    assert closest_km_bin(0)=="lt250km"
    assert closest_km_bin(249.999)=="lt250km"
    assert closest_km_bin(250)=="250_750km"
    assert closest_km_bin(750)=="750_2000km"
    assert closest_km_bin(2000)=="ge2000km"
    assert closest_km_bin(8000)=="ge2000km"
    with pytest.raises(ValueError):closest_km_bin(-1)
    with pytest.raises(ValueError):closest_km_bin(np.nan)
    with pytest.raises(ValueError):closest_km_bin(None)


def test_original_species_must_be_present_in_actual_ranked_source_region():
    byregion={"sourceA":{"A one","B one"},
              "sourceB":{"A one","C one"}}
    closest,ranks=source_opportunity_distances([
        ("sourceA",120.),("sourceB",350.)],byregion)
    assert closest=={"A one":120.,"B one":120.,"C one":350.}
    assert ranks[1]==2 and ranks[2]==1
    with pytest.raises(ValueError):
        source_opportunity_distances([
            ("sourceA",200.),("sourceB",100.)],byregion)
    with pytest.raises(ValueError):
        source_opportunity_distances([
            ("sourceA",float("nan"))],byregion)


def test_exact_taxonomic_genus_rank_distance_null_reconciles_uniform_case():
    focal={"A one","B one"}
    candidate={"A two","C one"}
    strata={s:("source_once","lt250km") for s in focal|candidate}
    v=exact_taxon_stratum_genus_reference(focal,candidate,strata,min_effective=1)
    assert v is not None
    result,u,s=v
    assert u==focal and s==candidate
    assert result["n_original_matched_effective_source_species"]==2
    assert result["n_original_joint_source_breadth_distance_strata"]==1
    assert result["original_actual_focal_shared_source_genera"]==1
    assert result["expected_shared_genera_under_breadth_distance_exact_source_null"]==pytest.approx(2/3)
    assert result["observed_minus_breadth_distance_expected_shared_genera"]==pytest.approx(1/3)
    assert result["original_observed_focal_genus_pair_concentration"]==0
    assert result["expected_focal_genus_pair_concentration_under_breadth_distance_exact_null"]==pytest.approx(1/6)
    z=exact_same_names_breadth_baseline(result,u,s,{x:1 for x in u|s})
    assert z["observed_minus_breadth_only_expected_shared_genera_on_SAME_names"]==pytest.approx(1/3)
    assert z["additional_genus_sharing_expected_from_not_fixing_source_distance"]==pytest.approx(0.)


def test_exact_two_distance_strata_conditioned_taxonomic_expectation():
    focal={"A near","A far"}
    candidate={"B near","C far"}
    strata={
        "A near":("source_once","lt250km"),
        "B near":("source_once","lt250km"),
        "A far":("source_once","750_2000km"),
        "C far":("source_once","750_2000km"),
    }
    v=exact_taxon_stratum_genus_reference(focal,candidate,strata,min_effective=1)
    assert v is not None
    q,u,s=v
    assert q["n_original_joint_source_breadth_distance_strata"]==2
    assert q["n_original_matched_effective_source_species"]==2
    assert q["original_observed_focal_genus_pair_concentration"]==1
    assert q["expected_focal_genus_pair_concentration_under_breadth_distance_exact_null"]==pytest.approx(.25)
    assert q["original_actual_focal_shared_source_genera"]==0
    assert q["expected_shared_genera_under_breadth_distance_exact_source_null"]==pytest.approx(.5)
    same=exact_same_names_breadth_baseline(
        q,u,s,{x:1 for x in u|s})
    assert same["original_breadth_only_expected_shared_genera_on_SAME_names"]>0
    assert same["original_source_breadth_N_named_species_still_same_between_nulls"]


def test_unmatched_distance_bin_is_not_implicitly_a_botanical_absence():
    focal={"A one","B two"}
    candidate={"C three","D four"}
    bins={
        "A one":("source_once","lt250km"),
        "B two":("source_once","lt250km"),
        "C three":("source_once","lt250km"),
        "D four":("source_once","ge2000km"),
    }
    v=exact_taxon_stratum_genus_reference(focal,candidate,bins)
    # One overlapping bin, but an effective source match count below n>=5.
    assert v is None
    with pytest.raises(ValueError):
        exact_taxon_stratum_genus_reference(
            {"A one"},{"A one"},bins)


def test_source_clustering_and_bootstrap_on_anonymous_exact_terms():
    plants={"A one","B one"}
    names={"A two","C one"}
    b={sp:("source_once","lt250km") for sp in plants|names}
    q,u,s=exact_taxon_stratum_genus_reference(plants,names,b,min_effective=1)
    q=exact_same_names_breadth_baseline(q,u,s,{s:1 for s in plants|names})
    r=[dict(**q,block="B"+str(i//2)) for i in range(60)]
    z=group_equal(r,"block",boots=99)
    assert z["status"]=="EXACT_GIFT_DISTANCE_AND_LISTING_BREADTH_SOURCE_GENUS_NULL"
    assert z["n_original_island_floras"]==60
    assert z["n_original_source_groups"]==30
    assert z["observed_minus_breadth_distance_expected_shared_genera_group_equal_mean"]==pytest.approx(1/3)
    assert z["observed_minus_breadth_distance_expected_shared_genera_group_conditional_95"]==pytest.approx([1/3,1/3])
    assert group_equal(r[:12],"block")["status"]=="INSUFFICIENT_ORIGINAL_SOURCE_GEO_BLOCK_COVERAGE"
