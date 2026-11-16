"""Synthetic WCVP exact accepted native L3 filters, source strata and sample controls."""
from __future__ import annotations
import zipfile
from math import comb
import itertools
import numpy as np
import pandas as pd
import pytest

from fd_wcvp_v16_native_breadth_null import (
    native_l3_bin,accepted_wcvp_native_l3_range,
    exact_genus_strata,geography_matched_one_island,grouped_summary
)


def test_predeclared_real_native_WCVP_L3_strata():
    assert native_l3_bin(1)=="nativeL3_1"
    assert native_l3_bin(2)=="nativeL3_2_3"
    assert native_l3_bin(3)=="nativeL3_2_3"
    assert native_l3_bin(4)=="nativeL3_4_10"
    assert native_l3_bin(10)=="nativeL3_4_10"
    assert native_l3_bin(11)=="nativeL3_11_30"
    assert native_l3_bin(31)=="nativeL3_31plus"
    with pytest.raises(ValueError):native_l3_bin(0)


def test_exact_official_WCVP_accepted_species_native_l3_filters(tmp_path):
    path=tmp_path/"source.zip"
    names=pd.DataFrame([
       dict(plant_name_id="1",taxon_rank="Species",taxon_status="Accepted",
            taxon_name="A one"),
       dict(plant_name_id="2",taxon_rank="Species",taxon_status="Accepted",
            taxon_name="A two"),
       dict(plant_name_id="3",taxon_rank="Species",taxon_status="Synonym",
            taxon_name="B discarded"),
       dict(plant_name_id="4",taxon_rank="Genus",taxon_status="Accepted",
            taxon_name="B"),
       dict(plant_name_id="5",taxon_rank="Species",taxon_status="Accepted",
            taxon_name="C noNative"),
    ])
    D=pd.DataFrame([
       dict(plant_name_id="1",area_code_l3="USA",introduced=0,extinct=0,location_doubtful=0),
       dict(plant_name_id="1",area_code_l3="USA",introduced=0,extinct=0,location_doubtful=0),
       dict(plant_name_id="1",area_code_l3="JPN",introduced=0,extinct=0,location_doubtful=0),
       dict(plant_name_id="2",area_code_l3="JPN",introduced=1,extinct=0,location_doubtful=0),
       dict(plant_name_id="2",area_code_l3="CHN",introduced=0,extinct=1,location_doubtful=0),
       dict(plant_name_id="2",area_code_l3="BRA",introduced=0,extinct=0,location_doubtful=1),
       dict(plant_name_id="2",area_code_l3="PER",introduced=0,extinct=0,location_doubtful=0),
       dict(plant_name_id="3",area_code_l3="ARG",introduced=0,extinct=0,location_doubtful=0),
       dict(plant_name_id="5",area_code_l3="FRA",introduced=1,extinct=0,location_doubtful=0),
    ])
    with zipfile.ZipFile(path,"w") as z:
        z.writestr("wcvp_names.csv",names.to_csv(index=False,sep="|"))
        z.writestr("wcvp_distribution.csv",D.to_csv(index=False,sep="|"))
    result,manifest=accepted_wcvp_native_l3_range(
       path,{"A one","A two","B discarded","C noNative","A missing"})
    assert result=={"A one":2,"A two":1}
    assert manifest["wcvp_accepted_exact_species_name_matches"]==3
    assert manifest["wcvp_accepted_exact_species_with_positive_native_TDWG_L3"]==2
    assert manifest["wcvp_exact_accepted_names_without_positive_native_L3"]==1
    assert manifest["n_input_original_GIFT_distinct_accepted_species"]==5


def test_exact_genus_sampling_same_species_names_across_source_nulls():
    # A genus two species, B/C each singleton; actual one shared genus.
    native={"A one","B one"}
    candidate={"A two","C one"}
    full=native|candidate
    strata={s:("source_once","lt250km","nativeL3_1") for s in full}
    test=exact_genus_strata(native,candidate,strata,min_effective=2)
    assert test is not None
    q,ns,cs=test
    assert ns==native and cs==candidate
    assert q["original_actual_shared_botanical_genera"]==1
    assert q["expected_shared_botanical_genera_exact_finite_source"]==pytest.approx(2/3)
    assert q["observed_minus_finite_expected_shared_genera"]==pytest.approx(1/3)
    baseline=exact_genus_strata(ns,cs,{s:strata[s][:2] for s in full},
                                min_effective=2)
    assert baseline[0]["observed_minus_finite_expected_shared_genera"]==pytest.approx(
        q["observed_minus_finite_expected_shared_genera"])
    with pytest.raises(ValueError):
        exact_genus_strata(native,candidate,{
            s:("source_once","lt250km","unlisted") for s in full})


def test_native_range_breadth_can_change_null_but_sample_name_set_remains_same():
    native={"A one","B one"}
    candidate={"A two","C one"}
    full=native|candidate
    three={
        "A one":("source_once","lt250km","nativeL3_1"),
        "B one":("source_once","lt250km","nativeL3_1"),
        "A two":("source_once","lt250km","nativeL3_1"),
        "C one":("source_once","lt250km","nativeL3_1"),
    }
    r=exact_genus_strata(native,candidate,three,min_effective=1)
    assert r[0]["n_original_botanical_focal_names_in_eligible_source_strata"]==2
    strict=three.copy()
    strict["C one"]=("source_once","lt250km","nativeL3_31plus")
    z=exact_genus_strata(native,candidate,strict,min_effective=1)
    assert "C one" not in z[2]
    assert z[0]["n_original_botanical_source_only_names_in_eligible_strata"]==1


def test_nearest_source_and_WCVP_count_missing_is_not_zero_native_area():
    a={"A first","B first","C first","D first","E first"}
    b={f"Z species{i}" for i in range(10)}
    alln=a|b
    gift={s:1 for s in alln}
    km={s:100. for s in alln}
    wcvp={s:1 for s in alln}
    x=geography_matched_one_island(a,b,gift,km,wcvp)
    assert x["n_original_focal_with_WCVP_positive_native_range"]==5
    assert x["n_original_source_only_with_WCVP_positive_native_range"]==10
    assert x["WCVP_source_focal_retention_after_triple_overlap"]==1.
    assert x["WCVP_source_candidate_retention_after_triple_overlap"]==1.
    assert x["observed_minus_range_and_source_expected_shared_genera"]==pytest.approx(
       x["observed_minus_source_km_breadth_only_expected_shared_genera_ON_SAME_SPECIES"])
    incomplete={s:1 for s in a}
    assert geography_matched_one_island(a,b,gift,km,incomplete) is None


def test_original_grouping_conditional_resampling_and_support_gate():
    rows=[]
    for i in range(70):
        rows.append(dict(
            spatial_block=f"Block{i//2}",
            observed_minus_range_and_source_expected_shared_genera=-2.,
            observed_minus_source_km_breadth_only_expected_shared_genera_ON_SAME_SPECIES=-1.,
            additional_expected_genus_sharing_after_WCVP_native_breadth_conditioning=1.,
        ))
    z=grouped_summary(rows,"spatial_block",boots=99)
    assert z["status"]=="ORIGINAL_WCVP_NATIVE_RANGE_CONDITIONAL_SOURCE_GENUS_NULL"
    assert z["n_original_matched_floras"]==70
    assert z["n_original_groups"]==35
    assert z["observed_minus_range_and_source_expected_shared_genera_conditional_group_95"]==pytest.approx([-2,-2])
    small=grouped_summary(rows[:10],"spatial_block",boots=99)
    assert small["status"]=="INSUFFICIENT_NATIVE_RANGE_SOURCE_SUPPORT"


def test_exact_cross_stratum_genus_pair_expectation_matches_full_enumeration():
    native={"A one","A two"}
    source={"B one","C two"}
    strata={
      "A one":("source_once","lt250km","nativeL3_1"),
      "B one":("source_once","lt250km","nativeL3_1"),
      "A two":("source_11_30","750_2000km","nativeL3_2_3"),
      "C two":("source_11_30","750_2000km","nativeL3_2_3")
    }
    z=exact_genus_strata(native,source,strata,min_effective=2)
    assert z is not None
    v=z[0]
    assert v["n_original_source_native_breadth_km_WCVP_strata"]==2
    assert v["observed_minus_expected_focal_genus_concentration"]==pytest.approx(.75)
    assert v["expected_shared_botanical_genera_exact_finite_source"]==pytest.approx(.5)
