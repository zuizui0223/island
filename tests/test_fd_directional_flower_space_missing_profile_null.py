"""Pure synthetic exact k/region/genus original categorical sector-absence null tests."""
import numpy as np
import pytest

from fd_directional_flower_space_missing_profile_null import (
    no_category_draw_probability,prepare_regional_reference,
    null_expectation_one_island,SECTORS,
)


def test_valid_original_unique_state_and_carrier_count_null():
    joint={
        "GenA one":("tubular","zygomorphic"),
        "GenA two":("open_radial","actinomorphic"),
        "GenB one":("bell_campanulate","zygomorphic"),
        "GenB two":("open_radial","actinomorphic"),
    }
    source=set(joint)
    ref=prepare_regional_reference(joint,source,"bilateral_tube_or_bell")
    z=null_expectation_one_island(
        {"GenA two","GenB two"},joint,source,"bilateral_tube_or_bell",
        reference=ref)
    assert z["n_region_reference_joint_profiles"]==3
    assert z["n_region_reference_special_joint_profiles"]==2
    assert z["n_observed_joint_profiles"]==1
    assert z["observed_special_region_occupied"]==0
    assert z["expected_occupied_special_uniform_given_k"]==pytest.approx(2/3)
    assert z["expected_occupied_special_category_prevalence_fixed_n"]==pytest.approx(1)
    assert z["expected_occupied_special_regional_genus_swap"]==pytest.approx(1)
    assert z["expected_absence_regional_genus_swap"]==pytest.approx(.25)
    assert z["excess_complete_sector_absence_genus_fixed_n"]==pytest.approx(.75)
    assert z["excess_special_missing_genus_fixed_n"]==pytest.approx(.5)


def test_sampling_without_replacement_checks_finite_profiles():
    assert no_category_draw_probability(4,1,2)==pytest.approx(.5)
    assert no_category_draw_probability(4,2,2)==pytest.approx(1/6)
    assert no_category_draw_probability(3,3,1)==0
    assert no_category_draw_probability(7,0,5)==1
    with pytest.raises(ValueError):no_category_draw_probability(4,2,7)
    with pytest.raises(ValueError):no_category_draw_probability(4,5,2)
    with pytest.raises(ValueError):no_category_draw_probability(4,-1,2)


def test_absence_under_independent_genera_not_infinite_replacement():
    joint={
        "A one":("open_radial","actinomorphic"),
        "A two":("tubular","zygomorphic"),
        "B one":("open_radial","actinomorphic"),
        "B two":("tubular","zygomorphic")
    }
    source=set(joint)
    x=null_expectation_one_island(
        {"A one","B one"},joint,source,"bilateral_tube_or_bell")
    assert x["expected_absence_regional_genus_swap"]==pytest.approx(.25)
    assert x["expected_absence_category_prevalence_fixed_n"]==pytest.approx(1/6)
    assert x["excess_special_missing_genus_fixed_n"]==pytest.approx(.75)


def test_biologically_unsupported_category_not_invented_by_cross_product():
    joint={
        "A one":("tubular","actinomorphic"),
        "B one":("open_radial","zygomorphic"),
    }
    z=prepare_regional_reference(joint,set(joint),"bilateral_tube_or_bell")
    assert z is None
    assert not any(SECTORS["bilateral_tube_or_bell"](v) for v in joint.values())
    with pytest.raises(ValueError):
        null_expectation_one_island(
            {"Alien absent"},joint,set(joint),"bilateral_tube_or_bell")
