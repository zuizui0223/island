"""Original source breadth matching on actual named native GIFT botanical species."""
import pytest
from fd_direct_source_breadth_matched import (
    listing_breadth_bin,breadth_matched_sector_difference,
)


def test_a_priori_original_GIFT_region_listing_breadth_categories():
    assert listing_breadth_bin(1)=="source_once"
    assert listing_breadth_bin(2)=="source_2_3"
    assert listing_breadth_bin(3)=="source_2_3"
    assert listing_breadth_bin(4)=="source_4_10"
    assert listing_breadth_bin(10)=="source_4_10"
    assert listing_breadth_bin(11)=="source_11_30"
    assert listing_breadth_bin(30)=="source_11_30"
    assert listing_breadth_bin(31)=="source_31plus"
    assert listing_breadth_bin(748)=="source_31plus"
    with pytest.raises(ValueError):listing_breadth_bin(0)
    with pytest.raises(ValueError):listing_breadth_bin(-2)


def test_breadth_matched_diff_removes_documented_source_frequency_imbalance():
    native={"A native1","A native2","B native3"}
    source={"C src1","C src2","D src3","E src4"}
    states={
        "A native1":1,"A native2":1,"B native3":0,
        "C src1":0,"C src2":0,"D src3":1,"E src4":1}
    breadth={
        "A native1":1,"A native2":1,"B native3":12,
        "C src1":1,"C src2":1,"D src3":12,"E src4":400}
    z=breadth_matched_sector_difference(native,source,states,breadth)
    assert z["n_original_breadth_overlapping_bin_support"]==2
    assert z["n_effectively_matched_source_listing_breadth_carriers"]==3
    assert z["source_breadth_matched_native_minus_candidate_sector_share"]==pytest.approx(1/3)
    assert z["same_island_original_unmatched_native_minus_candidate_sector_share"]==pytest.approx(2/3-1/2)
    assert z["fraction_native_carriers_from_overlapping_source_breadth_bins"]==1.
    assert z["fraction_mainland_candidate_carriers_from_overlap_breadth_bins"]==.75


def test_no_matching_source_breadth_is_not_biological_zero():
    states={"A one":1,"B two":0}
    breadth={"A one":1,"B two":400}
    z=breadth_matched_sector_difference(
        {"A one"},{"B two"},states,breadth)
    assert z["n_effectively_matched_source_listing_breadth_carriers"]==0
    assert z["source_breadth_matched_native_minus_candidate_sector_share"] is None
    assert z["same_island_original_unmatched_native_minus_candidate_sector_share"]==1.


def test_absent_direct_labels_and_duplicate_source_names_rejected():
    with pytest.raises(ValueError):
        breadth_matched_sector_difference(
            {"A one"},{"A one"},{"A one":1},{"A one":2})
    with pytest.raises(ValueError):
        breadth_matched_sector_difference(
            {"A one"},{"B two"},{"A one":1},{"A one":1,"B two":2})
    assert breadth_matched_sector_difference(
        set(),{"B two"},{"B two":0},{"B two":1}) is None
