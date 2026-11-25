"""Synthetic checks for fixed-absolute-distance island profile reference N1d."""
import numpy as np
import pytest

from fd.fixed_radius_ses import (
    RADII_KM, fixed_radius_support_summary, matched_radius_summary,
    within_radius_donors,
)


def test_same_absolute_radius_excludes_focal_and_out_of_regime():
    meta={name:{"analysis_regime":regime} for name,regime in (
        ("f","tropical"),("near1","tropical"),("near2","tropical"),
        ("just_out","tropical"),("far","tropical"),
        ("other_regime","northern_midlatitude"),
    )}
    geo={
        "f":{"island_latitude":"0","island_longitude":"0"},
        "near1":{"island_latitude":"0","island_longitude":"1"},
        "near2":{"island_latitude":"0","island_longitude":"8"},
        "just_out":{"island_latitude":"0","island_longitude":"10"},
        "far":{"island_latitude":"0","island_longitude":"30"},
        "other_regime":{"island_latitude":"0","island_longitude":"1"},
    }
    near,farthest=within_radius_donors("f","tropical",meta,geo,1000)
    assert near==["near1","near2"]
    assert farthest==pytest.approx(889.56,abs=.1)
    extended,distance=within_radius_donors("f","tropical",meta,geo,2000)
    assert extended==["near1","near2","just_out"]
    assert distance==pytest.approx(1111.95,abs=.1)
    assert RADII_KM==(1000,2000)
    with pytest.raises(ValueError):
        within_radius_donors("f","tropical",meta,geo,-1)
    with pytest.raises(ValueError):
        within_radius_donors("missing","tropical",meta,geo,1000)


def test_radius_support_reports_eligible_counts_and_suppresses_small_geo_cells():
    rows=[]
    for i in range(6):
        rows.append(dict(scope="all", geological_origin="oceanic",
                         radius_km=1000,n_donor_islands=4,
                         n_donor_profiles=40,farthest_km=400,
                         has_full=(i<3),has_equal30=True))
        rows.append(dict(scope="all", geological_origin="oceanic",
                         radius_km=2000,n_donor_islands=1,
                         n_donor_profiles=0,farthest_km=1500,
                         has_full=False,has_equal30=False))
    report=fixed_radius_support_summary(rows)
    assert len(report)==2
    score={x["radius_km"]:x for x in report}
    assert score[1000]["n_equal30_estimable"]==6
    assert score[1000]["median_donor_island_count"]==4.
    assert score[1000]["n_full_estimable"]==3
    assert score[2000]["status"]=="small_cell_suppressed"
    assert "median_donor_island_count" not in score[2000]
    assert all("island_id" not in row for row in report)


def test_fixed_radius_matched_same_island_and_original_observation():
    records=[]
    for j in range(7):
        for null, ses, ref in (
            ("N1b_regional_islandfreq_LOO",-.8,.32),
            ("N1d_radius1000km_geo_LOO",-.2,.28),
            ("N1d_radius2000km_geo_LOO",-.5,.30),
        ):
            if j>=4 and null=="N1d_radius2000km_geo_LOO":
                continue
            records.append(dict(scope="all",mode="equal30",
                endpoint="mean_nearest_profile_distance",region="tropical",
                geological_origin="oceanic",island_id=f"secret{j}",
                null=null,observed=.27,null_mean=ref,
                delta=.27-ref,ses=ses))
    got=matched_radius_summary(records)
    assert len(got)==4  # two radii x pooled/regime
    s={r:(e) for r,e in ((x["radius_km"],x) for x in got
                         if x["region"]=="all")}
    assert s[1000]["n_islands_both_nulls"]==7
    assert s[1000]["median_within_island_SES_change"]==pytest.approx(.6)
    assert s[1000]["median_N1d_minus_N1b_expected_D"]==pytest.approx(-.04)
    assert s[1000]["fraction_negative_N1d"]==1.
    assert s[2000]["n_islands_both_nulls"]==4
    assert s[2000]["status"]=="small_cell_suppressed"
    assert "median_N1d_SES" not in s[2000]
    assert all("island_id" not in row for row in got)

    corrupted=[dict(x) for x in records]
    corrupted[1]["observed"]=.25
    with pytest.raises(AssertionError):
        matched_radius_summary(corrupted)


def test_unscored_region_and_duplicate_reference_dont_hide_errors():
    records=[dict(scope="all",mode="full",endpoint="mean_pairwise_distance",
                  region="tropical",geological_origin="continental",
                  island_id="one",observed=.3,null="N1b_regional_islandfreq_LOO",
                  ses=.2,delta=.02,null_mean=.28)]
    assert matched_radius_summary(records)==[]
    with pytest.raises(ValueError):
        matched_radius_summary(records+records)
