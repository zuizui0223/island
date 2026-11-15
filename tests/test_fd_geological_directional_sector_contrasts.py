"""Source-geometry falsification tests for original oceanic/continental sectors."""
import numpy as np
import pandas as pd
import pytest
from fd_geological_directional_sector_contrasts import (
    SECTORS,valid_sector_state,sector_roster,genus_ref,genus_expected_share,
    source_group_cluster_fit,MIN_ISLANDS_PER_ORIGIN,
)


def test_ontologies_not_mistake_geological_distance_for_origin():
    assert len([r for r in SECTORS if r[3]])==7
    assert valid_sector_state("floral_symmetry","zygomorphic","zygomorphic")
    assert not valid_sector_state("floral_symmetry","zygomorphic","actinomorphic")
    assert valid_sector_state("floral_form","tube_or_bell","bell_campanulate")
    assert not valid_sector_state("floral_form","tube_or_bell","open_radial")
    assert valid_sector_state("tube_depth_class","deep","deep")
    assert valid_sector_state("floral_form__floral_symmetry","joint_special",
                              ("tubular","zygomorphic"))
    assert not valid_sector_state("floral_form__floral_symmetry","joint_special",
                              ("tubular","actinomorphic"))
    assert valid_sector_state("floral_form__floral_symmetry","joint_open",
                              ("open_radial","actinomorphic"))
    assert valid_sector_state("self_incompatibility","SC","SC")
    assert not valid_sector_state("self_incompatibility","SC","SI")
    assert valid_sector_state("autonomous_selfing_capacity","autonomous","autonomous")
    assert not valid_sector_state("autonomous_selfing_capacity","autonomous","absent")
    with pytest.raises(ValueError):
        valid_sector_state("floral_form","bee_guild","tubular")


def test_binary_source_ambiguous_and_joint_direct_exclusion():
    source={
      "self_incompatibility":{"A one":"SC","A two":"SI",
                              "B three":"mixed_or_variable"},
      "autonomous_selfing_capacity":{"A one":"absent","A two":"autonomous",
                                     "B three":"delayed"},
      "floral_form":{"A one":"open_radial","A two":"tubular",
                     "B three":"bell_campanulate"},
      "floral_symmetry":{"A one":"actinomorphic","B three":"zygomorphic"}}
    assert sector_roster(source,"self_incompatibility")=={
        "A one":"SC","A two":"SI"}
    assert sector_roster(source,"autonomous_selfing_capacity")=={
        "A one":"absent","A two":"autonomous"}
    assert sector_roster(source,"floral_form__floral_symmetry")=={
        "A one":("open_radial","actinomorphic"),
        "B three":("bell_campanulate","zygomorphic")}


def test_regional_genus_reference_exact_target_count_contributes():
    t={"A one":"SC","A two":"SI","B one":"SC","B two":"SC"}
    ref=genus_ref(t,set(t),"SC","self_incompatibility")
    assert ref["A"]==pytest.approx(.5)
    assert ref["B"]==pytest.approx(1.)
    assert genus_expected_share({"A one","B one"},ref)==pytest.approx(.75)
    assert genus_expected_share({"A two","B two"},ref)==pytest.approx(.75)
    with pytest.raises(AssertionError):
        genus_expected_share({"NotPresent case"},ref)


def make_data(seed=17,n=240):
    rng=np.random.default_rng(seed)
    idx=np.arange(n)
    clusters=idx//4
    origin=np.where(clusters%2==0,"oceanic","continental")
    oce=(origin=="oceanic").astype(float)
    xdist=rng.normal(4.0,.9,n)+.6*oce
    xarea=rng.normal(3.,.7,n)-.25*oce
    cov=rng.uniform(.2,.85,n)
    xn=np.log(rng.integers(5,35,n))
    lat=rng.uniform(15,50,n)
    y=(.15+.075*oce-.013*xdist+.018*xarea
       +.08*cov+.007*xn+.0005*lat+rng.normal(0,.014,n))
    return pd.DataFrame(dict(
        origin=origin,spatial_block=["B"+str(i) for i in clusters],
        first_source=["S"+str(i//2) for i in clusters],
        x_continent_distance=xdist,x_area=xarea,
        x_coverage=cov,x_direct_n=xn,x_latitude=lat,
        response=y,genus_relative_response=y-.13))


def test_oceanic_adjusted_effect_cluster_and_sampling_are_grounded():
    frame=make_data()
    q=source_group_cluster_fit(frame,"response","spatial_block",
                               n_boot=99,seed=20261115)
    assert q["status"]=="EXPLORATORY_ADJUSTED_GEOLOGICAL_ORIGIN_CATEGORY_CONTRAST"
    assert q["n_islands"]==240
    assert q["n_oceanic"]==120 and q["n_continental"]==120
    assert q["n_oceanic_clusters"]==30 and q["n_continental_clusters"]==30
    assert q["adjusted_oceanic_minus_continental_fraction"]==pytest.approx(.075,abs=.018)
    assert q["origin_cluster_999resample_quantile95"][0]>0
    assert q["n_valid_cluster_resamples"]>=80
    q2=source_group_cluster_fit(frame,"response","first_source",
                                n_boot=99,seed=7)
    assert q2["n_clusters"]==30
    assert q2["n_oceanic_clusters"]>=5 and q2["n_continental_clusters"]>=5
    assert q2["adjusted_oceanic_minus_continental_fraction"]==pytest.approx(.075,abs=.018)


def test_missing_source_cluster_does_not_censor_physical_effect():
    frame=make_data()
    frame.loc[frame.index[:30],"first_source"]=None
    a=source_group_cluster_fit(frame,"response","spatial_block",n_boot=99)
    b=source_group_cluster_fit(frame,"response","first_source",n_boot=99)
    assert a["n_islands"]==240
    assert b["n_islands"]==210


def test_insufficient_two_origin_or_spatial_blocks_fail_closed():
    frame=make_data(n=240)
    frame.loc[frame.origin.eq("continental"),"origin"]="oceanic"
    a=source_group_cluster_fit(frame,"response","spatial_block",n_boot=99)
    assert a["status"]=="INSUFFICIENT_GEOLOGY_ORIGIN_SOURCE_SUPPORT"
    assert a["n_continental"]==0
    frame=make_data(n=32)
    z=source_group_cluster_fit(frame,"response","spatial_block",n_boot=99)
    assert z["status"]=="INSUFFICIENT_GEOLOGY_ORIGIN_SOURCE_SUPPORT"
