"""Synthetic-only tests; real plant and animal identifiers are never committed."""
import numpy as np
import pandas as pd
import pytest
from fd.guild_sector_alignment import run


def make_data(n=180, seed=3, qualified=True, correlated=True):
    rng=np.random.default_rng(seed)
    i=np.arange(n)
    a=rng.gamma(2.,3.,size=n)
    b=(a+.2*rng.standard_normal(n)) if correlated else rng.gamma(2.,3.,size=n)
    island_n=rng.integers(25,51,size=n)
    source_n=rng.integers(20,61,size=n)
    island_pos=rng.binomial(island_n, np.clip(.3+.12*np.log1p(a)+.025*np.log1p(b),.05,.9))
    source_pos=rng.binomial(source_n,.32)
    base = pd.DataFrame(dict(
        island_id=[f"synthetic_{k}" for k in i],region="tropical",
        scope="all_recorded",source_mode="geo_k5",
        sector="actinomorphic_deep_tube",
        spatial_block=[f"b{k%18}" for k in i],
        first_source=[f"s{k%15}" for k in i],
        island_total=60,source_total=80,island_pair_n=island_n,
        source_pair_n=source_n,island_pair_positive=island_pos,
        source_pair_positive=source_pos,bird_island_richness=a+1,
        bird_source_richness=3.0,
        sphingid_island_richness=b+1,sphingid_source_richness=3.0,
        bird_range_qualified=qualified,sphingid_model_qualified=qualified,
        genus_reference_qualified=True,
        genus_reference_expected_delta=rng.normal(0,.04,size=n),
        log1p_distance_km=rng.uniform(0.4,7,size=n),
        log1p_area_km2=rng.uniform(1.,14.,size=n),
        absolute_latitude=rng.uniform(0,25,size=n),
        temperature_pc1=rng.normal(size=n),
        precipitation_pc1=rng.normal(size=n)))
    return base


def test_original_population_gate():
    with pytest.raises(ValueError,match="4505"):
        run(make_data(),bootstrap=100)


def test_independent_guild_map_required():
    ans=run(make_data(qualified=False),require_population=False,bootstrap=100)
    assert ans["status"]=="NOT_IDENTIFIABLE_GUILD_MAP_OR_CLUSTER_SUPPORT"
    assert ans["support"]["n_map_and_direct_qualified"]==0
    assert "fits" not in ans


def test_nested_direct_count_validation():
    d=make_data()
    d.loc[0,"island_pair_positive"]=d.loc[0,"island_pair_n"]+1
    with pytest.raises(ValueError,match="nested original recorded"):
        run(d,require_population=False,bootstrap=100)


def test_usable_map_yields_anonymous_model_and_missingness_bounds():
    d=make_data(correlated=False)
    ans=run(d,require_population=False,bootstrap=100)
    assert ans["status"]=="EXPLORATORY_GUILD_ASSOCIATION_ESTIMABLE"
    assert set(ans["fits"])=={"spatial_block","first_source"}
    assert "island_id" not in str(ans)
    low=ans["missingness_bounds_on_mean_delta"]["lower_unweighted"]
    high=ans["missingness_bounds_on_mean_delta"]["upper_unweighted"]
    assert low < 0 < high
    assert ans["missingness_bounds_on_mean_delta"]["sign_robust_to_arbitrary_missing_binary_states"] is False


def test_unqualified_genus_reference_is_not_usable():
    d=make_data(correlated=False)
    d["genus_reference_qualified"]=False
    ans=run(d,require_population=False,bootstrap=100)
    assert ans["status"]=="NOT_IDENTIFIABLE_GUILD_MAP_OR_CLUSTER_SUPPORT"
    assert ans["support"]["n_map_and_direct_qualified"]==0
    assert "fits" not in ans


def test_missing_genus_baseline_column_fails_closed():
    d=make_data(correlated=False).drop(columns=["genus_reference_expected_delta"])
    with pytest.raises(ValueError,match="genus_reference_expected_delta"):
        run(d,require_population=False,bootstrap=100)


def test_genus_baseline_fraction_must_be_valid():
    d=make_data(correlated=False)
    d.loc[0,"genus_reference_expected_delta"]=1.5
    with pytest.raises(ValueError,match="genus expected"):
        run(d,require_population=False,bootstrap=100)


def test_near_collinearity_is_rejected():
    d=make_data()
    d["sphingid_island_richness"]=d["bird_island_richness"]
    d["sphingid_source_richness"]=d["bird_source_richness"]
    ans=run(d,require_population=False,bootstrap=100)
    assert ans["status"]=="NOT_IDENTIFIABLE_CLUSTER_MODEL"
    assert ans["fits"]["spatial_block"]["status"]=="RANK_DEFICIENT"
