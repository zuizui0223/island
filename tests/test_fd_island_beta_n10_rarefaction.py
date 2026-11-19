"""Synthetic full-flora same-mask equal-10 source null tests (no private IDs)."""
import json
import pandas as pd
import pytest
from fd_island_beta_n10_rarefaction import analyze


def _fixture(n_species=12):
    states={"flower_primary_color":{},"floral_form":{}}
    flora=[];geo=[]
    for i in range(8):
        for k in range(n_species):
            sp=f"Gen{k%4} species{i}_{k}"
            states["flower_primary_color"][sp]=(
                "white" if (i+k)%2 else "yellow")
            states["floral_form"][sp]=(
                "open" if (i*3+k)%4 else "tube")
            flora.append(dict(island_id=f"I{i}",accepted_species=sp,
                              origin_status="native"))
        geo.append(dict(island_id=f"I{i}",
                        analysis_regime="northern_midlatitude",
                        spatial_block=f"blk{i//2}",
                        island_latitude=30, island_longitude=i*8))
    return states,pd.DataFrame(flora),pd.DataFrame(geo)


def test_same_pair_and_fixed_local_genus_pool_survive_thinning():
    s,f,g=_fixture(12)
    q=analyze(s,f,g,n_thinnings=3,n_matches=5,min_pairs=1)
    assert q["status"]=="ORIGINAL_N10_COVERAGE10_EQUAL10_SPATIAL_STRESS"
    assert q["n_full_botanical_islands"]==8
    assert q["n_original_region_n5_reference_islands"]==8
    assert q["n_original_region_genus_donor_species"]==96
    assert q["n_original_fixed_full_flora_species_disjoint_pairs"]>=1
    assert q["n_maximally_block_disjoint_pairs"]>=1
    assert q["n_deterministic_10species_thinnings"]==3
    assert q["n_outcome_blind_10degree_block_matchings"]==5
    assert q["claims"]["exact_same_full_flora_species_disjoint_pair_graph"]
    text=json.dumps(q)
    assert all(k not in text for k in (
        '"island_id"','"accepted_species"','"spatial_block"','"I0"'))


def test_identical_exact_ten_names_never_changes_observed_null():
    s,f,g=_fixture(10)
    q=analyze(s,f,g,n_thinnings=4,n_matches=5,min_pairs=1)
    assert q["original_full_carrier_pair_equal_mean_residual_D2"]==pytest.approx(
        q["equal10_pair_equal_median_residual_D2"])
    assert q["original_full_carrier_block_disjoint_median_residual_D2"]==pytest.approx(
        q["equal10_block_disjoint_median_residual_D2"])
    assert q["equal10_per_thinning_maximal_matching_median_min"]==pytest.approx(
        q["equal10_per_thinning_maximal_matching_median_max"])


def test_determinism_under_source_plant_row_order():
    s,f,g=_fixture(12)
    a=analyze(s,f,g,n_thinnings=3,n_matches=5,min_pairs=1)
    b=analyze(s,f.sample(frac=1,random_state=12),
              g.sample(frac=1,random_state=3),
              n_thinnings=3,n_matches=5,min_pairs=1)
    assert a==b


def test_source_coverage_and_spatial_block_fail_closed():
    s,f,g=_fixture(12)
    with pytest.raises(ValueError):
        analyze(s,pd.concat([f,f.iloc[:1]]),g,min_pairs=1)
    with pytest.raises(ValueError):
        analyze(s,f,g.drop(columns=["island_longitude"]),min_pairs=1)
    with pytest.raises(ValueError):
        analyze(s,f,g,n_thinnings=0,min_pairs=1)
    with pytest.raises(ValueError):
        analyze(s,f,g,n_matches=0,min_pairs=1)
