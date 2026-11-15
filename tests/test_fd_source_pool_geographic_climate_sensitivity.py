"""Original botanical source candidate geography-vs-climate robustness math."""
import numpy as np
import pandas as pd
import pytest
from fd_source_pool_geographic_climate_sensitivity import (
    source_species_min_km,paired_mode_group_summary
)


def test_climate_rank_is_NOT_geographical_distance_rank():
    original={
        "S1":{"A plant","B plant"},
        "S2":{"A plant","C plant"},
        "S3":{"D plant","C plant"},
    }
    # Climate model picks S1 as climate-rank1, but its distance is far.
    # Both A,C choose actual nearest among the regions that LIST them.
    z=source_species_min_km(
        [("S1",1500.),("S2",500.),("S3",700.)],original)
    assert z=={"A plant":500.,"B plant":1500.,
               "C plant":500.,"D plant":700.}
    assert "D plant" not in source_species_min_km(
        [("S1",1500.)],original)
    with pytest.raises(ValueError):
        source_species_min_km([("S1",-1.)],original)
    with pytest.raises(ValueError):
        source_species_min_km([("S1",float("nan"))],original)


def test_same_original_island_paired_source_mode_uses_NO_fake_evolution():
    z=[]
    for i in range(60):
        z.append(dict(
            block=f"B{i//2}",
            source_mode_climate_minus_geo5_original_genus_overlap_deficit=.12,
            source_mode_climate_minus_geo5_original_genus_concentration_residual=-.004,
            source_mode_original_gift_candidate_jaccard=.25,
            geo5_source_candidate_species_count=900,
            climate10_source_candidate_species_count=1400,
        ))
    s=paired_mode_group_summary(z,"block",n_boot=99,seed=3)
    assert s["status"]=="ORIGINAL_SAME_ISLAND_GEO_VS_CLIMATE_SOURCE_DEFINITION_SENSITIVITY"
    assert s["n_islands"]==60 and s["n_groups"]==30
    assert s["source_mode_climate_minus_geo5_original_genus_overlap_deficit_source_group_equal_mean"]==pytest.approx(.12)
    assert s["source_mode_climate_minus_geo5_original_genus_overlap_deficit_conditional_group_95"]==pytest.approx([.12,.12])
    assert s["source_mode_climate_minus_geo5_original_genus_concentration_residual_source_group_equal_mean"]==pytest.approx(-.004)
    assert s["mean_source_mode_union_candidate_name_jaccard"]==pytest.approx(.25)
    t=paired_mode_group_summary(z[:12],"block",n_boot=99)
    assert t["status"]=="INSUFFICIENT_SAME_ISLAND_SOURCE_PAIR_GROUPS"


def test_no_shared_source_mode_does_not_report_zero_source_effect():
    s=paired_mode_group_summary([],"block")
    assert s["status"]=="NO_BOTH_ORIGINAL_GIFT_SOURCE_MODES"
    assert s["n_islands"]==0
