"""Source-only, fixed-island GIFT geography/climate and missingness checks."""
import numpy as np
import pandas as pd
import pytest

from fd.source_mode_robustness import (
    MORPHOLOGY, SCOPES, same_island_source_mode_summary,
)


def sample(n=40):
    records=[]
    for i in range(n):
        for scope in SCOPES:
            for sector in MORPHOLOGY:
                records.append(dict(
                    island_id=f"test_island_{i}", regime="tropical",
                    scope=scope, sector=sector, spatial_block=f"b{i // 2}",
                    source_signature=f"r{(i // 2) % 6}",
                    n_island_total=16, n_island_direct=8,
                    island_sector_share=.3 + (i // 2) / 300,
                    n_geo5_total=24, n_geo5_direct=12,
                    geo5_sector_share=.20 + (i // 2) / 300,
                    n_climate10_total=24, n_climate10_direct=12,
                    climate10_sector_share=.25 + (i // 2) / 300,
                    climate10_source_assigned=True,
                ))
    return pd.DataFrame(records)


def test_paired_source_difference_exact_same_island_no_new_flower_response():
    obj=same_island_source_mode_summary(sample(), expected_islands=40, draws=299)
    assert obj["n_original_unique_botanical_islands"]==40
    assert obj["n_scope_sector_regional_panels"]==75
    x=next(q for q in obj["panels"] if q["region"]=="tropical"
           and q["scope"]=="all_recorded" and q["morphological_sector"]=="zygomorphic")
    assert x["n_same_islands_eligible_both_source_definitions"]==40
    assert x["n_same_islands_Direct_documentation_gap_both_le10pp"]==40
    y=x["same_island_paired_all"]
    assert y["status"]=="DESCRIPTIVE_SOURCE_DEFINITION_SENSITIVITY"
    assert np.isclose(y["delta_geo5_block_mean"], .10)
    assert np.isclose(y["delta_climate10_block_mean"], .05)
    assert np.isclose(y["paired_difference_climate_minus_geo5_block_mean"], -.05)
    assert np.isclose(y["same_island_mode_sign_concordance"], 1.0)


def test_native_and_all_are_separate_panels_not_mixed_independent_replicates():
    obj=same_island_source_mode_summary(sample(),expected_islands=40,draws=299)
    q=[v for v in obj["panels"] if v["region"]=="tropical" and
       v["morphological_sector"]=="zygomorphic"]
    assert {x["scope"] for x in q}==set(SCOPES)


def test_missing_climate_source_does_not_mean_zero_flower_proportion():
    df=sample()
    ix=df.island_id.eq("test_island_0")
    df.loc[ix,"n_climate10_total"]=0
    df.loc[ix,"n_climate10_direct"]=0
    df.loc[ix,"climate10_sector_share"]=np.nan
    df.loc[ix,"climate10_source_assigned"]=False
    x=same_island_source_mode_summary(df,expected_islands=40,draws=299)
    q=next(z for z in x["panels"] if z["scope"]=="all_recorded" and
           z["region"]=="tropical" and z["morphological_sector"]=="zygomorphic")
    assert q["n_same_islands_eligible_both_source_definitions"]==39
    assert x["n_original_unique_botanical_islands"]==40
    assert x["n_original_islands_with_climate10_source_assignment"]==39


def test_unbalanced_direct_documentation_does_not_remove_cohort_silently():
    df=sample()
    df["n_geo5_total"]=240  # only 5% Direct mainland coverage
    x=same_island_source_mode_summary(df,expected_islands=40,draws=299)
    q=next(z for z in x["panels"] if z["scope"]=="all_recorded" and
           z["region"]=="tropical" and z["morphological_sector"]=="zygomorphic")
    assert q["n_same_islands_eligible_both_source_definitions"]==40
    assert q["n_same_islands_Direct_documentation_gap_both_le10pp"]==0
    assert q["same_island_documentation_gap_both_le10pp"]["status"]=="NOT_IDENTIFIABLE_GEOGRAPHIC_SUPPORT"


def test_duplicate_and_missing_values_do_not_become_biological_absence():
    df=sample()
    with pytest.raises(ValueError,match="duplicate"):
        same_island_source_mode_summary(pd.concat([df,df.iloc[:1]]),draws=299)
    df=sample()
    df.loc[0,"geo5_sector_share"]=np.nan
    with pytest.raises(ValueError,match="missing trait share"):
        same_island_source_mode_summary(df,draws=299)


def test_source_carrier_zero_does_not_qualify():
    df=sample()
    ix=df.island_id.eq("test_island_1")
    df.loc[ix,"n_geo5_direct"]=0
    df.loc[ix,"geo5_sector_share"]=np.nan
    q=same_island_source_mode_summary(df,expected_islands=40,draws=299)
    z=next(
        x for x in q["panels"]
        if x["region"]=="tropical" and x["scope"]=="all_recorded"
        and x["morphological_sector"]=="zygomorphic"
    )
    assert z["n_same_islands_eligible_both_source_definitions"]==39
