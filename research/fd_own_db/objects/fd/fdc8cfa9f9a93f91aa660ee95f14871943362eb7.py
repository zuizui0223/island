"""Synthetic-only NHM original five-source-corpus heldout calibration audits."""
import numpy as np
import pandas as pd
import pytest
from fd.nhm_numeric_provider_holdout import (
    PROVIDER_NAMES,SCOPES,summarize_provider_disjoint
)


def sample(n=40):
    data=[]
    for i in range(n):
        for scope in SCOPES:
            for provider in PROVIDER_NAMES:
                available=provider in ("FloraChina_text","BHL_text")
                row=dict(
                    island_id=f"island_{i}",region="tropical",scope=scope,
                    provider=provider,spatial_block=f"block_{i//2}",
                    first_source=f"first_{(i//2)%8}",
                    calibration_supported=available,
                    n_target=(6 if available else 0),
                    n_target_genera=(2 if available else 0),
                    island_documentation_fraction=(.5 if available else np.nan),
                )
                for mode in ("geo5","climate10"):
                    row[mode+"__n_numeric_donors"]=12 if available else 0
                    row[mode+"__source_coverage"]=.5 if available else np.nan
                    row[mode+"__brier_genus"]=.10 if available else np.nan
                    row[mode+"__brier_pooled"]=.25 if available else np.nan
                data.append(row)
    return pd.DataFrame(data)


def test_real_provider_disjoint_numeric_support_and_five_corpus_panel():
    x=summarize_provider_disjoint(sample(),population=40,draws=299)
    assert x["n_original_island_scope_provider_rows"]==40*3*5
    assert x["n_regional_status_source_provider_panels"]==5*3*5
    a=next(q for q in x["panels"] if
           q["region"]=="tropical" and q["scope"]=="all_recorded" and
           q["source_numeric_heldout_provider"]=="FloraChina_text")
    assert a["n_qualified_source_provider_disjoint_islands"]==40
    assert a["n_direct_NHM_provider_coverage_balanced_both_le10pp"]==40
    z=a["all__spatial_block"]
    assert z["status"]=="EXPLORATORY_ORIGINAL_PROVIDER_DISJOINT_NUMERIC_SOURCE_PREDICTION"
    assert np.isclose(z["geo5__gain__block_equal"],.15)
    assert np.isclose(z["climate10__gain__block_equal"],.15)


def test_failed_source_provider_support_does_not_become_negative_result():
    x=summarize_provider_disjoint(sample(),population=40,draws=299)
    a=next(q for q in x["panels"] if
           q["region"]=="tropical" and q["scope"]=="all_recorded" and
           q["source_numeric_heldout_provider"]=="FloraPakistan_text")
    assert a["n_original_islands_without_numeric_model"]==40
    assert a["n_qualified_source_provider_disjoint_islands"]==0
    assert a["all__spatial_block"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PROVIDER_SUPPORT"


def test_documentation_gap_screen_uses_both_original_source_modes():
    v=sample()
    v.loc[v.provider.eq("FloraChina_text"),"geo5__source_coverage"]=.1
    x=summarize_provider_disjoint(v,population=40,draws=299)
    a=next(q for q in x["panels"] if
           q["region"]=="tropical" and q["scope"]=="strict_native" and
           q["source_numeric_heldout_provider"]=="FloraChina_text")
    assert a["n_qualified_source_provider_disjoint_islands"]==40
    assert a["n_direct_NHM_provider_coverage_balanced_both_le10pp"]==0
    assert a["documentation_balance10__spatial_block"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PROVIDER_SUPPORT"


def test_no_duplicate_islands_or_false_predicted_trait_values():
    v=sample()
    with pytest.raises(ValueError,match="duplicate"):
        summarize_provider_disjoint(pd.concat([v,v.iloc[:1]]),
                                    population=40,draws=299)
    v=sample()
    v.loc[v.provider.eq("FloraPakistan_text"),"n_target"]=6
    with pytest.raises(ValueError,match="failed NHM disjoint calibration"):
        summarize_provider_disjoint(v,population=40,draws=299)
    v=sample()
    v.loc[0,"geo5__brier_genus"]=np.nan
    with pytest.raises(ValueError,match="undefined is not a biological zero"):
        summarize_provider_disjoint(v,population=40,draws=299)


def test_nonrandom_native_and_no_full_population_contraction():
    v=sample()
    with pytest.raises(ValueError,match="original full"):
        summarize_provider_disjoint(v,population=41,draws=299)
    assert set(v.scope)==set(SCOPES)
