"""Synthetic-only same original island/accepted plant NHM corpus pair controls."""
import numpy as np
import pandas as pd
import pytest
from fd.nhm_provider_matched_source import (
    PAIR,SCOPES,summarize_paired_nhm_providers
)


def frame(n=40):
    rows=[]
    for i in range(n):
        for scope in SCOPES:
            q=dict(
                island_id=f"original_test_{i}",region="tropical",
                scope=scope,spatial_block=f"block_{i//2}",
                first_source=f"source_{(i//2)%8}",
                n_target=8,n_target_genera=2,
                island_documentation_fraction=.5,
            )
            for provider in PAIR:
                for mode in ("geo5","climate10"):
                    q[f"{provider}__{mode}__n_numeric"]=14
                    q[f"{provider}__{mode}__source_coverage"]=.5
                    q[f"{provider}__{mode}__brier_genus"]=(
                        .08 if provider=="BHL_text" else .06
                    )
                    q[f"{provider}__{mode}__brier_pooled"]=.20
                    q[f"{provider}__{mode}__n_accepted_donor_species_shared_with_other_provider"]=2
            rows.append(q)
    return pd.DataFrame(rows)


def test_same_plant_BHL_china_source_scores_and_region_cohort():
    x=summarize_paired_nhm_providers(frame(),population=40,draws=299)
    assert x["n_original_region_provenance_panels"]==15
    assert x["n_original_private_island_scope_matched_provider_rows"]==120
    p=next(z for z in x["panels"] if z["region"]=="tropical"
           and z["scope"]=="all_recorded")
    assert p["n_common_islands_same_direct_plants_with_both_NHM_donor_corpora"]==40
    assert p["n_original_islands_balanced_across_all_four_NHM_sources"]==40
    z=p["all__spatial_block"]
    assert z["status"]=="EXPLORATORY_SAME_PLANT_TWO_NHM_CORPUS_SOURCE_PREDICTION"
    assert np.isclose(z["BHL_text__geo5__gain__equal_group_mean"],.12)
    assert np.isclose(z["FloraChina_text__geo5__gain__equal_group_mean"],.14)
    assert np.isclose(z["geo5__FloraChina_minus_BHL_gain__equal_group_mean"],.02)


def test_provider_heldout_missingness_does_not_become_no_trait_difference():
    q=frame()
    cut=q.island_id.eq("original_test_0")
    q.loc[cut,"n_target"]=0
    q.loc[cut,"n_target_genera"]=0
    q.loc[cut,"island_documentation_fraction"]=np.nan
    for provider in PAIR:
        for mode in ("geo5","climate10"):
            q.loc[cut,f"{provider}__{mode}__n_numeric"]=0
            q.loc[cut,f"{provider}__{mode}__source_coverage"]=np.nan
            q.loc[cut,f"{provider}__{mode}__brier_genus"]=np.nan
            q.loc[cut,f"{provider}__{mode}__brier_pooled"]=np.nan
            q.loc[cut,f"{provider}__{mode}__n_accepted_donor_species_shared_with_other_provider"]=0
    out=summarize_paired_nhm_providers(q,population=40,draws=299)
    p=next(x for x in out["panels"] if x["region"]=="tropical"
           and x["scope"]=="all_recorded")
    assert p["n_original_region_islands"]==40
    assert p["n_common_islands_same_direct_plants_with_both_NHM_donor_corpora"]==39


def test_imbalanced_original_documentation_is_coverage_fail_not_effect_zero():
    a=frame()
    a["BHL_text__climate10__source_coverage"]=.10
    p=summarize_paired_nhm_providers(a,population=40,draws=299)
    z=next(q for q in p["panels"] if q["region"]=="tropical"
           and q["scope"]=="strict_native")
    assert z["n_original_islands_balanced_across_all_four_NHM_sources"]==0
    assert z["documentation_balance10__spatial_block"][
        "status"]=="NOT_IDENTIFIABLE_MATCHED_PROVIDER_GEO_SUPPORT"


def test_reject_duplicate_original_plant_and_brier_source_incoherence():
    q=frame()
    with pytest.raises(ValueError,match="duplicated"):
        summarize_paired_nhm_providers(
            pd.concat([q,q.iloc[:1]]),population=40,draws=299
        )
    q=frame()
    q.loc[0,"FloraChina_text__geo5__brier_genus"]=np.nan
    with pytest.raises(ValueError,match="missing original"):
        summarize_paired_nhm_providers(q,population=40,draws=299)
