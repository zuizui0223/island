"""Synthetic tests for signed genus identity and within-genus state accounting."""
import numpy as np
import pandas as pd
import pytest

from fd.source_genus_exact_partition_modes import (
    REGIMES, SCOPES, MORPHOLOGY, TERMS,
    summarize_exact_genus_partition,
)


def sample(n=40):
    rows=[]
    for i in range(n):
        for scope in SCOPES:
            for sector in MORPHOLOGY:
                row=dict(
                    island_id=f"i{i}",regime="tropical",scope=scope,sector=sector,
                    spatial_block=f"b{i//2}",source_signature=f"s{(i//2)%7}",
                    n_scope_total=20,n_scope_direct=10,
                )
                for mode in ("geo5","climate10"):
                    row[f"{mode}__n_source_total"]=24
                    row[f"{mode}__n_source_direct"]=12
                    row[f"{mode}__n_shared_genera"]=2
                    row[f"{mode}__n_native_carriers_from_shared_genera"]=6
                    comps=dict(
                        shared_genus_composition=-.05 if mode=="geo5" else -.04,
                        shared_genus_within_species_states=.01,
                        native_exclusive_genus=.20,
                        source_exclusive_genus_negative=(
                            -.26 if mode=="geo5" else -.27
                        ),
                    )
                    row.update({mode+"__"+k:v for k,v in comps.items()})
                    row[mode+"__observed_difference"]=sum(comps.values())
                rows.append(row)
    return pd.DataFrame(rows)


def test_exact_four_terms_reconcile_and_cluster_ranges():
    obj=summarize_exact_genus_partition(sample(),population=40,bootstrap=299)
    assert obj["n_region_scope_sector_panels"]==100
    x=next(r for r in obj["panels"] if r["region"]=="tropical"
           and r["scope"]=="all_recorded" and r["structural_sector"]=="zygomorphic")
    assert x["n_same_original_islands_with_both_sources"]==40
    assert x["n_documentation_gap_both_le10pp"]==40
    z=x["all__spatial_block"]
    assert z["status"]=="EXACT_SIGNED_GENUS_ACCOUNTING_DESCRIPTIVE_NOT_CAUSAL"
    assert np.isclose(z["geo5__observed_difference__block_equal"],-.10)
    assert np.isclose(z["climate10__observed_difference__block_equal"],-.10)
    for mode in ("geo5","climate10"):
        assert np.isclose(sum(z[f"{mode}__{t}__block_equal"] for t in TERMS),
                          z[f"{mode}__observed_difference__block_equal"])


def test_invalid_partition_fails_algebraic_equality():
    df=sample()
    df.loc[0,"geo5__shared_genus_composition"]+=.1
    with pytest.raises(ValueError, match="exact reconciliation"):
        summarize_exact_genus_partition(df,population=40,bootstrap=299)


def test_missing_mode_stays_nonidentified_not_null_zero():
    df=sample()
    missing=df.island_id.eq("i0")
    for field in ("observed_difference",*TERMS):
        df.loc[missing,"climate10__"+field]=np.nan
    df.loc[missing,"climate10__n_source_total"]=0
    df.loc[missing,"climate10__n_source_direct"]=0
    obj=summarize_exact_genus_partition(df,population=40,bootstrap=299)
    x=next(r for r in obj["panels"] if r["region"]=="tropical"
           and r["scope"]=="all_recorded" and r["structural_sector"]=="zygomorphic")
    assert x["n_same_original_islands_with_both_sources"]==39


def test_no_partial_terms_or_duplicate_taxon_rows():
    df=sample()
    df.loc[0,"geo5__shared_genus_composition"]=np.nan
    with pytest.raises(ValueError,match="partial"):
        summarize_exact_genus_partition(df,population=40,bootstrap=299)
    df=sample()
    with pytest.raises(ValueError,match="duplicate"):
        summarize_exact_genus_partition(pd.concat([df,df.iloc[:1]]),
                                        population=40,bootstrap=299)


def test_no_fake_source_opportunity_when_documentation_is_imbalanced():
    df=sample()
    df["geo5__n_source_total"]=300
    df["climate10__n_source_total"]=300
    o=summarize_exact_genus_partition(df,population=40,bootstrap=299)
    z=next(r for r in o["panels"] if r["region"]=="tropical"
           and r["scope"]=="all_recorded" and r["structural_sector"]=="zygomorphic")
    assert z["n_documentation_gap_both_le10pp"]==0
    assert z["documentation_gap_le10pp__spatial_block"][
        "status"]=="NOT_IDENTIFIABLE_INDEPENDENT_CLUSTER_SUPPORT"
