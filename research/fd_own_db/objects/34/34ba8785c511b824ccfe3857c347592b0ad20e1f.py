"""Pure synthetic true same-species pair genus/source accounting."""
import numpy as np
import pandas as pd
import pytest

from fd.real_joint_genus_accounting import (
    TARGETS,SCOPES,PARTS,summarize_real_joint_genus
)


def frame(n=40):
    rows=[]
    for i in range(n):
        for scope in SCOPES:
            for sector in TARGETS:
                q=dict(
                    island_id=f"synthetic{i}",region="tropical",scope=scope,sector=sector,
                    spatial_block=f"block{i//2}",first_source=f"source{(i//2)%6}",
                    n_island_total=20,n_island_direct=10,
                )
                for mode in ("geo5","climate10"):
                    pieces=dict(
                        shared_genus_composition=.10,
                        shared_genus_within_species_states=.04,
                        native_exclusive_genus=.20,
                        source_exclusive_genus_negative=-.24
                    )
                    q[f"{mode}__n_mainland_total"]=30
                    q[f"{mode}__n_mainland_direct"]=15
                    q[f"{mode}__n_shared_genera"]=3
                    q[f"{mode}__observed_difference"]=sum(pieces.values())
                    for p,v in pieces.items():q[f"{mode}__{p}"]=v
                rows.append(q)
    return pd.DataFrame(rows)


def test_correct_genus_partition_sum_and_full_original_population_gate():
    z=summarize_real_joint_genus(frame(),population=40,draws=299)
    assert z["n_regional_scope_true_pair_panels"]==5*3*3
    assert z["n_original_island_scope_true_pair_rows"]==40*3*3
    t=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert t["n_same_original_islands_both_mainland_sources"]==40
    assert t["n_original_islands_documentation_gap_both_le10pp"]==40
    a=t["all__spatial_block"]
    assert a["status"]=="EXACT_REAL_PAIRED_FLORA_GENUS_ACCOUNTING_NOT_CAUSAL"
    assert np.isclose(a["geo5__observed_difference__equal_group_mean"],.10)
    assert np.isclose(sum(a[f"geo5__{p}__equal_group_mean"] for p in PARTS),.10)


def test_source_missingness_not_biological_zero():
    x=frame()
    choose=x.island_id.eq("synthetic0")
    for key in ("observed_difference",*PARTS):
        x.loc[choose,"climate10__"+key]=np.nan
    x.loc[choose,"climate10__n_mainland_direct"]=0
    x.loc[choose,"climate10__n_mainland_total"]=0
    z=summarize_real_joint_genus(x,population=40,draws=299)
    q=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert q["n_same_original_islands_both_mainland_sources"]==39


def test_invalid_signed_exact_joint_and_original_scope_fails():
    x=frame()
    x.loc[0,"geo5__shared_genus_within_species_states"]+=.15
    with pytest.raises(AssertionError, match="do not sum"):
        summarize_real_joint_genus(x,population=40,draws=299)
    x=frame()
    with pytest.raises(ValueError,match="duplicate"):
        summarize_real_joint_genus(pd.concat([x,x.iloc[:1]]),population=40,draws=299)
    with pytest.raises(ValueError,match="full original"):
        summarize_real_joint_genus(x,population=41,draws=299)


def test_same_island_matched_documentation_filter_can_make_unidentifiable():
    x=frame()
    x["geo5__n_mainland_total"]=500
    x["climate10__n_mainland_total"]=500
    z=summarize_real_joint_genus(x,population=40,draws=299)
    q=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded"
           and p["sector"]=="actinomorphic_deep")
    assert q["n_same_original_islands_both_mainland_sources"]==40
    assert q["n_original_islands_documentation_gap_both_le10pp"]==0
    assert q["doc_gap_le10pp__spatial_block"]["status"]=="NOT_IDENTIFIABLE_BLOCK_OR_SOURCE_SUPPORT"
