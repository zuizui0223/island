"""Synthetic falsifications of real same-species floral genus-conditional pairing."""
import itertools
import numpy as np
import pandas as pd
import pytest

from fd.real_joint_within_genus_coupling import (
    PRIORITY,SCOPES, one_roster_analytic_genus_pair,
    summarize_original_genus_joint_coupling,
)


def test_within_genus_fixed_margins_is_analytically_exact():
    sp={f"Alpha flower{i}" for i in range(4)}
    names=sorted(sp)
    profile={
        "known":frozenset(sp),
        "A":frozenset(names[:2]),
        "B":frozenset(names[:2]),
        "joint":frozenset(names[:2]),
    }
    a=one_roster_analytic_genus_pair(sp,profile)
    assert a["n_direct"]==4
    assert a["n_genus_with_both_marginals_variable"]==1
    assert a["n_carriers_in_variable_genera"]==4
    assert np.isclose(a["coupling"],.25)
    assert np.isclose(a["genus_null_expected_coupling"],0)
    assert np.isclose(a["excess_coupling_beyond_genus_null"],.25)
    # Enumerate all C(4,2) assignments of B to original accepted species.
    co=[]
    for labels in itertools.combinations(names,2):
        co.append(len(set(labels)&set(names[:2]))/4-1/4)
    assert np.isclose(np.mean(co),a["genus_null_expected_coupling"])
    assert np.isclose(np.std(co,ddof=0),a["genus_null_conditional_sd"])


def test_between_genus_marginal_coupling_is_preserved_not_called_evolution():
    names=["Alpha flower1","Alpha flower2","Beta flower1","Beta flower2"]
    p={"known":frozenset(names),"A":frozenset(names[:2]),
       "B":frozenset(names[:2]),"joint":frozenset(names[:2])}
    out=one_roster_analytic_genus_pair(set(names),p)
    assert np.isclose(out["coupling"],.25)
    assert np.isclose(out["genus_null_expected_coupling"],.25)
    assert np.isclose(out["excess_coupling_beyond_genus_null"],0)
    assert out["n_genus_with_both_marginals_variable"]==0
    assert np.isclose(out["genus_null_conditional_sd"],0)


def test_only_one_direct_axis_is_unknown_not_false():
    x={"known":frozenset(["Alpha plant1"]),"A":frozenset(["Alpha plant1"]),
       "B":frozenset(),"joint":frozenset()}
    z=one_roster_analytic_genus_pair({"Beta absent","Alpha plant2"},x)
    assert z["n_direct"]==0
    assert np.isnan(z["coupling"]) and np.isnan(z["genus_null_expected_coupling"])


def _frame(n=40):
    out=[]
    for i in range(n):
        for scope in SCOPES:
            for sector in PRIORITY:
                row=dict(
                    island_id=f"unit{i}",regime="tropical",scope=scope,sector=sector,
                    spatial_block=f"b{i//2}",
                    source_signature=f"s{(i//2)%10}",
                )
                for g,parts in {
                    "island":(20,10,.2,.1,1,5),
                    "geo5":(40,20,.1,.08,1,8),
                    "climate10":(40,20,.15,.1,1,8),
                }.items():
                    total,ncar,C,expected,ng,inf=parts
                    row.update({
                        f"{g}_total":total,
                        f"{g}_n_direct":ncar,
                        f"{g}_coupling":C,
                        f"{g}_genus_null_expected_coupling":expected,
                        f"{g}_excess_coupling_beyond_genus_null":C-expected,
                        f"{g}_genus_null_conditional_sd":.03,
                        f"{g}_n_genus_with_both_marginals_variable":ng,
                        f"{g}_n_carriers_in_variable_genera":inf,
                    })
                out.append(row)
    return pd.DataFrame(out)


def test_genus_conditional_source_relative_coupling_adds_exactly():
    report=summarize_original_genus_joint_coupling(
        _frame(),population=40,boots=299
    )
    assert report["n_region_scope_structural_pair_panels"]==45
    r=next(x for x in report["panels"] if x["region"]=="tropical"
           and x["scope"]=="all_recorded" and
           x["genuine_joint_sector"]=="actinomorphic_deep")
    assert r["n_both_real_source_mode_qualified_islands"]==40
    assert r["n_both_informative_and_documentation_balanced"]==40
    assert r["n_islands_with_informative_genus_on_island"]==40
    assert r["n_islands_with_informative_genus_in_geo5_source"]==40
    assert r["n_islands_with_informative_genus_in_climate10_source"]==40
    q=r["both_balanced_and_informative__spatial_block"]
    assert q["status"]=="DESCRIPTIVE_GENUS_CONDITIONAL_COUPLING_NOT_CAUSAL"
    assert np.isclose(q["geo5__source_relative_original_coupling_group_equal_mean"],.1)
    assert np.isclose(q["geo5__source_relative_genus_expected_coupling_group_equal_mean"],.02)
    assert np.isclose(q["geo5__source_relative_beyond_genus_null_group_equal_mean"],.08)


def test_missing_source_cannot_be_turned_into_real_zero():
    x=_frame()
    z=x.island_id.eq("unit0")
    x.loc[z,"climate10_n_direct"]=0
    x.loc[z,"climate10_total"]=0
    for a in ("coupling","genus_null_expected_coupling",
              "excess_coupling_beyond_genus_null"):
        x.loc[z,"climate10_"+a]=np.nan
    obj=summarize_original_genus_joint_coupling(x,population=40,boots=299)
    r=next(z for z in obj["panels"] if z["region"]=="tropical" and
           z["scope"]=="all_recorded" and z["genuine_joint_sector"]=="actinomorphic_deep")
    assert r["n_both_real_source_mode_qualified_islands"]==39
    assert obj["n_original_GIFT_botanical_islands"]==40


def test_no_variable_genus_is_not_accepted_as_informative_support():
    x=_frame()
    x["island_n_genus_with_both_marginals_variable"]=0
    y=summarize_original_genus_joint_coupling(x,population=40,boots=299)
    z=next(q for q in y["panels"] if q["region"]=="tropical"
           and q["scope"]=="strict_native" and q["genuine_joint_sector"]=="zygomorphic_deep")
    assert z["n_with_both_source_and_island_informative_genera"]==0
    assert z["informative_genera_both_sources__spatial_block"][
        "status"]=="NOT_IDENTIFIABLE_GEO_OR_SOURCE_REPLICATION"


def test_fails_closed_for_bad_species_or_reference_accounting():
    with pytest.raises(ValueError,match="duplicate"):
        x=_frame()
        summarize_original_genus_joint_coupling(
            pd.concat([x,x.iloc[:1]]),population=40,boots=299
        )
    with pytest.raises(AssertionError,match="per-island"):
        x=_frame()
        x.loc[0,"island_excess_coupling_beyond_genus_null"]=.4
        summarize_original_genus_joint_coupling(x,population=40,boots=299)
