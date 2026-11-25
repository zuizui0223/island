import itertools
from collections import Counter

import numpy as np
import pandas as pd
import pytest

from fd.genus_flower_pair_covariance import (
    partition_genus_coupling,partition_island_source,
    aggregate_pair_genus_partition,
)
from fd.real_joint_syndrome_access import JOINT_SPECS


def profile(known,aa,bb):
    return dict(known=frozenset(known),A=frozenset(aa),B=frozenset(bb),
                joint=frozenset(set(aa)&set(bb)))


def test_exact_within_genus_pairing_and_finite_hypergeometric_variance():
    names={"Alpha a","Alpha b","Alpha c","Alpha d"}
    q=profile(names,{"Alpha a","Alpha b"},{"Alpha a","Alpha b"})
    z=partition_genus_coupling(names,q)
    assert z["n_genera"]==1
    assert z["n_variable_genera"]==1
    assert z["n_plants_in_variable_genera"]==4
    assert z["joint"]==pytest.approx(.5)
    assert z["genus_independent_joint"]==pytest.approx(.25)
    assert z["within_genus_coupling"]==pytest.approx(.25)
    assert z["between_genus_coupling"]==pytest.approx(0)
    assert z["coupling"]==pytest.approx(.25)
    # Enumerate all 6 distinct ways to pair B state on 4 accepted plants.
    vals=[]
    for b in itertools.combinations(sorted(names),2):
        vals.append(len(set(b)&set(q["A"]))/4)
    assert np.mean(vals)==pytest.approx(z["genus_independent_joint"])
    assert np.var(vals)==pytest.approx(z["finite_within_genus_joint_variance"])


def test_true_between_genus_covariance_without_any_within_genus_trait_variation():
    names={"Alpha a","Alpha b","Beta a","Beta b"}
    selected={"Alpha a","Alpha b"}
    q=profile(names,selected,selected)
    z=partition_genus_coupling(names,q)
    assert z["n_variable_genera"]==0
    assert z["coupling"]==pytest.approx(.25)
    assert z["within_genus_coupling"]==pytest.approx(0)
    assert z["between_genus_coupling"]==pytest.approx(.25)
    assert z["finite_within_genus_joint_variance"]==0
    # Raw pairing is positive although genus-fixed permutation CANNOT
    # produce any within-genus alternative. This is compositional only.


def test_historical_not_directly_observed_flower_axis_is_not_negative():
    names={"Alpha a","Alpha b"}
    q=profile({"Alpha a"},{"Alpha a"},{"Alpha a"})
    z=partition_genus_coupling(names,q)
    assert z["n_direct"]==1
    assert z["joint"]==1
    assert z["coupling"]==0
    empty=partition_genus_coupling({"Beta a"},q)
    assert empty["n_direct"]==0 and np.isnan(empty["joint"])
    assert np.isnan(empty["within_genus_coupling"])


def test_source_relative_coupling_exactly_equals_between_plus_within():
    island={"Alpha a","Alpha b","Beta a","Beta b"}
    donors={"Alpha c","Alpha d","Beta c","Beta d"}
    q=profile(island|donors,
              {"Alpha a","Alpha b","Alpha c","Beta c"},
              {"Alpha a","Alpha b","Alpha d","Beta c"})
    pairs={name:q for name in JOINT_SPECS}
    result=partition_island_source(
        full_flora=island,
        scope_rosters=dict(all_recorded=island,strict_native=island,
                           strict_native_nonendemic=island),
        source_rosters=dict(geo5=donors,climate10=donors),
        pairs=pairs)
    assert len(result)==3*len(JOINT_SPECS)
    one=result[0]
    assert one["island__n_direct"]==4
    assert one["geo5__n_direct"]==4
    for mode in ("geo5","climate10"):
        assert one[f"{mode}__delta_coupling"]==pytest.approx(
            one[f"{mode}__delta_between_genus_coupling"]+
            one[f"{mode}__delta_within_genus_coupling"])


def test_original_source_plant_leakage_is_rejected():
    island={"Alpha a"}
    donors={"Alpha a","Beta b"}
    q=profile(island|donors,{"Alpha a"},{"Alpha a"})
    with pytest.raises(ValueError,match="leaked"):
        partition_island_source(
            full_flora=island,scope_rosters={k:island for k in (
                "all_recorded","strict_native","strict_native_nonendemic")},
            source_rosters={"geo5":donors,"climate10":donors},
            pairs={x:q for x in JOINT_SPECS})


def fake_frame(n_islands=40,genus_var=True):
    rows=[]
    for i in range(n_islands):
        for scope in ("all_recorded","strict_native","strict_native_nonendemic"):
            for pair in JOINT_SPECS:
                q=dict(island_id=f"I{i}",regime="tropical",
                    spatial_block=f"B{i//2}",first_source=f"S{i//2}",
                    scope=scope,sector=pair,island_total=10,
                    island__n_direct=6,island__n_genera=3,
                    island__n_variable_genera=1 if genus_var else 0,
                    n_shared_informative_genera_both_source_modes=1 if genus_var else 0)
                for mode in ("geo5","climate10"):
                    q[f"{mode}__total"]=20
                    q[f"{mode}__n_direct"]=12
                    q[f"{mode}__n_variable_genera"]=1 if genus_var else 0
                    q[f"{mode}__delta_joint"]=.09
                    q[f"{mode}__delta_coupling"]=.07
                    q[f"{mode}__delta_between_genus_coupling"]=.04
                    q[f"{mode}__delta_within_genus_coupling"]=.03
                rows.append(q)
    return pd.DataFrame(rows)


def test_full_island_denominator_and_positive_genus_partition():
    z=aggregate_pair_genus_partition(fake_frame(),population=40,draws=199)
    assert z["original_botanical_island_population"]==40
    assert z["n_island_scope_pair_records"]==40*3*11
    assert z["n_region_scope_pair_panels"]==5*3*11
    p=next(x for x in z["panels"] if x["region"]=="tropical" and
           x["scope"]=="all_recorded" and x["sector"]=="actinomorphic_deep")
    assert p["n_both_mainland_source_supported"]==40
    assert p["n_islands_with_true_variation_on_all_three_rosters"]==40
    assert p["n_islands_with_same_variable_genus_on_island_and_both_sources"]==40
    assert p["same_variable_genus_across_island_and_sources__spatial_block"]["status"]==(
        "DESCRIPTIVE_EXACT_GENUS_WITHIN_BETWEEN_PAIR_PARTITION")
    for cohort in ("all","genuine_within_genus_var"):
        v=p[cohort+"__spatial_block"]
        assert v["status"]=="DESCRIPTIVE_EXACT_GENUS_WITHIN_BETWEEN_PAIR_PARTITION"
        assert v["geo5__delta_coupling__equal_group_mean"]==pytest.approx(.07)
        assert v["geo5__delta_between_genus_coupling__equal_group_mean"]==pytest.approx(.04)
        assert v["geo5__delta_within_genus_coupling__equal_group_mean"]==pytest.approx(.03)


def test_lacking_intrageneric_morphology_is_NOT_identifiable_not_a_true_zero():
    z=aggregate_pair_genus_partition(fake_frame(genus_var=False),
                                     population=40,draws=199)
    p=next(x for x in z["panels"] if x["region"]=="tropical" and
           x["scope"]=="all_recorded" and x["sector"]=="actinomorphic_deep")
    assert p["all__spatial_block"]["status"]=="DESCRIPTIVE_EXACT_GENUS_WITHIN_BETWEEN_PAIR_PARTITION"
    assert p["genuine_within_genus_var__spatial_block"]["status"]==(
        "NOT_IDENTIFIABLE_TRUE_WITHIN_GENUS_SOURCE_SUPPORT")
    assert "geo5__delta_within_genus_coupling__equal_group_mean" not in (
        p["genuine_within_genus_var__spatial_block"])
    assert p["n_islands_with_true_variation_on_all_three_rosters"]==0
    assert p["same_variable_genus_across_island_and_sources__spatial_block"]["status"]==(
        "NOT_IDENTIFIABLE_TRUE_WITHIN_GENUS_SOURCE_SUPPORT")


def test_duplicate_floras_rejected_not_pseudoreplicated():
    x=fake_frame()
    x=pd.concat([x,x.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="complete|duplicate|population"):
        aggregate_pair_genus_partition(x,population=40,draws=199)



def test_variable_genus_on_both_sides_is_not_enough_if_genera_differ():
    island={"Alpha a","Alpha b","Alpha c","Alpha d",
            "Beta a","Beta b"}
    donors={"Alpha e","Alpha f","Beta c","Beta d","Beta e","Beta f"}
    # Original island: Alpha genus has both traits varying.
    # Mainland: Beta genus has both traits varying, Alpha is fixed.
    A={"Alpha a","Alpha b","Beta a","Beta b",
       "Alpha e","Alpha f","Beta c","Beta d"}
    B={"Alpha a","Alpha c","Beta a","Beta b",
       "Alpha e","Alpha f","Beta c","Beta e"}
    prof=profile(island|donors,A,B)
    out=partition_island_source(
        full_flora=island,
        scope_rosters={k:island for k in (
            "all_recorded","strict_native","strict_native_nonendemic")},
        source_rosters={"geo5":donors,"climate10":donors},
        pairs={name:prof for name in JOINT_SPECS})
    q=out[0]
    assert q["island__n_variable_genera"]==1
    assert q["geo5__n_variable_genera"]==1
    assert q["climate10__n_variable_genera"]==1
    assert q["n_shared_informative_genera_both_source_modes"]==0


def test_strict_same_variable_genus_all_three_source_rosters():
    island={"Alpha a","Alpha b","Alpha c","Alpha d"}
    donors={"Alpha e","Alpha f","Alpha g","Alpha h"}
    A={"Alpha a","Alpha b","Alpha e","Alpha f"}
    B={"Alpha a","Alpha c","Alpha e","Alpha g"}
    prof=profile(island|donors,A,B)
    out=partition_island_source(
        full_flora=island,
        scope_rosters={k:island for k in (
            "all_recorded","strict_native","strict_native_nonendemic")},
        source_rosters={"geo5":donors,"climate10":donors},
        pairs={name:prof for name in JOINT_SPECS})
    q=out[0]
    assert q["n_shared_informative_genera_both_source_modes"]==1
