import numpy as np
import pandas as pd
import pytest

from fd.source_loo_flower_placebo import (
    heldout_source_placebo, summarize_source_loo_placebo,
)


def simple():
    targets=["Alpha a","Alpha b","Beta a","Beta b"]
    donor={"Alpha":{"Alpha c","Alpha d"},"Beta":{"Beta c","Beta d"}}
    known=set(targets)|set().union(*donor.values())
    rad={"Alpha a","Alpha b","Alpha c","Alpha d"}
    tub={"Beta a","Beta b","Beta c","Beta d"}
    profiles={
        "radial":{"known":known,"joint":rad},
        "tubular":{"known":known,"joint":tub},
    }
    return dict(targets=targets,donors_by_genus=donor,
                radial=profiles["radial"],tubular=profiles["tubular"])


def test_source_loo_matching_controls_no_island_excess_when_identical_genus_states():
    z=heldout_source_placebo(**simple())
    for name in ("radial","tubular"):
        assert z[f"island_gain_{name}"]==pytest.approx(4/9)
        assert z[f"source_loo_gain_{name}"]==pytest.approx(4/9)
    assert z["island_excess_over_source_loo"]==pytest.approx(0)


def test_perturbed_island_species_changes_excess_but_not_source_control():
    x=simple()
    before=heldout_source_placebo(**x)
    x["radial"]=dict(x["radial"],joint=set(x["radial"]["joint"])|{"Beta a","Beta b"})
    after=heldout_source_placebo(**x)
    assert after["source_loo_gain_radial"]==pytest.approx(before["source_loo_gain_radial"])
    assert after["island_excess_over_source_loo"]!=pytest.approx(0)


def test_one_source_donor_is_not_allowable_loo():
    x=simple()
    x["donors_by_genus"]["Beta"]={"Beta c"}
    with pytest.raises(ValueError,match=">=2"):
        heldout_source_placebo(**x)


def test_island_species_must_not_enter_source_donors():
    x=simple()
    x["donors_by_genus"]["Beta"].add("Beta a")
    with pytest.raises(ValueError,match="leaked"):
        heldout_source_placebo(**x)


def test_independent_three_trait_coverage_required():
    x=simple()
    x["tubular"]["known"].remove("Beta a")
    with pytest.raises(ValueError,match="three Direct"):
        heldout_source_placebo(**x)


def frame(n=40,eligible=True):
    rows=[]
    for i in range(n):
        for scope in ("all_recorded","strict_native","strict_native_nonendemic"):
            for gate in ("any_congener","at_least_two_congeners"):
                for stratum in ("triple_direct","triple_direct_deep_only"):
                    row=dict(island_id=f"I{i}",regime="tropical",
                        spatial_block=f"B{i//2}",source_signature=f"S{i//2}",
                        scope=scope,donor_gate=gate,stratum=stratum,
                        n_target=6 if eligible else 0,
                        n_target_genera=2 if eligible else 0,
                        island_direct_coverage=.5 if eligible else np.nan)
                    for mode in ("geo5","climate10"):
                        row[f"{mode}__n_donor"]=12 if eligible else 0
                        row[f"{mode}__donor_direct_coverage"]=.5 if eligible else np.nan
                        for fld,val in (("island_gain_radial",.2),
                                        ("island_gain_tubular",-.1),
                                        ("source_loo_gain_radial",.12),
                                        ("source_loo_gain_tubular",-.08),
                                        ("island_gain_difference",.3),
                                        ("source_loo_gain_difference",.2),
                                        ("island_excess_over_source_loo",.1)):
                            row[f"{mode}__{fld}"]=val if (eligible and gate=="at_least_two_congeners") else np.nan
                    rows.append(row)
    return pd.DataFrame(rows)


def test_global_population_and_source_placebo_excess():
    out=summarize_source_loo_placebo(frame(),population=40,draws=199)
    assert out["n_reporting_panels"]==5*3*2*2*2
    z=[r for r in out["panels"] if r["region"]=="tropical" and
       r["scope"]=="all_recorded" and r["stratum"]=="triple_direct"
       and r["source_mode"]=="geo5" and r["support"]=="all"][0]
    assert z["n_supported_islands"]==40
    assert z["status"]=="DESCRIPTIVE_REAL_ISLAND_VS_SOURCE_LOO"
    assert z["block__island_gain_difference"]==pytest.approx(.3)
    assert z["block__source_loo_gain_difference"]==pytest.approx(.2)
    assert z["block__island_excess_over_source_loo"]==pytest.approx(.1)


def test_unsupported_source_placebo_is_not_zero():
    out=summarize_source_loo_placebo(frame(eligible=False),population=40,draws=199)
    assert all(r["status"]=="NOT_IDENTIFIABLE_SOURCE_LOO_SUPPORT" for r in out["panels"])
    assert all("block__island_excess_over_source_loo" not in r for r in out["panels"])
