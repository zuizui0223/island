import numpy as np
import pandas as pd
import pytest

from fd.triple_matched_congener import (
    aggregate_triple, triple_matched_island,
    GATES, SCOPES, STRATA, MODES,
)


def fixture():
    target = {"Alpha a", "Alpha b", "Beta a", "Beta b", "Gamma missing_form"}
    donors = {"Alpha c", "Alpha d", "Beta c", "Beta d"}
    allnames = target | donors
    deep = set(allnames)
    radial_pos = {"Alpha a", "Alpha b", "Alpha c"}
    tubular_pos = {"Alpha a", "Alpha b", "Beta a", "Beta b",
                   "Alpha c", "Alpha d", "Beta d"}
    profiles = {
        "actinomorphic_deep": dict(known=set(allnames), A=radial_pos,
                                   B=deep, joint=radial_pos & deep),
        "tubular_bell_deep": dict(known=set(allnames)-
                                {"Gamma missing_form"}, A=tubular_pos,
                                B=deep, joint=tubular_pos & deep),
    }
    scopes = {k:target for k in SCOPES}
    return dict(full_island=target, scope_rosters=scopes,
                geo5_source=donors, climate10_source=donors,
                profiles=profiles)


def test_same_triply_documented_target_and_identical_source_donors():
    z = triple_matched_island(**fixture())
    assert len(z) == 12
    t = next(v for v in z if v["scope"]=="all_recorded"
             and v["donor_gate"]=="any_congener"
             and v["stratum"]=="triple_direct")
    assert t["n_target"] == 4   # not five: Gamma lacks form
    assert t["n_target_genera"] == 2
    assert t["geo5__n_donor"] == 4
    assert t["geo5__gain_radial"] == pytest.approx(.1875)
    assert t["geo5__gain_tubular"] == pytest.approx(-.0625)
    assert t["geo5__gain_delta"] == pytest.approx(.25)
    assert t["climate10__gain_delta"] == pytest.approx(.25)


def test_true_deep_only_is_a_separate_stratum():
    d = fixture()
    for profile in d["profiles"].values():
        profile["B"] = set(profile["B"]) - {"Beta b", "Beta d"}
        profile["joint"] = profile["A"] & profile["B"]
    rows = triple_matched_island(**d)
    q = next(x for x in rows if x["scope"]=="all_recorded"
             and x["stratum"]=="triple_direct_deep_only"
             and x["donor_gate"]=="any_congener")
    assert q["n_target"] == 3
    assert q["geo5__n_donor"] == 3


def test_source_leakage_rejected():
    d=fixture()
    d["geo5_source"]=set(d["geo5_source"])|{"Alpha a"}
    with pytest.raises(ValueError,match="leaked"):
        triple_matched_island(**d)


def test_inconsistent_depth_sources_rejected():
    d=fixture()
    d["profiles"]["tubular_bell_deep"]["B"]=(
        set(d["profiles"]["tubular_bell_deep"]["B"])-{"Alpha a"})
    with pytest.raises(AssertionError,match="disagree"):
        triple_matched_island(**d)


def _agg_frame(n=40,eligible=True):
    rows=[]
    for i in range(n):
        for scope in SCOPES:
            for gate in GATES:
                for stratum in STRATA:
                    row=dict(island_id=f"I{i}",regime="tropical",scope=scope,
                             donor_gate=gate,stratum=stratum,
                             spatial_block=f"B{i//2}",source_signature=f"S{i//2}",
                             n_target=6 if eligible else 0,
                             n_target_genera=2 if eligible else 0,
                             island_direct_coverage=.5 if eligible else np.nan)
                    for mode in MODES:
                        row.update({
                            f"{mode}__n_donor":12 if eligible else 0,
                            f"{mode}__donor_direct_coverage":.5 if eligible else np.nan,
                            f"{mode}__gain_radial":.2 if eligible else np.nan,
                            f"{mode}__gain_tubular":-.1 if eligible else np.nan,
                            f"{mode}__gain_delta":.3 if eligible else np.nan,
                            f"{mode}__relative_gain_delta":.1 if eligible else np.nan,
                            f"{mode}__prevalence_radial":.5 if eligible else np.nan,
                            f"{mode}__prevalence_tubular":.6 if eligible else np.nan,
                            f"{mode}__brier_pooled_radial":.4 if eligible else np.nan,
                            f"{mode}__brier_pooled_tubular":.2 if eligible else np.nan,
                        })
                    rows.append(row)
    return pd.DataFrame(rows)


def test_aggregate_preserves_original_population_and_common_source_support():
    z=aggregate_triple(_agg_frame(),population=40,draws=199)
    assert z["n_original_botanical_islands"]==40
    assert z["n_original_scope_gate_stratum_rows"]==480
    assert z["n_panels"]==5*3*2*2*2*2
    eligible=[x for x in z["panels"] if x["region"]=="tropical" and
              x["scope"]=="all_recorded" and
              x["donor_gate"]=="any_congener" and
              x["stratum"]=="triple_direct" and
              x["documentation"]=="all" and x["mode"]=="geo5"][0]
    assert eligible["n_supported"]==40
    assert eligible["status"]=="DESCRIPTIVE_TRIPLE_DIRECT_MATCHED_SOURCE"
    assert eligible["block__gain_delta"]==pytest.approx(.3)
    assert eligible["block__normalized_status"]=="SAME_COHORT_BASELINE_RISK_NORMALIZED"
    assert eligible["block__skill_radial"]==pytest.approx(.5)
    assert eligible["block__skill_tubular"]==pytest.approx(-.5)
    assert eligible["block__skill_delta"]==pytest.approx(1.0)


def test_unsupported_never_imputed_as_no_effect():
    z=aggregate_triple(_agg_frame(eligible=False),population=40,draws=199)
    assert all(x["status"]=="NOT_IDENTIFIABLE_TRIPLE_DIRECT_SUPPORT"
               for x in z["panels"])
    assert all("block__gain_delta" not in x for x in z["panels"])


def test_source_loo_negative_control_is_scored_only_when_donors_sufficient():
    d=fixture()
    new_donors={"Alpha e","Alpha f","Alpha g",
                "Beta e","Beta f","Beta g"}
    d["geo5_source"]=set(d["geo5_source"]) | new_donors
    d["climate10_source"]=set(d["climate10_source"]) | new_donors
    for prof in d["profiles"].values():
        prof["known"] = set(prof["known"]) | new_donors
        prof["B"] = set(prof["B"]) | new_donors
        prof["A"] = set(prof["A"]) | {"Alpha e","Alpha f","Beta e"}
        prof["joint"] = set(prof["A"]) & set(prof["B"])
    z=triple_matched_island(**d)
    q=next(x for x in z if x["scope"]=="all_recorded"
           and x["donor_gate"]=="at_least_two_congeners"
           and x["stratum"]=="triple_direct")
    assert q["geo5__n_donor"]==10
    assert np.isfinite(q["geo5__island_excess_over_source_loo"])
    assert np.isfinite(q["climate10__source_loo_gain_difference"])
    inadequate=triple_matched_island(**fixture())
    q0=next(x for x in inadequate if x["scope"]=="all_recorded"
            and x["donor_gate"]=="at_least_two_congeners"
            and x["stratum"]=="triple_direct")
    assert np.isnan(q0["geo5__island_excess_over_source_loo"])
