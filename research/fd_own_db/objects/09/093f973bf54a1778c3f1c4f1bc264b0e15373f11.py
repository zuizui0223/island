import numpy as np
import pandas as pd
import pytest

from fd.syndrome_ensemble_source import (
    MOTIFS,SCOPES,make_motifs,make_island_motif_rows,
    aggregate_motifs,fourpart,
)


def states():
    return {
      "flower_primary_color":{
        "Alpha a":"red_pink","Alpha b":"white",
        "Beta a":"white","Beta b":"white",
        "Gamma onlycolor":"red_pink"},
      "floral_symmetry":{
        "Alpha a":"actinomorphic","Alpha b":"actinomorphic",
        "Beta a":"actinomorphic","Beta b":"zygomorphic",
        "Gamma onlycolor":"actinomorphic"},
      "floral_form":{
        "Alpha a":"tubular","Alpha b":"open_radial",
        "Beta a":"open_radial","Beta b":"bell_campanulate"},
      "tube_depth_class":{
        "Alpha a":"deep","Alpha b":"shallow",
        "Beta a":"shallow","Beta b":"deep",
        "Gamma onlycolor":"deep"},
    }


def args():
    return dict(island_id="I",region="tropical",spatial_block="B",
       first_source="G",scope_rosters={
        "all_recorded":{"Alpha a","Beta a"},
        "strict_native":{"Alpha a"},
        "strict_native_nonendemic":{"Alpha a"}},
       sources={"geo5":{"Alpha b","Beta b"},
                "climate10":{"Alpha b","Beta b"}},
       motifs=make_motifs(states()))


def test_syndrome_motifs_are_genuine_same_species_3_or_4_axis():
    m=make_motifs(states())
    assert len(m)==len(MOTIFS)==12
    t=m["red_pink_radial_tubular_deep"]
    assert t["n_axes"]==4
    assert "Alpha a" in t["joint"]
    assert "Gamma onlycolor" not in t["known"]
    assert "Alpha b" not in t["joint"]
    assert "Beta b" not in t["joint"]


def test_source_genus_fourpart_exact_even_with_shared_and_disjoint_genera():
    m=make_motifs(states())["red_pink_radial_deep"]
    p=fourpart({"Alpha a","Beta a"},{"Alpha b","Beta b"},m)
    assert sum(p)==pytest.approx(.5)
    assert p[1]==pytest.approx(.5)
    v=fourpart({"Alpha a"},{"Beta a"},m)
    assert sum(v)==pytest.approx(1)
    assert v[2]==pytest.approx(1)


def test_source_excludes_all_accepted_island_species():
    a=args()
    a["sources"]["geo5"].add("Alpha a")
    with pytest.raises(ValueError,match="focal island"):
        make_island_motif_rows(**a)


def test_source_and_island_status_do_not_impute_missing_axes():
    rows=make_island_motif_rows(**args())
    assert len(rows)==len(MOTIFS)*len(SCOPES)
    r=next(q for q in rows if q["scope"]=="all_recorded" and
           q["motif"]=="red_pink_radial_tubular_deep")
    assert r["n_island_direct"]==2
    assert r["n_island_positive"]==1
    assert r["geo5__n_source_direct"]==2
    assert r["geo5__difference"]==pytest.approx(.5)
    assert sum(r["geo5__"+k] for k in (
       "shared_genus_source_weight","shared_genus_different_species_state",
       "island_only_recorded_genus","source_only_recorded_genus_negative"
       ))==pytest.approx(.5)
    assert not np.isnan(r["geo5__coupling_difference"])


def frame(n=40,active=True):
    base=make_island_motif_rows(**args())
    out=[]
    for i in range(n):
        for q in base:
            a=dict(q,island_id="I"+str(i),
                   spatial_block="B"+str(i//2),
                   first_source="S"+str(i//2))
            if active:
                a.update(n_island_total=10,n_island_direct=8,
                         n_island_positive=4,
                         island_frequency=.5,island_coupling=.15)
                for mode in ("geo5","climate10"):
                    a[mode+"__n_source_total"]=20
                    a[mode+"__n_source_direct"]=15
                    a[mode+"__source_frequency"]=.4
                    a[mode+"__source_coupling"]=.1
                    a[mode+"__difference"]=.1
                    a[mode+"__coupling_difference"]=.05
                    a[mode+"__shared_genus_source_weight"]=.04
                    a[mode+"__shared_genus_different_species_state"]=.02
                    a[mode+"__island_only_recorded_genus"]=.06
                    a[mode+"__source_only_recorded_genus_negative"]=-.02
            out.append(a)
    return pd.DataFrame(out)


def test_aggregate_full_island_population_qualified_and_signed():
    z=aggregate_motifs(frame(),population=40,draws=199)
    assert z["n_original_GIFT_botanical_islands"]==40
    assert z["n_reporting_panels"]==5*3*len(MOTIFS)
    a=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded"
           and p["motif"]=="red_pink_radial_deep")
    assert a["n_islands_two_source_supported"]==40
    assert a["all__spatial_block"]["status"]=="EXPLORATORY_REAL_SAME_PLANT_MOTIF_SOURCE_COMPOSITION"
    assert a["all__spatial_block"]["geo5__difference__equal_group_mean"]==pytest.approx(.1)
    assert a["n_documentation_balanced"]==40


def test_no_supported_morphology_must_not_be_scored_as_zero_guild():
    z=aggregate_motifs(frame(active=False),population=40,draws=199)
    assert all(p["all__spatial_block"]["status"]==
               "NOT_IDENTIFIABLE_DIRECT_SOURCE_OR_BLOCK_SUPPORT"
               for p in z["panels"])


def test_missing_signed_genus_terms_fail_closed():
    d=frame()
    d.loc[0,"geo5__island_only_recorded_genus"]=np.nan
    with pytest.raises(ValueError,match="incomplete signed"):
        aggregate_motifs(d,population=40,draws=199)

def test_no_false_documentation_balance_from_small_overall_annotation_fractions():
    # A giant unannotated mainland source and island can both appear 0% documented;
    # that is not valid support for a shared-genus matched trait comparison.
    d=frame()
    d.loc[d.motif.eq("red_pink_radial_deep"),"geo5__n_shared_direct_genera"]=0
    d.loc[d.motif.eq("red_pink_radial_deep"),"geo5__genus_conditioned_documentation_gap"]=np.nan
    z=aggregate_motifs(d,population=40,draws=199)
    x=next(p for p in z["panels"] if p["region"]=="tropical"
           and p["scope"]=="all_recorded"
           and p["motif"]=="red_pink_radial_deep")
    assert x["n_islands_two_source_supported"]==40
    assert x["n_documentation_balanced"]==0
    assert x["doc_gap_le10pp__spatial_block"]["status"]==(
        "NOT_IDENTIFIABLE_DIRECT_SOURCE_OR_BLOCK_SUPPORT")
