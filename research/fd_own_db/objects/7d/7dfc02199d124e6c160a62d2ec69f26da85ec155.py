import numpy as np
import pytest

from fd.dopi_curveball_paired_degree_null import (
    curveball_trades,site_curveball_null,audit_curveball_vs_swap,
    _rows_to_edges,_same_publication_bootstrap,FOCAL_MOTIFS,
)
from fd.dopi_degree_fixed_visitor_null import (
    site_degree_preserving_null,_degrees,
)


def source():
    edges={
        ("Alpha pos","Bombus terrestris"),
        ("Alpha pos","Episyrphus balteatus"),
        ("Alpha neg","Bombus lapidarius"),
        ("Alpha neg","Syrphus ribesii"),
        ("Beta pos","Bombus lapidarius"),
        ("Beta pos","Syrphus ribesii"),
        ("Beta neg","Bombus terrestris"),
        ("Beta neg","Episyrphus balteatus"),
    }
    guild={
        "Bombus terrestris":"Bombus",
        "Bombus lapidarius":"Bombus",
        "Episyrphus balteatus":"Syrphidae_hoverfly",
        "Syrphus ribesii":"Syrphidae_hoverfly",
    }
    return edges,guild


def traits():
    known={"Alpha pos","Alpha neg","Beta pos","Beta neg"}
    return {m:dict(known=known,joint={"Alpha pos","Beta pos"},axis_count=1)
            for m in FOCAL_MOTIFS}


def test_independent_curveball_proposal_preserves_every_plant_and_insect_degree():
    edges,guild=source()
    original=_degrees(edges)
    rows={}
    for plant,v in edges:rows.setdefault(plant,set()).add(v)
    rng=np.random.default_rng(20261118)
    out,n=curveball_trades(rows,rng,500)
    assert n>0
    assert _degrees(_rows_to_edges(out))==original
    assert sum(len(x) for x in out.values())==8
    assert _rows_to_edges(rows)==edges # deep copy does not mutate original
    assert set(v for _,v in _rows_to_edges(out))==set(guild)


def test_curveball_and_swap_same_original_species_and_flower_cohort():
    edges,guild=source()
    x=site_curveball_null(edges,guild,traits(),seed="same_pond",draws=100)
    y=site_degree_preserving_null(edges,guild,traits(),seed="same_pond",draws=100)
    assert x["status"]=="ORIGINAL_SOURCE_CURVEBALL_SITE_PROFILES"
    assert y["status"]=="ORIGINAL_SOURCE_DEGREE_FIXED_SITE_PROFILES"
    for m in FOCAL_MOTIFS:
        a=x["motif_profiles"][m]
        b=y["motif_profiles"][m]
        assert a["status"]=="DESCRIPTIVE_CURVEBALL_BOTH_DEGREES_FIXED"
        assert b["status"]=="DESCRIPTIVE_FIXED_BOTH_SIDES_DEGREES_NETWORK_NULL"
        assert a["original_observed_guild_richness_difference"]==pytest.approx(
            b["original_observed_guild_richness_difference"])
        assert abs(a["excess_over_curveball_null"]-
                   b["excess_over_network_degree_null"])<.1


def test_saturated_bipartite_source_has_no_actual_randomized_opportunity():
    plants=["Alpha pos","Alpha neg","Beta pos","Beta neg"]
    guilds=["Bombus terrestris","Episyrphus balteatus","Syrphus ribesii"]
    guild={"Bombus terrestris":"Bombus",
           "Episyrphus balteatus":"Syrphidae_hoverfly",
           "Syrphus ribesii":"Syrphidae_hoverfly"}
    all_edges={(p,v) for p in plants for v in guilds}
    q=site_curveball_null(all_edges,guild,traits(),seed="complete",draws=40)
    assert q["status"]=="ORIGINAL_SOURCE_CURVEBALL_SITE_PROFILES"
    cell=q["motif_profiles"]["open_radial"]
    assert cell["status"]=="NOT_IDENTIFIABLE_CURVEBALL_NO_TRADABLE_OPPORTUNITY"
    assert "excess_over_curveball_null" not in cell


def test_paired_source_study_estimate_is_not_raw_network_or_site_pseudoreplication():
    edges,guild=source()
    index={
        "private_positive_edges":{},
        "documented_guild_absence_not_inferred":True,
        "n_original_source_studies_qualified":12,
        "source_row_counts":dict(original_DoPI_v2_source_rows=676900),
    }
    for i in range(12):
        index["private_positive_edges"][(
            f"original_paper_{i}",("original_site","53.1","-1.4"),"2025","06"
        )]={(p,v,guild[v]) for p,v in edges}
    out=audit_curveball_vs_swap(index,motifs=traits(),draws=199,site_draws=40)
    assert out["status"]=="GIFT4505_DOPI2026_IDENTICAL_SITE_CURVEBALL_VS_EDGE_SWAP"
    assert out["n_original_GIFT_botanical_islands"]==4505
    assert out["n_original_DoPI_v2_raw_rows"]==676900
    assert out["n_original_site_month_positive_networks"]==12
    assert out["chains_per_method"]==4
    assert out["snapshots_per_chain"]==40
    for item in out["panels"]:
        # Exact common-source conditioning can exclude papers whose one
        # randomized chain fails the strict mixing gate. This is real missing
        # support, not a reason to fabricate twelve numerical estimates.
        assert 10<=item["n_original_source_publications"]<=12
        assert item["n_original_study_site_months"]==item["n_original_source_publications"]
        assert item["status"]=="DESCRIPTIVE_IDENTICAL_SOURCE_SITE_TWO_DEGREE_NULLS"
        assert np.isfinite(item["edge_swap_excess"])
        assert np.isfinite(item["curveball_excess"])
        assert item["curveball_minus_edge_swap_excess"]==pytest.approx(
            item["curveball_excess"]-item["edge_swap_excess"])
        assert item["observed_guild_breadth_difference"]==pytest.approx(0)
        assert item["n_source_papers_methods_have_different_effect_sign"]<=12
    assert out["limitations"]["no_unobserved_visitors_or_guild_absences_inferred"]
    assert out["limitations"]["four_focal_morphologies_posthoc_selected_from_27"]


def test_source_missing_real_morphology_is_not_null_pollinator_absence():
    edges,guild=source()
    m=traits()
    m["open_radial"]["joint"]=set()
    z=site_curveball_null(edges,guild,m,seed="no_motif",draws=40)
    assert "open_radial" not in z["motif_profiles"]
    idx={"private_positive_edges":{
        ("paper1","site","2025","6"):{(p,v,guild[v]) for p,v in edges}
    },"documented_guild_absence_not_inferred":True,
       "n_original_source_studies_qualified":1,
       "source_row_counts":dict(original_DoPI_v2_source_rows=676900)}
    report=audit_curveball_vs_swap(idx,motifs=m,draws=199,site_draws=40)
    a=next(x for x in report["panels"] if x["motif"]=="open_radial")
    assert a["status"]=="NOT_IDENTIFIABLE_PAIRED_SOURCE_PUBLICATIONS"


def test_missing_source_positive_or_bad_seed_input_fails_closed():
    with pytest.raises(ValueError,match="negative"):
        curveball_trades({},np.random.default_rng(1),-1)
    edges,guild=source()
    with pytest.raises(ValueError,match="guild identity"):
        site_curveball_null(edges,{"Bombus terrestris":"Bombus"},
                            traits(),seed="wrong",draws=40)
    with pytest.raises(ValueError,match="draws"):
        site_curveball_null(edges,guild,traits(),seed="wrong",draws=1)
    with pytest.raises(ValueError,match="positive observer"):
        audit_curveball_vs_swap(
            {"documented_guild_absence_not_inferred":False},
            motifs=traits(),draws=199,site_draws=40)


def test_paired_paper_bootstrap_is_descriptive_and_never_fabricates_missing_papers():
    values={
        f"paper{i}":[dict(
            observed_guild_breadth_difference=.2,
            edge_swap_expected=.1+i/1000,
            curveball_expected=.101+i/1000,
            edge_swap_excess=.1-i/1000,
            curveball_excess=.099-i/1000,
            curveball_minus_edge_swap_excess=-.001
        )] for i in range(12)}
    result=_same_publication_bootstrap(values,seed="synthetic",draws=199)
    assert result["status"]=="DESCRIPTIVE_IDENTICAL_SOURCE_SITE_TWO_DEGREE_NULLS"
    assert result["n_original_source_publications"]==12
    assert result["curveball_minus_edge_swap_excess"]==pytest.approx(-.001)
    assert result["n_source_papers_methods_have_different_effect_sign"]==0
    assert result["source_pub_leave_one_out_curveball_min"]>0
    assert result["curveball_excess_ci95"][0] <= result["curveball_excess"]
