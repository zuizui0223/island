import pytest
from fd.dopi_site_period_guild import (
    _site,_time,site_time_positive_index,site_period_audit,
)


def guild(order,family,name):
    return "Bombus" if name=="Bombus terrestris" else "Syrphidae_hoverfly"


def rows(nstudies=12,*,different_year=False,different_month=False,
         same_location=True,pollen=False):
    q=[]
    for i in range(nstudies):
        for group,plants,animal,family in (
            ("bee",("Alpha a","Beta a"),"Bombus terrestris","Apidae"),
            ("fly",("Alpha b","Beta b"),"Episyrphus balteatus","Syrphidae"),
        ):
            for plant in plants:
                q.append(dict(doi=f"doi{i}",publication=f"Ref {i}",
                              plant_species=plant,pollinator_species=animal,
                              pollinator_order="Diptera" if group=="fly" else "Hymenoptera",
                              pollinator_family=family,
                              pollination="2" if pollen else "3",
                              latitude="51.523" if group=="bee" or same_location else "51.525",
                              longitude="-0.123",grid_letter="",grid_code="",
                              country="UK",
                              year="2025" if group=="bee" or not different_year else "2024",
                              month="6" if group=="bee" or not different_month else "7"))
    return q


def motif():
    return {"open_radial":dict(known={"Alpha a","Alpha b","Beta a","Beta b"},
                               joint={"Alpha a"},n_axes=1)}


def run(source):
    return site_time_positive_index(
        source,fields=set(source[0]),island_names=set(motif()["open_radial"]["known"]),
        motifs=motif(),classifier=guild)


def panel(result,resolution,tier):
    return next(q for q in result["panels"] if
                q["guild_a"]=="Bombus" and q["guild_b"]=="Syrphidae_hoverfly"
                and q["floral_motif"]=="open_radial" and
                q["resolution"]==resolution and q["evidence_tier"]==tier)


def test_same_true_site_and_year_month_guild_plus_genus():
    indexed=run(rows())
    z=site_period_audit(indexed,motifs=motif(),draws=299)
    assert z["n_original_global_GIFT_botanical_islands"]==4505
    assert z["n_reporting_panels"]==4*5
    for resolution in ("same_study_site_year","same_study_site_year_month"):
        q=panel(z,resolution,"positive_visitation_123")
        assert q["n_study_site_period_units_with_both_guilds"]==12
        assert q["n_study_site_period_units_with_two_shared_genera"]==12
        a=q["same_site_positive_visitor"]
        g=q["same_site_equal_genus"]
        assert a["status"]=="DESCRIPTIVE_POSITIVE_GUILD_COMMON_SITE_PERIOD"
        assert g["status"]=="DESCRIPTIVE_POSITIVE_GUILD_COMMON_SITE_PERIOD"
        assert a["guild_a_minus_guild_b_motif_frequency"]==pytest.approx(.5)
        assert g["guild_a_minus_guild_b_motif_frequency"]==pytest.approx(.5)


def test_different_observation_year_never_becomes_same_site_positive_guild():
    source=rows(different_year=True)
    z=site_period_audit(run(source),motifs=motif(),draws=299)
    q=panel(z,"same_study_site_year","positive_visitation_123")
    assert q["n_study_site_period_units_with_both_guilds"]==0
    assert q["same_site_positive_visitor"]["status"]=="NOT_IDENTIFIABLE_COMMON_SITE_INDEPENDENT_STUDIES"


def test_different_month_site_year_is_available_but_site_month_is_not():
    z=site_period_audit(run(rows(different_month=True)),motifs=motif(),draws=299)
    assert panel(z,"same_study_site_year","positive_visitation_123")[
        "same_site_positive_visitor"]["status"]=="DESCRIPTIVE_POSITIVE_GUILD_COMMON_SITE_PERIOD"
    assert panel(z,"same_study_site_year_month","positive_visitation_123")[
        "same_site_positive_visitor"]["status"]=="NOT_IDENTIFIABLE_COMMON_SITE_INDEPENDENT_STUDIES"


def test_nearby_but_different_published_sites_never_join():
    z=site_period_audit(run(rows(same_location=False)),motifs=motif(),draws=299)
    assert panel(z,"same_study_site_year","positive_visitation_123")[
        "n_study_site_period_units_with_both_guilds"]==0


def test_no_site_month_or_full_year_never_invented():
    assert _site(dict(country="UK",latitude="51.5",longitude="-0.1")) is None
    assert _site(dict(country="UK",grid_letter="TQ",grid_code="2345")) is not None
    assert _time(dict(year="NA",month="5")) is None
    assert _time(dict(year="2025",month="NA"))==("2025",None)
    assert _time(dict(year="2025",month="13"))==("2025",None)


def test_grade3_visits_not_evidence_of_pollen_transfer():
    z=site_period_audit(run(rows()),motifs=motif(),draws=299)
    q=panel(z,"same_study_site_year","pollen_supported_12")
    assert q["same_site_positive_visitor"]["status"]=="NOT_IDENTIFIABLE_COMMON_SITE_INDEPENDENT_STUDIES"
    zz=site_period_audit(run(rows(pollen=True)),motifs=motif(),draws=299)
    assert panel(zz,"same_study_site_year","pollen_supported_12")[
        "same_site_positive_visitor"]["status"]=="DESCRIPTIVE_POSITIVE_GUILD_COMMON_SITE_PERIOD"


def test_source_schema_and_missing_site_do_not_become_effect_zero():
    src=rows()
    with pytest.raises(ValueError,match="DoPI site/year"):
        site_time_positive_index(src,fields={"plant_species"},
            island_names={"Alpha a"},motifs=motif(),classifier=guild)
    for q in src:q["latitude"]="";q["grid_code"]=""
    idx=run(src)
    z=site_period_audit(idx,motifs=motif(),draws=299)
    assert all(p["same_site_positive_visitor"]["status"]==
               "NOT_IDENTIFIABLE_COMMON_SITE_INDEPENDENT_STUDIES"
               for p in z["panels"])
