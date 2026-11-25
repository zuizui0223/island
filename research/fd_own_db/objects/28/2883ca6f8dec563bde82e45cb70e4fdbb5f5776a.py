from collections import defaultdict
import pytest

from fd.dopi_site_period_guild import site_time_positive_index
from fd.dopi_positive_coaccess import summarize_positive_coaccess


def classifier(order,family,animal):
    return ("Bombus" if str(animal)=="Bombus terrestris"
            else "Syrphidae_hoverfly")


def fake_rows(nstudies=12,*,different_month=False,level="3",
              no_same_accepted_species=False):
    rows=[]
    for i in range(nstudies):
        a={"Alpha p","Alpha n","Beta p","Gamma n"}
        b={"Alpha p","Beta p","Beta n","Gamma p"}
        if no_same_accepted_species:
            b={"Beta n","Gamma p","Delta n"}
        for guild,animals in (("Bombus",a),("Syrphidae",b)):
            for sp in animals:
                rows.append(dict(
                    doi=f"original_study{i}",publication=f"Publication {i}",
                    plant_species=sp,
                    pollinator_species=(
                        "Bombus terrestris" if guild=="Bombus"
                        else "Episyrphus balteatus"),
                    pollinator_family=(
                        "Apidae" if guild=="Bombus" else "Syrphidae"),
                    pollinator_order=(
                        "Hymenoptera" if guild=="Bombus" else "Diptera"),
                    pollination=level,
                    year="2025",month="06" if guild=="Bombus" or
                                     not different_month else "07",
                    latitude="51.235",longitude="-0.123",
                    country="United Kingdom",grid_letter="",grid_code="",
                ))
    return rows


def profile(*,genus_fixed=False):
    known={"Alpha p","Alpha n","Beta p","Beta n","Gamma p","Gamma n"}
    if genus_fixed:
        positive={"Alpha p","Alpha n","Beta p","Beta n"}
    else:
        positive={"Alpha p","Beta p","Gamma p"}
    return {"open_radial":dict(known=known,joint=positive,axis_count=1)}


def indexed(rows,motifs):
    return site_time_positive_index(
        rows,fields=set(rows[0]),
        island_names=set(motifs["open_radial"]["known"]),
        motifs=motifs,classifier=classifier,
    )


def panel(z,*,resolution="same_study_site_year_month",
          evidence="positive_visitation_123"):
    return next(q for q in z["panels"]
         if q["published_location_scope"]==resolution
         and q["visitor_evidence"]==evidence
         and q["guild_a"]=="Bombus"
         and q["guild_b"]=="Syrphidae_hoverfly"
         and q["floral_motif"]=="open_radial")


def test_same_accepted_species_observed_by_two_guilds_preserves_study_period():
    motifs=profile()
    d=indexed(fake_rows(),motifs)
    z=summarize_positive_coaccess(d,motifs=motifs,draws=199)
    assert z["original_GIFT_botanical_flora_population"]==4505
    assert z["n_panels"]==2*2*5
    q=panel(z)
    assert q["n_source_site_periods_with_both_guilds_positive"]==12
    assert q["n_source_site_periods_with_at_least_one_same_plant_positive_both"]==12
    assert q["n_distinct_source_positive_co_documented_plant_names"]==2
    assert q["n_site_periods_with_informative_within_genus_morphology"]==12
    assert q["genus_fixed_co_documentation_excess"]["n_source_publications"]==12
    assert q["co_documented_vs_positive_union"]["status"]=="DESCRIPTIVE_REAL_GUILD_CO_DOCUMENTATION"
    assert q["co_documented_vs_positive_union"]["study_equal_difference"]==pytest.approx(.5)
    assert q["genus_fixed_co_documentation_excess"]["status"]=="DESCRIPTIVE_GENUS_CONDITIONED_CO_DOCUMENTATION_EXCESS"
    assert q["genus_fixed_co_documentation_excess"]["study_equal_difference"]==pytest.approx(.5)
    assert z["limitations"]["guild_not_reported_on_same_plant_is_not_true_absence"]


def test_genus_fixed_flower_motif_cannot_identify_within_genus_assoc():
    motifs=profile(genus_fixed=True)
    z=summarize_positive_coaccess(indexed(fake_rows(),motifs),
        motifs=motifs,draws=199)
    q=panel(z)
    assert q["co_documented_vs_positive_union"]["status"]=="DESCRIPTIVE_REAL_GUILD_CO_DOCUMENTATION"
    assert q["co_documented_vs_positive_union"]["study_equal_difference"]==pytest.approx(1/3)
    assert q["n_site_periods_with_informative_within_genus_morphology"]==0
    assert q["genus_fixed_co_documentation_excess"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PUBLICATIONS"
    assert "study_equal_difference" not in q["genus_fixed_co_documentation_excess"]


def test_different_months_not_same_month_co_documentation():
    motifs=profile()
    z=summarize_positive_coaccess(
        indexed(fake_rows(different_month=True),motifs),
        motifs=motifs,draws=199)
    assert panel(z)["co_documented_vs_positive_union"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PUBLICATIONS"
    assert panel(z,resolution="same_study_site_year")[
        "co_documented_vs_positive_union"]["status"]=="DESCRIPTIVE_REAL_GUILD_CO_DOCUMENTATION"


def test_different_accepted_species_never_infer_both_visitor_guilds():
    motifs=profile()
    z=summarize_positive_coaccess(
        indexed(fake_rows(no_same_accepted_species=True),motifs),
        motifs=motifs,draws=199)
    q=panel(z)
    assert q["n_source_site_periods_with_at_least_one_same_plant_positive_both"]==0
    assert q["co_documented_vs_positive_union"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PUBLICATIONS"


def test_visitation_grade3_not_pollen_transfer():
    motifs=profile()
    z=summarize_positive_coaccess(indexed(fake_rows(level="3"),motifs),
        motifs=motifs,draws=199)
    assert panel(z,evidence="pollen_supported_12")[
        "co_documented_vs_positive_union"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PUBLICATIONS"
    z2=summarize_positive_coaccess(indexed(fake_rows(level="2"),motifs),
        motifs=motifs,draws=199)
    assert panel(z2,evidence="pollen_supported_12")[
        "co_documented_vs_positive_union"]["status"]=="DESCRIPTIVE_REAL_GUILD_CO_DOCUMENTATION"


def test_not_documented_flower_motif_does_not_become_zero_effect():
    m=profile()
    m["open_radial"]["joint"]=set()
    z=summarize_positive_coaccess(indexed(fake_rows(),m),
        motifs=m,draws=199)
    q=panel(z)
    assert q["co_documented_vs_positive_union"]["status"]=="NOT_IDENTIFIABLE_MOTIF_POSITIVE_DOCUMENTATION"
    assert q["genus_fixed_co_documentation_excess"]["status"]=="NOT_IDENTIFIABLE_SOURCE_PUBLICATIONS"
