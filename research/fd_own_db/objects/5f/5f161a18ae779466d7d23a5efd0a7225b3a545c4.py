import pytest

from fd.dopi_guild_portfolio_effort_null import (
    accepted_visitor_binomial,expected_positive_guild_richness,
    positive_visitor_effort_index,
    guild_richness_motif_audit,
)


def classifier(order,family,animal):
    return ("Syrphidae_hoverfly" if "Syrphidae" in str(family)
            else "Bombus")


def original_visitors(nstudies=12, *,all_bees=False,motif_only=False,
                      different_month=False):
    out=[]
    for i in range(nstudies):
        for plant, visitors in (
            ("Alpha pos",("Bombus terrestris","Episyrphus balteatus")),
            ("Alpha neg",("Bombus terrestris","Bombus lapidarius")),
            ("Beta pos",("Bombus terrestris","Episyrphus balteatus")),
            ("Beta neg",("Bombus terrestris","Bombus lapidarius")),
        ):
            if all_bees:
                visitors=tuple("Bombus impatiens" if
                               s.startswith("Episyrphus") else s
                               for s in visitors)
            for animal in visitors:
                fly=animal.startswith("Episyrphus")
                out.append(dict(
                    doi=f"study{i}",publication=f"ref{i}",
                    plant_species=plant,pollinator_species=animal,
                    pollinator_family="Syrphidae" if fly else "Apidae",
                    pollinator_order="Diptera" if fly else "Hymenoptera",
                    pollination="3",country="UK",
                    grid_letter="",grid_code="",
                    latitude="51.235",longitude="-0.123",
                    year="2025",
                    month="7" if different_month and fly else "6"))
    return out


def traits(*,single_genus=False,all_same=False):
    known={"Alpha pos","Alpha neg","Beta pos","Beta neg"}
    positive=known if all_same else {"Alpha pos","Beta pos"}
    return {"open_radial":dict(known=known,joint=positive,n_axes=1)}


def run(rows,profiles=None):
    profiles=profiles or traits()
    return positive_visitor_effort_index(rows,fields=set(rows[0]),
        GIFT_accepted_names=set(profiles["open_radial"]["known"]),
        motifs=profiles,visitor_classifier=classifier)


def single(out):
    return out["panels"][0]


def test_exact_finite_guild_margins_for_one_positive_source_plant():
    # Eight total positive plant×visitor-species records: six bees, two flies.
    # A plant with TWO observed visitor species has E[K]=1+12/28=10/7.
    assert expected_positive_guild_richness(2,{"Bombus":6,"Syrphidae_hoverfly":2})==pytest.approx(10/7)
    assert expected_positive_guild_richness(1,{"Bombus":6,"Syrphidae_hoverfly":2})==pytest.approx(1)
    assert expected_positive_guild_richness(8,{"Bombus":6,"Syrphidae_hoverfly":2})==pytest.approx(2)
    assert expected_positive_guild_richness(0,{"Bombus":6,"Syrphidae_hoverfly":2})==0
    with pytest.raises(ValueError):
        expected_positive_guild_richness(9,{"Bombus":6,"Syrphidae_hoverfly":2})


def test_true_positive_guild_portfolio_remains_after_fixed_exposure_null():
    d=run(original_visitors())
    out=guild_richness_motif_audit(d,motifs=traits(),draws=299)
    p=single(out)
    assert out["original_GIFT_botanical_island_population"]==4505
    assert out["n_original_DOPI_source_rows"]==12*8
    assert out["n_scored_source_sites_with_two_or_more_visitor_taxa"]==12
    assert p["n_original_site_months_with_both_motif_states_and_two_visitor_taxa"]==12
    assert p["n_original_site_months_with_same_genus_both_flower_states"]==12
    for key in ("study_site","same_genus_study_site"):
        z=p[key]
        assert z["status"]=="DESCRIPTIVE_EFFORT_CONDITIONED_POSITIVE_VISITOR_GUILD_PORTFOLIO"
        assert z["n_original_source_publications"]==12
        assert z["observed_positive_guild_richness_difference"]==pytest.approx(1)
        assert z["expected_difference_fixed_visitor_opportunity"]==pytest.approx(0)
        assert z["excess_over_fixed_source_effort_null"]==pytest.approx(1)
    assert out["limitations"]["guild_nondetection_not_biological_absence"]


def test_duplicate_animal_edges_never_create_fake_positive_visitor_species_opportunity():
    rows=original_visitors()
    d1=run(rows)
    d2=run(rows+rows*3)
    out1=guild_richness_motif_audit(d1,motifs=traits(),draws=199)
    out2=guild_richness_motif_audit(d2,motifs=traits(),draws=199)
    assert single(out1)["study_site"]["excess_over_fixed_source_effort_null"]==pytest.approx(
        single(out2)["study_site"]["excess_over_fixed_source_effort_null"])
    assert d2["source_row_counts"]["original_DoPI_v2_source_rows"]==len(rows)*4


def test_all_same_guild_never_establishes_actual_multiguild_breadth():
    d=run(original_visitors(all_bees=True))
    out=guild_richness_motif_audit(d,motifs=traits(),draws=199)
    p=single(out)
    assert p["n_real_positive_multiguild_plant_site_month_units"]==0
    assert p["study_site"]["status"]=="NOT_IDENTIFIABLE_REAL_MULTIGUILD_DOCUMENTATION"


def test_no_within_genus_flower_variation_is_not_no_effect():
    d=run(original_visitors())
    out=guild_richness_motif_audit(d,motifs=traits(all_same=True),draws=199)
    p=single(out)
    assert p["n_original_site_months_with_both_motif_states_and_two_visitor_taxa"]==0
    assert p["same_genus_study_site"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_PUBLICATIONS"


def test_different_published_months_never_pooled_into_shared_observation():
    d=run(original_visitors(different_month=True))
    out=guild_richness_motif_audit(d,motifs=traits(),draws=199)
    p=single(out)
    assert p["study_site"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_PUBLICATIONS"
    assert p["n_original_site_months_with_both_motif_states_and_two_visitor_taxa"]==0


def test_nonpollination_grade4_and_missing_source_time_not_inferred():
    rows=original_visitors()
    for row in rows:
        if row["pollinator_family"]=="Syrphidae":
            row["pollination"]="4"
    index=run(rows)
    out=guild_richness_motif_audit(index,motifs=traits(),draws=199)
    # Once grade 4 non-pollinating interactions are removed, focal plants
    # lack enough comparable observed visitor species: no study controls.
    assert single(out)["study_site"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_PUBLICATIONS"

    bad=original_visitors()
    for row in bad:
        row["month"]="NA"
    with pytest.raises(ValueError,match="no original study-site"):
        guild_richness_motif_audit(run(bad),motifs=traits(),draws=199)


def test_true_visitor_species_not_bombus_sp_unresolved_or_author_variants():
    assert accepted_visitor_binomial("Bombus terrestris")=="bombus terrestris"
    assert accepted_visitor_binomial("Bombus terrestris (Linnaeus)")=="bombus terrestris"
    assert accepted_visitor_binomial("Bombus_terrestris")=="bombus terrestris"
    assert accepted_visitor_binomial("Episyrphus balteatus")=="episyrphus balteatus"
    for label in ("Bombus sp.","Bombus spp.","Bombus cf. terrestris",
                  "Bombus aff. terrestris","Bombus", "unidentified Syrphidae"):
        assert accepted_visitor_binomial(label) is None
    rows=original_visitors()
    row=dict(rows[0]);row["pollinator_species"]="Bombus sp."
    rows.append(row)
    row2=dict(rows[0]);row2["pollinator_species"]="Bombus terrestris (Linnaeus)"
    rows.append(row2)
    idx=run(rows)
    assert idx["source_row_counts"]["original_visitor_taxon_not_resolved_to_species"]==1
    # The original and the author-suffixed accepted visitor name must be
    # identical positive-plant/source-site degree records.
    assert len(idx["private_positive_edges"])==12
    assert sum(len(v) for v in idx["private_positive_edges"].values())==12*8
