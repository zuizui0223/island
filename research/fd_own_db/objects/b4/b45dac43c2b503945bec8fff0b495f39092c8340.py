import pytest
from fd.dopi_gift_guild_motif_bridge import (
    screen_positive_dopi_rows, summarize_study_guild_motifs,
    SCOPES, PAIRS, MIN_STUDIES,
)


def fake_guild(order,family,species):
    return "Bombus" if str(species).startswith("Bombus") else (
        "butterfly_family" if str(species).startswith("Papilio") else "other_bee"
    )


def source_rows(nstudies=12,*,all_all=True):
    rows=[]
    for i in range(nstudies):
        survey="all" if all_all else "single species"
        for plant,animal,family,grade in (
            ("Alpha a","Papilio machaon","Papilionidae","3"),
            ("Beta a","Papilio machaon","Papilionidae","2"),
            ("Alpha b","Bombus terrestris","Apidae","3"),
            ("Beta b","Bombus terrestris","Apidae","1"),
        ):
            rows.append(dict(doi=f"ref{i}",publication=f"Paper {i}",
                plant_species=plant,pollinator_species=animal,
                pollinator_family=family,pollinator_order="Lepidoptera",
                methodology="do",pollinator_survey=survey,
                plant_survey=survey,pollination=grade))
    return rows


def profiles():
    known={"Alpha a","Alpha b","Beta a","Beta b"}
    return {
        "yellow_orange_radial_deep":dict(
            known=known, joint={"Alpha a","Beta a"}
        ),
        "open_radial_shallow":dict(
            known=known, joint={"Alpha b","Beta b"}
        ),
    }


def test_study_control_avoids_false_guild_absence_and_preserves_original_roster():
    source=source_rows()
    idx=screen_positive_dopi_rows(source,columns=set(source[0]),
        island_plant_names={"Alpha a","Alpha b","Beta a","Beta b","Other x"},
        motifs=profiles(),visitor_classifier=fake_guild)
    assert idx["original_source_row_counts"]["source_rows"]==48
    assert idx["n_distinct_publication_keys"]==12
    assert idx["n_original_island_accepted_plants_matching_DoPI"]==4
    assert idx["no_true_absence_inferred"]
    report=summarize_study_guild_motifs(idx,motifs=profiles(),draws=199)
    a=next(p for p in report["panels"] if p["guild_a"]=="butterfly_family" and
        p["guild_b"]=="Bombus" and
        p["motif"]=="yellow_orange_radial_deep" and
        p["evidence_scope"]=="positive_visitor_or_pollen_123")
    assert a["n_within_publication_trait_scored_contrasts"]==12
    assert a["all_study"]["status"]=="DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION"
    assert a["all_study"]["mean_guildA_minus_guildB"]==pytest.approx(1)
    assert a["shared_genus_study"]["status"]=="DESCRIPTIVE_POSITIVE_VISITORS_WITHIN_PUBLICATION"
    assert a["shared_genus_study"]["mean_guildA_minus_guildB"]==pytest.approx(1)
    assert report["limits"]["no_true_guild_absence_inferred"]


def test_pollen_grade_only_requires_both_eligible_guilds_each_study():
    source=source_rows()
    idx=screen_positive_dopi_rows(source,columns=set(source[0]),
        island_plant_names={"Alpha a","Alpha b","Beta a","Beta b"},
        motifs=profiles(),visitor_classifier=fake_guild)
    report=summarize_study_guild_motifs(idx,motifs=profiles(),draws=199)
    x=next(p for p in report["panels"] if p["guild_a"]=="butterfly_family" and
       p["guild_b"]=="Bombus" and p["motif"]=="yellow_orange_radial_deep" and
       p["evidence_scope"]=="pollen_evidence_12")
    assert x["n_within_publication_trait_scored_contrasts"]==0
    assert x["all_study"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_STUDIES"


def test_grade4_robber_cannot_be_visitor_or_pollination():
    rows=source_rows()
    rows[0]["pollination"]="4"
    idx=screen_positive_dopi_rows(rows,columns=set(rows[0]),
        island_plant_names={"Alpha a","Alpha b","Beta a","Beta b"},
        motifs=profiles(),visitor_classifier=fake_guild)
    assert idx["original_source_row_counts"]["unsupported_evidence_or_nonpollinating_rows"]==1
    assert "Alpha a" not in idx["indexed_units"][("ref0","butterfly_family",SCOPES[0])]


def test_guild_record_missing_is_not_controlled_negative():
    rows=[r for r in source_rows() if not r["pollinator_species"].startswith("Bombus")]
    idx=screen_positive_dopi_rows(rows,columns=set(rows[0]),
        island_plant_names={"Alpha a","Alpha b","Beta a","Beta b"},
        motifs=profiles(),visitor_classifier=fake_guild)
    output=summarize_study_guild_motifs(idx,motifs=profiles(),draws=199)
    assert all(p["all_study"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_STUDIES"
       for p in output["panels"] if p["guild_b"]=="Bombus" or p["guild_a"]=="Bombus")


def test_inconsistent_source_columns_or_empty_original_roster_fail_closed():
    with pytest.raises(ValueError,match="missing DoPI"):
        screen_positive_dopi_rows([],columns={"plant_species"},
            island_plant_names={"Alpha a"},motifs=profiles(),
            visitor_classifier=fake_guild)
    rows=source_rows()
    with pytest.raises(ValueError,match="original GIFT"):
        screen_positive_dopi_rows(rows,columns=set(rows[0]),
            island_plant_names=set(),motifs=profiles(),
            visitor_classifier=fake_guild)


def test_publication_bootstrap_support_not_interaction_row_pseudoreplication():
    rows=source_rows(nstudies=1)*100
    idx=screen_positive_dopi_rows(rows,columns=set(rows[0]),
        island_plant_names={"Alpha a","Alpha b","Beta a","Beta b"},
        motifs=profiles(),visitor_classifier=fake_guild)
    out=summarize_study_guild_motifs(idx,motifs=profiles(),draws=199)
    p=next(x for x in out["panels"] if x["guild_a"]=="butterfly_family" and
       x["guild_b"]=="Bombus" and
       x["motif"]=="yellow_orange_radial_deep" and
       x["evidence_scope"]=="positive_visitor_or_pollen_123")
    assert p["n_within_publication_trait_scored_contrasts"]==1
    assert p["all_study"]["status"]=="NOT_IDENTIFIABLE_INDEPENDENT_STUDIES"



def test_genuinely_unseen_flower_motif_cannot_look_like_precise_null():
    rows=source_rows()
    motif={"white_radial_tubular_deep":{
        "known":{"Alpha a","Alpha b","Beta a","Beta b"},
        "joint":set(),
    }}
    idx=screen_positive_dopi_rows(rows,columns=set(rows[0]),
        island_plant_names={"Alpha a","Alpha b","Beta a","Beta b"},
        motifs=motif,visitor_classifier=fake_guild)
    out=summarize_study_guild_motifs(idx,motifs=motif,draws=199)
    q=next(p for p in out["panels"] if p["guild_a"]=="butterfly_family"
        and p["guild_b"]=="Bombus" and
        p["evidence_scope"]=="positive_visitor_or_pollen_123")
    assert q["n_within_publication_trait_scored_contrasts"]==12
    assert q["n_studies_with_positive_motif_either"]==0
    assert q["all_study"]["status"]=="NOT_IDENTIFIABLE_MOTIF_POSITIVE_STUDY_SUPPORT"
    assert "mean_guildA_minus_guildB" not in q["all_study"]
