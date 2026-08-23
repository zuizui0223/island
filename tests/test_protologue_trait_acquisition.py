import importlib
from pathlib import Path

import pandas as pd
import pytest
import typer

protologue = importlib.import_module("island_v2.protologue_trait_acquisition")

CONFIG = protologue.load_config(Path("config/protologue_acquisition.yml"))


def _extract(text):
    return protologue.extract_floral_colour(text, CONFIG)


# --- folding -------------------------------------------------------------


def test_folding_maps_every_character_back_to_the_original():
    folded, origin = protologue.fold("Blüten weiß")
    assert folded == "bluten weiss"
    assert len(folded) == len(origin)
    # both characters of the expanded eszett point at the one original character
    assert origin[folded.index("ss")] == origin[folded.index("ss") + 1]
    assert "Blüten weiß"[origin[-1]] == "ß"


def test_diacritics_do_not_shift_offsets():
    folded, origin = protologue.fold("corolle blanchâtre")
    assert folded == "corolle blanchatre"
    assert origin == list(range(len("corolle blanchâtre")))


# --- multilingual extraction --------------------------------------------


@pytest.mark.parametrize(
    ("description", "value"),
    [
        ("Corolla alba, tubo brevi.", "white"),
        ("Flores lutei, in racemis dispositi.", "yellow_orange"),
        ("Corolla purpurea, limbo patente.", "blue_purple"),
        ("Fleurs blanches, odorantes.", "white"),
        ("Corolle rouge vif, tube court.", "red_pink"),
        ("Blüten gelb, in Trauben.", "yellow_orange"),
        ("Blumenkrone weiß, aussen bräunlich am Kelch.", "white"),
        ("Flores blancos, perfumados.", "white"),
        ("Corola azul, tubo corto.", "blue_purple"),
    ],
)
def test_colour_is_read_in_every_supported_language(description, value):
    outcome, normalized, _, quote = _extract(description)
    assert outcome == protologue.ACCEPTED, description
    assert normalized == value
    assert quote


def test_english_only_rules_would_have_missed_these():
    # The measured failure this module exists to avoid: the recoverable slice of
    # the WFO rejection ledger was non-English statements nothing bound.
    english_terms = {"white", "yellow", "red", "blue", "green", "purple", "pink"}
    for description in ("Corolla alba.", "Fleurs jaunes.", "Blüten rötlich."):
        assert not english_terms & set(description.lower().split())
        assert _extract(description)[0] == protologue.ACCEPTED


# --- word boundaries -----------------------------------------------------


def test_an_epithet_is_not_read_as_a_colour():
    # "rotundifolia" contains German "rot"; "albacea" contains Latin "alba".
    # Substring matching would turn both into flower colour.
    assert _extract("Folia rotundifolia, caule glabro.")[0] != protologue.ACCEPTED
    outcome, _, _, _ = _extract("Corolla infundibuliformis, rotundifolia lamina.")
    assert outcome != protologue.ACCEPTED


def test_a_compound_organ_word_is_not_split():
    # "Blütenblätter" must read as a floral organ, never as "Blätter".
    outcome, value, _, _ = _extract("Blütenblätter weiß.")
    assert outcome == protologue.ACCEPTED
    assert value == "white"


# --- organ attachment ----------------------------------------------------


def test_fruit_colour_never_becomes_flower_colour():
    assert _extract("Fructus ruber, maturus.")[0] == protologue.REJECT_COMPETING
    assert _extract("Fruit rouge à maturité.")[0] == protologue.REJECT_COMPETING
    assert _extract("Frucht rot.")[0] == protologue.REJECT_COMPETING


def test_calyx_and_style_are_competing_organs():
    assert _extract("Calyx viridis, corolla decidua.")[0] == protologue.REJECT_COMPETING
    assert _extract("Stylus purpureus, exsertus.")[0] == protologue.REJECT_COMPETING


def test_a_tie_between_organs_resolves_to_the_competing_one():
    # "rubri" sits one character from each organ, so nothing decides it on
    # distance. The conservative reading wins and the row is dropped.
    outcome, _, _, _ = _extract("Flores rubri fructus.")
    assert outcome == protologue.REJECT_COMPETING


def test_the_nearer_organ_still_decides_when_both_are_present():
    # Not a tie: "flores" is adjacent and "fructus" is eleven characters away.
    # The sentence says both organs are red, and flower colour red is correct.
    outcome, value, _, _ = _extract("Fructus et flores rubri.")
    assert outcome == protologue.ACCEPTED
    assert value == "red_pink"


def test_a_distant_colour_word_is_not_anchored_to_the_flower():
    far = "Flores in racemis terminalibus " + ("dispositi " * 12) + "solo rubro."
    assert _extract(far)[0] in {
        protologue.REJECT_NO_ORGAN,
        protologue.REJECT_COMPETING,
    }


# --- negation, variability, absence --------------------------------------


def test_declining_to_state_the_colour_is_not_an_observation():
    for description in (
        "Flores ignoti.",
        "Corolla non visa.",
        "Fleurs inconnues.",
        "Blüten unbekannt.",
        "Flores ? albi.",
    ):
        assert _extract(description)[0] == protologue.REJECT_NEGATED, description


def test_two_ontology_values_become_variable_rather_than_the_first_hit():
    outcome, value, _, _ = _extract("Corolla alba demum rosea.")
    assert outcome == protologue.ACCEPTED
    assert value == "multicolored_variable"


def test_two_terms_in_the_same_group_stay_one_value():
    outcome, value, matched, _ = _extract("Corolla rosea, rubro-striata.")
    assert outcome == protologue.ACCEPTED
    assert value == "red_pink"
    assert "+" in matched


def test_text_without_colour_is_distinguished_from_text_without_text():
    assert _extract("Frutex 2 m altus, ramis glabris.")[0] == protologue.REJECT_NO_COLOUR
    assert _extract("   ")[0] == protologue.REJECT_NO_TEXT


def test_the_quote_is_verbatim_not_folded():
    _, _, _, quote = _extract("Arbor 5 m. Blumenkrone weiß, duftend. Frucht rot.")
    assert "weiß" in quote
    assert "Frucht" not in quote


# --- target queue --------------------------------------------------------


def _master():
    return pd.DataFrame(
        {
            "accepted_species": ["Narrow one", "Narrow two", "Widespread sp"],
            "genus": ["Narrow", "Narrow", "Widespread"],
            "family": ["Fam", "Fam", "Fam"],
            "n_islands": [1, 3, 40],
            "n_records": [4, 9, 5000],
        }
    )


def test_the_queue_is_the_tail_not_the_poorly_recorded():
    queue = protologue.build_target_queue(_master(), CONFIG)
    assert list(queue["accepted_species"]) == ["Narrow one", "Narrow two"]


def test_already_resolved_species_are_skipped():
    queue = protologue.build_target_queue(_master(), CONFIG, {"Narrow one"})
    assert list(queue["accepted_species"]) == ["Narrow two"]


# --- citation parsing ----------------------------------------------------


def test_a_citation_for_a_different_name_is_not_used():
    payload = {"results": [{"name": "Other species", "publication": "Fl.", "page": "1", "publicationYear": "1900"}]}
    assert protologue.parse_ipni_citation(payload, "Cyanea aspera", CONFIG)["outcome"] == (
        protologue.REJECT_NO_CITATION
    )


def test_a_citation_without_a_page_cannot_locate_a_scan():
    payload = {"results": [{"name": "Cyanea aspera", "publication": "Fl. Hawaii", "publicationYear": "1888"}]}
    result = protologue.parse_ipni_citation(payload, "Cyanea aspera", CONFIG)
    assert result["outcome"] == protologue.REJECT_CITATION_INCOMPLETE


def test_a_complete_citation_is_carried_forward():
    payload = {
        "results": [
            {
                "name": "Cyanea aspera",
                "publication": "Fl. Hawaii",
                "referenceCollation": "45",
                "publicationYear": "1888",
                "id": "123-1",
            }
        ]
    }
    citation = protologue.parse_ipni_citation(payload, "Cyanea aspera", CONFIG)
    assert citation["outcome"] == ""
    assert citation["ipni_id"] == "123-1"
    assert "Fl. Hawaii" in citation["source_citation"]


# --- copyright gate ------------------------------------------------------


CITATION = {"publication": "Fl. Hawaii", "page": "45", "year": "1888"}


def _item(rights, pub_date, page_numbers="45", ocr="Corolla alba."):
    return {
        "Result": [
            {
                "ItemID": "9",
                "Rights": rights,
                "PublicationDate": pub_date,
                "Pages": [
                    {
                        "PageID": "77",
                        "PageNumbers": page_numbers,
                        "PageUrl": "https://www.biodiversitylibrary.org/page/77",
                        "OcrText": ocr,
                    }
                ],
            }
        ]
    }


def test_an_item_in_copyright_is_skipped_not_scraped():
    payload = _item("In copyright", "1975")
    assert protologue.select_public_domain_page(payload, CITATION, CONFIG)["outcome"] == (
        protologue.REJECT_IN_COPYRIGHT
    )


def test_a_recent_item_is_skipped_even_when_labelled_public_domain():
    payload = _item("Public domain", "1975")
    assert protologue.select_public_domain_page(payload, CITATION, CONFIG)["outcome"] == (
        protologue.REJECT_IN_COPYRIGHT
    )


def test_a_public_domain_page_carries_its_licence_and_url():
    page = protologue.select_public_domain_page(_item("Public domain", "1888"), CITATION, CONFIG)
    assert page["outcome"] == ""
    assert page["license"] == "Public domain"
    assert page["source_url"].endswith("/page/77")
    assert page["page_id"] == "77"


def test_a_page_that_is_not_the_cited_page_is_not_used():
    payload = _item("Public domain", "1888", page_numbers="12")
    assert protologue.select_public_domain_page(payload, CITATION, CONFIG)["outcome"] == (
        protologue.REJECT_NO_SCAN
    )


# --- end to end over a fake network --------------------------------------


def _fetch(url, params):
    if "ipni" in url:
        return {
            "results": [
                {
                    "name": params["q"],
                    "publication": "Fl. Hawaii",
                    "referenceCollation": "45",
                    "publicationYear": "1888",
                    "id": "123-1",
                }
            ]
        }
    return _item("Public domain", "1888")


def test_acquire_produces_a_candidate_with_full_lineage():
    queue = pd.DataFrame({"accepted_species": ["Cyanea aspera"]})
    audit = protologue.acquire(queue, _fetch, CONFIG, "key", retrieval_date="2026-08-23")
    row = audit.iloc[0]

    assert row["outcome"] == protologue.ACCEPTED
    assert row["normalized_value"] == "white"
    assert row["evidence_scope"] == "species_direct"
    assert row["source_type"] == "protologue"
    assert row["review_status"] == "unreviewed"
    assert not protologue.require_lineage(row.to_dict(), CONFIG)


def test_a_record_missing_lineage_is_rejected_not_downgraded():
    def fetch_without_url(url, params):
        if "ipni" in url:
            return _fetch(url, params)
        payload = _item("Public domain", "1888")
        payload["Result"][0]["Pages"][0]["PageUrl"] = ""
        payload["Result"][0]["Pages"][0].pop("ItemUrl", None)
        return payload

    queue = pd.DataFrame({"accepted_species": ["Cyanea aspera"]})
    row = protologue.acquire(queue, fetch_without_url, CONFIG, "key").iloc[0]
    assert row["outcome"] == protologue.REJECT_MISSING_LINEAGE
    assert row["normalized_value"] == ""


def test_a_species_with_no_citation_is_unresolved_not_colourless():
    def fetch_nothing(url, params):
        return {"results": []}

    queue = pd.DataFrame({"accepted_species": ["Cyanea aspera"]})
    row = protologue.acquire(queue, fetch_nothing, CONFIG, "key").iloc[0]
    assert row["outcome"] == protologue.REJECT_NO_CITATION
    assert row["normalized_value"] == ""


# --- configuration -------------------------------------------------------


def test_config_maps_only_onto_ontology_values():
    ontology = {
        "white",
        "green_brown_inconspicuous",
        "yellow_orange",
        "red_pink",
        "blue_purple",
        "multicolored_variable",
    }
    assert set(protologue.expand_colour_terms(CONFIG).values()) <= ontology
    assert CONFIG["multiple_value_result"] in ontology


def test_no_organ_is_both_floral_and_competing():
    assert not set(CONFIG["floral_organ_terms"]) & set(CONFIG["competing_organ_terms"])


def test_every_supported_language_contributes_vocabulary():
    terms = set(protologue.expand_colour_terms(CONFIG))
    assert {"white", "alba", "blanc", "weiss", "blanco"} <= terms


def test_latin_stems_decline_to_the_forms_protologues_actually_use():
    # Agreement with the organ noun decides the form, so the plural and the
    # ablative are not optional extras -- "flores albi" and "floribus albis"
    # are as common as the dictionary "albus".
    terms = protologue.expand_colour_terms(CONFIG)
    for form in ("albus", "alba", "album", "albi", "albae", "albo", "albis"):
        assert terms[form] == "white", form
    assert terms["lutei"] == "yellow_orange"
    assert terms["rubri"] == "red_pink"


def test_a_hand_listed_irregular_is_not_overwritten_by_a_generated_form():
    # "ruber" cannot be produced from the stem rubr-, so it is listed literally.
    assert protologue.expand_colour_terms(CONFIG)["ruber"] == "red_pink"


def test_third_declension_colours_decline_separately():
    terms = protologue.expand_colour_terms(CONFIG)
    for form in ("viridis", "viride", "virides", "viridia", "viridibus"):
        assert terms[form] == "green_brown_inconspicuous", form


def test_permitted_routes_only():
    assert CONFIG["citation_index"]["access"] == "official_api"
    assert CONFIG["scan_provider"]["access"] == "official_api"


def test_malformed_config_fails_closed():
    with pytest.raises(typer.BadParameter):
        protologue.load_config(Path("config/colour_schema.yml"))


# --- sentence boundaries -------------------------------------------------


def test_leaf_colour_does_not_leak_across_a_sentence_into_the_flower():
    # The 38.6% failure class: "fusca" describes the leaf, and the next
    # sentence's "Flores" is nearer than its own "Folia". Without a sentence
    # wall the record reads as variable rather than white.
    outcome, value, matched, _ = _extract(
        "Folia coriacea, subtus fusca. Flores albi, corolla tubulosa."
    )
    assert outcome == protologue.ACCEPTED
    assert value == "white"
    assert "fusca" not in matched


def test_a_semicolon_separates_organ_blocks_too():
    outcome, value, _, _ = _extract("Fructus ruber; flores albi.")
    assert outcome == protologue.ACCEPTED
    assert value == "white"


def test_a_colour_alone_in_its_sentence_has_no_organ_to_anchor_to():
    assert _extract("Flores numerosi. Rubri.")[0] == protologue.REJECT_NO_ORGAN


def test_sentence_spans_cover_the_delimited_segments():
    spans = protologue.sentence_spans("ab. cd; ef")
    assert [("ab. cd; ef")[s:e] for s, e in spans] == ["ab", " cd", " ef"]


def test_a_nearer_trailing_organ_costs_a_nuance_rather_than_causing_an_error():
    # "Corolla alba demum rosea, calyce viridi" -- the corolla is white then
    # pink, the calyx green. Distance puts "calyce" two characters from "rosea"
    # and "corolla" twelve, so the second colour is dropped. Proximity cannot
    # resolve this without modelling clause structure, and the module is built
    # to lose a nuance rather than invent an attribution: the surviving answer,
    # that the corolla is white, is true. Pinned so the trade is deliberate.
    outcome, value, matched, _ = _extract("Corolla alba demum rosea, calyce viridi.")
    assert outcome == protologue.ACCEPTED
    assert value == "white"
    assert matched == "alba"
