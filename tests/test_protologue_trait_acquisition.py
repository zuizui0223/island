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


def test_japanese_voicing_survives_folding():
    """Dakuten is not an accent, so stripping it is not folding but corruption.

    Discarding the mark the way an acute is discarded turns ピンク into ヒンク and
    collapses ハラ, バラ and パラ onto one form. Nothing was mis-valued by it --
    the mangled terms still matched each other -- but the audit ledger recorded
    `matched_terms` that do not appear on the page the reviewer is reading, and
    any future kana pair separated only by voicing would silently merge.
    """
    for text in ("ピンク色", "オレンジ色", "バラ", "パラ", "ハラ"):
        assert protologue.fold(text)[0] == text

    # Halfwidth kana, which is how a good deal of scanned Japanese arrives,
    # folds onto the precomposed form rather than staying a separate spelling.
    assert protologue.fold("ﾋﾟﾝｸ")[0] == "ピンク"


def test_the_matched_terms_appear_in_the_quote_they_came_from():
    outcome, _, matched, quote = _extract("花はピンク色")
    assert outcome == protologue.ACCEPTED
    for term in matched.split("+"):
        assert term in quote, (term, quote)


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


# --- scripts beyond Latin -----------------------------------------------


@pytest.mark.parametrize(
    ("description", "value"),
    [
        ("花は白色", "white"),
        ("花は淡黄色。", "yellow_orange"),
        ("初夏に紫色の花をたくさん枝先に付ける。", "blue_purple"),
        ("Flores brancas, perfumadas.", "white"),
        ("Corola amarela, tubo curto.", "yellow_orange"),
        ("Цветки белые, душистые.", "white"),
        ("Венчик синий.", "blue_purple"),
    ],
)
def test_colour_is_read_in_the_added_languages(description, value):
    outcome, normalized, _, _ = _extract(description)
    assert outcome == protologue.ACCEPTED, description
    assert normalized == value


def test_a_cyrillic_term_does_not_match_inside_a_longer_word():
    # The Latin-only boundary class left Cyrillic unguarded, so "белый" matched
    # inside "белыйцветок" -- the same defect as reading "rot" in
    # "rotundifolia", which the boundary rule exists to prevent.
    folded, _ = protologue.fold("белыйцветок")
    assert protologue._boundary_spans(folded, protologue.fold("белый")[0]) == []
    folded, _ = protologue.fold("цветки белые")
    assert protologue._boundary_spans(folded, protologue.fold("белые")[0])


def test_a_spaceless_script_still_matches_flush_against_its_neighbours():
    # Japanese and Chinese do not separate words, so a boundary guard would
    # reject every real match.
    folded, _ = protologue.fold("花は白色")
    assert protologue._boundary_spans(folded, "白色") == [(2, 4)]
    assert protologue._boundary_spans(folded, "花") == [(0, 1)]


def test_a_japanese_full_stop_walls_the_fruit_off_from_the_flower():
    outcome, value, matched, _ = _extract("花は白色。果実は赤い。")
    assert outcome == protologue.ACCEPTED
    assert value == "white"
    assert "赤" not in matched


# --- elision and abbreviation -------------------------------------------


def test_an_ellipsis_is_elision_not_a_full_stop():
    # A reviewer quoting around omitted text still means one statement. Split on
    # the dots, "albida" loses "petala" and the row is rejected as unanchored.
    outcome, value, _, _ = _extract("petala ... valde conspicue glandulosi-punctata albida")
    assert outcome == protologue.ACCEPTED
    assert value == "white"


def test_an_abbreviations_period_does_not_end_the_sentence():
    outcome, value, _, _ = _extract(
        "Origanum Vetteri ... Karpathos ... fl. pink, 21 July 1950, Davis 18005"
    )
    assert outcome == protologue.ACCEPTED
    assert value == "red_pink"


def test_masking_a_period_does_not_move_any_offset():
    # The spans feed back into the caller's text, so both substitutions must
    # preserve length exactly.
    for text in ("petala ... albida", "fl. pink", "a … b", "sp. nov. fl. alba"):
        masked = protologue._mask_non_terminal_periods(
            protologue.fold(text)[0], protologue.DEFAULT_ABBREVIATIONS
        )
        assert len(masked) == len(protologue.fold(text)[0]), text


# --- compound hue versus genuine variability -----------------------------


def test_a_compound_hue_is_not_a_corolla_that_changes_colour():
    # "pale reddish purple" is one hue written with a modifier. Calling it
    # variable asserts a change the description does not describe.
    outcome, value, _, _ = _extract("Corolla pale reddish purple, 17-21 mm long")
    assert outcome == protologue.ACCEPTED
    assert value == "ontology_unresolved"


def test_a_japanese_modifier_between_the_words_still_reads_as_one_hue():
    # Japanese puts its intensity word between the colours, where English puts
    # it in front: 紫がかった濃いピンク色 is "deep pink with a purple tinge".
    outcome, value, _, _ = _extract("花は紫がかった濃いピンク色")
    assert outcome == protologue.ACCEPTED
    assert value == "ontology_unresolved"


def test_a_range_marker_still_means_genuine_variability():
    for description in (
        "Corolla alba demum rosea.",
        "Perianth greenish yellow or yellowish green",
        # A bare "to" is a range across a population, not a gradient on one
        # flower: "white shading to pink" is one hue, "white to pink" is two.
        "Corolla white to pink.",
    ):
        outcome, value, _, _ = _extract(description)
        assert outcome == protologue.ACCEPTED
        assert value == "multicolored_variable", description


# --- the binary class the model consumes ---------------------------------


def test_the_plain_class_matches_the_analysis_definition():
    # Declared in config so acquisition stays free of the analysis dependency
    # chain; this is what stops the two drifting apart.
    import importlib

    analysis = importlib.import_module("island_v2.status_category_decomposition")
    assert protologue.plain_colour_values(CONFIG) == set(
        analysis.CATEGORY_SPECS["colour"]["categories"]["plain"]
    )


@pytest.mark.parametrize(
    ("description", "binary"),
    [
        # The six judgements the reviewer recorded in the wave-24 ledger.
        ("flowers yellowish-cream with greenish tint and pink-brown central veins", "unresolved"),
        ("花は淡黄色。", "nonplain"),
        ("purely yellowish-green flowers", "unresolved"),
        ("Corolla pale reddish purple", "nonplain"),
        ("Flowers bluish-violet", "nonplain"),
        ("Perianth greenish yellow or yellowish green", "unresolved"),
    ],
)
def test_the_binary_class_reproduces_the_reviewers_judgements(description, binary):
    _, _, matched, _ = _extract(description)
    assert protologue.binary_plain_class(matched, CONFIG) == binary


def test_a_compound_hue_keeps_a_usable_binary_class():
    # The whole point of the second level: the ontology cannot decide red-pink
    # versus blue-purple, but both sit on the non-plain side of the line, so the
    # value the model consumes survives.
    _, value, matched, _ = _extract("Corolla pale reddish purple")
    assert value == "ontology_unresolved"
    assert protologue.binary_plain_class(matched, CONFIG) == "nonplain"


def test_nothing_matched_leaves_the_binary_class_unresolved():
    assert protologue.binary_plain_class("", CONFIG) == "unresolved"


# --- scored against the frozen review ledger -----------------------------

REVIEW_LEDGER = Path(
    "data/v2/curation/minimum_endemic_colour_manual_review_20260823.csv"
)

# Two rows disagree with this ledger by design rather than by defect. Both are
# compound hues straddling red-pink and blue-purple, and the reviewer's later
# frozen call on exactly that construction -- "Corolla pale reddish purple" in
# the wave-24 ledger -- was that the five-value ontology cannot decide it while
# the binary class certainly can. Resolving them would mean picking a hue the
# reviewer declined to pick.
COMPOUND_HUE_ROWS = {"Cirsium umezawanum", "Cirsium tanegashimense"}

# The score this file reached when the added languages, elision handling and
# compound-hue rule landed. It went 5/12 -> 10/12; the floor guards the gain.
LEDGER_FLOOR = 10


def _score_ledger():
    import csv

    rows = list(csv.DictReader(REVIEW_LEDGER.open(encoding="utf-8")))
    correct, wrong = [], []
    for row in rows:
        outcome, value, _, _ = _extract(row["source_excerpt"])
        got = value if outcome == protologue.ACCEPTED else outcome
        (correct if got == row["normalized_value"].strip() else wrong).append(
            (row["accepted_species"], row["normalized_value"].strip(), got)
        )
    return rows, correct, wrong


@pytest.mark.skipif(not REVIEW_LEDGER.exists(), reason="review ledger not present")
def test_the_extractor_does_not_regress_against_reviewed_evidence():
    rows, correct, wrong = _score_ledger()
    assert len(correct) >= LEDGER_FLOOR, (
        f"scored {len(correct)}/{len(rows)}, below the {LEDGER_FLOOR} floor. "
        f"Missed: {wrong}"
    )


@pytest.mark.skipif(not REVIEW_LEDGER.exists(), reason="review ledger not present")
def test_every_remaining_disagreement_is_a_compound_hue():
    # Pins why the score is not perfect, so a new kind of failure cannot hide
    # behind the two rows that are deliberate.
    _, _, wrong = _score_ledger()
    for species, _, got in wrong:
        assert species in COMPOUND_HUE_ROWS, f"unexpected miss: {species} -> {got}"
        assert got == "ontology_unresolved", got


@pytest.mark.skipif(not REVIEW_LEDGER.exists(), reason="review ledger not present")
def test_the_binary_class_resolves_for_every_reviewed_row_with_a_colour():
    # The ontology fails on two rows; the level the model consumes fails on none.
    import csv

    for row in csv.DictReader(REVIEW_LEDGER.open(encoding="utf-8")):
        outcome, _, matched, _ = _extract(row["source_excerpt"])
        if outcome != protologue.ACCEPTED:
            continue
        assert protologue.binary_plain_class(matched, CONFIG) in {"plain", "nonplain"}, (
            row["accepted_species"]
        )
