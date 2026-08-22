import importlib
from pathlib import Path

import pandas as pd
import pytest
import typer

recovery = importlib.import_module("island_v2.bulk_recovery_audit")

CONFIG = recovery.load_config(Path("config/bulk_recovery.yml"))


def _audit(rows):
    return pd.DataFrame(
        rows,
        columns=[
            "source_scientific_name",
            "match_status",
            "gbif_match_type",
            "gbif_status",
            "reason",
            "gbif_accepted_usage_key",
        ],
    )


def test_confident_gbif_rejection_with_a_key_is_recoverable():
    audit = _audit(
        [("Aus bus L.", "unmatched", "EXACT", "ACCEPTED", "reject_not_exact_synonym", "123")]
    )
    result = recovery.classify_name_match(audit, CONFIG)

    assert result.iloc[0]["recovery_class"] == recovery.RECOVERABLE_KEY
    assert result.iloc[0]["gbif_key"] == "123"
    assert result.iloc[0]["recovered_tier"] == "gbif_usage_key_resolved"


def test_matched_rows_are_not_part_of_the_recovery_population():
    audit = _audit(
        [("Aus bus", "matched", "EXACT", "ACCEPTED", "accepted_name_exact", "123")]
    )
    assert recovery.classify_name_match(audit, CONFIG).empty


def test_family_conflict_and_non_species_rank_stay_finally_rejected():
    audit = _audit(
        [
            ("Aus bus", "unmatched", "EXACT", "ACCEPTED", "reject_family_conflict", "1"),
            ("Aus bus", "unmatched", "EXACT", "ACCEPTED", "accepted_name_family_conflict", "2"),
            ("Aus sp.", "unmatched", "EXACT", "ACCEPTED", "reject_non_species_rank", "3"),
        ]
    )
    result = recovery.classify_name_match(audit, CONFIG)
    assert set(result["recovery_class"]) == {recovery.BLOCKED_REASON}


def test_fuzzy_and_higherrank_matches_are_never_recovered():
    audit = _audit(
        [
            ("Aus bus", "unmatched", "FUZZY", "ACCEPTED", "reject_non_exact_match", "1"),
            ("Aus", "unmatched", "HIGHERRANK", "ACCEPTED", "reject_non_exact_match", "2"),
            ("Aus bus", "unmatched", "EXACT", "DOUBTFUL", "reject_not_exact_synonym", "3"),
        ]
    )
    result = recovery.classify_name_match(audit, CONFIG)
    assert set(result["recovery_class"]) == {recovery.BLOCKED_GBIF}


def test_confident_rejection_without_a_key_cannot_be_resolved():
    audit = _audit(
        [("Aus bus", "unmatched", "EXACT", "ACCEPTED", "reject_not_exact_synonym", "")]
    )
    result = recovery.classify_name_match(audit, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.BLOCKED_NO_KEY


def test_queue_deduplicates_keys_and_counts_rows():
    audit = _audit(
        [
            ("Aus bus", "unmatched", "EXACT", "ACCEPTED", "reject_not_exact_synonym", "7"),
            ("Aus bus var. x", "unmatched", "EXACT", "ACCEPTED", "reject_not_exact_synonym", "7"),
            ("Cus dus", "unmatched", "EXACT", "ACCEPTED", "reject_not_exact_synonym", "9"),
        ]
    )
    queue = recovery.name_match_queue(recovery.classify_name_match(audit, CONFIG))

    assert queue["gbif_key"].tolist() == ["7", "9"]
    assert queue.loc[queue["gbif_key"] == "7", "n_rows"].item() == 2


def _rejected(rows):
    return pd.DataFrame(
        rows,
        columns=["trait_name", "exact_supporting_quote", "context_gate_reason", "language"],
    )


def test_measurement_of_a_non_target_organ_stays_rejected():
    rejected = _rejected(
        [
            (
                "flower_size_class",
                "the style 3.7-7 mm long, glabrous, the stigma branches 2.4-3.4 mm",
                "measurement_not_directly_attached_to_flower_subject",
                "en",
            )
        ]
    )
    result = recovery.classify_extraction(rejected, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.CORRECT_ORGAN


def test_enumeration_template_stays_rejected():
    rejected = _rejected(
        [
            (
                "inflorescence_display",
                "Inflorescencia: cima(s) número de flor(es) pauciflora(s)/multiflora(s);",
                "pauciflora_taxonomic_name_or_citation_only",
                "es",
            )
        ]
    )
    result = recovery.classify_extraction(rejected, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.CORRECT_TEMPLATE


def test_explicit_non_english_statement_is_recoverable():
    rejected = _rejected(
        [
            (
                "inflorescence_display",
                "Blüten klein, weiss, in dichter, end- oder blattwinkelständiger Traube .",
                "inflorescence_state_not_explicit_in_description",
                "de",
            )
        ]
    )
    result = recovery.classify_extraction(rejected, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.RECOVERABLE_LANG
    assert result.iloc[0]["matched_multilingual_term"] == "traube->raceme"


def test_recovery_needs_the_term_present_not_just_a_language_tag():
    rejected = _rejected(
        [("inflorescence_display", "eine unklare Beschreibung ohne Fachbegriff", "x", "de")]
    )
    result = recovery.classify_extraction(rejected, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.NEEDS_REVIEW


def test_a_french_term_is_recovered_even_when_the_language_tag_is_missing():
    # 37 of the real rejected rows carry an empty language tag; the quote still
    # states the trait plainly and must not be lost to the missing tag.
    rejected = _rejected(
        [
            (
                "inflorescence_display",
                "les supérieurs forment une panicule;",
                "inflorescence_state_not_explicit_in_description",
                "",
            )
        ]
    )
    result = recovery.classify_extraction(rejected, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.RECOVERABLE_LANG
    assert result.iloc[0]["language_code"] == "unknown"


def test_non_target_organ_wins_over_a_multilingual_term():
    # "grappe" is a French raceme term, but the quote attributes it to the fruit
    rejected = _rejected(
        [("inflorescence_display", "fruit en grappe dense", "x", "fr")]
    )
    result = recovery.classify_extraction(rejected, CONFIG)
    assert result.iloc[0]["recovery_class"] == recovery.CORRECT_ORGAN


def test_language_normalization_never_guesses_english():
    aliases = CONFIG["language_aliases"]
    assert recovery.normalize_language("English", aliases) == "en"
    assert recovery.normalize_language("EN", aliases) == "en"
    assert recovery.normalize_language("fre", aliases) == "fr"
    assert recovery.normalize_language("Spanish", aliases) == "es"
    assert recovery.normalize_language("", aliases) == "unknown"
    assert recovery.normalize_language("gibberish", aliases) == "unknown"


def test_malformed_inputs_fail_closed():
    with pytest.raises(typer.BadParameter):
        recovery.classify_name_match(pd.DataFrame({"match_status": ["unmatched"]}), CONFIG)
    with pytest.raises(typer.BadParameter):
        recovery.classify_extraction(pd.DataFrame({"trait_name": ["x"]}), CONFIG)
    with pytest.raises(typer.BadParameter):
        recovery.load_config(Path("config/island_weighted_acquisition.yml"))


MASTER = {"Aus bus", "Cus dus"}


def _queue(rows):
    return pd.DataFrame(rows, columns=["gbif_key", "n_rows", "example_source_name"])


def _record(species, rank="SPECIES", status="ACCEPTED"):
    return {"species": species, "rank": rank, "taxonomicStatus": status}


def test_resolver_recovers_only_accepted_species_present_in_the_master():
    queue = _queue([("1", 3, "Aus bus L.")])
    result = recovery.resolve_name_match_queue(
        queue, MASTER, lambda key: _record("Aus bus")
    )
    row = result.iloc[0]
    assert row["outcome"] == "recovered"
    assert row["resolved_accepted_species"] == "Aus bus"
    assert bool(row["on_island_master"])
    assert row["n_rows"] == 3


def test_resolver_drops_names_outside_the_island_universe():
    queue = _queue([("1", 1, "Xus yus")])
    result = recovery.resolve_name_match_queue(
        queue, MASTER, lambda key: _record("Xus yus")
    )
    assert result.iloc[0]["outcome"] == "rejected_off_island_master"
    assert not bool(result.iloc[0]["on_island_master"])


def test_resolver_fails_closed_on_unresolved_non_species_and_non_accepted():
    queue = _queue([("1", 1, "a"), ("2", 1, "b"), ("3", 1, "c")])
    lookups = {
        "1": None,
        "2": _record("Aus bus", rank="GENUS"),
        "3": _record("Aus bus", status="SYNONYM"),
    }
    result = recovery.resolve_name_match_queue(queue, MASTER, lambda key: lookups[key])
    assert result["outcome"].tolist() == [
        "unresolved_key",
        "rejected_non_species_rank",
        "rejected_not_accepted",
    ]
    assert not result["on_island_master"].any()


def test_resolver_falls_back_to_canonical_name_when_species_field_is_absent():
    queue = _queue([("1", 1, "a")])
    result = recovery.resolve_name_match_queue(
        queue,
        MASTER,
        lambda key: {"canonicalName": "Cus dus", "rank": "SPECIES", "taxonomicStatus": "ACCEPTED"},
    )
    assert result.iloc[0]["outcome"] == "recovered"
    assert result.iloc[0]["resolved_accepted_species"] == "Cus dus"


def test_resolver_on_an_empty_queue_returns_an_empty_frame_with_the_contract_columns():
    result = recovery.resolve_name_match_queue(_queue([]), MASTER, lambda key: None)
    assert result.empty
    assert "outcome" in result.columns
