import importlib
from pathlib import Path

import pandas as pd
import pytest
import typer

mining = importlib.import_module("island_v2.gbif_label_trait_mining")

CONFIG = mining.load_config(Path("config/gbif_label_mining.yml"))
MASTER = {"Cyanea aspera", "Melicope wawraeana"}

OCCURRENCE_COLUMNS = [
    "gbifID",
    "species",
    "taxonRank",
    "basisOfRecord",
    "occurrenceRemarks",
    "fieldNotes",
    "dynamicProperties",
    "habitat",
    "recordedBy",
    "eventDate",
    "locality",
    "datasetKey",
    "institutionCode",
    "occurrenceID",
]


def _occ(**overrides):
    row = {c: "" for c in OCCURRENCE_COLUMNS}
    row.update(
        {
            "gbifID": "1",
            "species": "Cyanea aspera",
            "taxonRank": "SPECIES",
            "basisOfRecord": "PRESERVED_SPECIMEN",
        }
    )
    row.update(overrides)
    return pd.DataFrame([row], columns=OCCURRENCE_COLUMNS)


def _outcome(**overrides):
    return mining.mine_occurrences(_occ(**overrides), MASTER, CONFIG).iloc[0]


def test_plain_label_statement_becomes_a_candidate():
    row = _outcome(occurrenceRemarks="Shrub 2 m; flowers white, fragrant.")
    assert row["outcome"] == mining.ACCEPTED
    assert row["normalized_value"] == "white"
    assert row["evidence_scope"] == "species_direct"
    assert row["review_status"] == "unreviewed"
    assert "flowers white" in row["exact_supporting_quote"]


def test_colour_on_a_competing_organ_is_rejected():
    row = _outcome(occurrenceRemarks="Fruits red, ripe; collected on ridge.")
    assert row["outcome"] == mining.REJECT_COMPETING


def test_leaf_colour_never_becomes_flower_colour():
    row = _outcome(fieldNotes="Leaves purple beneath, glabrous above.")
    assert row["outcome"] == mining.REJECT_COMPETING


def test_calyx_colour_is_a_competing_organ_not_the_flower():
    row = _outcome(occurrenceRemarks="Calyx green, corolla lobes spreading.")
    assert row["outcome"] == mining.REJECT_COMPETING


def test_a_distant_colour_word_is_not_anchored_to_the_flower():
    far = "flowers borne on a long peduncle " + ("x" * 80) + " soil red clay"
    row = _outcome(occurrenceRemarks=far)
    assert row["outcome"] in {mining.REJECT_NO_ORGAN, mining.REJECT_COMPETING}


def test_negated_statements_are_voided():
    for note in (
        "flowers not seen",
        "Sterile, no flowers",
        "flowers ? white",
        "not in flower at time of collection",
    ):
        assert _outcome(occurrenceRemarks=note)["outcome"] == mining.REJECT_NEGATED


def test_two_colour_groups_become_variable_rather_than_the_first_hit():
    row = _outcome(occurrenceRemarks="Corolla white to pink.")
    assert row["outcome"] == mining.ACCEPTED
    assert row["normalized_value"] == "multicolored_variable"


def test_two_terms_in_the_same_group_stay_one_value():
    row = _outcome(occurrenceRemarks="Flowers pink, rose-tinged.")
    assert row["normalized_value"] == "red_pink"


def test_non_vouchered_records_are_excluded():
    row = _outcome(basisOfRecord="HUMAN_OBSERVATION", occurrenceRemarks="flowers white")
    assert row["outcome"] == mining.REJECT_BASIS


def test_records_not_determined_to_species_are_excluded():
    row = _outcome(taxonRank="GENUS", occurrenceRemarks="flowers white")
    assert row["outcome"] == mining.REJECT_RANK


def test_species_outside_the_island_master_are_excluded():
    row = _outcome(species="Zea mays", occurrenceRemarks="flowers white")
    assert row["outcome"] == mining.REJECT_OFF_MASTER


def test_label_without_colour_or_text_is_reported_distinctly():
    assert _outcome(occurrenceRemarks="Shrub 2 m tall on ridge")["outcome"] == (
        mining.REJECT_NO_COLOUR
    )
    assert _outcome()["outcome"] == mining.REJECT_NO_LABEL


def test_all_label_fields_are_read_not_only_occurrence_remarks():
    for field in ("occurrenceRemarks", "fieldNotes", "dynamicProperties", "habitat"):
        row = _outcome(**{field: "flowers yellow"})
        assert row["outcome"] == mining.ACCEPTED, field
        assert row["normalized_value"] == "yellow_orange"


def test_duplicate_sheets_of_one_gathering_collapse_to_one_event():
    rows = []
    for sheet in range(4):
        rows.append(
            {
                **{c: "" for c in OCCURRENCE_COLUMNS},
                "gbifID": str(sheet),
                "species": "Cyanea aspera",
                "taxonRank": "SPECIES",
                "basisOfRecord": "PRESERVED_SPECIMEN",
                "occurrenceRemarks": "flowers white",
                "recordedBy": "Forbes, C.N.",
                "eventDate": "1912-06-01",
                "locality": "Kauai, Wahiawa",
                "institutionCode": f"HERB{sheet}",
            }
        )
    audit = mining.mine_occurrences(pd.DataFrame(rows, columns=OCCURRENCE_COLUMNS), MASTER, CONFIG)
    events = mining.collapse_to_independent_events(audit, CONFIG)

    assert (audit["outcome"] == mining.ACCEPTED).sum() == 4
    assert len(events) == 1
    assert events.iloc[0]["n_duplicate_sheets"] == 4


def test_separate_gatherings_stay_separate():
    base = {
        **{c: "" for c in OCCURRENCE_COLUMNS},
        "species": "Cyanea aspera",
        "taxonRank": "SPECIES",
        "basisOfRecord": "PRESERVED_SPECIMEN",
        "occurrenceRemarks": "flowers white",
        "locality": "Kauai",
    }
    rows = [
        {**base, "gbifID": "1", "recordedBy": "Forbes", "eventDate": "1912-06-01"},
        {**base, "gbifID": "2", "recordedBy": "Rock", "eventDate": "1919-03-02"},
    ]
    audit = mining.mine_occurrences(pd.DataFrame(rows, columns=OCCURRENCE_COLUMNS), MASTER, CONFIG)
    assert len(mining.collapse_to_independent_events(audit, CONFIG)) == 2


def test_config_maps_only_onto_ontology_values():
    ontology = {
        "white",
        "green_brown_inconspicuous",
        "yellow_orange",
        "red_pink",
        "blue_purple",
        "multicolored_variable",
    }
    assert set(CONFIG["colour_terms"].values()) <= ontology
    assert CONFIG["multiple_value_result"] in ontology
    # a fruit must never be listed as a floral organ
    assert not set(CONFIG["floral_organ_terms"]) & set(CONFIG["competing_organ_terms"])


def test_malformed_config_fails_closed():
    with pytest.raises(typer.BadParameter):
        mining.load_config(Path("config/colour_schema.yml"))
