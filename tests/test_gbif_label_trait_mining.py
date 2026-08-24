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


# --- the rules this lane used to be missing ------------------------------
#
# The four cases below are measured, not imagined. Each was run through both
# lanes at the head where the label miner still had its own matching code, and
# each is a sentence the two lanes answered differently -- with the label lane
# giving the wrong answer in the first two.


@pytest.mark.parametrize(
    ("label", "was", "reason"),
    [
        (
            "Flowers on reduced peduncles",
            "red_pink",
            'matched "red" inside "reduced"',
        ),
        (
            "Leaves coriaceous, beneath brown. Flowers seen.",
            "green_brown_inconspicuous",
            "a leaf clause reached across a full stop",
        ),
    ],
)
def test_the_measured_false_positives_no_longer_produce_a_candidate(label, was, reason):
    row = _outcome(occurrenceRemarks=label)
    assert row["outcome"] != mining.ACCEPTED, reason
    assert row["normalized_value"] != was


def test_a_colour_word_inside_a_longer_word_is_not_a_colour():
    # The whole class the boundary guard exists for, not just the one instance
    # that was found: a substring match turns ordinary botanical prose into
    # evidence.
    for label in (
        "Flowers on reduced peduncles",
        "Leaves rotundifolia, flowers seen",
        "Corolla lobes greenishly veined nowhere",
    ):
        assert _outcome(occurrenceRemarks=label)["outcome"] != mining.ACCEPTED, label


def test_an_organ_may_not_reach_across_a_full_stop():
    row = _outcome(occurrenceRemarks="Leaves brown beneath. Flowers seen, fragrant.")
    assert row["outcome"] == mining.REJECT_COMPETING

    # ...but the wall must not fall inside a collector's abbreviation, which is
    # most of what a label is made of.
    row = _outcome(occurrenceRemarks="Shrub 2 m; infl. white, alt. 300 m")
    assert row["outcome"] == mining.ACCEPTED
    assert row["normalized_value"] == "white"


def test_labels_that_are_not_in_english_are_read():
    # A third of the first frozen review ledger was Japanese sources an
    # English-only vocabulary read as having no colour at all. Labels are not
    # exempt from that: a Kyoto or Paris sheet carries its label in its own
    # language.
    for label, value in (
        ("花は白色", "white"),
        ("Fleurs blanchâtres", "white"),
        ("Blüten weiß", "white"),
        ("цветки белые", "white"),
    ):
        row = _outcome(occurrenceRemarks=label)
        assert row["outcome"] == mining.ACCEPTED, label
        assert row["normalized_value"] == value, label


# --- the binary class the model consumes ---------------------------------


def test_the_ledger_carries_the_binary_class_the_model_reads():
    row = _outcome(occurrenceRemarks="Flowers white, fragrant.")
    assert row["binary_plain_class"] == "plain"
    assert _outcome(occurrenceRemarks="Flowers purple.")["binary_plain_class"] == (
        "nonplain"
    )
    # No colour statement means no class, not a plain flower.
    assert _outcome(occurrenceRemarks="Shrub 2 m tall")["binary_plain_class"] == (
        "unresolved"
    )


def test_a_compound_hue_the_ontology_cannot_name_still_has_a_binary_class():
    row = _outcome(occurrenceRemarks="Corolla pale reddish purple")
    assert row["outcome"] == mining.ACCEPTED
    assert row["normalized_value"] == "ontology_unresolved"
    assert row["binary_plain_class"] == "nonplain"


def test_a_range_across_the_line_is_variable_not_one_unnameable_hue():
    # "white shading to pink" is one gradient; "white to pink" is a population
    # producing both, which is variability rather than a hue the five-value
    # ontology merely cannot name.
    assert _outcome(occurrenceRemarks="Corolla white shading to pink.")[
        "normalized_value"
    ] == "ontology_unresolved"
    assert _outcome(occurrenceRemarks="Corolla white to pink.")[
        "normalized_value"
    ] == "multicolored_variable"


def test_the_statement_and_the_whole_label_are_both_retained():
    # A reviewer needs the sentence to judge the claim and the label to judge
    # the specimen, so neither replaces the other.
    row = _outcome(
        occurrenceRemarks="Shrub 2 m on ridge; flowers white, fragrant.",
        habitat="Wet forest",
    )
    assert row["exact_supporting_quote"] == "flowers white, fragrant."
    assert "Wet forest" in row["label_text"]


# --- the two lanes may not drift apart again -----------------------------


def test_both_lanes_answer_the_same_sentence_the_same_way():
    """The guard against the drift this module was rewritten to end.

    Protologue prose and specimen labels are different sources, but they are
    read by the same reviewer into the same ledger, so a sentence appearing in
    both must not become evidence in one lane and nothing in the other.
    """
    protologue = importlib.import_module("island_v2.protologue_trait_acquisition")
    protologue_config = protologue.load_config(
        Path("config/protologue_acquisition.yml")
    )

    for text in (
        "Flowers white, fragrant.",
        "Flowers on reduced peduncles",
        "Leaves coriaceous, beneath brown. Flowers seen.",
        "Corolla pale reddish purple",
        "Corolla white to pink.",
        "花は白色",
        "Fleurs blanchâtres",
        "Blüten weiß",
        "fl. pink",
        "Fruits red, ripe.",
    ):
        label = mining.extract_colour(text, CONFIG)
        prose = protologue.extract_floral_colour(text, protologue_config)
        assert label[:3] == prose[:3], text
        assert mining.binary_plain_class(label[2], CONFIG) == (
            protologue.binary_plain_class(prose[2], protologue_config)
        ), text


def test_the_lanes_differ_only_where_the_source_genuinely_differs():
    """What each lane is still allowed to decide for itself.

    Everything else -- the colour vocabulary, the organ lists, the abbreviations
    and the ontology results -- comes from the shared file, so this test fails
    the moment a rule is quietly forked back into one lane.
    """
    protologue = importlib.import_module("island_v2.protologue_trait_acquisition")
    other = protologue.load_config(Path("config/protologue_acquisition.yml"))

    differing = {key for key in set(CONFIG) | set(other) if CONFIG.get(key) != other.get(key)}
    assert differing == {
        # How the source is reached, which is not a matching rule at all.
        "label_fields",
        "accepted_basis_of_record",
        "required_taxon_rank",
        "independence_key",
        "citation_index",
        "scan_provider",
        "target_selection",
        "required_lineage_fields",
        # Protologue prose is denser than a label written in the field.
        "organ_proximity_chars",
        # A collector's "sterile" means there was nothing to observe; a
        # protologue's "flores steriles" describes the floret.
        "negation_markers",
    }


def test_both_lanes_resolve_the_same_vocabulary():
    protologue = importlib.import_module("island_v2.protologue_trait_acquisition")
    other = protologue.load_config(Path("config/protologue_acquisition.yml"))
    for key in (
        "colour_terms",
        "floral_organ_terms",
        "competing_organ_terms",
        "latin_adjective_stems",
        "abbreviations",
        "compound_gap_terms",
        "plain_colour_values",
        "multiple_value_result",
        "compound_hue_result",
    ):
        assert CONFIG[key] == other[key], key


def test_an_organ_is_never_both_floral_and_competing_after_inheritance():
    assert not set(CONFIG["floral_organ_terms"]) & set(CONFIG["competing_organ_terms"])
