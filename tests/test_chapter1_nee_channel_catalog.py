from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.chapter1_nee_channel_catalog import (
    build_catalog,
    filter_interaction_chunk,
    interaction_fingerprint,
    load_config,
)


CONFIG = Path("config/chapter1_nee_channel_taxon_catalog.yml")
VISITS = "http://purl.obolibrary.org/obo/RO_0002622"
VISITED_BY = "http://purl.obolibrary.org/obo/RO_0002623"
POLLINATES = "http://purl.obolibrary.org/obo/RO_0002455"


def _config() -> dict[str, object]:
    return load_config(CONFIG)


def _row(**overrides: object) -> dict[str, object]:
    base = {column: "" for column in _config()["required_columns"]}
    base.update(
        {
            "sourceTaxonName": "Bombus terrestris",
            "sourceTaxonSpeciesName": "Bombus terrestris",
            "sourceTaxonGenusName": "Bombus",
            "sourceTaxonFamilyName": "Apidae",
            "sourceTaxonOrderName": "Hymenoptera",
            "sourceTaxonClassName": "Insecta",
            "sourceTaxonKingdomName": "Animalia",
            "interactionTypeId": VISITS,
            "interactionTypeName": "visitsFlowersOf",
            "targetTaxonName": "Campanula rotundifolia",
            "targetTaxonSpeciesName": "Campanula rotundifolia",
            "targetTaxonGenusName": "Campanula",
            "targetTaxonFamilyName": "Campanulaceae",
            "targetTaxonOrderName": "Asterales",
            "targetTaxonClassName": "Magnoliopsida",
            "targetTaxonKingdomName": "Plantae",
            "referenceCitation": "Example study A",
            "referenceUrl": "https://example.org/a",
            "sourceNamespace": "mock",
        }
    )
    base.update(overrides)
    return base


def test_bombus_precedes_non_bombus_apidae_classification() -> None:
    evidence, holdouts, removed = filter_interaction_chunk(pd.DataFrame([_row()]), _config())
    assert removed == 0
    assert holdouts.empty
    assert evidence.loc[0, "channel_id"] == "bombus"


def test_non_bombus_bee_is_classified_after_explicit_flower_visit() -> None:
    row = _row(
        sourceTaxonName="Xylocopa virginica",
        sourceTaxonSpeciesName="Xylocopa virginica",
        sourceTaxonGenusName="Xylocopa",
    )
    evidence, _, _ = filter_interaction_chunk(pd.DataFrame([row]), _config())
    assert evidence.loc[0, "channel_id"] == "non_bombus_bees"


def test_reverse_flowervisitedby_direction_extracts_pollinator_from_target() -> None:
    row = _row(
        sourceTaxonName="Salvia elegans",
        sourceTaxonSpeciesName="Salvia elegans",
        sourceTaxonGenusName="Salvia",
        sourceTaxonFamilyName="Lamiaceae",
        sourceTaxonOrderName="Lamiales",
        sourceTaxonClassName="Magnoliopsida",
        sourceTaxonKingdomName="Plantae",
        interactionTypeId=VISITED_BY,
        interactionTypeName="flowersVisitedBy",
        targetTaxonName="Selasphorus rufus",
        targetTaxonSpeciesName="Selasphorus rufus",
        targetTaxonGenusName="Selasphorus",
        targetTaxonFamilyName="Trochilidae",
        targetTaxonOrderName="Apodiformes",
        targetTaxonClassName="Aves",
        targetTaxonKingdomName="Animalia",
    )
    evidence, holdouts, _ = filter_interaction_chunk(pd.DataFrame([row]), _config())
    assert holdouts.empty
    assert evidence.loc[0, "channel_id"] == "flower_visiting_birds"
    assert evidence.loc[0, "pollinator_species"] == "Selasphorus rufus"
    assert evidence.loc[0, "plant_taxon"] == "Salvia elegans"


def test_lepidoptera_and_diptera_require_explicit_flower_interaction_first() -> None:
    butterfly = _row(
        sourceTaxonName="Danaus plexippus",
        sourceTaxonSpeciesName="Danaus plexippus",
        sourceTaxonGenusName="Danaus",
        sourceTaxonFamilyName="Nymphalidae",
        sourceTaxonOrderName="Lepidoptera",
    )
    fly = _row(
        sourceTaxonName="Eristalis tenax",
        sourceTaxonSpeciesName="Eristalis tenax",
        sourceTaxonGenusName="Eristalis",
        sourceTaxonFamilyName="Syrphidae",
        sourceTaxonOrderName="Diptera",
        referenceCitation="Example study B",
    )
    evidence, _, _ = filter_interaction_chunk(pd.DataFrame([butterfly, fly]), _config())
    assert set(evidence["channel_id"]) == {"lepidoptera", "diptera"}


def test_nonflower_interaction_is_never_used_to_define_channel_catalog() -> None:
    row = _row(
        interactionTypeId="http://purl.obolibrary.org/obo/RO_0002434",
        interactionTypeName="interactsWith",
    )
    evidence, holdouts, removed = filter_interaction_chunk(pd.DataFrame([row]), _config())
    assert evidence.empty
    assert holdouts.empty
    assert removed == 0


def test_nonplant_flower_claim_is_held_out() -> None:
    row = _row(targetTaxonKingdomName="Fungi")
    evidence, holdouts, _ = filter_interaction_chunk(pd.DataFrame([row]), _config())
    assert evidence.empty
    assert holdouts.loc[0, "reason"] == "plant_side_not_resolved_as_Plantae"


def test_species_unresolved_interaction_is_holdout_not_genus_imputation() -> None:
    row = _row(sourceTaxonSpeciesName="", sourceTaxonName="Bombus", sourceTaxonGenusName="Bombus")
    evidence, holdouts, _ = filter_interaction_chunk(pd.DataFrame([row]), _config())
    assert evidence.empty
    assert holdouts.loc[0, "reason"] == "pollinator_species_unresolved"


def test_refuted_interaction_is_subtracted_before_cataloguing() -> None:
    row = _row()
    fingerprint = interaction_fingerprint(row)
    evidence, holdouts, removed = filter_interaction_chunk(
        pd.DataFrame([row]), _config(), {fingerprint}
    )
    assert evidence.empty
    assert holdouts.empty
    assert removed == 1


def test_single_pollination_claim_qualifies_confirmatory_catalog() -> None:
    row = _row(
        interactionTypeId=POLLINATES,
        interactionTypeName="pollinates",
        referenceCitation="Pollination effectiveness study",
    )
    evidence, _, _ = filter_interaction_chunk(pd.DataFrame([row]), _config())
    catalog = build_catalog(
        evidence,
        _config(),
        source_version_doi="10.5281/zenodo.9999999",
        source_sha256="abc",
    )
    assert catalog.loc[0, "n_pollination_rows"] == 1
    assert catalog.loc[0, "catalog_tier"] == "confirmatory"


def test_single_flower_visit_reference_is_sensitivity_only() -> None:
    evidence, _, _ = filter_interaction_chunk(pd.DataFrame([_row()]), _config())
    catalog = build_catalog(
        evidence,
        _config(),
        source_version_doi="10.5281/zenodo.9999999",
        source_sha256="abc",
    )
    assert catalog.loc[0, "catalog_tier"] == "sensitivity"


def test_two_independent_visit_references_across_two_plants_qualify_confirmatory() -> None:
    rows = [
        _row(referenceCitation="Study A", targetTaxonSpeciesName="Plant alpha", targetTaxonName="Plant alpha"),
        _row(referenceCitation="Study B", targetTaxonSpeciesName="Plant beta", targetTaxonName="Plant beta"),
    ]
    evidence, _, _ = filter_interaction_chunk(pd.DataFrame(rows), _config())
    catalog = build_catalog(
        evidence,
        _config(),
        source_version_doi="10.5281/zenodo.9999999",
        source_sha256="abc",
    )
    assert catalog.loc[0, "n_independent_references"] == 2
    assert catalog.loc[0, "n_distinct_plant_taxa"] == 2
    assert catalog.loc[0, "catalog_tier"] == "confirmatory"
