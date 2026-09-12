from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_nee_channel_catalog_supported import run_supported_catalog


def _row(config: dict) -> dict[str, str]:
    row = {column: "" for column in config["required_columns"]}
    row.update(
        {
            "sourceTaxonName": "Bombus terrestris",
            "sourceTaxonSpeciesName": "Bombus terrestris",
            "sourceTaxonGenusName": "Bombus",
            "sourceTaxonFamilyName": "Apidae",
            "sourceTaxonOrderName": "Hymenoptera",
            "sourceTaxonClassName": "Insecta",
            "sourceTaxonKingdomName": "Animalia",
            "interactionTypeId": "http://purl.obolibrary.org/obo/RO_0002455",
            "interactionTypeName": "pollinates",
            "targetTaxonName": "Campanula rotundifolia",
            "targetTaxonSpeciesName": "Campanula rotundifolia",
            "targetTaxonGenusName": "Campanula",
            "targetTaxonFamilyName": "Campanulaceae",
            "targetTaxonOrderName": "Asterales",
            "targetTaxonClassName": "Magnoliopsida",
            "targetTaxonKingdomName": "Plantae",
            "referenceCitation": "Synthetic supported claim",
            "sourceNamespace": "synthetic",
        }
    )
    return row


def test_provider_semantics_are_support_refute_split() -> None:
    config = yaml.safe_load(Path("config/chapter1_nee_channel_taxon_catalog.yml").read_text())
    policy = config["source_policy"]
    assert policy["interactions_export_argument_type"] == "SUPPORTS"
    assert policy["refuted_export_argument_type"] == "REFUTES"
    assert policy["refuted_interactions_must_be_subtracted"] is False
    assert policy["refuted_interactions_role"] == "contradiction_audit_only"


def test_same_pair_in_refutes_does_not_delete_supported_argument(tmp_path: Path) -> None:
    config = yaml.safe_load(Path("config/chapter1_nee_channel_taxon_catalog.yml").read_text())
    supported = tmp_path / "interactions.tsv.gz"
    refuted = tmp_path / "refuted-interactions.tsv.gz"
    pd.DataFrame([_row(config)]).to_csv(supported, sep="\t", index=False, compression="gzip")
    # Deliberately put the identical synthetic row in the REFUTES product. Provider
    # argument semantics say this is contradictory evidence, not a row to subtract
    # from the SUPPORTS export.
    pd.DataFrame([_row(config)]).to_csv(refuted, sep="\t", index=False, compression="gzip")

    out = tmp_path / "out"
    receipt = run_supported_catalog(
        interactions_tsv_gz=supported,
        refuted_tsv_gz=refuted,
        source_version_doi="synthetic",
        output_dir=out,
        chunksize=10,
        config_path=Path("config/chapter1_nee_channel_taxon_catalog.yml"),
        skip_digest_check=True,
    )

    catalog = pd.read_csv(out / "channel_taxon_catalog.csv")
    confirmatory = pd.read_csv(out / "channel_taxon_catalog_confirmatory.csv")
    assert len(catalog) == 1
    assert len(confirmatory) == 1
    assert catalog.iloc[0]["pollinator_species"] == "Bombus terrestris"
    assert receipt["support_refute_cross_subtraction_applied"] is False
    assert receipt["n_refuted_flower_interactions_removed"] == 0
    assert receipt["n_confirmatory_catalog_taxa"] == 1
