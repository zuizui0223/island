from __future__ import annotations

import pandas as pd
import yaml

from island_v2.trait_candidate_logic import audit_table


def _row(**updates: str) -> dict[str, str]:
    row = {
        "trait_name": "self_incompatibility",
        "standardized_value": "SC",
        "candidate_kind": "source_backed",
        "evidence_scope": "species_direct",
        "source_type": "institutional_web",
        "source_reliability": "D_unvetted_web",
        "source_citation": "Example species page",
        "source_url": "https://example.org/species",
        "raw_value": "self-compatible",
        "evidence_excerpt": "The focal species is self-compatible.",
        "supporting_taxa": "[]",
        "inference_rule": "",
        "inference_note": "Species page; unreviewed source tier.",
        "confidence": "low",
    }
    row.update(updates)
    return row


def _ontology() -> dict:
    return yaml.safe_load(open("config/trait_ontology.yml", encoding="utf-8"))


def test_low_reliability_source_passes_logic_audit() -> None:
    audited = audit_table(pd.DataFrame([_row()]), _ontology())
    assert audited.loc[0, "logic_status"] == "pass"
    assert audited.loc[0, "source_reliability"] == "D_unvetted_web"


def test_autonomous_selfing_cannot_be_coded_from_sc_alone() -> None:
    audited = audit_table(
        pd.DataFrame(
            [
                _row(
                    trait_name="autonomous_selfing_capacity",
                    standardized_value="autonomous",
                    raw_value="self-compatible",
                    evidence_excerpt="The focal species is self-compatible.",
                )
            ]
        ),
        _ontology(),
    )
    assert audited.loc[0, "logic_status"] == "quarantine"
    assert "autonomous_selfing_inferred_from_compatibility_only" in audited.loc[0, "logic_errors_json"]


def test_hierarchical_inference_requires_scope_support_and_rule() -> None:
    audited = audit_table(
        pd.DataFrame(
            [
                _row(
                    trait_name="flower_primary_color",
                    standardized_value="red_pink",
                    candidate_kind="hierarchical_inference",
                    evidence_scope="genus_inference",
                    source_type="flora_or_monograph",
                    source_reliability="A_primary_or_monograph",
                    source_citation="Genus flora",
                    supporting_taxa='["Species alpha", "Species beta"]',
                    inference_rule="Two documented congeners share red flowers.",
                    confidence="low",
                )
            ]
        ),
        _ontology(),
    )
    assert audited.loc[0, "logic_status"] == "pass"


def test_pollination_guild_requires_pollination_evidence_not_flower_shape() -> None:
    audited = audit_table(
        pd.DataFrame(
            [
                _row(
                    trait_name="pollination_functional_guild",
                    standardized_value="birds",
                    raw_value="red tubular flowers",
                    evidence_excerpt="Flowers are red and tubular.",
                )
            ]
        ),
        _ontology(),
    )
    assert audited.loc[0, "logic_status"] == "quarantine"
    assert "pollination_guild_without_pollination_evidence" in audited.loc[0, "logic_errors_json"]
