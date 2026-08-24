from __future__ import annotations

import pandas as pd

from island_v2.bsdb_morphology_checkpoint import (
    _trait_source_lineage_key,
    build_checkpoint,
)


def _row(species: str, **updates: str) -> dict[str, str]:
    row = {
        "Infrasp": "",
        "FlrSymCheck": "radial",
        "trait.source": "FNA",
        "genus_trait_source": "",
        "family_trait_source": "",
        "tnrs_status": "Accepted",
        "tnrs_family": "Exampleaceae",
        "tnrs_Sciname": species,
        "tnrs_infrasp": "",
    }
    row.update(updates)
    return row


def test_build_checkpoint_keeps_only_species_level_symmetry_sources() -> None:
    source = pd.DataFrame(
        [
            _row("Alpha one"),
            _row(
                "Beta two",
                FlrSymCheck="bilateral",
                tnrs_status="Synonym",
                **{"trait.source": "TRY|BIEN"},
            ),
            _row("Gamma three", FlrSymCheck="slightly zygomorphic"),
            _row(
                "Delta four",
                **{
                    "trait.source": "",
                    "genus_trait_source": "FNA",
                    "family_trait_source": "Kubitzki",
                },
            ),
            _row("Epsilon five", FlrSymCheck="wind"),
            _row("Zeta six", tnrs_family="Otheraceae"),
            _row("Eta seven", Infrasp="subsp. minor"),
            _row("Theta eight"),
            _row("Outside nine"),
        ]
    )
    master = pd.DataFrame(
        {
            "accepted_species": [
                "Alpha one",
                "Beta two",
                "Gamma three",
                "Delta four",
                "Epsilon five",
                "Zeta six",
                "Eta seven",
                "Theta eight",
            ],
            "family": ["Exampleaceae"] * 8,
        }
    )
    baseline = pd.DataFrame(
        {"accepted_species": ["Theta eight"], "trait_name": ["floral_symmetry"]}
    )

    evidence, audit, selection = build_checkpoint(
        source,
        baseline,
        master,
        audit_size=3,
    )

    assert evidence["accepted_species"].tolist() == [
        "Alpha one",
        "Beta two",
        "Gamma three",
    ]
    assert evidence["normalized_value"].tolist() == [
        "actinomorphic",
        "zygomorphic",
        "zygomorphic",
    ]
    assert evidence["trait_name"].eq("floral_symmetry").all()
    assert evidence["axis"].eq("floral_structural_complexity").all()
    assert evidence["quality"].eq("high").all()
    assert evidence["evidence_scope"].tolist() == [
        "species_direct",
        "synonym_direct",
        "species_direct",
    ]
    assert evidence["source_excerpt"].str.contains("trait.source=").all()
    assert evidence["lineage_method"].eq(
        "individually_collected_species_trait_source_not_family_or_genus"
    ).all()
    assert len(audit) == 3
    assert set(audit["candidate_id"]) == set(evidence["source_record_id"])

    reasons = selection.set_index("tnrs_Sciname")["selection_reason"].to_dict()
    assert reasons["Delta four"] == "missing_individually_collected_species_source"
    assert reasons["Epsilon five"] == "unsupported_floral_symmetry_state"
    assert reasons["Zeta six"] == "family_conflict"
    assert reasons["Eta seven"] == "infraspecific_record"
    assert reasons["Theta eight"] == "already_resolved_direct_pair"
    assert reasons["Outside nine"] == "not_in_strict_island_universe"


def test_build_checkpoint_rejects_missing_required_columns() -> None:
    source = pd.DataFrame([_row("Alpha one")]).drop(columns=["trait.source"])
    master = pd.DataFrame(
        {"accepted_species": ["Alpha one"], "family": ["Exampleaceae"]}
    )
    baseline = pd.DataFrame(columns=["accepted_species", "trait_name"])

    try:
        build_checkpoint(source, baseline, master, audit_size=1)
    except ValueError as error:
        assert "trait.source" in str(error)
    else:
        raise AssertionError("missing provenance column must fail closed")


def test_trait_source_lineage_is_order_insensitive_for_provider_lists() -> None:
    assert _trait_source_lineage_key("GBIF|KEW") == _trait_source_lineage_key(
        "KEW | GBIF"
    )
    assert _trait_source_lineage_key("FNA, GBIF") == "fna__gbif"
