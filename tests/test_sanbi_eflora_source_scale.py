from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from analysis.acquire_sanbi_eflora_traits_20260809 import (
    EXPECTED_ARCHIVE_SHA256,
    build_evidence,
    canonical_species_name,
    has_ancillary_colour_conflict,
    is_clean_primary_colour,
)
from island_v2.open_web_common import reviewed_source_package_evidence
from island_v2.trait_measurements import load_rules


def test_canonical_species_name_rejects_hybrids_and_infraspecific_records() -> None:
    assert canonical_species_name("Erica testacea L.") == "Erica testacea"
    assert canonical_species_name("Erica testacea (L.) Smith") == "Erica testacea"
    assert canonical_species_name("Erica testacea F.Muell.") == "Erica testacea"
    assert canonical_species_name("Erica testacea subsp. minor Smith") == ""
    assert canonical_species_name("Erica testacea var. minor Smith") == ""
    assert canonical_species_name("Erica × testacea") == ""


def test_build_evidence_keeps_exact_quote_citation_and_trait_identity() -> None:
    taxon = pd.DataFrame(
        [
            {
                "id": "wfo-1",
                "scientificName": "Erica testacea L.",
            },
            {
                "id": "wfo-2",
                "scientificName": "Erica testacea subsp. minor Smith",
            },
            {
                "id": "wfo-3",
                "scientificName": "Outside absentia L.",
            },
        ]
    )
    sentence = (
        "Flowers white, campanulate and actinomorphic, 12 mm long, in racemes; "
        "corolla tube 8 mm long. Leaves red in autumn."
    )
    descriptions = pd.DataFrame(
        [
            {
                "id": "wfo-1",
                "description": sentence,
                "type": "Morphology",
                "source": "ref-1",
                "language": "English",
            },
            {
                "id": "wfo-1",
                "description": "Flowers pink in cultivated forms.",
                "type": "Morphology",
                "source": "ref-1",
                "language": "English",
            },
            {
                "id": "wfo-1",
                "description": "Flowers blue.",
                "type": "Habitat",
                "source": "ref-1",
                "language": "English",
            },
            {
                "id": "wfo-2",
                "description": "Flowers yellow.",
                "type": "Morphology",
                "source": "ref-2",
                "language": "English",
            },
        ]
    )
    citation = "Botanist, A. 2012. Species treatment. [CC BY]"
    references = pd.DataFrame(
        [
            {
                "id": "wfo-1",
                "identifier": "ref-1",
                "bibliographicCitation": citation,
            }
        ]
    )
    evidence, inventory = build_evidence(
        taxon=taxon,
        descriptions=descriptions,
        references=references,
        accepted_species={"Erica testacea"},
        current_direct_pairs={("Erica testacea", "flower_primary_color")},
        rules=load_rules(Path("config/measurement_classification.yml")),
    )

    assert inventory["exact_fixed_species"] == 1
    assert not evidence.empty
    assert "flower_primary_color" not in set(evidence["trait_name"])
    assert {
        "floral_form",
        "floral_symmetry",
        "flower_size_class",
        "inflorescence_display",
        "tube_depth_class",
    }.issubset(set(evidence["trait_name"]))
    assert set(evidence["accepted_species"]) == {"Erica testacea"}
    assert set(evidence["source_citation"]) == {citation}
    assert all(quote in sentence for quote in evidence["source_excerpt"])
    assert not evidence["source_excerpt"].str.contains("cultivated", case=False).any()
    assert evidence["source_lineage"].str.startswith("citation:").all()
    assert evidence["source_record_id"].is_unique


def test_sanbi_archive_is_pinned_to_observed_release_hash() -> None:
    assert EXPECTED_ARCHIVE_SHA256 == (
        "157fe16ca6a5698df24fca4796aec2d69a8da20fd38a3bd36633eb0d5cdedfe4"
    )


def test_ancillary_organ_colours_are_not_primary_flower_colours() -> None:
    assert has_ancillary_colour_conflict("Flowers yellow; staminodes black.")
    assert has_ancillary_colour_conflict("Buds reddish. Petals white.")
    assert has_ancillary_colour_conflict("Sepals green or white.")
    assert has_ancillary_colour_conflict("Involucral bracts whitish; florets yellow.")
    assert not has_ancillary_colour_conflict("Flowers white with a red star in the centre.")
    assert not has_ancillary_colour_conflict("Corolla greenish yellow, petals purple.")
    assert is_clean_primary_colour("Flowers white with a red star.", "red_pink|white")
    assert not is_clean_primary_colour("Flowers yellow; staminodes black.", "yellow_orange")
    assert not is_clean_primary_colour("Flowers probably yellow.", "yellow_orange")
    assert not is_clean_primary_colour("Corolla red; lobes green.", "red_pink")


def test_cultivated_form_measurements_are_not_source_evidence() -> None:
    taxon = pd.DataFrame(
        [{"id": "wfo-1", "scientificName": "Cucumis testacea L."}]
    )
    descriptions = pd.DataFrame(
        [
            {
                "id": "wfo-1",
                "description": (
                    "Corolla in cultivated forms 20 mm long, but in wild forms "
                    "usually 8 mm long."
                ),
                "type": "Morphology",
                "source": "ref-1",
            }
        ]
    )
    references = pd.DataFrame(
        [
            {
                "id": "wfo-1",
                "identifier": "ref-1",
                "bibliographicCitation": "Flora treatment.",
            }
        ]
    )

    evidence, _ = build_evidence(
        taxon=taxon,
        descriptions=descriptions,
        references=references,
        accepted_species={"Cucumis testacea"},
        current_direct_pairs=set(),
        rules=load_rules(Path("config/measurement_classification.yml")),
    )

    assert evidence.empty


def test_committed_sanbi_package_passes_trait_gate_and_manifest() -> None:
    root = Path(
        "data/v2/staging/traits/direct_llm_pilot/"
        "20260809_sanbi_eflora_source_acquisition"
    )
    evidence = pd.read_csv(
        root / "sanbi_eflora_south_africa_evidence.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "sanbi_eflora_south_africa_manual_audit_200.csv", dtype=str
    ).fillna("")
    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)
    review_summary = json.loads(
        (root / "sanbi_eflora_south_africa_review_summary.json").read_text(
            encoding="utf-8"
        )
    )

    assert summary["package_gate_passed"]
    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 198
    assert summary["precision"] == 0.99
    assert review_summary["precision_contract"] == "accepted_correct/all_reviewed"
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "floral_form",
        "floral_symmetry",
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
        "tube_depth_class",
    }
    assert len(evidence) == 4_069
    assert len(selected) == 4_067
    assert selected.drop_duplicates(["accepted_species", "trait_name"]).shape[0] == 2_548
    assert selected.drop_duplicates(["accepted_species", "axis"]).shape[0] == 1_947
    rejected = set(
        audit.loc[audit["accepted_correct"].str.casefold().eq("false"), "candidate_id"]
    )
    assert len(rejected) == 2
    assert not rejected.intersection(selected["source_record_id"])

    manifest = {}
    for line in (root / "file_manifest.sha256").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", maxsplit=1)
        manifest[name] = digest
    assert set(manifest) == {
        path.name for path in root.iterdir() if path.is_file() and path.name != "file_manifest.sha256"
    }
    for name, expected in manifest.items():
        assert hashlib.sha256((root / name).read_bytes()).hexdigest() == expected
