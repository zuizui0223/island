from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.llm_evidence_extraction import prepare_packets, validate_results


def _ontology(path: Path) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "version": "test",
                "traits": {
                    "flower_primary_color": {
                        "domain": "floral_signal",
                        "allowed_values": ["white", "red_pink", "unresolved"],
                    },
                    "floral_form": {
                        "domain": "floral_architecture",
                        "allowed_values": ["tubular", "open_radial", "unresolved"],
                    },
                    "self_incompatibility": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["SI", "SC", "unresolved"],
                    },
                    "autonomous_selfing_capacity": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["autonomous", "absent", "unresolved"],
                    },
                    "mating_system": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["mixed_mating", "unresolved"],
                    },
                    "cleistogamy": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["facultative", "unresolved"],
                    },
                    "pollen_vector_mode": {
                        "domain": "pollen_vector",
                        "allowed_values": ["biotic", "abiotic_wind", "unresolved"],
                    },
                    "pollination_functional_guild": {
                        "domain": "pollination",
                        "allowed_values": ["other_bees", "birds", "unresolved"],
                    },
                    "floral_symmetry": {
                        "domain": "floral_architecture",
                        "allowed_values": ["actinomorphic", "unresolved"],
                    },
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def _prompt(path: Path) -> None:
    path.write_text("Use only packet_input.json and return JSONL.\n", encoding="utf-8")


def _sources(path: Path) -> None:
    pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "source_text": "Plantus alba has flowers white and tubular. It is self-incompatible.",
                "source_url": "https://example.org/alba",
                "source_citation": "Example flora",
                "source_type": "flora_or_monograph",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Plantus indirecta",
                "source_text": "Flowers are red.",
                "source_url": "https://example.org/indirecta",
                "source_citation": "Genus overview",
                "source_type": "encyclopedia",
                "evidence_scope": "species_indirect",
            },
        ]
    ).to_csv(path, index=False)


def test_prepare_packets_is_gap_only_and_hash_locked(tmp_path: Path) -> None:
    ontology = tmp_path / "ontology.yml"
    prompt = tmp_path / "prompt.txt"
    sources = tmp_path / "sources.csv"
    existing = tmp_path / "existing.csv"
    _ontology(ontology)
    _prompt(prompt)
    _sources(sources)
    pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "white",
            }
        ]
    ).to_csv(existing, index=False)

    report = prepare_packets(
        source_csv=sources,
        output_dir=tmp_path / "packets",
        prompt_path=prompt,
        ontology_path=ontology,
        candidate_csvs=[existing],
        batch_size=10,
    )

    assert report["n_species_with_unresolved_target_traits"] == 1
    assert report["n_packets"] == 1
    packet_dir = tmp_path / "packets/batch_00001"
    packet = json.loads((packet_dir / "packet_input.json").read_text())
    task = packet["tasks"][0]
    assert task["accepted_species"] == "Plantus alba"
    assert "flower_primary_color" not in task["target_traits"]
    assert "unresolved" not in task["allowed_values"]["floral_form"]
    assert {source["accepted_species"] for source in packet["sources"]} == {"Plantus alba"}
    manifest = json.loads((packet_dir / "packet_manifest.json").read_text())
    assert len(manifest["packet_input_sha256"]) == 64
    assert len(manifest["prompt_sha256"]) == 64


def test_ingest_accepts_only_quote_backed_ontology_value(tmp_path: Path) -> None:
    ontology = tmp_path / "ontology.yml"
    prompt = tmp_path / "prompt.txt"
    sources = tmp_path / "sources.csv"
    _ontology(ontology)
    _prompt(prompt)
    _sources(sources)
    prepare_packets(
        source_csv=sources,
        output_dir=tmp_path / "packets",
        prompt_path=prompt,
        ontology_path=ontology,
        batch_size=10,
    )
    packet_dir = tmp_path / "packets/batch_00001"
    packet = json.loads((packet_dir / "packet_input.json").read_text())
    chunk_id = packet["tasks"][0]["source_chunk_ids"][0]
    result = tmp_path / "result.jsonl"
    result.write_text(
        json.dumps(
            {
                "accepted_species": "Plantus alba",
                "trait_name": "floral_form",
                "standardized_value": "tubular",
                "source_chunk_id": chunk_id,
                "evidence_quote": "flowers white and tubular",
                "confidence": "high",
            }
        )
        + "\n",
        encoding="utf-8",
    )

    report = validate_results(
        packet_dir=packet_dir,
        result_jsonl=result,
        output_dir=tmp_path / "ingest",
        model_id="test-model-1",
        model_provider="test-provider",
    )

    assert report["n_validated_claims"] == 1
    row = pd.read_csv(tmp_path / "ingest/v2_llm_reported_candidates.csv").iloc[0]
    assert row["trait_name"] == "floral_form"
    assert row["provisional_candidate_value"] == "tubular"
    assert row["review_status"] == "unreviewed_llm_extraction"
    assert row["model_id"] == "test-model-1"
    assert row["source_text_sha256"]
    assert row["extraction_run_id"].startswith("llm_")


def test_ingest_rejects_hallucinated_quote_and_invalid_value(tmp_path: Path) -> None:
    ontology = tmp_path / "ontology.yml"
    prompt = tmp_path / "prompt.txt"
    sources = tmp_path / "sources.csv"
    _ontology(ontology)
    _prompt(prompt)
    _sources(sources)
    prepare_packets(
        source_csv=sources,
        output_dir=tmp_path / "packets",
        prompt_path=prompt,
        ontology_path=ontology,
    )
    packet_dir = tmp_path / "packets/batch_00001"
    packet = json.loads((packet_dir / "packet_input.json").read_text())
    chunk_id = packet["tasks"][0]["source_chunk_ids"][0]
    result = tmp_path / "bad.jsonl"
    result.write_text(
        "\n".join(
            [
                json.dumps(
                    {
                        "accepted_species": "Plantus alba",
                        "trait_name": "floral_form",
                        "standardized_value": "tubular",
                        "source_chunk_id": chunk_id,
                        "evidence_quote": "flowers are visited by hummingbirds",
                        "confidence": "high",
                    }
                ),
                json.dumps(
                    {
                        "accepted_species": "Plantus alba",
                        "trait_name": "floral_form",
                        "standardized_value": "bell_campanulate",
                        "source_chunk_id": chunk_id,
                        "evidence_quote": "flowers white and tubular",
                        "confidence": "high",
                    }
                ),
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    with pytest.raises(Exception, match="failed validation"):
        validate_results(
            packet_dir=packet_dir,
            result_jsonl=result,
            output_dir=tmp_path / "ingest",
            model_id="test-model-1",
            strict=True,
        )

    rejected = pd.read_csv(tmp_path / "ingest/llm_rejected_claims.csv")
    assert set(rejected["reason"]) == {
        "evidence_quote_not_in_source_chunk",
        "value_outside_packet_ontology",
    }


def test_ingest_detects_modified_packet_input(tmp_path: Path) -> None:
    ontology = tmp_path / "ontology.yml"
    prompt = tmp_path / "prompt.txt"
    sources = tmp_path / "sources.csv"
    _ontology(ontology)
    _prompt(prompt)
    _sources(sources)
    prepare_packets(
        source_csv=sources,
        output_dir=tmp_path / "packets",
        prompt_path=prompt,
        ontology_path=ontology,
    )
    packet_dir = tmp_path / "packets/batch_00001"
    input_path = packet_dir / "packet_input.json"
    input_path.write_text(input_path.read_text() + " ", encoding="utf-8")
    result = tmp_path / "result.jsonl"
    result.write_text("", encoding="utf-8")

    with pytest.raises(Exception, match="hash does not match"):
        validate_results(
            packet_dir=packet_dir,
            result_jsonl=result,
            output_dir=tmp_path / "ingest",
            model_id="test-model-1",
        )
