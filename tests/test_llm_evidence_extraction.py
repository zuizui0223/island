from __future__ import annotations

import hashlib
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
                    "tube_depth_class": {
                        "domain": "floral_architecture",
                        "allowed_values": ["shallow", "deep", "unresolved"],
                    },
                    "flower_size_class": {
                        "domain": "floral_architecture",
                        "allowed_values": ["small", "large", "unresolved"],
                    },
                    "inflorescence_display": {
                        "domain": "floral_architecture",
                        "allowed_values": ["solitary", "few_flowered", "unresolved"],
                    },
                    "reward_type": {
                        "domain": "floral_architecture",
                        "allowed_values": ["nectar", "pollen", "unresolved"],
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
                    "sex_system": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["hermaphroditic", "dioecious", "unresolved"],
                    },
                    "herkogamy": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["pronounced", "absent", "unresolved"],
                    },
                    "dichogamy": {
                        "domain": "reproductive_assurance",
                        "allowed_values": ["protandry", "protogyny", "unresolved"],
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
    path.write_bytes(b"Use only packet_input.json and return JSONL.\r\nReturn no prose.\r\n")


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
    assert "sex_system" in task["target_traits"]
    assert "pollination_functional_guild" not in task["target_traits"]
    assert "unresolved" not in task["allowed_values"]["floral_form"]
    assert {source["accepted_species"] for source in packet["sources"]} == {"Plantus alba"}
    manifest = json.loads((packet_dir / "packet_manifest.json").read_text())
    written_prompt = (packet_dir / "prompt.txt").read_bytes()
    assert written_prompt == (b"Use only packet_input.json and return JSONL.\nReturn no prose.\n")
    assert manifest["prompt_sha256"] == hashlib.sha256(written_prompt).hexdigest()
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
    manifest = json.loads((tmp_path / "ingest/llm_ingest_manifest.json").read_text())
    accepted_path = tmp_path / "ingest/v2_llm_reported_candidates.csv"
    assert manifest["validation_status"] == "success"
    assert manifest["accepted_csv_filename"] == accepted_path.name
    assert manifest["accepted_csv_row_count"] == 1
    assert manifest["accepted_csv_sha256"] == hashlib.sha256(accepted_path.read_bytes()).hexdigest()
    assert manifest["result_jsonl_sha256"] == hashlib.sha256(result.read_bytes()).hexdigest()
    assert (
        manifest["packet_input_sha256"]
        == hashlib.sha256((packet_dir / "packet_input.json").read_bytes()).hexdigest()
    )
    assert manifest["model_provider"] == "test-provider"
    assert manifest["model_id"] == "test-model-1"
    assert manifest["extraction_run_id"] == row["extraction_run_id"]


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

    ingest_dir = tmp_path / "ingest"
    ingest_dir.mkdir()
    accepted_path = ingest_dir / "v2_llm_reported_candidates.csv"
    accepted_path.write_text("stale,accepted\n", encoding="utf-8")
    (ingest_dir / "llm_ingest_manifest.json").write_text(
        json.dumps({"validation_status": "success"}), encoding="utf-8"
    )

    with pytest.raises(Exception, match="failed validation"):
        validate_results(
            packet_dir=packet_dir,
            result_jsonl=result,
            output_dir=ingest_dir,
            model_id="test-model-1",
            strict=True,
        )

    rejected = pd.read_csv(ingest_dir / "llm_rejected_claims.csv")
    assert set(rejected["reason"]) == {
        "evidence_quote_not_in_source_chunk",
        "value_outside_packet_ontology",
    }
    assert not accepted_path.exists()
    failure_manifest = json.loads(
        (ingest_dir / "llm_ingest_manifest.json").read_text(encoding="utf-8")
    )
    assert failure_manifest["validation_status"] == "failed"
    assert failure_manifest["n_validated_claims"] == 0
    assert failure_manifest["n_rejected_claims"] == 2
    assert "accepted_csv_sha256" not in failure_manifest


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
