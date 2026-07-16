from __future__ import annotations

import json
import hashlib
from pathlib import Path

import pandas as pd
import pytest
import typer

from island_v2.all_species_trait_export import OUTPUT_COLUMNS, build_export, map_candidate_value
from island_v2.llm_evidence_extraction import prepare_packets, validate_results
from island_v2.search_enabled_llm_campaign import _job_record, canonical_result_sha


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_inference_bundle(
    root: Path,
    *,
    species: str,
    genus: str,
    family: str,
    values: dict[str, str],
) -> Path:
    root.mkdir(parents=True)
    target_fields = "|".join(
        field
        for field in (
            "flower_color",
            "flower_shape",
            "pollination_guild",
            "mating_system",
            "self_incompatibility",
        )
        if values.get(field, "unknown") != "unknown"
    )
    job = _job_record(
        {
            "species": species,
            "genus": genus,
            "family": family,
            "target_fields": target_fields,
        },
        0,
        1,
    )
    row = {
        "species": species,
        "flower_color": values.get("flower_color", "unknown"),
        "flower_shape": values.get("flower_shape", "unknown"),
        "pollination_guild": values.get("pollination_guild", "unknown"),
        "pollination_notes": values.get("pollination_notes", "unknown"),
        "mating_system": values.get("mating_system", "unknown"),
        "self_incompatibility": values.get("self_incompatibility", "unknown"),
        "evidence_type": values.get("evidence_type", "inference"),
        "confidence": values.get("confidence", "low"),
        "job_id": job["job_id"],
        "job_sha256": job["job_sha256"],
        "target_fields": target_fields,
        "model_provider": "test",
        "model_id": "test-model",
        "prompt_version": job["prompt_version"],
        "prompt_sha256": job["prompt_sha256"],
        "result_sha256": "",
    }
    row["result_sha256"] = canonical_result_sha(row)
    results_path = root / "validated_trait_results.csv"
    pd.DataFrame([row]).to_csv(results_path, index=False)
    history_path = root / "job_history.csv"
    pd.DataFrame(
        [
            {
                "job_id": row["job_id"],
                "job_sha256": row["job_sha256"],
                "species": species,
                "target_fields": target_fields,
                "prompt_version": row["prompt_version"],
                "prompt_sha256": row["prompt_sha256"],
                "result_sha256": row["result_sha256"],
                "status": "completed",
            }
        ]
    ).to_csv(history_path, index=False)
    (root / "validated_llm_manifest.json").write_text(
        json.dumps(
            {
                "contract_version": "validated_llm_traits_v1",
                "validation_status": "success",
                "n_validated_results": 1,
                "validated_trait_results_sha256": _sha256(results_path),
                "job_history_sha256": _sha256(history_path),
            }
        ),
        encoding="utf-8",
    )
    return root


def test_candidate_mapping_is_direct_and_conservative() -> None:
    assert map_candidate_value("flower_primary_color", "red|white") == [
        ("flower_color", "red, white")
    ]
    assert map_candidate_value("floral_form", "campanulate|tubular") == [
        ("flower_shape", "bell-shaped"),
        ("flower_shape", "tubular"),
    ]
    assert map_candidate_value("pollen_vector_mode", "abiotic_wind") == [
        ("pollination_guild", "wind")
    ]
    assert map_candidate_value("pollination_functional_guild", "bats") == []
    assert map_candidate_value("pollen_vector_mode", "biotic") == []
    assert map_candidate_value("mating_system", "predominantly_outcrossing") == []


def test_genus_si_candidates_are_likely_low_confidence_and_never_source_backed(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Alpha target", "Beta target", "Gamma conflict"],
            "genus": ["Alpha", "Beta", "Gamma"],
            "family": ["FamA", "FamB", "FamC"],
        }
    ).to_csv(master, index=False)
    candidates = tmp_path / "candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha target",
                "trait_name": "self_incompatibility",
                "standardized_value": "SC",
                "candidate_kind": "source_backed",
                "evidence_scope": "species_direct",
                "source_type": "curated_trait_database",
                "source_name": "direct_db",
                "source_url": "https://example.test/direct",
                "source_record_id": "direct:1",
                "evidence_excerpt": "Alpha target is self-compatible.",
                "source_reliability": "B_curated_database_or_institution",
                "confidence": "medium",
            },
            {
                "accepted_species": "Alpha target",
                "trait_name": "self_incompatibility",
                "standardized_value": "SI",
                "candidate_kind": "hierarchical_inference",
                "evidence_scope": "genus_inference",
                "source_type": "curated_trait_database",
                "source_name": "genus_db",
                "source_url": "https://example.test/genus",
                "source_record_id": "genus:alpha",
                "evidence_excerpt": "Five directly reported congeners are SI.",
                "source_reliability": "B_curated_database_or_institution",
                "confidence": "low",
            },
            {
                "accepted_species": "Beta target",
                "trait_name": "self_incompatibility",
                "standardized_value": "SI",
                "candidate_kind": "hierarchical_inference",
                "evidence_scope": "genus_inference",
                "source_type": "curated_trait_database",
                "source_name": "genus_db",
                "source_url": "https://example.test/genus",
                "source_record_id": "genus:beta",
                "evidence_excerpt": "Six directly reported congeners are SI.",
                "source_reliability": "B_curated_database_or_institution",
                "confidence": "low",
            },
            *[
                {
                    "accepted_species": "Gamma conflict",
                    "trait_name": "self_incompatibility",
                    "standardized_value": value,
                    "candidate_kind": "source_backed",
                    "evidence_scope": "species_direct",
                    "source_type": "curated_trait_database",
                    "source_name": f"conflict_{value}",
                    "source_url": f"https://example.test/{value}",
                    "source_record_id": f"conflict:{value}",
                    "evidence_excerpt": f"Gamma conflict is {value}.",
                    "source_reliability": "B_curated_database_or_institution",
                    "confidence": "medium",
                }
                for value in ("SI", "SC")
            ],
        ]
    ).to_csv(candidates, index=False)

    output = tmp_path / "output"
    report = build_export(master_csv=master, output_dir=output, candidate_csvs=[candidates])
    result = pd.read_csv(output / "all_species_traits.csv", dtype=str).set_index("species")
    assert result.loc["Alpha target", "self_incompatibility"] == "SC"
    assert result.loc["Alpha target", "confidence"] == "medium"
    assert result.loc["Beta target", "self_incompatibility"] == "likely_SI"
    assert result.loc["Beta target", "evidence_type"] == "inference"
    assert result.loc["Beta target", "confidence"] == "low"
    assert result.loc["Gamma conflict", "self_incompatibility"] == "unknown"
    assert result.loc["Gamma conflict", "confidence"] == "low"
    assert report["n_source_backed_evidence_rows"] == 3


def test_build_export_retains_every_species_and_keeps_llm_separate(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two", "Gamma three"],
            "genus": ["Alpha", "Beta", "Gamma"],
            "family": ["Aaceae", "Baceae", "Gaceae"],
        }
    ).to_csv(master, index=False)
    applicability = tmp_path / "applicability.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Beta two", "Gamma three"],
            "angiosperm_analysis_eligible": [True, True, False],
        }
    ).to_csv(applicability, index=False)

    efloras = tmp_path / "efloras_trait_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "flower_primary_color",
                "candidate_value": "red|white",
                "flora_id": "1",
                "flora_name": "Example Flora",
                "source_url": "https://flora.test/alpha",
                "source_record_id": "flora:1",
                "source_excerpt": "Flowers red or white.",
            },
            {
                "accepted_species": "Alpha one",
                "trait_name": "self_incompatibility",
                "candidate_value": "self_incompatible",
                "flora_id": "1",
                "flora_name": "Example Flora",
                "source_url": "https://flora.test/alpha",
                "source_record_id": "flora:2",
                "source_excerpt": "Plants self-incompatible.",
            },
            {
                "accepted_species": "Gamma three",
                "trait_name": "flower_primary_color",
                "candidate_value": "yellow",
                "flora_id": "1",
                "flora_name": "Example Flora",
                "source_url": "https://flora.test/gamma",
                "source_record_id": "flora:3",
                "source_excerpt": "Flowers yellow.",
            },
        ]
    ).to_csv(efloras, index=False)

    database = tmp_path / "database_trait_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "pollen_vector_mode",
                "standardized_value": "abiotic_wind",
                "source_type": "curated_trait_database",
                "source_name": "Example DB",
                "source_url": "https://db.test/alpha",
                "source_record_id": "db:1",
                "evidence_excerpt": "Pollen vector: wind.",
                "source_reliability": "B_curated_database_or_institution",
            }
        ]
    ).to_csv(database, index=False)

    llm = _write_inference_bundle(
        tmp_path / "validated_llm_bundle",
        species="Beta two",
        genus="Beta",
        family="Baceae",
        values={
            "flower_color": "blue",
            "flower_shape": "tubular",
            "pollination_guild": "bees",
            "pollination_notes": "bee pollination inferred",
            "mating_system": "mixed_mating",
            "self_incompatibility": "likely_SC",
            "evidence_type": "horticulture",
            "confidence": "high",
        },
    )

    output = tmp_path / "output"
    report = build_export(
        master_csv=master,
        output_dir=output,
        applicability_csv=applicability,
        candidate_csvs=[efloras, database],
        llm_inference_bundle_dirs=[llm],
    )

    result = pd.read_csv(output / "all_species_traits.csv", dtype=str).fillna("")
    assert list(result.columns) == OUTPUT_COLUMNS
    assert result["species"].tolist() == ["Alpha one", "Beta two", "Gamma three"]
    alpha = result.set_index("species").loc["Alpha one"]
    assert alpha["flower_color"] == "red, white"
    assert alpha["pollination_guild"] == "wind"
    assert alpha["self_incompatibility"] == "SI"
    assert set(alpha["evidence_type"].split("|")) == {"review", "flora"}
    assert alpha["confidence"] == "medium"
    beta = result.set_index("species").loc["Beta two"]
    assert beta["flower_color"] == "blue"
    assert beta["self_incompatibility"] == "likely_SC"
    assert beta["evidence_type"] == "inference"
    assert beta["confidence"] == "low"
    gamma = result.set_index("species").loc["Gamma three"]
    assert set(gamma[[*OUTPUT_COLUMNS[1:7]]]) == {"unknown"}
    assert gamma["evidence_type"] == "inference"
    assert gamma["confidence"] == "low"
    evidence = pd.read_csv(output / "all_species_trait_evidence.csv.gz", dtype=str).fillna("")
    assert evidence.loc[evidence["source_kind"].eq("llm_only"), "source_backed"].eq("False").all()
    channels = pd.read_csv(
        output / "all_species_trait_value_channels.csv.gz", dtype=str
    ).fillna("")
    assert len(channels) == 3 * 5
    alpha_color = channels.loc[
        channels["species"].eq("Alpha one") & channels["field"].eq("flower_color")
    ].iloc[0]
    assert alpha_color["reported_value"] == "red, white"
    assert alpha_color["inferred_value"] == "unknown"
    assert alpha_color["selected_origin"] == "reported"
    beta_color = channels.loc[
        channels["species"].eq("Beta two") & channels["field"].eq("flower_color")
    ].iloc[0]
    assert beta_color["reported_value"] == "unknown"
    assert beta_color["inferred_value"] == "blue"
    assert beta_color["selected_origin"] == "inferred"
    assert channels.loc[channels["species"].eq("Gamma three"), "selected_origin"].eq(
        "not_applicable"
    ).all()
    assert report["n_master_species"] == 3
    assert report["n_output_rows"] == 3
    assert report["n_trait_applicable_species"] == 2
    assert report["n_trait_not_applicable_species"] == 1
    assert report["n_species_all_unknown"] == 1
    assert report["n_llm_only_evidence_rows"] == 5
    assert report["n_reported_cells"] == 3
    assert report["n_inferred_cells"] == 5
    assert report["n_not_applicable_cells"] == 5
    assert report["n_unresolved_applicable_cells"] == 2


def test_source_backed_field_outranks_llm_for_same_field(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame({"accepted_species": ["Alpha one"]}).to_csv(master, index=False)
    candidates = tmp_path / "trait_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "flower_primary_color",
                "standardized_value": "white",
                "source_type": "curated_trait_database",
                "source_url": "https://db.test/alpha",
                "source_record_id": "db:1",
                "evidence_excerpt": "Flower color = white.",
                "source_reliability": "B_curated_database_or_institution",
            }
        ]
    ).to_csv(candidates, index=False)
    llm = _write_inference_bundle(
        tmp_path / "validated_llm_bundle",
        species="Alpha one",
        genus="Alpha",
        family="",
        values={"flower_color": "red"},
    )
    build_export(
        master_csv=master,
        output_dir=tmp_path / "out",
        candidate_csvs=[candidates],
        llm_inference_bundle_dirs=[llm],
    )
    result = pd.read_csv(tmp_path / "out" / "all_species_traits.csv", dtype=str)
    assert result.loc[0, "flower_color"] == "white"
    assert result.loc[0, "evidence_type"] == "review"
    assert result.loc[0, "confidence"] == "medium"
    channels = pd.read_csv(
        tmp_path / "out" / "all_species_trait_value_channels.csv.gz", dtype=str
    )
    color = channels.loc[channels["field"].eq("flower_color")].iloc[0]
    assert color["reported_value"] == "white"
    assert color["inferred_value"] == "red"
    assert color["selected_value"] == "white"
    assert color["selected_origin"] == "reported"


def test_inference_bundle_rejects_jointly_rehashed_fake_result_hash(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Alpha one"],
            "genus": ["Alpha"],
            "family": ["Aaceae"],
        }
    ).to_csv(master, index=False)
    bundle = _write_inference_bundle(
        tmp_path / "bundle",
        species="Alpha one",
        genus="Alpha",
        family="Aaceae",
        values={"flower_color": "red"},
    )
    results_path = bundle / "validated_trait_results.csv"
    results = pd.read_csv(results_path, dtype=str).fillna("")
    results.loc[0, "result_sha256"] = "b" * 64
    results.to_csv(results_path, index=False)
    manifest_path = bundle / "validated_llm_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["validated_trait_results_sha256"] = _sha256(results_path)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    with pytest.raises(typer.BadParameter, match="canonical result_sha256 mismatch"):
        build_export(
            master_csv=master,
            output_dir=tmp_path / "out",
            llm_inference_bundle_dirs=[bundle],
        )


def test_evidence_locked_llm_candidates_cannot_bypass_bundle_validation(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame({"accepted_species": ["Alpha one"]}).to_csv(master, index=False)
    candidates = tmp_path / "v2_llm_reported_candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "provisional_candidate_value": "campanulate",
                "source_type": "flora_or_monograph",
                "source_url": "https://flora.test/alpha",
                "evidence_quote": "Corolla campanulate.",
                "source_chunk_id": "src_123",
                "extraction_method": "provider_neutral_llm_from_frozen_source_text",
                "confidence": "high",
            }
        ]
    ).to_csv(candidates, index=False)

    with pytest.raises(typer.BadParameter, match="validated extraction bundle"):
        build_export(
            master_csv=master,
            output_dir=tmp_path / "out",
            candidate_csvs=[candidates],
        )


def test_evidence_locked_llm_bundle_is_revalidated_against_frozen_text(
    tmp_path: Path,
) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Plantus alba"],
            "genus": ["Plantus"],
            "family": ["Plantaceae"],
        }
    ).to_csv(master, index=False)
    sources = tmp_path / "sources.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alba",
                "source_text": "Flowers are white and tubular.",
                "source_url": "https://flora.test/plantus-alba",
                "source_citation": "Example Flora",
                "source_type": "flora_or_monograph",
                "evidence_scope": "species_direct",
            }
        ]
    ).to_csv(sources, index=False)
    packets = tmp_path / "packets"
    prepare_packets(
        source_csv=sources,
        output_dir=packets,
        prompt_path=Path("config/v2_llm_evidence_prompt.txt"),
        ontology_path=Path("config/trait_ontology.yml"),
    )
    packet_dir = packets / "batch_00001"
    packet = json.loads((packet_dir / "packet_input.json").read_text(encoding="utf-8"))
    result_jsonl = tmp_path / "result.jsonl"
    result_jsonl.write_text(
        json.dumps(
            {
                "accepted_species": "Plantus alba",
                "trait_name": "floral_form",
                "standardized_value": "tubular",
                "source_chunk_id": packet["tasks"][0]["source_chunk_ids"][0],
                "evidence_quote": "Flowers are white and tubular.",
                "confidence": "high",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    ingest_dir = tmp_path / "ingest"
    validate_results(
        packet_dir=packet_dir,
        result_jsonl=result_jsonl,
        output_dir=ingest_dir,
        model_id="test-model",
        model_provider="test-provider",
    )
    bundle = tmp_path / "llm_extracted_bundle"
    bundle.mkdir()
    for source, target_name in (
        (ingest_dir / "v2_llm_reported_candidates.csv", "v2_llm_reported_candidates.csv"),
        (ingest_dir / "llm_ingest_manifest.json", "llm_ingest_manifest.json"),
        (packet_dir / "packet_input.json", "packet_input.json"),
        (packet_dir / "packet_manifest.json", "packet_manifest.json"),
        (packet_dir / "prompt.txt", "prompt.txt"),
        (result_jsonl, "result.jsonl"),
    ):
        (bundle / target_name).write_bytes(source.read_bytes())

    build_export(
        master_csv=master,
        output_dir=tmp_path / "out",
        llm_extracted_bundle_dirs=[bundle],
    )
    result = pd.read_csv(tmp_path / "out" / "all_species_traits.csv", dtype=str)
    assert result.loc[0, "flower_shape"] == "tubular"
    evidence = pd.read_csv(
        tmp_path / "out" / "all_species_trait_evidence.csv.gz", dtype=str
    ).fillna("")
    assert evidence.loc[0, "source_kind"] == "llm_extracted_flora_or_monograph"
    assert evidence.loc[0, "source_backed"] == "True"

    packet_payload = json.loads((bundle / "packet_input.json").read_text(encoding="utf-8"))
    packet_payload["sources"][0]["chunk_text"] = "Tampered source text."
    (bundle / "packet_input.json").write_text(json.dumps(packet_payload), encoding="utf-8")
    with pytest.raises(typer.BadParameter, match="packet_input_sha256 mismatch"):
        build_export(
            master_csv=master,
            output_dir=tmp_path / "tampered",
            llm_extracted_bundle_dirs=[bundle],
        )


def test_raw_llm_jsonl_is_not_an_accepted_input(tmp_path: Path) -> None:
    master = tmp_path / "master.csv"
    pd.DataFrame({"accepted_species": ["Alpha one"]}).to_csv(master, index=False)
    raw = tmp_path / "raw.jsonl"
    raw.write_text(
        json.dumps({"job_id": "job-alpha", "result": "Alpha one,red"}) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(typer.BadParameter):
        build_export(
            master_csv=master,
            output_dir=tmp_path / "out",
            llm_inference_bundle_dirs=[raw],
        )
