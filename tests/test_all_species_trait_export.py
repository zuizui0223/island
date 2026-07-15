from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.all_species_trait_export import OUTPUT_COLUMNS, build_export, map_candidate_value


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

    llm = tmp_path / "batch_0000_results.jsonl"
    llm.write_text(
        json.dumps(
            {
                "job_id": "job-beta",
                "result": (
                    "Beta two,blue,tubular,bees,bee pollination inferred,"
                    "mixed_mating,likely_SC,horticulture,high"
                ),
            }
        )
        + "\n",
        encoding="utf-8",
    )

    output = tmp_path / "output"
    report = build_export(
        master_csv=master,
        output_dir=output,
        applicability_csv=applicability,
        candidate_csvs=[efloras, database],
        llm_results_jsonls=[llm],
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
    assert evidence.loc[evidence["source_kind"].eq("llm_only"), "source_backed"].eq(
        "False"
    ).all()
    assert report["n_master_species"] == 3
    assert report["n_output_rows"] == 3
    assert report["n_trait_applicable_species"] == 2
    assert report["n_trait_not_applicable_species"] == 1
    assert report["n_species_all_unknown"] == 1
    assert report["n_llm_only_evidence_rows"] == 5


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
    llm = tmp_path / "llm.jsonl"
    llm.write_text(
        json.dumps(
            {
                "job_id": "job-alpha",
                "result": "Alpha one,red,unknown,unknown,unknown,unknown,unknown,inference,low",
            }
        )
        + "\n",
        encoding="utf-8",
    )
    build_export(
        master_csv=master,
        output_dir=tmp_path / "out",
        candidate_csvs=[candidates],
        llm_results_jsonls=[llm],
    )
    result = pd.read_csv(tmp_path / "out" / "all_species_traits.csv", dtype=str)
    assert result.loc[0, "flower_color"] == "white"
    assert result.loc[0, "evidence_type"] == "review"
    assert result.loc[0, "confidence"] == "medium"
