import json

import pandas as pd

from island_v2 import global_candidate_index as indexer


def _candidate(
    species: str,
    value: str,
    *,
    trait: str = "pollination_functional_guild",
    task: str = "alternative_guild_openalex",
    phase: str = "alternative_pollinator_function",
    url: str = "https://example.test/source",
    citation: str = "Example",
    confidence: str = "machine_extracted",
) -> dict:
    return {
        "accepted_species": species,
        "trait_name": trait,
        "candidate_value": value,
        "source_url": url,
        "source_citation": citation,
        "source_excerpt": f"{species} is associated with {value}.",
        "source_type": "paper_abstract",
        "evidence_scope": "species_direct",
        "matched_term": value,
        "confidence": confidence,
        "campaign_task": task,
        "campaign_phase": phase,
        "target_for_task": True,
    }


def _write_wave(campaign_dir, wave_id: str, rows: list[dict]) -> None:
    wave_dir = campaign_dir / "waves" / wave_id
    wave_dir.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(
        wave_dir / "machine_candidates.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )


def _write_ledger(campaign_dir) -> None:
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "family": "FamA",
                "n_islands": "3",
                "n_records": "12",
                "machine_biotic_candidate": "True",
                "reproductive_wikimedia_status": "processed",
                "reproductive_openalex_status": "processed",
                "floral_access_wikimedia_status": "processed",
                "alternative_guild_openalex_status": "processed",
            },
            {
                "accepted_species": "Beta two",
                "family": "FamB",
                "n_islands": "1",
                "n_records": "4",
                "machine_biotic_candidate": "True",
                "reproductive_wikimedia_status": "processed",
                "reproductive_openalex_status": "processed",
                "floral_access_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
            },
            {
                "accepted_species": "Gamma three",
                "family": "FamC",
                "n_islands": "2",
                "n_records": "8",
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "processed",
                "reproductive_openalex_status": "pending",
                "floral_access_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
            },
        ]
    )
    ledger.to_csv(
        campaign_dir / "campaign_ledger.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )


def test_candidate_index_deduplicates_source_rows_and_retains_wave_history(tmp_path):
    campaign_dir = tmp_path / "campaign"
    duplicate = _candidate("Alpha one", "other_bees")
    _write_wave(campaign_dir, "wave_00001_alternative", [duplicate])
    _write_wave(
        campaign_dir,
        "wave_00002_alternative",
        [duplicate, _candidate("Alpha one", "bumblebees", url="https://example.test/bee")],
    )

    raw, wave_ids = indexer.load_wave_candidates(campaign_dir)
    index = indexer.build_candidate_index(raw)

    assert wave_ids == ["wave_00001_alternative", "wave_00002_alternative"]
    assert len(index) == 2
    other = index.loc[index["candidate_value"].eq("other_bees")].iloc[0]
    assert other["n_wave_observations"] == 2
    assert other["first_wave_id"] == "wave_00001_alternative"
    assert other["last_wave_id"] == "wave_00002_alternative"
    assert other["source_citation"] == "Example"
    assert other["confidence"] == "machine_extracted"
    assert bool(other["target_for_task"])


def test_index_does_not_collapse_distinct_citations_or_confidence():
    raw = pd.DataFrame(
        [
            {
                **_candidate("Alpha one", "flies", citation="Paper A"),
                "wave_id": "wave_1",
            },
            {
                **_candidate("Alpha one", "flies", citation="Paper B"),
                "wave_id": "wave_2",
            },
            {
                **_candidate(
                    "Alpha one",
                    "flies",
                    citation="Paper A",
                    confidence="explicit_term_match",
                ),
                "wave_id": "wave_3",
            },
        ]
    )

    index = indexer.build_candidate_index(raw)

    assert len(index) == 3
    assert set(index["source_citation"]) == {"Paper A", "Paper B"}
    assert set(index["confidence"]) == {"machine_extracted", "explicit_term_match"}


def test_counterfactual_screening_is_machine_only_and_keeps_unknown_guilds(tmp_path):
    campaign_dir = tmp_path / "campaign"
    campaign_dir.mkdir()
    _write_ledger(campaign_dir)
    raw = pd.DataFrame(
        [
            {**_candidate("Alpha one", "bumblebees"), "wave_id": "wave_1"},
            {**_candidate("Alpha one", "flies"), "wave_id": "wave_2"},
            {
                **_candidate(
                    "Alpha one",
                    "tubular",
                    trait="floral_form",
                    task="floral_access_wikimedia",
                    phase="floral_access",
                    url="https://example.test/floral",
                ),
                "wave_id": "wave_3",
            },
        ]
    )
    index = indexer.build_candidate_index(raw)
    ledger = indexer.load_campaign_ledger(campaign_dir)

    screening = indexer.build_counterfactual_screening(index, ledger)

    assert set(screening["accepted_species"]) == {"Alpha one", "Beta two"}
    alpha = screening.loc[screening["accepted_species"].eq("Alpha one")].iloc[0]
    assert bool(alpha["machine_has_bumblebee_candidate"])
    assert bool(alpha["machine_has_nonbee_candidate"])
    assert bool(alpha["machine_has_alternative_guild_candidate"])
    assert alpha["machine_screening_label"] == "alternative_or_generalist_candidate"
    assert alpha["n_floral_access_candidates"] == 1
    beta = screening.loc[screening["accepted_species"].eq("Beta two")].iloc[0]
    assert beta["machine_screening_label"] == "biotic_vector_candidate_without_guild_yet"
    assert beta["alternative_guild_openalex_status"] == "pending"
    assert set(screening["evidence_use_tier"]) == {"machine_unreviewed_screening_only"}


def test_build_command_writes_reproducible_outputs_and_status(tmp_path):
    campaign_dir = tmp_path / "campaign"
    campaign_dir.mkdir()
    _write_ledger(campaign_dir)
    _write_wave(
        campaign_dir,
        "wave_00001_alternative",
        [_candidate("Alpha one", "mixed_or_generalist")],
    )

    indexer.build(campaign_dir=campaign_dir, output_dir=None)

    index = pd.read_csv(campaign_dir / "machine_candidate_index.csv.gz")
    screening = pd.read_csv(campaign_dir / "counterfactual_screening.csv.gz")
    status = json.loads((campaign_dir / "candidate_index_status.json").read_text())
    first_bytes = (campaign_dir / "machine_candidate_index.csv.gz").read_bytes()
    indexer.build(campaign_dir=campaign_dir, output_dir=None)
    second_bytes = (campaign_dir / "machine_candidate_index.csv.gz").read_bytes()

    assert len(index) == 1
    assert len(screening) == 2
    assert status["n_waves_indexed"] == 1
    assert status["n_unique_candidate_rows"] == 1
    assert status["n_machine_biotic_screening_species"] == 2
    assert status["n_by_source_type"] == {"paper_abstract": 1}
    assert status["n_by_confidence"] == {"machine_extracted": 1}
    assert "not accepted trait values" in status["interpretation"]
    assert first_bytes == second_bytes
