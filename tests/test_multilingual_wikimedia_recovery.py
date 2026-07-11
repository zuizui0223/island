from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.multilingual_wikimedia_recovery import (
    RECOVERY_TASK,
    extract_recovery_candidates,
    fetch_language_extracts,
    run_recovery_wave,
)


def test_fetch_language_extracts_batches_redirected_scientific_names() -> None:
    calls: list[tuple[str, dict]] = []

    def getter(url: str, params: dict) -> dict:
        calls.append((url, params))
        assert params["titles"] == "Campanula microdonta|Poa annua"
        return {
            "query": {
                "normalized": [{"from": "Poa_annua", "to": "Poa annua"}],
                "redirects": [
                    {"from": "Campanula microdonta", "to": "シマホタルブクロ"}
                ],
                "pages": {
                    "1": {
                        "title": "シマホタルブクロ",
                        "extract": "送粉者はマルハナバチである。",
                    },
                    "2": {
                        "title": "Poa annua",
                        "extract": "イチゴツナギ属の風媒花。",
                    },
                },
            }
        }

    sources, completed, errors = fetch_language_extracts(
        getter,
        ["Campanula microdonta", "Poa annua"],
        "ja",
    )

    assert errors == []
    assert completed == {"Campanula microdonta", "Poa annua"}
    assert len(calls) == 1
    assert set(sources["accepted_species"]) == {"Campanula microdonta", "Poa annua"}
    assert sources.loc[
        sources["accepted_species"].eq("Campanula microdonta"), "source_url"
    ].item().endswith("/%E3%82%B7%E3%83%9E%E3%83%9B%E3%82%BF%E3%83%AB%E3%83%96%E3%82%AF%E3%83%AD") is False
    # URL provenance keeps the resolved Wikipedia title; percent encoding is left to the client.
    assert "シマホタルブクロ" in sources.loc[
        sources["accepted_species"].eq("Campanula microdonta"), "source_url"
    ].item()


def test_extract_recovery_candidates_maps_japanese_chinese_and_russian() -> None:
    sources = pd.DataFrame(
        [
            {
                "accepted_species": "Campanula microdonta",
                "language": "ja",
                "source_text": "本種は両性花をもち、送粉者はマルハナバチである。",
                "source_url": "https://ja.wikipedia.org/wiki/example",
                "source_citation": "Japanese source",
                "source_type": "wikipedia_extract_ja",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Examplea sinensis",
                "language": "zh",
                "source_text": "该植物为雌雄异株，主要由熊蜂授粉。",
                "source_url": "https://zh.wikipedia.org/wiki/example",
                "source_citation": "Chinese source",
                "source_type": "wikipedia_extract_zh",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Examplea rossica",
                "language": "ru",
                "source_text": "Это двудомное растение, опыляемое шмелями.",
                "source_url": "https://ru.wikipedia.org/wiki/example",
                "source_citation": "Russian source",
                "source_type": "wikipedia_extract_ru",
                "evidence_scope": "species_direct",
            },
        ]
    )

    candidates = extract_recovery_candidates(sources)
    observed = set(
        candidates[["accepted_species", "trait_name", "candidate_value"]].itertuples(
            index=False, name=None
        )
    )

    assert ("Campanula microdonta", "sex_system", "hermaphroditic") in observed
    assert (
        "Campanula microdonta",
        "pollination_functional_guild",
        "bumblebees",
    ) in observed
    assert ("Campanula microdonta", "pollen_vector_mode", "biotic") in observed
    assert ("Examplea sinensis", "sex_system", "dioecious") in observed
    assert (
        "Examplea sinensis",
        "pollination_functional_guild",
        "bumblebees",
    ) in observed
    assert ("Examplea rossica", "sex_system", "dioecious") in observed
    assert (
        "Examplea rossica",
        "pollination_functional_guild",
        "bumblebees",
    ) in observed
    assert set(candidates["campaign_task"]) == {RECOVERY_TASK}
    assert candidates["target_for_task"].astype(bool).all()


def test_run_recovery_wave_reopens_downstream_tasks(tmp_path: Path) -> None:
    campaign_dir = tmp_path / "campaign"
    campaign_dir.mkdir()
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Campanula microdonta",
                "genus": "Campanula",
                "family": "Campanulaceae",
                "n_islands": "2",
                "n_records": "10",
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "processed",
                "floral_access_wikimedia_status": "not_applicable",
                "alternative_guild_wikimedia_status": "not_applicable",
                "alternative_guild_openalex_status": "not_applicable",
            },
            {
                "accepted_species": "Poa annua",
                "genus": "Poa",
                "family": "Poaceae",
                "n_islands": "3",
                "n_records": "20",
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "pending",
                "floral_access_wikimedia_status": "pending",
                "alternative_guild_wikimedia_status": "pending",
                "alternative_guild_openalex_status": "pending",
            },
        ]
    )
    ledger.to_csv(campaign_dir / "campaign_ledger.csv.gz", index=False, compression="gzip")
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps({"n_machine_biotic_candidates": 0}) + "\n",
        encoding="utf-8",
    )

    def getter(url: str, params: dict) -> dict:
        title = params["titles"]
        if "ja.wikipedia.org" in url:
            return {
                "query": {
                    "pages": {
                        "1": {
                            "title": title,
                            "extract": "送粉者はマルハナバチである。",
                        }
                    }
                }
            }
        return {"query": {"pages": {"-1": {"title": title, "missing": ""}}}}

    summary = run_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=100,
        languages=("ja", "zh", "ru"),
        getter=getter,
    )

    assert summary["n_species_attempted"] == 1
    assert summary["n_new_biotic_species"] == 1
    updated = pd.read_csv(campaign_dir / "campaign_ledger.csv.gz", dtype=str).fillna("")
    campanula = updated.loc[updated["accepted_species"].eq("Campanula microdonta")].iloc[0]
    poa = updated.loc[updated["accepted_species"].eq("Poa annua")].iloc[0]
    assert campanula["machine_biotic_candidate"].lower() == "true"
    assert campanula["floral_access_wikimedia_status"] == "pending"
    assert campanula["alternative_guild_wikimedia_status"] == "pending"
    assert campanula[f"{RECOVERY_TASK}_status"] == "processed"
    assert poa[f"{RECOVERY_TASK}_status"] == "pending"

    waves = list((campaign_dir / "waves").glob("wave_*_reproductive_multilingual_recovery"))
    assert len(waves) == 1
    candidates = pd.read_csv(waves[0] / "machine_candidates.csv.gz", dtype=str).fillna("")
    assert {
        ("pollination_functional_guild", "bumblebees"),
        ("pollen_vector_mode", "biotic"),
    } <= set(candidates[["trait_name", "candidate_value"]].itertuples(index=False, name=None))
