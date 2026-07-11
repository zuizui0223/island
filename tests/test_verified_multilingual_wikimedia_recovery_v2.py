from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.verified_multilingual_wikimedia_recovery_v2 import (
    WAVE_TASK,
    fetch_lead_language_extracts,
    run_lead_full_recovery_wave,
)


def test_lead_extracts_use_twenty_page_batches() -> None:
    species_to_title = {
        f"Examplea species{index}": f"現地名{index}"
        for index in range(25)
    }
    batch_sizes: list[int] = []

    def getter(url: str, params: dict) -> dict:
        assert "ja.wikipedia.org" in url
        assert params["exintro"] == 1
        assert params["exlimit"] == 20
        requested = params["titles"].split("|")
        batch_sizes.append(len(requested))
        return {
            "query": {
                "pages": {
                    str(index): {
                        "title": title,
                        "extract": "雌雄異株である。",
                    }
                    for index, title in enumerate(requested)
                }
            }
        }

    sources, completed, errors = fetch_lead_language_extracts(
        getter,
        species_to_title,
        "ja",
    )

    assert batch_sizes == [20, 5]
    assert len(sources) == 25
    assert len(completed) == 25
    assert errors == []


def test_lead_then_full_fallback_recovers_and_rejects_uncertainty(
    tmp_path: Path,
) -> None:
    campaign_dir = tmp_path / "campaign"
    campaign_dir.mkdir()
    species = [
        "Examplea dioica",
        "Examplea incompatibilis",
        "Examplea incerta",
    ]
    pd.DataFrame(
        [
            {
                "accepted_species": name,
                "genus": "Examplea",
                "family": "Exampleaceae",
                "n_islands": "1",
                "n_records": str(10 - index),
                "machine_biotic_candidate": "False",
                "reproductive_wikimedia_status": "processed",
                "floral_access_wikimedia_status": "not_applicable",
                "alternative_guild_wikimedia_status": "not_applicable",
                "alternative_guild_openalex_status": "not_applicable",
            }
            for index, name in enumerate(species)
        ]
    ).to_csv(
        campaign_dir / "campaign_ledger.csv.gz",
        index=False,
        compression="gzip",
    )
    (campaign_dir / "campaign_status.json").write_text(
        json.dumps({"n_machine_biotic_candidates": 0}) + "\n",
        encoding="utf-8",
    )

    full_requests: list[str] = []

    def getter(url: str, params: dict) -> dict:
        if "en.wikipedia.org" in url:
            return {
                "query": {
                    "pages": {
                        str(index): {
                            "title": name,
                            "pageprops": {"wikibase_item": f"Q{index}"},
                        }
                        for index, name in enumerate(species)
                    }
                }
            }
        if "wikidata.org" in url:
            entities = {}
            for qid in params["ids"].split("|"):
                index = int(qid[1:])
                entities[qid] = {
                    "claims": {
                        "P225": [
                            {
                                "mainsnak": {
                                    "datavalue": {"value": species[index]}
                                }
                            }
                        ]
                    },
                    "sitelinks": {
                        "jawiki": {"title": f"現地名{index}"},
                    },
                }
            return {"entities": entities}
        if "ja.wikipedia.org" in url and params.get("exintro") == 1:
            texts = {
                "現地名0": "雌雄異株である。",
                "現地名1": "多年草である。",
                "現地名2": "化石植物である。",
            }
            requested = params["titles"].split("|")
            return {
                "query": {
                    "pages": {
                        str(index): {
                            "title": title,
                            "extract": texts[title],
                        }
                        for index, title in enumerate(requested)
                    }
                }
            }
        if "ja.wikipedia.org" in url:
            title = params["titles"]
            full_requests.append(title)
            texts = {
                "現地名1": "繁殖様式として自家不和合性が報告される。",
                "現地名2": "雌雄同株であるかは不明である。",
            }
            return {
                "query": {
                    "pages": {
                        "1": {"title": title, "extract": texts[title]}
                    }
                }
            }
        raise AssertionError((url, params))

    summary = run_lead_full_recovery_wave(
        campaign_dir=campaign_dir,
        batch_size=3,
        languages=("ja",),
        getter=getter,
    )

    assert summary["task"] == WAVE_TASK
    assert summary["n_species_with_lead_text"] == 3
    assert summary["n_species_with_lead_candidates"] == 1
    assert summary["n_species_full_fallback_attempted"] == 2
    assert summary["n_species_with_full_text"] == 2
    assert summary["n_candidates_total"] == 2
    assert summary["n_species_with_candidates"] == 2
    assert summary["n_uncertain_candidates_rejected"] == 1
    assert full_requests == ["現地名1", "現地名2"]

    wave_dir = next(
        (campaign_dir / "waves").glob(
            "wave_*_reproductive_multilingual_verified_lead_full_v2"
        )
    )
    candidates = pd.read_csv(
        wave_dir / "machine_candidates.csv.gz",
        dtype=str,
    ).fillna("")
    observed = set(
        candidates[["accepted_species", "trait_name", "candidate_value"]]
        .itertuples(index=False, name=None)
    )
    assert (
        "Examplea dioica",
        "sex_system",
        "dioecious",
    ) in observed
    assert (
        "Examplea incompatibilis",
        "self_incompatibility",
        "SI",
    ) in observed
    assert not candidates["accepted_species"].eq("Examplea incerta").any()
