from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.v1_category_traits import OUTPUT_COLUMNS
from island_v2.web_trait_shard_campaign import (
    _retryable_errors_by_species,
    build_shard_plan,
    reconcile_checkpoint,
    run_shard,
    stable_shard,
)

CONFIG = Path("config/angiosperm_scope.yml")


def _master(path: Path, names: list[str]) -> Path:
    pd.DataFrame(
        {
            "accepted_species": names,
            "family": ["Exampleaceae"] * len(names),
            "n_islands": [1] * len(names),
            "n_records": [1] * len(names),
        }
    ).to_csv(path, index=False)
    return path


def _successful_collector(species: list[str], **_: object) -> tuple[pd.DataFrame, pd.DataFrame]:
    sources = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "source_text": (
                    "Flowers are red and tubular. The species is pollinated by bees "
                    "and shows autonomous selfing."
                ),
                "source_url": f"https://example.org/{name.replace(' ', '-')}",
                "source_citation": "Example flora",
                "source_type": "flora_or_monograph",
                "evidence_scope": "species_direct",
            }
            for name in species
        ]
    )
    return sources, pd.DataFrame(columns=["species", "source", "error"])


def test_stable_shards_are_deterministic_and_plan_preserves_denominator() -> None:
    names = [f"Plantus species{i}" for i in range(100)]
    master = pd.DataFrame(
        {
            "accepted_species": names,
            "genus": ["Plantus"] * len(names),
            "family": ["Exampleaceae"] * len(names),
        }
    )

    assert stable_shard("Campanula punctata", 128) == stable_shard("Campanula punctata", 128)
    plan = build_shard_plan(master, 8)

    assert len(plan) == 8
    assert plan["n_species"].sum() == 100
    assert set(plan["shard_index"]) == set(range(8))


def test_run_shard_freezes_packet_validates_rows_and_resumes(tmp_path: Path) -> None:
    master = _master(
        tmp_path / "master.csv",
        ["Plantus alpha", "Plantus beta", "Plantus gamma"],
    )
    campaign = tmp_path / "campaign"

    first = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=2,
        expected_species=0,
        collector=_successful_collector,
    )

    assert first["n_species_attempted_this_run"] == 2
    assert first["status_counts"] == {"completed": 2, "pending": 1}
    packet = campaign / "packets" / "batch_000001"
    results = pd.read_csv(packet / "v1_category_traits.csv", dtype=str)
    assert list(results.columns) == OUTPUT_COLUMNS
    assert set(results["flower_color"]) == {"red"}
    assert set(results["flower_shape"]) == {"tubular"}
    assert (packet / "source_texts.csv").exists()
    manifest = json.loads((packet / "packet_manifest.json").read_text(encoding="utf-8"))
    assert manifest["n_species"] == 2
    assert manifest["files"]["source_texts.csv"]["sha256"]

    second = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=2,
        expected_species=0,
        collector=_successful_collector,
    )

    assert second["n_species_attempted_this_run"] == 1
    assert second["complete"] is True
    assert second["status_counts"] == {"completed": 3}
    cumulative = pd.read_csv(campaign / "cumulative" / "trait_results.csv", dtype=str)
    assert list(cumulative.columns) == OUTPUT_COLUMNS
    assert cumulative["species"].nunique() == 3
    assert (campaign / "packets" / "batch_000002").exists()

    complete = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=2,
        expected_species=0,
        collector=_successful_collector,
    )
    assert complete["n_species_attempted_this_run"] == 0


def test_per_species_lookup_error_retries_then_exhausts(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus retryii"])
    campaign = tmp_path / "campaign"

    def failing_collector(species: list[str], **_: object) -> tuple[pd.DataFrame, pd.DataFrame]:
        return (
            pd.DataFrame(
                columns=[
                    "accepted_species",
                    "source_text",
                    "source_url",
                    "source_citation",
                    "source_type",
                    "evidence_scope",
                ]
            ),
            pd.DataFrame([{"species": species[0], "source": "world_flora", "error": "timeout"}]),
        )

    first = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=1,
        max_attempts=2,
        expected_species=0,
        collector=failing_collector,
    )
    assert first["status_counts"] == {"retry": 1}

    second = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=1,
        max_attempts=2,
        expected_species=0,
        collector=failing_collector,
    )
    assert second["status_counts"] == {"exhausted": 1}
    assert second["complete"] is True
    checkpoint = pd.read_csv(campaign / "checkpoint.csv", dtype=str)
    assert checkpoint.loc[0, "attempts"] == "2"
    assert checkpoint.loc[0, "last_error"] == "timeout"


def test_interrupted_rows_are_recovered_for_resume() -> None:
    shard_species = pd.DataFrame(
        {
            "accepted_species": ["Plantus interrupted"],
            "genus": ["Plantus"],
            "family": ["Exampleaceae"],
        }
    )
    existing = pd.DataFrame(
        [
            {
                "species": "Plantus interrupted",
                "status": "running",
                "attempts": "1",
                "last_packet_id": "batch_000001",
            }
        ]
    )

    checkpoint = reconcile_checkpoint(shard_species, existing, 7, 3)

    assert checkpoint.loc[0, "status"] == "retry"
    assert checkpoint.loc[0, "last_error"] == "interrupted_previous_run"
    assert checkpoint.loc[0, "shard_index"] == 7


def test_wikimedia_error_strings_retry_transient_failures_not_missing_pages() -> None:
    errors = pd.DataFrame(
        [
            {
                "species": "",
                "source": "wikimedia",
                "error": "wikidata:Plantus retryii:HTTP 503",
            },
            {
                "species": "",
                "source": "wikimedia",
                "error": "wikipedia:Plantus absentia:no_sitelink",
            },
            {
                "species": "Plantus certificateii",
                "source": "world_flora_online",
                "error": "[SSL: CERTIFICATE_VERIFY_FAILED] unable to get local issuer certificate",
            },
        ]
    )

    observed = _retryable_errors_by_species(
        errors,
        ["Plantus retryii", "Plantus absentia", "Plantus certificateii"],
    )

    assert observed == {"Plantus retryii": ["wikidata:Plantus retryii:HTTP 503"]}


def test_v4_checkpoint_migrates_and_queues_newly_eligible_species(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha", "Plantus beta"])
    campaign = tmp_path / "campaign"
    cumulative = campaign / "cumulative"
    cumulative.mkdir(parents=True)
    (campaign / "campaign_manifest.json").write_text(
        json.dumps(
            {
                "contract_version": "public_web_9col_shards_v4",
                "master_fingerprint": "legacy-fingerprint",
                "n_global_species": 1,
                "shard_index": 0,
                "shard_count": 1,
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame(
        [
            {
                "species": "Plantus alpha",
                "genus": "Plantus",
                "family": "Exampleaceae",
                "shard_index": 0,
                "status": "completed",
                "attempts": 1,
                "last_packet_id": "batch_000001",
                "last_error": "",
                "result_sha256": "legacy-result",
                "updated_at": "2026-07-12T00:00:00Z",
            }
        ]
    ).to_csv(campaign / "checkpoint.csv", index=False)
    pd.DataFrame(
        [
            {
                "species": "Plantus alpha",
                "flower_color": "red",
                "flower_shape": "tubular",
                "pollination_guild": "bees",
                "pollination_notes": "legacy source-backed row",
                "mating_system": "unknown",
                "self_incompatibility": "unknown",
                "evidence_type": "flora",
                "confidence": "medium",
            }
        ],
        columns=OUTPUT_COLUMNS,
    ).to_csv(cumulative / "trait_results.csv", index=False)

    report = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=1,
        expected_species=0,
        migrate_v4=True,
        collector=_successful_collector,
    )

    assert report["status_counts"] == {"completed": 2}
    results = pd.read_csv(cumulative / "trait_results.csv", dtype=str)
    assert set(results["species"]) == {"Plantus alpha", "Plantus beta"}
    manifest = json.loads((campaign / "campaign_manifest.json").read_text(encoding="utf-8"))
    assert manifest["contract_version"] == "public_web_9col_shards_v5"
    assert manifest["migrated_from"]["contract_version"] == "public_web_9col_shards_v4"
