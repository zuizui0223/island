from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest
import typer

from island_v2 import trait_acquisition_priority as acquisition_priority
from island_v2.v1_category_traits import OUTPUT_COLUMNS
from island_v2.web_trait_shard_campaign import (
    PROVIDERS,
    _atomic_csv_gzip,
    _provider_policy_versions,
    _retryable_errors_by_species,
    _retryable_provider_errors,
    _select_checkpoint_batch,
    _validate_priority_audit,
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


def _priority_audit(
    names: list[str],
    regions: list[str],
    *,
    unresolved: list[int] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    unresolved = unresolved or [15] * len(names)
    for name, region, unresolved_cells in zip(names, regions, unresolved, strict=True):
        row: dict[str, object] = {
            column: "unknown" for column in acquisition_priority.PRIORITY_AUDIT_COLUMNS
        }
        row.update(
            {
                "accepted_species": name,
                "scheduling_region_category": region,
                "any_reviewed_NH_temperate": str(
                    region == "reviewed_nh_temperate"
                ).lower(),
                "any_NH_temperate_unreviewed": str(
                    region == "nh_temperate_unreviewed"
                ).lower(),
                "unresolved_trait_cells": unresolved_cells,
                "potential_genus_rescue_cells": 0,
                "frozen_species_direct_source": "false",
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


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


def _empty_collector(
    species: list[str],
    **_: object,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    del species
    return pd.DataFrame(), pd.DataFrame(columns=["species", "source", "error"])


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


def test_atomic_priority_csv_is_content_deterministic(tmp_path: Path) -> None:
    frame = pd.DataFrame(
        {"accepted_species": ["Plantus alpha", "Plantus beta"], "priority": [1, 2]}
    )
    first = tmp_path / "first.csv.gz"
    second = tmp_path / "second.csv.gz"
    _atomic_csv_gzip(frame, first)
    _atomic_csv_gzip(frame, second)

    assert hashlib.sha256(first.read_bytes()).digest() == hashlib.sha256(
        second.read_bytes()
    ).digest()


def test_priority_selection_puts_geography_before_yield_and_checkpoint_status() -> None:
    names = ["Plantus southern", "Plantus reviewed", "Plantus temperate"]
    checkpoint = pd.DataFrame(
        {
            "species": names,
            "status": ["pending", "retry", "pending"],
            "attempts": [0, 2, 0],
            "updated_at": ["", "2026-07-16T00:00:00Z", ""],
        }
    )
    audit = _priority_audit(
        names,
        ["southern_hemisphere", "reviewed_nh_temperate", "nh_temperate_unreviewed"],
        unresolved=[15, 1, 2],
    )

    selected = _select_checkpoint_batch(
        checkpoint,
        batch_size=3,
        retry_exhausted=False,
        priority_audit=audit,
    )

    assert selected["species"].tolist() == [
        "Plantus reviewed",
        "Plantus temperate",
        "Plantus southern",
    ]


def test_priority_audit_requires_the_exact_master_species_set() -> None:
    master = pd.DataFrame(
        {
            "accepted_species": ["Plantus alpha", "Plantus beta"],
            "genus": ["Plantus", "Plantus"],
            "family": ["Exampleaceae", "Exampleaceae"],
        }
    )
    incomplete = _priority_audit(
        ["Plantus alpha"],
        ["reviewed_nh_temperate"],
    )

    with pytest.raises(typer.BadParameter, match="exact campaign denominator"):
        _validate_priority_audit(incomplete, master, label="test audit")


def test_run_shard_consumes_frozen_priority_audit(tmp_path: Path) -> None:
    names = ["Plantus southern", "Plantus reviewed", "Plantus temperate"]
    master = _master(tmp_path / "master.csv", names)
    audit_path = tmp_path / "acquisition_priority.csv.gz"
    _priority_audit(
        names,
        ["southern_hemisphere", "reviewed_nh_temperate", "nh_temperate_unreviewed"],
        unresolved=[15, 1, 2],
    ).to_csv(audit_path, index=False, compression="gzip")
    seen: list[str] = []

    def recording_collector(
        species: list[str],
        **_: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        seen.extend(species)
        return _successful_collector(species)

    campaign = tmp_path / "campaign"
    report = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=2,
        expected_species=0,
        priority_audit_csv=audit_path,
        include_wikimedia=False,
        include_web_descriptions=False,
        collector=recording_collector,
    )

    assert seen == ["Plantus reviewed", "Plantus temperate"]
    assert report["acquisition_priority"]["enabled"] is True
    packet_priority = pd.read_csv(
        campaign / "packets" / "batch_000001" / "acquisition_priority.csv",
        dtype=str,
    )
    assert packet_priority["accepted_species"].tolist() == seen
    manifest = json.loads((campaign / "campaign_manifest.json").read_text(encoding="utf-8"))
    assert manifest["acquisition_priority"]["sha256"]


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
    provider_checkpoint = pd.read_csv(campaign / "provider_checkpoint.csv", dtype=str)
    assert len(provider_checkpoint) == 3 * 6
    assert not provider_checkpoint.duplicated(["species", "provider"]).any()
    assert set(provider_checkpoint["policy_version"]) == set(_provider_policy_versions().values())

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


def test_blank_provider_failure_fans_out_to_selected_species() -> None:
    errors = pd.DataFrame([{"species": "", "source": "openalex", "error": "HTTP 503"}])

    observed = _retryable_errors_by_species(
        errors,
        ["Plantus alpha", "Plantus beta"],
    )

    assert observed == {
        "Plantus alpha": ["HTTP 503"],
        "Plantus beta": ["HTTP 503"],
    }


def test_wikipedia_language_error_only_retries_wikimedia_provider() -> None:
    errors = pd.DataFrame(
        [{"species": "Plantus alpha", "source": "wikipedia_en", "error": "timeout"}]
    )
    running_pairs = {
        ("Plantus alpha", "gbif"),
        ("Plantus alpha", "europe_pmc"),
        ("Plantus alpha", "wikimedia"),
    }

    observed = _retryable_provider_errors(
        errors,
        ["Plantus alpha"],
        running_pairs,
    )

    assert observed == {("Plantus alpha", "wikimedia"): ["timeout"]}


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
        include_world_flora=False,
        collector=_successful_collector,
    )

    assert report["status_counts"] == {"completed": 1, "pending": 1}
    report = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=1,
        expected_species=0,
        migrate_v4=True,
        include_world_flora=False,
        collector=_successful_collector,
    )
    assert report["status_counts"] == {"completed": 2}
    results = pd.read_csv(cumulative / "trait_results.csv", dtype=str)
    assert set(results["species"]) == {"Plantus alpha", "Plantus beta"}
    manifest = json.loads((campaign / "campaign_manifest.json").read_text(encoding="utf-8"))
    assert manifest["contract_version"] == "public_web_9col_shards_v6"
    assert manifest["migrated_from"]["contract_version"] == "public_web_9col_shards_v4"


def test_v5_completed_species_requeues_new_provider_and_preserves_evidence(
    tmp_path: Path,
) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    campaign = tmp_path / "campaign"
    cumulative = campaign / "cumulative"
    cumulative.mkdir(parents=True)
    (campaign / "campaign_manifest.json").write_text(
        json.dumps(
            {
                "contract_version": "public_web_9col_shards_v5",
                "master_fingerprint": "legacy-fingerprint",
                "n_global_species": 1,
                "shard_index": 0,
                "shard_count": 1,
                "enabled_providers": ["gbif", "wikimedia"],
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
    pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alpha",
                "source_text": "Legacy flowers are red and tubular.",
                "source_url": "https://legacy.example/alpha",
                "source_citation": "Legacy flora",
                "source_type": "flora_or_monograph",
                "evidence_scope": "species_direct",
            }
        ]
    ).to_csv(cumulative / "search_evidence.csv", index=False)

    report = run_shard(
        master_csv=master,
        campaign_dir=campaign,
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=1,
        expected_species=0,
        include_gbif=False,
        include_wikimedia=False,
        include_openalex=True,
        include_web_descriptions=False,
        include_world_flora=False,
        collector=_empty_collector,
    )

    assert report["n_species_attempted_this_run"] == 1
    assert report["status_counts"] == {"completed": 1}
    assert report["provider_status_counts"] == {"openalex": {"no_hit": 1}}
    results = pd.read_csv(cumulative / "trait_results.csv", dtype=str)
    assert results.loc[0, "flower_color"] == "red"
    assert results.loc[0, "flower_shape"] == "tubular"
    assert results.loc[0, "pollination_notes"] == "legacy source-backed row"
    evidence = pd.read_csv(cumulative / "search_evidence.csv", dtype=str)
    assert evidence.loc[0, "source_url"] == "https://legacy.example/alpha"
    provider_checkpoint = pd.read_csv(campaign / "provider_checkpoint.csv", dtype=str)
    openalex = provider_checkpoint.loc[provider_checkpoint["provider"].eq("openalex")].iloc[0]
    assert openalex["status"] == "no_hit"
    assert openalex["policy_version"] == _provider_policy_versions()["openalex"]
    packet_provider = pd.read_csv(
        campaign / "packets" / "batch_000002" / "provider_status.csv",
        dtype=str,
    )
    assert set(packet_provider["provider"]) == {"openalex"}
    manifest = json.loads((campaign / "campaign_manifest.json").read_text(encoding="utf-8"))
    assert manifest["contract_version"] == "public_web_9col_shards_v6"
    assert manifest["migrated_from"]["contract_version"] == "public_web_9col_shards_v5"
    assert manifest["provider_policy"]["providers"]["openalex"] is True


def test_provider_wide_failure_retries_then_exhausts_every_selected_species(
    tmp_path: Path,
) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha", "Plantus beta"])
    campaign = tmp_path / "campaign"

    def provider_failure(
        species: list[str],
        **_: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        del species
        return (
            pd.DataFrame(),
            pd.DataFrame([{"species": "", "source": "openalex", "error": "HTTP 503"}]),
        )

    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 2,
        "max_attempts": 2,
        "expected_species": 0,
        "include_gbif": False,
        "include_wikimedia": False,
        "include_openalex": True,
        "include_web_descriptions": False,
        "include_world_flora": False,
        "collector": provider_failure,
    }
    first = run_shard(**arguments)
    assert first["status_counts"] == {"retry": 2}
    assert first["provider_status_counts"] == {"openalex": {"retry": 2}}
    assert first["complete"] is False

    second = run_shard(**arguments)
    assert second["status_counts"] == {"exhausted": 2}
    assert second["provider_status_counts"] == {"openalex": {"exhausted": 2}}
    assert second["complete"] is True
    checkpoint = pd.read_csv(campaign / "checkpoint.csv", dtype=str)
    assert not checkpoint["status"].eq("completed").any()


def test_exhausted_repair_rotates_to_the_least_attempted_species(tmp_path: Path) -> None:
    names = ["Plantus alpha", "Plantus beta", "Plantus gamma"]
    master = _master(tmp_path / "master.csv", names)
    campaign = tmp_path / "campaign"
    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 3,
        "max_attempts": 1,
        "expected_species": 0,
        "include_gbif": False,
        "include_wikimedia": False,
        "include_openalex": True,
        "include_web_descriptions": False,
        "include_world_flora": False,
    }

    def fail_all(
        species: list[str],
        **_: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        return (
            pd.DataFrame(),
            pd.DataFrame(
                [
                    {"species": name, "source": "openalex", "error": "HTTP 503"}
                    for name in species
                ]
            ),
        )

    run_shard(**arguments, collector=fail_all)
    checkpoint_path = campaign / "checkpoint.csv"
    checkpoint = pd.read_csv(checkpoint_path, dtype=str).fillna("")
    checkpoint.loc[checkpoint["species"].eq("Plantus alpha"), "attempts"] = "5"
    checkpoint.loc[checkpoint["species"].eq("Plantus beta"), "attempts"] = "2"
    checkpoint.loc[checkpoint["species"].eq("Plantus gamma"), "attempts"] = "3"
    checkpoint.to_csv(checkpoint_path, index=False)

    seen: list[str] = []

    def record_one(
        species: list[str],
        **_: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        seen.extend(species)
        return pd.DataFrame(), pd.DataFrame()

    run_shard(
        **{**arguments, "batch_size": 1},
        retry_exhausted=True,
        collector=record_one,
    )
    assert seen == ["Plantus beta"]


def test_each_provider_receives_only_species_with_pending_work(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha", "Plantus beta"])
    campaign = tmp_path / "campaign"
    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 2,
        "expected_species": 0,
        "include_gbif": True,
        "include_wikimedia": False,
        "include_openalex": True,
        "include_web_descriptions": False,
        "include_world_flora": False,
    }
    run_shard(**arguments, collector=_empty_collector)

    provider_path = campaign / "provider_checkpoint.csv"
    providers = pd.read_csv(provider_path, dtype=str)
    providers.loc[
        providers["species"].eq("Plantus alpha") & providers["provider"].eq("gbif"),
        "status",
    ] = "pending"
    providers.loc[
        providers["species"].eq("Plantus beta") & providers["provider"].eq("openalex"),
        "status",
    ] = "pending"
    providers.to_csv(provider_path, index=False)
    checkpoint_path = campaign / "checkpoint.csv"
    checkpoint = pd.read_csv(checkpoint_path, dtype=str)
    checkpoint["status"] = "pending"
    checkpoint.to_csv(checkpoint_path, index=False)

    calls: list[tuple[str, tuple[str, ...]]] = []

    def recording_collector(
        species: list[str],
        **flags: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        provider = next(
            name.removeprefix("include_")
            for name, enabled in flags.items()
            if name.startswith("include_") and enabled is True
        )
        calls.append((provider, tuple(species)))
        return pd.DataFrame(), pd.DataFrame(columns=["species", "source", "error"])

    run_shard(**arguments, collector=recording_collector)

    assert calls == [
        ("gbif", ("Plantus alpha",)),
        ("openalex", ("Plantus beta",)),
    ]


def test_changed_provider_policy_version_requeues_terminal_provider(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    campaign = tmp_path / "campaign"
    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 1,
        "expected_species": 0,
        "include_gbif": False,
        "include_wikimedia": False,
        "include_openalex": True,
        "include_web_descriptions": False,
        "include_world_flora": False,
        "collector": _empty_collector,
    }
    first = run_shard(**arguments)
    assert first["provider_status_counts"] == {"openalex": {"no_hit": 1}}
    provider_path = campaign / "provider_checkpoint.csv"
    providers = pd.read_csv(provider_path, dtype=str)
    providers.loc[providers["provider"].eq("openalex"), "policy_version"] = "stale_policy"
    providers.to_csv(provider_path, index=False)

    second = run_shard(**arguments)

    assert second["n_species_attempted_this_run"] == 1
    providers = pd.read_csv(provider_path, dtype=str)
    openalex = providers.loc[providers["provider"].eq("openalex")].iloc[0]
    assert openalex["status"] == "no_hit"
    assert openalex["policy_version"] == _provider_policy_versions()["openalex"]
    assert openalex["attempts"] == "1"


def test_legacy_sitemap_digest_policy_suffix_is_normalized_without_requeue(
    tmp_path: Path,
) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    campaign = tmp_path / "campaign"
    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 1,
        "expected_species": 0,
        "include_gbif": False,
        "include_wikimedia": False,
        "include_openalex": False,
        "include_europe_pmc": False,
        "include_web_descriptions": True,
        "include_world_flora": False,
        "collector": _empty_collector,
    }
    run_shard(**arguments)
    provider_path = campaign / "provider_checkpoint.csv"
    providers = pd.read_csv(provider_path, dtype=str)
    mask = providers["provider"].eq("web_descriptions")
    providers.loc[mask, "policy_version"] += ":0123456789abcdef"
    providers.to_csv(provider_path, index=False)

    second = run_shard(**arguments)

    assert second["n_species_attempted_this_run"] == 0
    normalized = pd.read_csv(provider_path, dtype=str).loc[mask, "policy_version"].iloc[0]
    assert normalized == _provider_policy_versions()["web_descriptions"]


def test_invalid_restored_provider_status_fails_closed(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    campaign = tmp_path / "campaign"
    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 1,
        "expected_species": 0,
        "include_gbif": False,
        "include_wikimedia": False,
        "include_openalex": True,
        "include_web_descriptions": False,
        "include_world_flora": False,
        "collector": _empty_collector,
    }
    run_shard(**arguments)
    provider_path = campaign / "provider_checkpoint.csv"
    providers = pd.read_csv(provider_path, dtype=str)
    providers.loc[providers["provider"].eq("openalex"), "status"] = "corrupt"
    providers.to_csv(provider_path, index=False)

    with pytest.raises(typer.BadParameter, match="invalid statuses"):
        run_shard(**arguments)


def test_sitemap_snapshot_change_is_provenance_without_global_requeue(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    campaign = tmp_path / "campaign"
    first_index = tmp_path / "web-index-1.json"
    second_index = tmp_path / "web-index-2.json"
    first_index.write_text('{"snapshot": 1}', encoding="utf-8")
    second_index.write_text('{"snapshot": 2}', encoding="utf-8")
    calls: list[tuple[str, Path | None]] = []

    def recording_collector(
        species: list[str],
        **flags: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        del species
        provider = next(
            name.removeprefix("include_")
            for name, enabled in flags.items()
            if name.startswith("include_") and enabled is True
        )
        calls.append((provider, flags.get("web_description_index_path")))
        return pd.DataFrame(), pd.DataFrame(columns=["species", "source", "error"])

    arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 1,
        "expected_species": 0,
        "include_gbif": True,
        "include_wikimedia": False,
        "include_openalex": False,
        "include_europe_pmc": False,
        "include_web_descriptions": True,
        "include_world_flora": False,
        "collector": recording_collector,
    }
    run_shard(**arguments, web_description_index_path=first_index)
    assert calls == [("gbif", None), ("web_descriptions", first_index)]

    calls.clear()
    run_shard(**arguments, web_description_index_path=second_index)

    assert calls == []
    providers = pd.read_csv(campaign / "provider_checkpoint.csv", dtype=str)
    versions = providers.set_index("provider")["policy_version"].to_dict()
    assert set(versions) == set(PROVIDERS)
    assert versions["gbif"] == _provider_policy_versions()["gbif"]
    assert versions["web_descriptions"] == _provider_policy_versions()["web_descriptions"]
    manifest = json.loads((campaign / "campaign_manifest.json").read_text(encoding="utf-8"))
    assert manifest["provider_policy"]["web_description_index_sha256"]


def test_seed_complete_provider_is_skipped_without_becoming_no_hit(tmp_path: Path) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    candidates = tmp_path / "candidates.csv"
    pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alpha",
                "trait_name": "flower_primary_color",
                "candidate_value": "red_pink",
            },
            {
                "accepted_species": "Plantus alpha",
                "trait_name": "floral_form",
                "candidate_value": "tubular",
            },
        ]
    ).to_csv(candidates, index=False)

    def unexpected_collector(
        species: list[str],
        **_: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        raise AssertionError(f"covered species should not be queried: {species}")

    report = run_shard(
        master_csv=master,
        campaign_dir=tmp_path / "campaign",
        scope_config_path=CONFIG,
        shard_index=0,
        shard_count=1,
        batch_size=1,
        expected_species=0,
        candidate_csv=candidates,
        include_gbif=True,
        include_wikimedia=False,
        include_openalex=False,
        include_europe_pmc=False,
        include_web_descriptions=False,
        include_world_flora=False,
        collector=unexpected_collector,
    )

    assert report["provider_status_counts"] == {"gbif": {"skipped_covered": 1}}
    result = pd.read_csv(
        tmp_path / "campaign" / "cumulative" / "trait_results.csv",
        dtype=str,
    )
    assert result.loc[0, "flower_color"] == "red/pink"
    assert result.loc[0, "flower_shape"] == "tubular"


def test_cumulative_source_ranking_prevents_weaker_provider_regression(
    tmp_path: Path,
) -> None:
    master = _master(tmp_path / "master.csv", ["Plantus alpha"])
    campaign = tmp_path / "campaign"

    def ranked_collector(
        species: list[str],
        **flags: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        if flags.get("include_gbif"):
            source_type = "gbif_species_description"
            text = "Flowers are red."
            citation = "GBIF species description"
        else:
            source_type = "wikipedia_extract"
            text = "Plantus alpha has blue flowers."
            citation = "Plantus alpha Wikipedia extract (en)"
        return (
            pd.DataFrame(
                [
                    {
                        "accepted_species": species[0],
                        "source_text": text,
                        "source_url": "https://example.test/source",
                        "source_citation": citation,
                        "source_type": source_type,
                        "evidence_scope": "species_direct",
                    }
                ]
            ),
            pd.DataFrame(columns=["species", "source", "error"]),
        )

    base_arguments = {
        "master_csv": master,
        "campaign_dir": campaign,
        "scope_config_path": CONFIG,
        "shard_index": 0,
        "shard_count": 1,
        "batch_size": 1,
        "expected_species": 0,
        "include_gbif": True,
        "include_openalex": False,
        "include_europe_pmc": False,
        "include_web_descriptions": False,
        "include_world_flora": False,
        "collector": ranked_collector,
    }
    run_shard(**base_arguments, include_wikimedia=False)
    run_shard(**base_arguments, include_wikimedia=True)

    result = pd.read_csv(
        campaign / "cumulative" / "trait_results.csv",
        dtype=str,
    )
    assert result.loc[0, "flower_color"] == "red"
    assert result.loc[0, "evidence_type"] == "flora"


def test_valid_shard_checkpoints_are_uploaded_even_after_later_step_failure() -> None:
    workflow = Path(".github/workflows/run-public-web-trait-shards.yml").read_text(encoding="utf-8")
    assert (
        "- name: Upload current search evidence packet and cumulative snapshot\n"
        "        if: always() && steps.status.outputs.ready == 'true'"
    ) in workflow
