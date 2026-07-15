"""Run resumable public-web trait extraction over deterministic species shards.

Each invocation advances one bounded batch in one shard. Public source text is
frozen before deterministic nine-column extraction, and transient per-species
lookup failures remain retryable. A successful zero-hit search is completed as
``unknown``; it is never interpreted as biological absence.
"""

from __future__ import annotations

import hashlib
import json
import os
from collections.abc import Callable, Mapping
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config
from island_v2.v1_category_search import (
    EVIDENCE_COLUMNS,
    SOURCE_COLUMNS,
    collapse_v1_rows,
    compact_source_texts,
    collect_free_sources,
    evidence_from_bulk_candidates,
    write_search_outputs,
)
from island_v2.v1_category_traits import OUTPUT_COLUMNS, validate_result_table

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT_VERSION = "public_web_9col_shards_v6"
LEGACY_CONTRACT_VERSION = "public_web_9col_shards_v5"
V4_CONTRACT_VERSION = "public_web_9col_shards_v4"
PROVIDER_POLICY_VERSION = "public_web_provider_accounting_v1"
PROVIDERS = (
    "gbif",
    "world_flora",
    "web_descriptions",
    "wikimedia",
    "openalex",
    "europe_pmc",
)
PROVIDER_IMPLEMENTATION_VERSIONS = {
    "gbif": "gbif_exact_species_v1",
    "world_flora": "world_flora_v1",
    "web_descriptions": "web_descriptions_indexed_v1",
    "wikimedia": "wikimedia_exact_title_v1",
    "openalex": "openalex_title_abstract_v1",
    "europe_pmc": "europe_pmc_title_abstract_reproduction_v2",
}
PROVIDER_TARGET_FIELDS = {
    "gbif": frozenset({"flower_color", "flower_shape"}),
    "world_flora": frozenset(
        {
            "flower_color",
            "flower_shape",
            "mating_system",
            "self_incompatibility",
        }
    ),
    "web_descriptions": frozenset(
        {
            "flower_color",
            "flower_shape",
            "mating_system",
            "self_incompatibility",
        }
    ),
    "wikimedia": frozenset(
        {
            "flower_color",
            "flower_shape",
            "mating_system",
            "self_incompatibility",
        }
    ),
    "openalex": frozenset({"mating_system", "self_incompatibility"}),
    "europe_pmc": frozenset({"mating_system", "self_incompatibility"}),
}
ACTIVE_TARGET_FIELDS = frozenset(
    {"flower_color", "flower_shape", "mating_system", "self_incompatibility"}
)
# v4/v5 did not persist provider attempts. These were the production campaign's
# enabled providers; any provider outside this set remains pending on migration.
LEGACY_DEFAULT_ENABLED_PROVIDERS = {
    "gbif",
    "web_descriptions",
    "wikimedia",
}
# The legacy GBIF and web-description implementations changed materially in v6
# (strict taxon identity plus repaired/complete plant-site extraction). Preserve
# their old evidence, but do not treat their old species-level completion flag as
# proof that the repaired provider was attempted.
LEGACY_REUSABLE_PROVIDERS = {"wikimedia"}
CHECKPOINT_COLUMNS = [
    "species",
    "genus",
    "family",
    "shard_index",
    "status",
    "attempts",
    "last_packet_id",
    "last_error",
    "result_sha256",
    "updated_at",
]
PROVIDER_CHECKPOINT_COLUMNS = [
    "species",
    "provider",
    "policy_version",
    "enabled",
    "status",
    "attempts",
    "last_packet_id",
    "last_error",
    "updated_at",
]
RETRYABLE_STATUSES = {"pending", "retry"}
TERMINAL_STATUSES = {"completed", "exhausted"}
PROVIDER_STATUSES = {
    "pending",
    "running",
    "retry",
    "hit",
    "no_hit",
    "skipped_covered",
    "legacy_completed",
    "exhausted",
}
PROVIDER_TERMINAL_STATUSES = {
    "hit",
    "no_hit",
    "skipped_covered",
    "legacy_completed",
    "exhausted",
}

Collector = Callable[..., tuple[pd.DataFrame, pd.DataFrame]]


@app.callback()
def main() -> None:
    """Build, validate, checkpoint, and resume public-web trait shards."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_shard(species: str, shard_count: int) -> int:
    """Return a stable shard independent of input row order."""
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    return int(_sha_text(_text(species))[:16], 16) % shard_count


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _atomic_csv(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    table.to_csv(temporary, index=False)
    temporary.replace(path)


def _atomic_json(payload: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)


def load_campaign_master(master_csv: Path, scope_config_path: Path) -> pd.DataFrame:
    """Retain every family-classified angiosperm from the current 115,328-name master."""
    master = pd.read_csv(master_csv, dtype=str).fillna("")
    if "accepted_species" not in master.columns:
        raise typer.BadParameter("master CSV must contain accepted_species")
    if "genus" not in master.columns:
        master["genus"] = master["accepted_species"].map(
            lambda value: _text(value).split()[0] if _text(value) else ""
        )
    scope = classify_scope(master, load_scope_config(scope_config_path))
    eligible = set(scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"].map(_text))
    for column in ("accepted_species", "genus", "family"):
        if column not in master.columns:
            master[column] = ""
        master[column] = master[column].map(_text)
    return (
        master.loc[master["accepted_species"].isin(eligible)]
        .drop_duplicates("accepted_species", keep="first")
        .sort_values("accepted_species")
        .reset_index(drop=True)
    )


def build_shard_plan(master: pd.DataFrame, shard_count: int) -> pd.DataFrame:
    assignments = master["accepted_species"].map(lambda species: stable_shard(species, shard_count))
    counts = assignments.value_counts().reindex(range(shard_count), fill_value=0)
    return pd.DataFrame(
        {
            "shard_index": range(shard_count),
            "n_species": counts.astype(int).tolist(),
        }
    )


def _master_fingerprint(master: pd.DataFrame) -> str:
    text = "\n".join(
        "\t".join(
            (
                _text(row.accepted_species),
                _text(row.genus),
                _text(row.family),
            )
        )
        for row in master[["accepted_species", "genus", "family"]].itertuples(index=False)
    )
    return _sha_text(text)


def _shard_species(master: pd.DataFrame, shard_index: int, shard_count: int) -> pd.DataFrame:
    assigned = master["accepted_species"].map(lambda species: stable_shard(species, shard_count))
    return master.loc[assigned.eq(shard_index), ["accepted_species", "genus", "family"]].copy()


def _provider_policy(
    *,
    include_gbif: bool,
    include_wikimedia: bool,
    include_openalex: bool,
    include_europe_pmc: bool,
    include_web_descriptions: bool,
    include_world_flora: bool,
) -> dict[str, bool]:
    return {
        "gbif": include_gbif,
        "world_flora": include_world_flora,
        "web_descriptions": include_web_descriptions,
        "wikimedia": include_wikimedia,
        "openalex": include_openalex,
        "europe_pmc": include_europe_pmc,
    }


def _provider_policy_versions(web_description_index_sha256: str = "") -> dict[str, str]:
    """Return independently bumpable policy versions for every provider."""
    versions = {
        provider: f"{PROVIDER_POLICY_VERSION}:{PROVIDER_IMPLEMENTATION_VERSIONS[provider]}"
        for provider in PROVIDERS
    }
    snapshot = web_description_index_sha256[:16] or "live-sitemaps"
    versions["web_descriptions"] = f"{versions['web_descriptions']}:{snapshot}"
    return versions


def _manifest_enabled_providers(manifest: dict[str, Any]) -> set[str]:
    policy = manifest.get("provider_policy", {})
    providers = policy.get("providers", {}) if isinstance(policy, dict) else {}
    if isinstance(policy, dict) and "providers" in policy and isinstance(providers, dict):
        return {
            provider
            for provider, value in providers.items()
            if provider in PROVIDERS and value is True
        }
    listed = manifest.get("enabled_providers", [])
    if "enabled_providers" in manifest and isinstance(listed, list):
        return {str(value) for value in listed if str(value) in PROVIDERS}
    return set(LEGACY_DEFAULT_ENABLED_PROVIDERS)


def reconcile_checkpoint(
    shard_species: pd.DataFrame,
    existing: pd.DataFrame | None,
    shard_index: int,
    max_attempts: int,
) -> pd.DataFrame:
    """Merge a restored checkpoint and recover interrupted running rows."""
    base = shard_species.rename(columns={"accepted_species": "species"}).copy()
    if existing is not None and not existing.empty:
        old = existing.copy().fillna("")
        preserved = [
            column
            for column in CHECKPOINT_COLUMNS
            if column in old.columns and column not in {"genus", "family", "shard_index"}
        ]
        base = base.merge(
            old[preserved],
            on="species",
            how="left",
            validate="one_to_one",
        )

    defaults: dict[str, object] = {
        "status": "pending",
        "attempts": 0,
        "last_packet_id": "",
        "last_error": "",
        "result_sha256": "",
        "updated_at": "",
    }
    for column, default in defaults.items():
        if column not in base.columns:
            base[column] = default
        base[column] = base[column].fillna("")
        if column != "attempts":
            base[column] = base[column].replace("", default)
    base["attempts"] = pd.to_numeric(base["attempts"], errors="coerce").fillna(0).astype(int)
    base["shard_index"] = shard_index

    interrupted = base["status"].eq("running")
    if interrupted.any():
        retryable = interrupted & base["attempts"].lt(max_attempts)
        base.loc[retryable, "status"] = "retry"
        base.loc[interrupted & ~retryable, "status"] = "exhausted"
        base.loc[interrupted, "last_error"] = "interrupted_previous_run"
    return base[CHECKPOINT_COLUMNS].sort_values("species").reset_index(drop=True)


def reconcile_provider_checkpoint(
    shard_species: pd.DataFrame,
    existing: pd.DataFrame | None,
    species_checkpoint: pd.DataFrame,
    enabled_providers: set[str],
    max_attempts: int,
    *,
    legacy_enabled_providers: set[str] | None = None,
    provider_policy_versions: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Build the species-by-provider ledger and migrate legacy species states."""
    policy_versions = dict(provider_policy_versions or _provider_policy_versions())
    missing_versions = set(PROVIDERS).difference(policy_versions)
    blank_versions = {
        provider for provider, version in policy_versions.items() if not _text(version)
    }
    if missing_versions or blank_versions:
        raise typer.BadParameter(
            "provider policy versions must be nonblank for every provider; "
            f"missing={sorted(missing_versions)}, blank={sorted(blank_versions)}"
        )
    species = [_text(value) for value in shard_species["accepted_species"]]
    base = pd.DataFrame(
        [{"species": name, "provider": provider} for name in species for provider in PROVIDERS],
        columns=["species", "provider"],
    )
    has_existing = existing is not None and not existing.empty
    if has_existing:
        old = existing.copy().fillna("")
        duplicates = old.duplicated(["species", "provider"], keep=False)
        if duplicates.any():
            raise typer.BadParameter("provider checkpoint has duplicate species-provider rows")
        preserved = [
            column
            for column in PROVIDER_CHECKPOINT_COLUMNS
            if column in old.columns and column not in {"species", "provider", "enabled"}
        ]
        base = base.merge(
            old[["species", "provider", *preserved]],
            on=["species", "provider"],
            how="left",
            validate="one_to_one",
        )

    restored_policy_version = (
        base["policy_version"].fillna("").map(_text)
        if "policy_version" in base.columns
        else pd.Series("", index=base.index, dtype=str)
    )

    defaults: dict[str, object] = {
        "policy_version": "",
        "status": "pending",
        "attempts": 0,
        "last_packet_id": "",
        "last_error": "",
        "updated_at": "",
    }
    for column, default in defaults.items():
        if column not in base.columns:
            base[column] = default
        base[column] = base[column].fillna("")
        if column not in {
            "policy_version",
            "attempts",
            "last_packet_id",
            "last_error",
            "updated_at",
        }:
            base[column] = base[column].replace("", default)
    base["attempts"] = pd.to_numeric(base["attempts"], errors="coerce").fillna(0).astype(int)
    invalid_statuses = sorted(set(base["status"]).difference(PROVIDER_STATUSES))
    if invalid_statuses:
        raise typer.BadParameter(
            f"provider checkpoint contains invalid statuses: {invalid_statuses}"
        )

    if not has_existing and legacy_enabled_providers:
        legacy = set(legacy_enabled_providers)
        high_level = species_checkpoint.set_index("species", drop=False)
        for index, row in base.iterrows():
            provider = _text(row["provider"])
            if provider not in legacy or provider not in LEGACY_REUSABLE_PROVIDERS:
                continue
            name = _text(row["species"])
            if name not in high_level.index:
                continue
            species_status = _text(high_level.loc[name, "status"])
            if species_status == "completed":
                base.loc[index, "status"] = "legacy_completed"
            elif species_status == "exhausted":
                base.loc[index, "status"] = "exhausted"
                base.loc[index, "attempts"] = max_attempts
                base.loc[index, "last_error"] = _text(high_level.loc[name, "last_error"])
            elif species_status in {"retry", "running"}:
                base.loc[index, "status"] = "retry"

    expected_policy_version = base["provider"].map(policy_versions).map(_text)
    stale_policy = has_existing & restored_policy_version.ne(expected_policy_version)
    if stale_policy.any():
        base.loc[stale_policy, "status"] = "pending"
        base.loc[stale_policy, "attempts"] = 0
        base.loc[stale_policy, "last_packet_id"] = ""
        base.loc[stale_policy, "last_error"] = "policy_version_changed"
        base.loc[stale_policy, "updated_at"] = _now()
    base["policy_version"] = expected_policy_version
    base["enabled"] = base["provider"].map(
        lambda provider: "true" if provider in enabled_providers else "false"
    )

    interrupted = base["status"].eq("running")
    if interrupted.any():
        retryable = interrupted & base["attempts"].lt(max_attempts)
        base.loc[retryable, "status"] = "retry"
        base.loc[interrupted & ~retryable, "status"] = "exhausted"
        base.loc[interrupted, "last_error"] = "interrupted_previous_run"
        base.loc[interrupted, "updated_at"] = _now()
    return (
        base[PROVIDER_CHECKPOINT_COLUMNS]
        .sort_values(["species", "provider"])
        .reset_index(drop=True)
    )


def _sync_species_checkpoint_from_providers(
    checkpoint: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    enabled_providers: set[str],
) -> None:
    """Maintain the legacy species-level status as a view of provider states."""
    for index, row in checkpoint.iterrows():
        name = _text(row["species"])
        provider_rows = provider_checkpoint.loc[
            provider_checkpoint["species"].eq(name)
            & provider_checkpoint["provider"].isin(enabled_providers)
        ]
        statuses = set(provider_rows["status"].map(_text))
        if "running" in statuses:
            status = "running"
        elif "pending" in statuses:
            status = "pending"
        elif "retry" in statuses:
            status = "retry"
        elif "exhausted" in statuses:
            status = "exhausted"
        else:
            status = "completed"
        errors = [
            _text(value)
            for value in provider_rows.loc[
                provider_rows["status"].isin({"retry", "exhausted"}), "last_error"
            ]
            if _text(value)
        ]
        error = " | ".join(dict.fromkeys(errors))[:2000]
        if (
            _text(checkpoint.loc[index, "status"]) != status
            or _text(checkpoint.loc[index, "last_error"]) != error
        ):
            checkpoint.loc[index, "updated_at"] = _now()
        checkpoint.loc[index, "status"] = status
        checkpoint.loc[index, "last_error"] = error


def _next_packet_id(checkpoint: pd.DataFrame) -> str:
    sequences: list[int] = []
    for value in checkpoint["last_packet_id"]:
        text = _text(value)
        if text.startswith("batch_") and text[6:].isdigit():
            sequences.append(int(text[6:]))
    return f"batch_{max(sequences, default=0) + 1:06d}"


def _result_sha(row: pd.Series) -> str:
    return _sha_text("\x1f".join(_text(row[column]) for column in OUTPUT_COLUMNS))


def _merge_cumulative(
    destination: Path,
    incoming: pd.DataFrame,
    *,
    dedupe_subset: list[str] | None = None,
    keep: str = "last",
) -> pd.DataFrame:
    if destination.exists():
        existing = pd.read_csv(destination, dtype=str).fillna("")
        combined = pd.concat([existing, incoming], ignore_index=True, sort=False).fillna("")
    else:
        combined = incoming.copy().fillna("")
    if dedupe_subset:
        combined = combined.drop_duplicates(dedupe_subset, keep=keep)
    else:
        combined = combined.drop_duplicates(keep=keep)
    _atomic_csv(combined, destination)
    return combined


def _merge_trait_results_preserving_known(
    prior: pd.DataFrame,
    incoming: pd.DataFrame,
) -> pd.DataFrame:
    """Merge a new provider pass without erasing known values from older evidence."""
    trait_fields = [
        "flower_color",
        "flower_shape",
        "pollination_guild",
        "mating_system",
        "self_incompatibility",
    ]
    confidence_rank = {"low": 0, "medium": 1, "high": 2}
    records = {
        _text(record["species"]): dict(record)
        for record in prior[OUTPUT_COLUMNS].fillna("").to_dict("records")
    }
    for incoming_record in incoming[OUTPUT_COLUMNS].fillna("").to_dict("records"):
        name = _text(incoming_record["species"])
        old = records.get(name)
        merged = dict(incoming_record)
        if old is not None:
            incoming_all_unknown = all(
                _text(incoming_record[field]).casefold() in {"", "unknown"}
                for field in trait_fields
            )
            used_old = False
            for field in trait_fields:
                if _text(merged[field]).casefold() in {"", "unknown"} and _text(
                    old[field]
                ).casefold() not in {"", "unknown"}:
                    merged[field] = old[field]
                    used_old = True
            if used_old:
                evidence = {
                    value
                    for source in (old["evidence_type"], merged["evidence_type"])
                    for value in _text(source).split("|")
                    if value
                }
                merged["evidence_type"] = "|".join(sorted(evidence))
                merged["confidence"] = max(
                    (_text(old["confidence"]), _text(merged["confidence"])),
                    key=lambda value: confidence_rank.get(value, -1),
                )
                if incoming_all_unknown and _text(old["pollination_notes"]):
                    merged["pollination_notes"] = old["pollination_notes"]
        records[name] = merged
    combined = pd.DataFrame(
        [records[name] for name in sorted(records)],
        columns=OUTPUT_COLUMNS,
    )
    return validate_result_table(combined, expected_species=sorted(records))


def _status_report(
    checkpoint: pd.DataFrame,
    provider_checkpoint: pd.DataFrame,
    results: pd.DataFrame,
    *,
    enabled_providers: set[str],
    global_species: int,
    shard_index: int,
    shard_count: int,
    packet_id: str,
    attempted: int,
) -> dict[str, Any]:
    counts = {
        str(status): int(count)
        for status, count in checkpoint["status"].value_counts().sort_index().items()
    }
    coverage = {
        field: int(results[field].ne("unknown").sum()) if field in results else 0
        for field in (
            "flower_color",
            "flower_shape",
            "pollination_guild",
            "mating_system",
            "self_incompatibility",
        )
    }
    n_terminal = int(checkpoint["status"].isin(TERMINAL_STATUSES).sum())
    active_provider_rows = provider_checkpoint.loc[
        provider_checkpoint["provider"].isin(enabled_providers)
    ]
    provider_status_counts = {
        provider: {
            str(status): int(count)
            for status, count in group["status"].value_counts().sort_index().items()
        }
        for provider, group in active_provider_rows.groupby("provider", sort=True)
    }
    n_provider_terminal = int(active_provider_rows["status"].isin(PROVIDER_TERMINAL_STATUSES).sum())
    return {
        "contract_version": CONTRACT_VERSION,
        "provider_policy_version": PROVIDER_POLICY_VERSION,
        "enabled_providers": sorted(enabled_providers),
        "active_target_fields": sorted(ACTIVE_TARGET_FIELDS),
        "updated_at": _now(),
        "n_global_species": global_species,
        "shard_index": shard_index,
        "shard_count": shard_count,
        "n_species_in_shard": int(len(checkpoint)),
        "n_species_attempted_this_run": attempted,
        "last_packet_id": packet_id,
        "status_counts": counts,
        "provider_status_counts": provider_status_counts,
        "n_provider_terminal": n_provider_terminal,
        "n_provider_remaining": int(len(active_provider_rows) - n_provider_terminal),
        "n_terminal": n_terminal,
        "n_remaining": int(len(checkpoint) - n_terminal),
        "complete": n_terminal == len(checkpoint),
        "coverage": coverage,
        "interpretation": (
            "provider no_hit is a terminal successful zero-result lookup, while "
            "skipped_covered means seed evidence already covered every provider target "
            "field. Pollination guild remains in the compatibility output but is not an "
            "active acquisition target; none of these statuses is biological absence."
        ),
    }


def _packet_manifest(packet_dir: Path, packet_id: str, species: list[str]) -> dict[str, Any]:
    files: dict[str, dict[str, Any]] = {}
    for path in sorted(packet_dir.iterdir()):
        if not path.is_file() or path.name == "packet_manifest.json":
            continue
        rows: int | None = None
        if path.suffix == ".csv":
            rows = int(len(pd.read_csv(path, dtype=str)))
        files[path.name] = {
            "sha256": _file_sha256(path),
            "bytes": path.stat().st_size,
            "rows": rows,
        }
    return {
        "contract_version": CONTRACT_VERSION,
        "provider_policy_version": PROVIDER_POLICY_VERSION,
        "active_target_fields": sorted(ACTIVE_TARGET_FIELDS),
        "packet_id": packet_id,
        "created_at": _now(),
        "species": species,
        "n_species": len(species),
        "output_columns": OUTPUT_COLUMNS,
        "validation": "exact species set, exact nine-column order, and controlled enums",
        "files": files,
    }


def _candidate_subset(candidate_csv: Path | None, species: list[str]) -> pd.DataFrame | None:
    if candidate_csv is None:
        return None
    candidates = pd.read_csv(candidate_csv, dtype=str).fillna("")
    if "accepted_species" not in candidates.columns:
        raise typer.BadParameter("candidate CSV must contain accepted_species")
    return candidates.loc[candidates["accepted_species"].isin(species)].copy()


def _seed_covered_provider_pairs(
    seed_evidence: pd.DataFrame,
    running_pairs: set[tuple[str, str]],
) -> set[tuple[str, str]]:
    """Identify provider passes that need no lookup because seed fields are complete."""
    covered_fields = {
        (_text(row.species), _text(row.field))
        for row in seed_evidence[["species", "field"]].itertuples(index=False)
    }
    return {
        (species, provider)
        for species, provider in running_pairs
        if all((species, field) in covered_fields for field in PROVIDER_TARGET_FIELDS[provider])
    }


def _terminal_no_hit_error(error: str) -> bool:
    terminal_markers = (
        "no_sitelink",
        "no_results",
        "no_species_direct_reproductive_source",
        "certificate_verify_failed",
        "404 not found",
        "response exceeded",
    )
    return any(marker in error.casefold() for marker in terminal_markers)


def _provider_terminal_no_hit_error(error: str) -> bool:
    """Only true lookup misses are provider-terminal no-hit outcomes in v6."""
    return any(
        marker in error.casefold()
        for marker in (
            "no_sitelink",
            "no_results",
            "no_species_direct_reproductive_source",
            "404 not found",
        )
    )


def _error_species(record: dict[str, object], selected_species: set[str]) -> str:
    species = _text(record.get("species"))
    error = _text(record.get("error"))
    if species or not error:
        return species
    parts = error.split(":", 2)
    if len(parts) == 3 and parts[1] in selected_species:
        return parts[1]
    return ""


def _retryable_errors_by_species(
    errors: pd.DataFrame,
    selected_species: list[str],
) -> dict[str, list[str]]:
    """Fan provider-wide failures out instead of falsely completing the packet."""
    wanted = set(selected_species)
    result: dict[str, list[str]] = {}
    for record in errors.fillna("").to_dict("records"):
        error = _text(record.get("error"))
        if not error or _terminal_no_hit_error(error):
            continue
        species = _error_species(record, wanted)
        targets = [species] if species in wanted else selected_species if not species else []
        for target in targets:
            result.setdefault(target, []).append(error)
    return result


def _provider_from_error_source(source: object) -> str:
    value = _text(source).casefold()
    if value in PROVIDERS:
        return value
    if value == "gbif":
        return "gbif"
    if value == "world_flora_online":
        return "world_flora"
    if value == "openalex":
        return "openalex"
    if value == "europe_pmc":
        return "europe_pmc"
    if value in {"wikimedia", "wikidata", "wikipedia"} or value.startswith("wikipedia_"):
        return "wikimedia"
    if value in {"ncsu_plant_toolbox", "nzpcn_flora", "pfaf", "useful_tropical_plants"}:
        return "web_descriptions"
    return ""


def _provider_from_source_record(record: dict[str, object]) -> str:
    source_type = _text(record.get("source_type")).casefold()
    url = _text(record.get("source_url")).casefold()
    citation = _text(record.get("source_citation")).casefold()
    if source_type == "gbif_species_description":
        return "gbif"
    if source_type == "openalex_title_abstract":
        return "openalex"
    if source_type == "europe_pmc_title_abstract":
        return "europe_pmc"
    if source_type in {"wikipedia_extract", "wikipedia_sitelink_extract_ja"}:
        return "wikimedia"
    if "worldfloraonline.org" in url or "world flora online" in citation:
        return "world_flora"
    if source_type in {"horticulture_site", "flora_or_monograph"}:
        return "web_descriptions"
    return ""


def _retryable_provider_errors(
    errors: pd.DataFrame,
    selected_species: list[str],
    running_pairs: set[tuple[str, str]],
) -> dict[tuple[str, str], list[str]]:
    wanted = set(selected_species)
    result: dict[tuple[str, str], list[str]] = {}
    for record in errors.fillna("").to_dict("records"):
        error = _text(record.get("error"))
        if not error or _provider_terminal_no_hit_error(error):
            continue
        species = _error_species(record, wanted)
        targets = [species] if species in wanted else selected_species if not species else []
        provider = _provider_from_error_source(record.get("source"))
        for target in targets:
            target_providers = (
                [provider]
                if provider and (target, provider) in running_pairs
                else [candidate for name, candidate in running_pairs if name == target]
            )
            for target_provider in target_providers:
                result.setdefault((target, target_provider), []).append(error)
    return result


def _provider_hits(
    sources: pd.DataFrame,
    running_pairs: set[tuple[str, str]],
) -> set[tuple[str, str]]:
    hits: set[tuple[str, str]] = set()
    for record in sources.fillna("").to_dict("records"):
        species = _text(record.get("accepted_species"))
        provider = _provider_from_source_record(record)
        if (species, provider) in running_pairs:
            hits.add((species, provider))
            continue
        candidates = [pair for pair in running_pairs if pair[0] == species]
        if not provider and len(candidates) == 1:
            hits.add(candidates[0])
    return hits


def run_shard(
    *,
    master_csv: Path,
    campaign_dir: Path,
    scope_config_path: Path,
    shard_index: int,
    shard_count: int = 128,
    batch_size: int = 50,
    max_attempts: int = 3,
    expected_species: int = 106295,
    candidate_csv: Path | None = None,
    retry_exhausted: bool = False,
    migrate_v4: bool = False,
    include_gbif: bool = True,
    include_wikimedia: bool = True,
    include_openalex: bool = False,
    include_europe_pmc: bool = False,
    include_web_descriptions: bool = True,
    include_world_flora: bool = True,
    web_description_index_path: Path | None = None,
    pause_seconds: float = 0.1,
    collector: Collector | None = None,
) -> dict[str, Any]:
    """Advance one shard batch and persist everything needed for resume."""
    if not 0 <= shard_index < shard_count:
        raise typer.BadParameter("shard_index must be within shard_count")
    if batch_size < 1:
        raise typer.BadParameter("batch_size must be positive")
    if max_attempts < 1:
        raise typer.BadParameter("max_attempts must be positive")

    master = load_campaign_master(master_csv, scope_config_path)
    if expected_species and len(master) != expected_species:
        raise typer.BadParameter(
            f"global campaign denominator changed: expected {expected_species}, got {len(master)}"
        )
    fingerprint = _master_fingerprint(master)
    shard_species = _shard_species(master, shard_index, shard_count)
    campaign_dir.mkdir(parents=True, exist_ok=True)
    cumulative_dir = campaign_dir / "cumulative"
    packet_dir_root = campaign_dir / "packets"
    checkpoint_path = campaign_dir / "checkpoint.csv"
    provider_checkpoint_path = campaign_dir / "provider_checkpoint.csv"
    campaign_manifest_path = campaign_dir / "campaign_manifest.json"
    provider_policy = _provider_policy(
        include_gbif=include_gbif,
        include_wikimedia=include_wikimedia,
        include_openalex=include_openalex,
        include_europe_pmc=include_europe_pmc,
        include_web_descriptions=include_web_descriptions,
        include_world_flora=include_world_flora,
    )
    web_description_index_sha256 = (
        _file_sha256(web_description_index_path) if web_description_index_path is not None else ""
    )
    provider_policy_versions = _provider_policy_versions(web_description_index_sha256)
    enabled_providers = {provider for provider, enabled in provider_policy.items() if enabled}
    current_policy_manifest = {
        "version": PROVIDER_POLICY_VERSION,
        "providers": provider_policy,
        "provider_versions": provider_policy_versions,
        "web_description_index_sha256": web_description_index_sha256,
    }
    previous: dict[str, Any] = {}
    legacy_contract_version = ""
    legacy_enabled_providers: set[str] = set()
    if campaign_manifest_path.exists():
        previous = json.loads(campaign_manifest_path.read_text(encoding="utf-8"))
        previous_version = _text(previous.get("contract_version"))
        shard_mismatch = {
            key: (previous.get(key), value)
            for key, value in {
                "shard_index": shard_index,
                "shard_count": shard_count,
            }.items()
            if int(previous.get(key, -1)) != value
        }
        if shard_mismatch:
            raise typer.BadParameter(f"restored campaign manifest mismatch: {shard_mismatch}")
        if previous_version == CONTRACT_VERSION:
            immutable = {"master_fingerprint": fingerprint}
            mismatched = {
                key: (previous.get(key), value)
                for key, value in immutable.items()
                if previous.get(key) != value
            }
            if mismatched:
                raise typer.BadParameter(f"restored campaign manifest mismatch: {mismatched}")
        elif previous_version == LEGACY_CONTRACT_VERSION or (
            migrate_v4 and previous_version == V4_CONTRACT_VERSION
        ):
            legacy_contract_version = previous_version
            legacy_enabled_providers = _manifest_enabled_providers(previous)
        else:
            raise typer.BadParameter(
                "restored campaign manifest requires an explicit supported migration; "
                f"got {previous_version!r}"
            )

    policy_history = previous.get("provider_policy_history", [])
    if not isinstance(policy_history, list):
        policy_history = []
    previous_policy = previous.get("provider_policy")
    if (
        previous.get("contract_version") == CONTRACT_VERSION
        and isinstance(previous_policy, dict)
        and previous_policy != current_policy_manifest
        and previous_policy not in policy_history
    ):
        policy_history.append(previous_policy)
    manifest_payload: dict[str, Any] = {
        "contract_version": CONTRACT_VERSION,
        "created_at": previous.get("created_at") or _now(),
        "updated_at": _now(),
        "master_csv": str(master_csv),
        "master_fingerprint": fingerprint,
        "n_global_species": len(master),
        "shard_index": shard_index,
        "shard_count": shard_count,
        "n_species_in_shard": len(shard_species),
        "provider_policy": current_policy_manifest,
        "provider_policy_history": policy_history,
        "enabled_providers": sorted(enabled_providers),
        "source_policy": (
            "free public sources; bulk candidates seed gaps; provider attempts are "
            "versioned per species; no paid API; island Bombus state is never an input"
        ),
    }
    if legacy_contract_version:
        manifest_payload["migrated_from"] = {
            "contract_version": legacy_contract_version,
            "master_fingerprint": previous.get("master_fingerprint"),
            "n_global_species": previous.get("n_global_species"),
            "assumed_enabled_providers": sorted(legacy_enabled_providers),
        }
    elif isinstance(previous.get("migrated_from"), dict):
        manifest_payload["migrated_from"] = previous["migrated_from"]
    _atomic_json(manifest_payload, campaign_manifest_path)

    existing = (
        pd.read_csv(checkpoint_path, dtype=str).fillna("") if checkpoint_path.exists() else None
    )
    checkpoint = reconcile_checkpoint(shard_species, existing, shard_index, max_attempts)
    existing_provider_checkpoint = (
        pd.read_csv(provider_checkpoint_path, dtype=str).fillna("")
        if provider_checkpoint_path.exists()
        else None
    )
    provider_checkpoint = reconcile_provider_checkpoint(
        shard_species,
        existing_provider_checkpoint,
        checkpoint,
        enabled_providers,
        max_attempts,
        legacy_enabled_providers=(legacy_enabled_providers if legacy_contract_version else None),
        provider_policy_versions=provider_policy_versions,
    )
    _sync_species_checkpoint_from_providers(
        checkpoint,
        provider_checkpoint,
        enabled_providers,
    )
    allowed = set(RETRYABLE_STATUSES)
    if retry_exhausted:
        allowed.add("exhausted")
    eligible = checkpoint.loc[checkpoint["status"].isin(allowed)].copy()
    eligible["_priority"] = eligible["status"].map({"pending": 0, "retry": 1, "exhausted": 2})
    selected = eligible.sort_values(["_priority", "species"]).head(batch_size)
    packet_id = _next_packet_id(checkpoint)
    selected_indexes = selected.index.tolist()

    prior_results_path = cumulative_dir / "trait_results.csv"
    prior_results = (
        pd.read_csv(prior_results_path, dtype=str).fillna("")
        if prior_results_path.exists()
        else pd.DataFrame(columns=OUTPUT_COLUMNS)
    )
    prior_results = prior_results.loc[
        prior_results["species"].isin(set(shard_species["accepted_species"]))
    ].copy()
    if selected.empty:
        report = _status_report(
            checkpoint,
            provider_checkpoint,
            prior_results,
            enabled_providers=enabled_providers,
            global_species=len(master),
            shard_index=shard_index,
            shard_count=shard_count,
            packet_id="",
            attempted=0,
        )
        _atomic_csv(checkpoint, checkpoint_path)
        _atomic_csv(provider_checkpoint, provider_checkpoint_path)
        _atomic_json(report, campaign_dir / "campaign_status.json")
        return report

    species = [_text(value) for value in selected["species"]]
    for index in selected_indexes:
        checkpoint.loc[index, "status"] = "running"
        checkpoint.loc[index, "attempts"] = int(checkpoint.loc[index, "attempts"]) + 1
        checkpoint.loc[index, "last_packet_id"] = packet_id
        checkpoint.loc[index, "last_error"] = ""
        checkpoint.loc[index, "updated_at"] = _now()
    provider_allowed = set(RETRYABLE_STATUSES)
    if retry_exhausted:
        provider_allowed.add("exhausted")
    running_mask = (
        provider_checkpoint["species"].isin(species)
        & provider_checkpoint["provider"].isin(enabled_providers)
        & provider_checkpoint["status"].isin(provider_allowed)
    )
    if not running_mask.any():
        raise RuntimeError("selected species have no provider work under the active policy")
    provider_checkpoint.loc[running_mask, "status"] = "running"
    provider_checkpoint.loc[running_mask, "attempts"] = (
        provider_checkpoint.loc[running_mask, "attempts"].astype(int) + 1
    )
    provider_checkpoint.loc[running_mask, "last_packet_id"] = packet_id
    provider_checkpoint.loc[running_mask, "last_error"] = ""
    provider_checkpoint.loc[running_mask, "updated_at"] = _now()
    running_pairs = {
        (_text(row.species), _text(row.provider))
        for row in provider_checkpoint.loc[running_mask, ["species", "provider"]].itertuples(
            index=False
        )
    }
    attempted_providers = {provider for _, provider in running_pairs}
    _sync_species_checkpoint_from_providers(
        checkpoint,
        provider_checkpoint,
        enabled_providers,
    )
    _atomic_csv(checkpoint, checkpoint_path)
    _atomic_csv(provider_checkpoint, provider_checkpoint_path)

    packet_dir = packet_dir_root / packet_id
    packet_dir.mkdir(parents=True, exist_ok=False)
    pd.DataFrame({"species": species}).to_csv(packet_dir / "species.csv", index=False)
    candidates = _candidate_subset(candidate_csv, species)
    seed_evidence = (
        evidence_from_bulk_candidates(candidates)
        if candidates is not None
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )
    skipped_covered_pairs = _seed_covered_provider_pairs(seed_evidence, running_pairs)

    search = collector or collect_free_sources
    try:
        source_parts: list[pd.DataFrame] = []
        error_parts: list[pd.DataFrame] = []
        for provider in sorted(attempted_providers):
            provider_species = sorted(
                name
                for name, candidate_provider in running_pairs
                if candidate_provider == provider
                and (name, candidate_provider) not in skipped_covered_pairs
            )
            if not provider_species:
                continue
            provider_seed = seed_evidence.loc[
                seed_evidence["species"].isin(provider_species)
            ].copy()
            try:
                provider_sources, provider_errors = search(
                    provider_species,
                    include_gbif=provider == "gbif",
                    include_wikimedia=provider == "wikimedia",
                    include_openalex=provider == "openalex",
                    include_europe_pmc=provider == "europe_pmc",
                    include_web_descriptions=provider == "web_descriptions",
                    include_world_flora=provider == "world_flora",
                    web_description_index_path=(
                        web_description_index_path if provider == "web_descriptions" else None
                    ),
                    pause_seconds=pause_seconds,
                    seed_evidence=provider_seed,
                )
                source_parts.append(
                    pd.DataFrame(provider_sources, columns=SOURCE_COLUMNS).fillna("")
                )
                error_parts.append(
                    pd.DataFrame(
                        provider_errors,
                        columns=["species", "source", "error"],
                    ).fillna("")
                )
            except Exception as exc:  # noqa: BLE001
                error_parts.append(
                    pd.DataFrame(
                        [
                            {
                                "species": "",
                                "source": provider,
                                "error": f"{type(exc).__name__}: {exc}",
                            }
                        ],
                        columns=["species", "source", "error"],
                    )
                )
        sources = (
            pd.concat(source_parts, ignore_index=True).drop_duplicates()
            if source_parts
            else pd.DataFrame(columns=SOURCE_COLUMNS)
        )
        errors = (
            pd.concat(error_parts, ignore_index=True).drop_duplicates()
            if error_parts
            else pd.DataFrame(columns=["species", "source", "error"])
        )
        sources = compact_source_texts(sources)
        sources.to_csv(packet_dir / "source_texts.csv", index=False)
        errors.to_csv(packet_dir / "lookup_errors.csv", index=False)
        frozen_sources = pd.read_csv(
            packet_dir / "source_texts.csv",
            dtype=str,
        ).fillna("")
        write_search_outputs(
            packet_dir,
            species,
            frozen_sources,
            errors,
            bulk_candidates=candidates,
        )
        extracted = pd.read_csv(packet_dir / "v1_category_traits.csv", dtype=str).fillna("")
        extracted = validate_result_table(extracted, expected_species=species)
        if list(extracted.columns) != OUTPUT_COLUMNS:
            raise ValueError("validated result lost the exact nine-column contract")
    except Exception as exc:
        error = f"{type(exc).__name__}: {exc}"
        failed_mask = provider_checkpoint["status"].eq("running") & provider_checkpoint[
            "species"
        ].isin(species)
        for index in provider_checkpoint.index[failed_mask]:
            provider_checkpoint.loc[index, "status"] = (
                "exhausted"
                if int(provider_checkpoint.loc[index, "attempts"]) >= max_attempts
                else "retry"
            )
            provider_checkpoint.loc[index, "last_error"] = error
            provider_checkpoint.loc[index, "updated_at"] = _now()
        _sync_species_checkpoint_from_providers(
            checkpoint,
            provider_checkpoint,
            enabled_providers,
        )
        _atomic_csv(provider_checkpoint, provider_checkpoint_path)
        _atomic_csv(checkpoint, checkpoint_path)
        _atomic_json(
            {
                "contract_version": CONTRACT_VERSION,
                "packet_id": packet_id,
                "error": error,
                "species": species,
                "providers": sorted(attempted_providers),
                "retry_policy": "retry until max_attempts, then retain as exhausted",
            },
            packet_dir / "fatal_error.json",
        )
        raise

    errors = errors.copy().fillna("")
    retryable_provider_errors = _retryable_provider_errors(
        errors,
        species,
        running_pairs,
    )
    provider_hits = _provider_hits(sources, running_pairs)
    for index in provider_checkpoint.index[
        provider_checkpoint["status"].eq("running") & provider_checkpoint["species"].isin(species)
    ]:
        key = (
            _text(provider_checkpoint.loc[index, "species"]),
            _text(provider_checkpoint.loc[index, "provider"]),
        )
        provider_checkpoint.loc[index, "updated_at"] = _now()
        if key in skipped_covered_pairs:
            provider_checkpoint.loc[index, "status"] = "skipped_covered"
            provider_checkpoint.loc[index, "last_error"] = ""
        elif key in retryable_provider_errors:
            provider_checkpoint.loc[index, "last_error"] = " | ".join(
                dict.fromkeys(value for value in retryable_provider_errors[key] if value)
            )[:2000]
            provider_checkpoint.loc[index, "status"] = (
                "exhausted"
                if int(provider_checkpoint.loc[index, "attempts"]) >= max_attempts
                else "retry"
            )
        elif key in provider_hits:
            provider_checkpoint.loc[index, "status"] = "hit"
            provider_checkpoint.loc[index, "last_error"] = ""
        else:
            provider_checkpoint.loc[index, "status"] = "no_hit"
            provider_checkpoint.loc[index, "last_error"] = ""

    cumulative_dir.mkdir(parents=True, exist_ok=True)
    packet_evidence = pd.read_csv(
        packet_dir / "v1_category_evidence.csv",
        dtype=str,
    ).fillna("")
    cumulative_evidence = _merge_cumulative(
        cumulative_dir / "trait_evidence.csv",
        packet_evidence,
    )
    ranked_evidence = cumulative_evidence.loc[cumulative_evidence["species"].isin(species)].copy()
    extracted = collapse_v1_rows(species, ranked_evidence)
    extracted = validate_result_table(extracted, expected_species=species)
    cumulative_results = _merge_trait_results_preserving_known(
        prior_results,
        extracted,
    )
    _atomic_csv(cumulative_results, prior_results_path)
    for index in selected_indexes:
        name = _text(checkpoint.loc[index, "species"])
        row = cumulative_results.loc[cumulative_results["species"].eq(name)].iloc[0]
        checkpoint.loc[index, "result_sha256"] = _result_sha(row)
        checkpoint.loc[index, "updated_at"] = _now()
    _sync_species_checkpoint_from_providers(
        checkpoint,
        provider_checkpoint,
        enabled_providers,
    )
    file_specs = (
        ("source_texts.csv", "search_evidence.csv", SOURCE_COLUMNS, None),
        (
            "v2_reported_candidates.csv",
            "v2_reported_candidates.csv",
            None,
            None,
        ),
    )
    for packet_name, cumulative_name, columns, subset in file_specs:
        incoming_path = packet_dir / packet_name
        if incoming_path.exists():
            incoming = pd.read_csv(incoming_path, dtype=str).fillna("")
        else:
            incoming = pd.DataFrame(columns=columns or [])
        _merge_cumulative(
            cumulative_dir / cumulative_name,
            incoming,
            dedupe_subset=subset,
        )

    packet_errors = errors.copy()
    packet_errors["packet_id"] = packet_id
    _merge_cumulative(
        cumulative_dir / "lookup_errors.csv",
        packet_errors,
    )
    packet_provider_status = provider_checkpoint.loc[
        provider_checkpoint["species"].isin(species)
        & provider_checkpoint["last_packet_id"].eq(packet_id)
    ].copy()
    packet_provider_status.to_csv(packet_dir / "provider_status.csv", index=False)
    _atomic_csv(provider_checkpoint, provider_checkpoint_path)
    _atomic_csv(checkpoint, checkpoint_path)
    manifest = _packet_manifest(packet_dir, packet_id, species)
    _atomic_json(manifest, packet_dir / "packet_manifest.json")
    report = _status_report(
        checkpoint,
        provider_checkpoint,
        cumulative_results,
        enabled_providers=enabled_providers,
        global_species=len(master),
        shard_index=shard_index,
        shard_count=shard_count,
        packet_id=packet_id,
        attempted=len(species),
    )
    report["packet_manifest"] = str(packet_dir / "packet_manifest.json")
    report["run_id"] = os.environ.get("GITHUB_RUN_ID", "")
    _atomic_json(report, campaign_dir / "campaign_status.json")
    return report


@app.command("plan")
def plan_command(
    master_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    scope_config_path: Path = typer.Option(Path("config/angiosperm_scope.yml"), exists=True),
    shard_count: int = typer.Option(128, min=1, max=4096),
    expected_species: int = typer.Option(106295, min=0),
) -> None:
    """Validate the denominator and write deterministic shard counts."""
    master = load_campaign_master(master_csv, scope_config_path)
    if expected_species and len(master) != expected_species:
        raise typer.BadParameter(
            f"global campaign denominator changed: expected {expected_species}, got {len(master)}"
        )
    shard_plan = build_shard_plan(master, shard_count)
    output_dir.mkdir(parents=True, exist_ok=True)
    shard_plan.to_csv(output_dir / "shard_plan.csv", index=False)
    report = {
        "contract_version": CONTRACT_VERSION,
        "provider_policy_version": PROVIDER_POLICY_VERSION,
        "providers": list(PROVIDERS),
        "n_global_species": len(master),
        "shard_count": shard_count,
        "min_species_per_shard": int(shard_plan["n_species"].min()),
        "max_species_per_shard": int(shard_plan["n_species"].max()),
        "master_fingerprint": _master_fingerprint(master),
    }
    _atomic_json(report, output_dir / "shard_plan.json")
    typer.echo(json.dumps(report, ensure_ascii=False))


@app.command("run")
def run_command(
    master_csv: Path = typer.Option(..., exists=True),
    campaign_dir: Path = typer.Option(...),
    scope_config_path: Path = typer.Option(Path("config/angiosperm_scope.yml"), exists=True),
    shard_index: int = typer.Option(..., min=0),
    shard_count: int = typer.Option(128, min=1, max=4096),
    batch_size: int = typer.Option(50, min=1, max=500),
    max_attempts: int = typer.Option(3, min=1, max=10),
    expected_species: int = typer.Option(106295, min=0),
    candidate_csv: Path | None = typer.Option(None, exists=True),
    retry_exhausted: bool = typer.Option(False),
    migrate_v4: bool = typer.Option(False),
    gbif: bool = typer.Option(True, "--gbif/--no-gbif"),
    wikimedia: bool = typer.Option(True, "--wikimedia/--no-wikimedia"),
    openalex: bool = typer.Option(False, "--openalex/--no-openalex"),
    europe_pmc: bool = typer.Option(False, "--europe-pmc/--no-europe-pmc"),
    web_descriptions: bool = typer.Option(
        True,
        "--web-descriptions/--no-web-descriptions",
    ),
    world_flora: bool = typer.Option(True, "--world-flora/--no-world-flora"),
    web_description_index: Path | None = typer.Option(None, exists=True),
    pause_seconds: float = typer.Option(0.1, min=0.0),
) -> None:
    """Advance one shard batch; restored checkpoints resume automatically."""
    report = run_shard(
        master_csv=master_csv,
        campaign_dir=campaign_dir,
        scope_config_path=scope_config_path,
        shard_index=shard_index,
        shard_count=shard_count,
        batch_size=batch_size,
        max_attempts=max_attempts,
        expected_species=expected_species,
        candidate_csv=candidate_csv,
        retry_exhausted=retry_exhausted,
        migrate_v4=migrate_v4,
        include_gbif=gbif,
        include_wikimedia=wikimedia,
        include_openalex=openalex,
        include_europe_pmc=europe_pmc,
        include_web_descriptions=web_descriptions,
        include_world_flora=world_flora,
        web_description_index_path=web_description_index,
        pause_seconds=pause_seconds,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
