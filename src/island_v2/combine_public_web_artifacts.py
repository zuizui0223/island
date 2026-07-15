"""Promote a complete public-web shard run into one immutable combined artifact."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.v1_category_traits import OUTPUT_COLUMNS, validate_result_table
from island_v2.web_trait_shard_campaign import (
    PROVIDERS,
    PROVIDER_CHECKPOINT_COLUMNS,
    PROVIDER_STATUSES,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

TABLE_SPECS: dict[str, tuple[str, str | None]] = {
    "trait_results": ("trait_results.csv", "species"),
    "trait_evidence": ("trait_evidence.csv", None),
    "search_evidence": ("search_evidence.csv", None),
    "v2_reported_candidates": ("v2_reported_candidates.csv", None),
    "lookup_errors": ("lookup_errors.csv", None),
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_metadata(path: Path | None) -> list[dict[str, Any]]:
    if path is None:
        return []
    records: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            records.append(json.loads(line))
    return records


def _discover_shards(root: Path) -> list[tuple[int, Path, dict[str, Any]]]:
    shards: list[tuple[int, Path, dict[str, Any]]] = []
    for manifest_path in root.rglob("campaign_manifest.json"):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if "shard_index" not in manifest or "shard_count" not in manifest:
            continue
        shards.append((int(manifest["shard_index"]), manifest_path.parent, manifest))
    return sorted(shards, key=lambda item: item[0])


def _combine_optional_tables(
    shards: list[tuple[int, Path, dict[str, Any]]],
) -> dict[str, pd.DataFrame]:
    combined: dict[str, pd.DataFrame] = {}
    for key, (name, dedupe_column) in TABLE_SPECS.items():
        parts: list[pd.DataFrame] = []
        for _, shard_dir, _ in shards:
            path = shard_dir / "cumulative" / name
            if path.exists():
                parts.append(pd.read_csv(path, dtype=str).fillna(""))
        table = (
            pd.concat(parts, ignore_index=True, sort=False).fillna("") if parts else pd.DataFrame()
        )
        if dedupe_column and not table.empty:
            if table[dedupe_column].duplicated().any():
                duplicates = table.loc[
                    table[dedupe_column].duplicated(keep=False), dedupe_column
                ].unique()
                raise ValueError(f"duplicate {dedupe_column} across shards: {duplicates[:10]}")
        elif not table.empty:
            table = table.drop_duplicates().reset_index(drop=True)
        combined[key] = table
    return combined


def combine_public_web_run(
    *,
    artifact_root: Path,
    output_dir: Path,
    source_run_id: int,
    expected_shards: int = 128,
    expected_species: int = 106295,
    expected_contract: str = "public_web_9col_shards_v6",
    artifact_metadata_jsonl: Path | None = None,
) -> dict[str, Any]:
    shards = _discover_shards(artifact_root)
    indexes = [index for index, _, _ in shards]
    if indexes != list(range(expected_shards)):
        raise ValueError(f"expected shard indexes 0..{expected_shards - 1}, got {indexes[:10]}...")
    contracts = {_text(manifest.get("contract_version")) for _, _, manifest in shards}
    fingerprints = {_text(manifest.get("master_fingerprint")) for _, _, manifest in shards}
    shard_counts = {int(manifest.get("shard_count", -1)) for _, _, manifest in shards}
    denominators = {int(manifest.get("n_global_species", -1)) for _, _, manifest in shards}
    if contracts != {expected_contract}:
        raise ValueError(f"shard contract mismatch: {contracts}")
    if len(fingerprints) != 1 or "" in fingerprints:
        raise ValueError(f"shard master fingerprint mismatch: {fingerprints}")
    if shard_counts != {expected_shards}:
        raise ValueError(f"shard-count mismatch: {shard_counts}")
    if denominators != {expected_species}:
        raise ValueError(f"source denominator mismatch: {denominators}")

    source_provider_policy: dict[str, Any] | None = None
    if expected_contract == "public_web_9col_shards_v6":
        policies: list[dict[str, Any]] = []
        for _, shard_dir, manifest in shards:
            policy = manifest.get("provider_policy")
            if not isinstance(policy, dict):
                raise ValueError(f"v6 shard has no provider policy: {shard_dir}")
            if not _text(policy.get("version")):
                raise ValueError(f"provider policy has no version: {shard_dir}")
            providers = policy.get("providers")
            versions = policy.get("provider_versions")
            if not isinstance(providers, dict) or set(providers) != set(PROVIDERS):
                raise ValueError(f"provider policy has an invalid provider set: {shard_dir}")
            if any(type(value) is not bool for value in providers.values()):
                raise ValueError(f"provider policy enabled flags are not booleans: {shard_dir}")
            if not isinstance(versions, dict) or set(versions) != set(PROVIDERS):
                raise ValueError(
                    f"provider policy versions have an invalid provider set: {shard_dir}"
                )
            if any(not _text(value) for value in versions.values()):
                raise ValueError(f"provider policy contains a blank version: {shard_dir}")
            enabled_providers = manifest.get("enabled_providers")
            expected_enabled_providers = sorted(
                provider for provider, enabled in providers.items() if enabled
            )
            if (
                not isinstance(enabled_providers, list)
                or sorted(_text(provider) for provider in enabled_providers)
                != expected_enabled_providers
            ):
                raise ValueError(
                    f"manifest enabled providers disagree with provider policy: {shard_dir}"
                )
            policies.append(policy)
        canonical_policies = {
            json.dumps(policy, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
            for policy in policies
        }
        if len(canonical_policies) != 1:
            raise ValueError("v6 shards do not share one provider policy/index snapshot")
        source_provider_policy = policies[0]

    checkpoint_parts: list[pd.DataFrame] = []
    for index, shard_dir, _ in shards:
        checkpoint = pd.read_csv(shard_dir / "checkpoint.csv", dtype=str).fillna("")
        if not checkpoint["shard_index"].astype(int).eq(index).all():
            raise ValueError(f"checkpoint rows assigned to the wrong shard: {index}")
        if checkpoint["species"].duplicated().any():
            raise ValueError(f"duplicate species inside checkpoint shard {index}")
        checkpoint_parts.append(checkpoint)
    checkpoint = pd.concat(checkpoint_parts, ignore_index=True)
    if len(checkpoint) != expected_species or checkpoint["species"].nunique() != expected_species:
        raise ValueError(
            f"checkpoint population mismatch: rows={len(checkpoint)}, "
            f"unique={checkpoint['species'].nunique()}"
        )

    provider_checkpoint = pd.DataFrame(columns=PROVIDER_CHECKPOINT_COLUMNS)
    if expected_contract == "public_web_9col_shards_v6":
        provider_parts: list[pd.DataFrame] = []
        assert source_provider_policy is not None
        expected_versions = source_provider_policy["provider_versions"]
        expected_enabled = source_provider_policy["providers"]
        for _, shard_dir, _ in shards:
            path = shard_dir / "provider_checkpoint.csv"
            if not path.exists():
                raise ValueError(f"v6 shard has no provider checkpoint: {shard_dir}")
            part = pd.read_csv(path, dtype=str).fillna("")
            missing = set(PROVIDER_CHECKPOINT_COLUMNS).difference(part.columns)
            if missing:
                raise ValueError(f"provider checkpoint missing columns: {sorted(missing)}")
            part = part[PROVIDER_CHECKPOINT_COLUMNS]
            if not set(part["provider"]).issubset(PROVIDERS):
                raise ValueError(f"provider checkpoint contains an unknown provider: {shard_dir}")
            expected_part_versions = part["provider"].map(expected_versions).map(_text)
            if not part["policy_version"].map(_text).eq(expected_part_versions).all():
                raise ValueError(
                    f"provider checkpoint policy versions disagree with manifest: {shard_dir}"
                )
            expected_part_enabled = part["provider"].map(
                lambda provider: "true" if expected_enabled[provider] else "false"
            )
            if not part["enabled"].map(_text).str.casefold().eq(expected_part_enabled).all():
                raise ValueError(
                    f"provider checkpoint enabled flags disagree with manifest: {shard_dir}"
                )
            invalid_statuses = sorted(set(part["status"]).difference(PROVIDER_STATUSES))
            if invalid_statuses:
                raise ValueError(
                    f"provider checkpoint contains invalid statuses: {invalid_statuses}"
                )
            provider_parts.append(part)
        provider_checkpoint = pd.concat(provider_parts, ignore_index=True)
        if provider_checkpoint.duplicated(["species", "provider"]).any():
            raise ValueError("duplicate species-provider rows across shards")
        if set(provider_checkpoint["species"]) != set(checkpoint["species"]):
            raise ValueError("provider checkpoint species do not match checkpoint")
        if set(provider_checkpoint["provider"]) != set(PROVIDERS):
            raise ValueError("provider checkpoint provider set does not match v6 contract")
        expected_provider_rows = expected_species * len(PROVIDERS)
        if len(provider_checkpoint) != expected_provider_rows:
            raise ValueError(
                f"provider checkpoint rows={len(provider_checkpoint)}, "
                f"expected={expected_provider_rows}"
            )
        if provider_checkpoint["policy_version"].eq("").any():
            raise ValueError("provider checkpoint contains blank policy versions")

    tables = _combine_optional_tables(shards)
    results = tables["trait_results"]
    if not results.empty:
        if list(results.columns) != OUTPUT_COLUMNS:
            raise ValueError("combined trait_results lost the exact nine-column order")
        foreign_species = sorted(set(results["species"]).difference(checkpoint["species"]))
        if foreign_species:
            raise ValueError(
                f"combined trait_results contain species outside checkpoint: {foreign_species[:10]}"
            )
        results = validate_result_table(results, expected_species=results["species"].tolist())
        tables["trait_results"] = results

    metadata = _read_metadata(artifact_metadata_jsonl)
    source_artifacts: list[dict[str, Any]] = []
    if metadata:
        pattern = re.compile(rf"^public-web-trait-shard-(\d+)-{source_run_id}$")
        for record in metadata:
            match = pattern.match(_text(record.get("name")))
            if not match:
                continue
            if bool(record.get("expired")):
                raise ValueError(f"source artifact expired: {record.get('name')}")
            digest = _text(record.get("digest"))
            if not digest:
                raise ValueError(f"source artifact has no digest: {record.get('name')}")
            run = record.get("workflow_run") or {}
            head_sha = _text(run.get("head_sha"))
            if not head_sha:
                raise ValueError(f"source artifact has no workflow head SHA: {record.get('name')}")
            source_artifacts.append(
                {
                    "shard_index": int(match.group(1)),
                    "artifact_id": int(record["id"]),
                    "name": record["name"],
                    "digest": digest,
                    "head_sha": head_sha,
                }
            )
        source_artifacts.sort(key=lambda item: item["shard_index"])
        if [item["shard_index"] for item in source_artifacts] != list(range(expected_shards)):
            raise ValueError("artifact metadata does not account for all shard artifacts")
        if len({item["head_sha"] for item in source_artifacts}) != 1:
            raise ValueError("source artifacts do not share one workflow head SHA")

    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output_dir / "combined_checkpoint.csv.gz"
    checkpoint.to_csv(checkpoint_path, index=False, compression="gzip")
    output_paths: dict[str, Path] = {"checkpoint": checkpoint_path}
    table_lengths: dict[str, int] = {"checkpoint": len(checkpoint)}
    if not provider_checkpoint.empty:
        provider_path = output_dir / "combined_provider_checkpoint.csv.gz"
        provider_checkpoint.to_csv(provider_path, index=False, compression="gzip")
        output_paths["provider_checkpoint"] = provider_path
        table_lengths["provider_checkpoint"] = len(provider_checkpoint)
    for key, table in tables.items():
        path = output_dir / f"combined_{key}.csv"
        table.to_csv(path, index=False)
        output_paths[key] = path
        table_lengths[key] = len(table)

    status_counts = {
        str(status): int(count)
        for status, count in checkpoint["status"].value_counts().sort_index().items()
    }
    fields = [
        "flower_color",
        "flower_shape",
        "pollination_guild",
        "mating_system",
        "self_incompatibility",
    ]
    coverage = {
        field: int(results[field].ne("unknown").sum()) if not results.empty else 0
        for field in fields
    }
    any_known = (
        results[fields].ne("unknown").any(axis=1) if not results.empty else pd.Series(dtype=bool)
    )
    provider_status_counts = (
        {
            provider: {
                str(status): int(count)
                for status, count in group["status"].value_counts().sort_index().items()
            }
            for provider, group in provider_checkpoint.groupby("provider", sort=True)
        }
        if not provider_checkpoint.empty
        else {}
    )
    files = {
        key: {
            "path": path.name,
            "rows": int(table_lengths[key]),
            "sha256": _file_sha256(path),
        }
        for key, path in output_paths.items()
    }
    manifest = {
        "version": "1.0",
        "source_run_id": source_run_id,
        "source_contract": expected_contract,
        "source_master_fingerprint": next(iter(fingerprints)),
        "source_provider_policy": source_provider_policy,
        "source_denominator": expected_species,
        "shard_count": expected_shards,
        "checkpoint_status_counts": status_counts,
        "provider_status_counts": provider_status_counts,
        "n_provider_checkpoint_rows": len(provider_checkpoint),
        "n_result_species": len(results),
        "n_species_any_requested_value": int(any_known.sum()) if len(any_known) else 0,
        "coverage": coverage,
        "n_trait_evidence_rows": len(tables["trait_evidence"]),
        "n_search_evidence_rows": len(tables["search_evidence"]),
        "n_v2_candidate_rows": len(tables["v2_reported_candidates"]),
        "n_lookup_error_rows": len(tables["lookup_errors"]),
        "source_artifact_count": len(source_artifacts),
        "source_artifacts": source_artifacts,
        "files": files,
        "interpretation": (
            "Provider no_hit rows are successful zero-result lookups; skipped_covered "
            "rows used complete seed coverage without a lookup. Neither is biological "
            "absence. URLs and excerpts remain in the combined evidence tables."
        ),
    }
    (output_dir / "combined_public_web_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("run")
def run_command(
    artifact_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
    source_run_id: int = typer.Option(..., min=1),
    expected_shards: int = typer.Option(128, min=1),
    expected_species: int = typer.Option(106295, min=1),
    expected_contract: str = typer.Option("public_web_9col_shards_v6"),
    artifact_metadata_jsonl: Path | None = typer.Option(None, exists=True),
) -> None:
    manifest = combine_public_web_run(
        artifact_root=artifact_root,
        output_dir=output_dir,
        source_run_id=source_run_id,
        expected_shards=expected_shards,
        expected_species=expected_species,
        expected_contract=expected_contract,
        artifact_metadata_jsonl=artifact_metadata_jsonl,
    )
    typer.echo(json.dumps(manifest, ensure_ascii=False))


if __name__ == "__main__":
    app()
