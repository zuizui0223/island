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
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.angiosperm_scope import classify_scope, load_config as load_scope_config
from island_v2.v1_category_search import (
    EVIDENCE_COLUMNS,
    SOURCE_COLUMNS,
    collect_free_sources,
    evidence_from_bulk_candidates,
    write_search_outputs,
)
from island_v2.v1_category_traits import OUTPUT_COLUMNS, validate_result_table

app = typer.Typer(add_completion=False, no_args_is_help=True)

CONTRACT_VERSION = "public_web_9col_shards_v5"
LEGACY_CONTRACT_VERSION = "public_web_9col_shards_v4"
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
RETRYABLE_STATUSES = {"pending", "retry"}
TERMINAL_STATUSES = {"completed", "exhausted"}

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
    eligible = set(
        scope.loc[scope["angiosperm_analysis_eligible"], "accepted_species"].map(_text)
    )
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


def _status_report(
    checkpoint: pd.DataFrame,
    results: pd.DataFrame,
    *,
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
    return {
        "contract_version": CONTRACT_VERSION,
        "updated_at": _now(),
        "n_global_species": global_species,
        "shard_index": shard_index,
        "shard_count": shard_count,
        "n_species_in_shard": int(len(checkpoint)),
        "n_species_attempted_this_run": attempted,
        "last_packet_id": packet_id,
        "status_counts": counts,
        "n_terminal": n_terminal,
        "n_remaining": int(len(checkpoint) - n_terminal),
        "complete": n_terminal == len(checkpoint),
        "coverage": coverage,
        "interpretation": (
            "completed unknown rows are successful zero-hit searches, not biological "
            "absences; exhausted rows retain failed lookup diagnostics."
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


def _write_failure_checkpoint(
    checkpoint: pd.DataFrame,
    selected_indexes: list[int],
    max_attempts: int,
    error: str,
    checkpoint_path: Path,
) -> None:
    for index in selected_indexes:
        checkpoint.loc[index, "status"] = (
            "exhausted" if int(checkpoint.loc[index, "attempts"]) >= max_attempts else "retry"
        )
        checkpoint.loc[index, "last_error"] = error
        checkpoint.loc[index, "updated_at"] = _now()
    _atomic_csv(checkpoint, checkpoint_path)


def _retryable_errors_by_species(
    errors: pd.DataFrame,
    selected_species: list[str],
) -> dict[str, list[str]]:
    """Recover species names from legacy Wikimedia error strings when possible."""
    wanted = set(selected_species)
    result: dict[str, list[str]] = {}
    for record in errors.fillna("").to_dict("records"):
        species = _text(record.get("species"))
        error = _text(record.get("error"))
        terminal_markers = (
            "no_sitelink",
            "certificate_verify_failed",
            "404 not found",
            "response exceeded",
        )
        if any(marker in error.casefold() for marker in terminal_markers):
            continue
        if not species and error and not error.endswith(":no_sitelink"):
            parts = error.split(":", 2)
            if len(parts) == 3 and parts[1] in wanted:
                species = parts[1]
        if species in wanted:
            result.setdefault(species, []).append(error)
    return result


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
    include_web_descriptions: bool = True,
    include_world_flora: bool = True,
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
    campaign_manifest_path = campaign_dir / "campaign_manifest.json"

    if campaign_manifest_path.exists():
        previous = json.loads(campaign_manifest_path.read_text(encoding="utf-8"))
        immutable = {
            "contract_version": CONTRACT_VERSION,
            "master_fingerprint": fingerprint,
            "shard_index": shard_index,
            "shard_count": shard_count,
        }
        mismatched = {
            key: (previous.get(key), value)
            for key, value in immutable.items()
            if previous.get(key) != value
        }
        if mismatched:
            legacy_migration = (
                migrate_v4
                and previous.get("contract_version") == LEGACY_CONTRACT_VERSION
                and int(previous.get("shard_index", -1)) == shard_index
                and int(previous.get("shard_count", -1)) == shard_count
            )
            if not legacy_migration:
                raise typer.BadParameter(f"restored campaign manifest mismatch: {mismatched}")
            _atomic_json(
                {
                    "contract_version": CONTRACT_VERSION,
                    "created_at": _now(),
                    "master_csv": str(master_csv),
                    "master_fingerprint": fingerprint,
                    "n_global_species": len(master),
                    "shard_index": shard_index,
                    "shard_count": shard_count,
                    "n_species_in_shard": len(shard_species),
                    "migrated_from": {
                        "contract_version": previous.get("contract_version"),
                        "master_fingerprint": previous.get("master_fingerprint"),
                        "n_global_species": previous.get("n_global_species"),
                    },
                    "source_policy": (
                        "free public sources; v4 evidence is retained only for current "
                        "angiosperm-eligible names; no paid API"
                    ),
                },
                campaign_manifest_path,
            )
    else:
        _atomic_json(
            {
                "contract_version": CONTRACT_VERSION,
                "created_at": _now(),
                "master_csv": str(master_csv),
                "master_fingerprint": fingerprint,
                "n_global_species": len(master),
                "shard_index": shard_index,
                "shard_count": shard_count,
                "n_species_in_shard": len(shard_species),
                "source_policy": (
                    "free public sources; bulk candidates seed gaps; no paid API; "
                    "island Bombus state is never a retrieval or inference input"
                ),
            },
            campaign_manifest_path,
        )

    existing = (
        pd.read_csv(checkpoint_path, dtype=str).fillna("") if checkpoint_path.exists() else None
    )
    checkpoint = reconcile_checkpoint(shard_species, existing, shard_index, max_attempts)
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
            prior_results,
            global_species=len(master),
            shard_index=shard_index,
            shard_count=shard_count,
            packet_id="",
            attempted=0,
        )
        _atomic_csv(checkpoint, checkpoint_path)
        _atomic_json(report, campaign_dir / "campaign_status.json")
        return report

    for index in selected_indexes:
        checkpoint.loc[index, "status"] = "running"
        checkpoint.loc[index, "attempts"] = int(checkpoint.loc[index, "attempts"]) + 1
        checkpoint.loc[index, "last_packet_id"] = packet_id
        checkpoint.loc[index, "last_error"] = ""
        checkpoint.loc[index, "updated_at"] = _now()
    _atomic_csv(checkpoint, checkpoint_path)

    species = [_text(value) for value in selected["species"]]
    packet_dir = packet_dir_root / packet_id
    packet_dir.mkdir(parents=True, exist_ok=False)
    pd.DataFrame({"species": species}).to_csv(packet_dir / "species.csv", index=False)
    candidates = _candidate_subset(candidate_csv, species)
    seed_evidence = (
        evidence_from_bulk_candidates(candidates)
        if candidates is not None
        else pd.DataFrame(columns=EVIDENCE_COLUMNS)
    )

    search = collector or collect_free_sources
    try:
        sources, errors = search(
            species,
            include_gbif=include_gbif,
            include_wikimedia=include_wikimedia,
            include_openalex=include_openalex,
            include_web_descriptions=include_web_descriptions,
            include_world_flora=include_world_flora,
            pause_seconds=pause_seconds,
            seed_evidence=seed_evidence,
        )
        sources = pd.DataFrame(sources, columns=SOURCE_COLUMNS).fillna("")
        errors = pd.DataFrame(errors, columns=["species", "source", "error"]).fillna("")
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
        _write_failure_checkpoint(
            checkpoint,
            selected_indexes,
            max_attempts,
            error,
            checkpoint_path,
        )
        _atomic_json(
            {
                "contract_version": CONTRACT_VERSION,
                "packet_id": packet_id,
                "error": error,
                "species": species,
                "retry_policy": "retry until max_attempts, then retain as exhausted",
            },
            packet_dir / "fatal_error.json",
        )
        raise

    errors = errors.copy().fillna("")
    retryable_errors = _retryable_errors_by_species(errors, species)
    for index in selected_indexes:
        name = _text(checkpoint.loc[index, "species"])
        row = extracted.loc[extracted["species"].eq(name)].iloc[0]
        checkpoint.loc[index, "result_sha256"] = _result_sha(row)
        checkpoint.loc[index, "updated_at"] = _now()
        if name in retryable_errors:
            checkpoint.loc[index, "last_error"] = " | ".join(
                dict.fromkeys(value for value in retryable_errors[name] if value)
            )[:2000]
            checkpoint.loc[index, "status"] = (
                "exhausted" if int(checkpoint.loc[index, "attempts"]) >= max_attempts else "retry"
            )
        else:
            checkpoint.loc[index, "status"] = "completed"
            checkpoint.loc[index, "last_error"] = ""

    cumulative_dir.mkdir(parents=True, exist_ok=True)
    cumulative_results = _merge_cumulative(
        prior_results_path,
        extracted,
        dedupe_subset=["species"],
    )
    file_specs = (
        ("source_texts.csv", "search_evidence.csv", SOURCE_COLUMNS, None),
        ("v1_category_evidence.csv", "trait_evidence.csv", EVIDENCE_COLUMNS, None),
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
    _atomic_csv(checkpoint, checkpoint_path)
    manifest = _packet_manifest(packet_dir, packet_id, species)
    _atomic_json(manifest, packet_dir / "packet_manifest.json")
    report = _status_report(
        checkpoint,
        cumulative_results,
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
    web_descriptions: bool = typer.Option(
        True,
        "--web-descriptions/--no-web-descriptions",
    ),
    world_flora: bool = typer.Option(True, "--world-flora/--no-world-flora"),
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
        include_web_descriptions=web_descriptions,
        include_world_flora=world_flora,
        pause_seconds=pause_seconds,
    )
    typer.echo(json.dumps(report, ensure_ascii=False))


if __name__ == "__main__":
    app()
