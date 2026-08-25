"""Build a durable full-coverage species-by-axis acquisition checkpoint.

An axis is complete only when a non-empty species- or synonym-direct candidate
exists for one of its source traits. Provider state is retained independently,
so every enabled species/provider pair is attempted once and terminal work is
never repeated. Missing values remain unresolved; no global fallback is used.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXIS_TRAITS = {
    "colour_vividness": {"flower_primary_color", "flower_color"},
    "floral_complexity": {
        "floral_form",
        "floral_symmetry",
        "tube_depth_class",
        "flower_size_class",
        "inflorescence_display",
        "flower_shape",
    },
    "reproductive_assurance": {
        "self_incompatibility",
        "autonomous_selfing_capacity",
        "mating_system",
        "cleistogamy",
    },
}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}
EMPTY_VALUES = {"", "unknown", "unresolved", "na", "nan", "none"}
PROVIDER_TERMINAL_STATUSES = {"hit", "no_hit", "skipped_covered", "legacy_completed", "exhausted"}
PROVIDER_SUCCESS_STATUSES = {"hit", "no_hit", "skipped_covered", "legacy_completed"}
PROVIDER_ERROR_STATUSES = {"retry", "exhausted"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _read_csv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, dtype=str).fillna("")
    except (OSError, ValueError, pd.errors.ParserError):
        return pd.DataFrame()


def _candidate_value_column(frame: pd.DataFrame) -> str | None:
    for column in (
        "provisional_candidate_value",
        "candidate_value",
        "standardized_value",
        "filled_value",
        "reported_value",
        "raw_candidate_value",
        "value",
    ):
        if column in frame.columns:
            return column
    return None


def _species_column(frame: pd.DataFrame) -> str | None:
    for column in ("accepted_species", "species"):
        if column in frame.columns:
            return column
    return None


def _iter_direct_evidence(paths: Iterable[Path]) -> Iterable[tuple[str, str, str, str]]:
    """Yield accepted species, axis, trait, and provider/source for direct records."""
    for path in paths:
        frame = _read_csv(path)
        if frame.empty:
            continue
        species_column = _species_column(frame)
        if species_column is None:
            continue

        if "trait_name" in frame.columns:
            value_column = _candidate_value_column(frame)
            if value_column is None:
                continue
            for row in frame.to_dict("records"):
                species = _text(row.get(species_column))
                trait = _text(row.get("trait_name"))
                value = _text(row.get(value_column))
                scope = _text(row.get("evidence_scope"))
                if not species or value.casefold() in EMPTY_VALUES:
                    continue
                if scope and scope not in DIRECT_SCOPES:
                    continue
                for axis, traits in AXIS_TRAITS.items():
                    if trait in traits:
                        source = (
                            _text(row.get("provider"))
                            or _text(row.get("source"))
                            or _text(row.get("source_type"))
                            or path.stem
                        )
                        yield species, axis, trait, source
            continue

        # Legacy public-web tables are wide. A populated source-backed field is
        # direct only when the row is not explicitly marked as inference.
        for row in frame.to_dict("records"):
            species = _text(row.get(species_column))
            if not species or _text(row.get("evidence_type")).casefold() == "inference":
                continue
            for axis, traits in AXIS_TRAITS.items():
                for trait in traits:
                    if trait not in frame.columns:
                        continue
                    value = _text(row.get(trait))
                    if value.casefold() in EMPTY_VALUES:
                        continue
                    source = (
                        _text(row.get("provider"))
                        or _text(row.get("source_kind"))
                        or _text(row.get("source"))
                        or path.stem
                    )
                    yield species, axis, trait, source


def build_checkpoint(campaign_dir: Path, output_dir: Path) -> dict[str, object]:
    species_checkpoint = _read_csv(campaign_dir / "checkpoint.csv")
    provider_checkpoint = _read_csv(campaign_dir / "provider_checkpoint.csv")
    if species_checkpoint.empty or "species" not in species_checkpoint.columns:
        raise typer.BadParameter(f"missing usable species checkpoint: {campaign_dir / 'checkpoint.csv'}")

    evidence_paths = sorted((campaign_dir / "cumulative").rglob("*.csv"))
    evidence_paths += sorted((campaign_dir / "cumulative").rglob("*.csv.gz"))
    evidence: dict[tuple[str, str], dict[str, set[str]]] = {}
    direct_record_keys: set[tuple[str, str, str, str]] = set()
    for species, axis, trait, source in _iter_direct_evidence(evidence_paths):
        item = evidence.setdefault((species, axis), {"traits": set(), "sources": set()})
        item["traits"].add(trait)
        item["sources"].add(source)
        direct_record_keys.add((species, axis, trait, source))

    provider_terminal: dict[str, int] = {}
    provider_remaining: dict[str, int] = {}
    provider_errors: dict[str, int] = {}
    successful_species: set[str] = set()
    attempted_species: set[str] = set()
    enabled_provider_tasks = 0
    terminal_provider_tasks = 0
    if not provider_checkpoint.empty and {"species", "provider", "status"}.issubset(provider_checkpoint):
        enabled = (
            provider_checkpoint["enabled"].map(_text).str.casefold().isin({"true", "1", "yes"})
            if "enabled" in provider_checkpoint.columns
            else pd.Series(True, index=provider_checkpoint.index)
        )
        enabled_table = provider_checkpoint.loc[enabled].copy()
        enabled_provider_tasks = int(len(enabled_table))
        for record in enabled_table.to_dict("records"):
            species = _text(record.get("species"))
            status = _text(record.get("status"))
            if status not in {"", "pending"}:
                attempted_species.add(species)
            if status in PROVIDER_SUCCESS_STATUSES:
                successful_species.add(species)
        for provider, group in enabled_table.groupby("provider"):
            statuses = group["status"].map(_text)
            provider_name = _text(provider)
            terminal_count = int(statuses.isin(PROVIDER_TERMINAL_STATUSES).sum())
            provider_terminal[provider_name] = terminal_count
            provider_remaining[provider_name] = int((~statuses.isin(PROVIDER_TERMINAL_STATUSES)).sum())
            provider_errors[provider_name] = int(statuses.isin(PROVIDER_ERROR_STATUSES).sum())
            terminal_provider_tasks += terminal_count

    rows: list[dict[str, object]] = []
    for record in species_checkpoint.to_dict("records"):
        species = _text(record.get("species"))
        row: dict[str, object] = {
            "accepted_species": species,
            "genus": _text(record.get("genus")),
            "family": _text(record.get("family")),
            "shard_index": _text(record.get("shard_index")),
            "species_campaign_status": _text(record.get("status")),
            "species_attempts": _text(record.get("attempts")),
        }
        missing: list[str] = []
        for axis in AXIS_TRAITS:
            item = evidence.get((species, axis), {"traits": set(), "sources": set()})
            complete = bool(item["traits"])
            row[f"{axis}_direct"] = complete
            row[f"{axis}_traits"] = "|".join(sorted(item["traits"]))
            row[f"{axis}_sources"] = "|".join(sorted(item["sources"]))
            if not complete:
                missing.append(axis)
        row["all_three_axes_direct"] = not missing
        row["missing_direct_axes"] = "|".join(missing)
        row["all_enabled_providers_terminal"] = _text(record.get("status")) in {"completed", "exhausted"}
        row["continue_direct_acquisition"] = not row["all_enabled_providers_terminal"]
        # No species subset or trait axis receives campaign priority.
        row["collection_priority"] = 0
        rows.append(row)

    status = pd.DataFrame(rows).sort_values("accepted_species").reset_index(drop=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    status.to_csv(output_dir / "three_axis_checkpoint.csv.gz", index=False)
    status.loc[status["continue_direct_acquisition"]].to_csv(
        output_dir / "three_axis_next_queue.csv.gz", index=False
    )

    axis_counts = {axis: int(status[f"{axis}_direct"].sum()) for axis in AXIS_TRAITS}
    unresolved_species = int((~status["all_three_axes_direct"]).sum())
    remaining_tasks = int(sum(provider_remaining.values()))
    summary = {
        "contract_version": "three_axis_full_coverage_checkpoint_v2",
        "total_species_processed": int(len(status)),
        "species_attempted": int(len(attempted_species)),
        "species_successfully_queried": int(len(successful_species)),
        "direct_evidence_records": int(len(direct_record_keys)),
        "direct_coverage_counts": axis_counts,
        "direct_coverage_fraction": {
            axis: (count / len(status) if len(status) else 0.0) for axis, count in axis_counts.items()
        },
        "species_all_three_axes_direct": int(status["all_three_axes_direct"].sum()),
        "unresolved_species_count": unresolved_species,
        "enabled_provider_tasks": enabled_provider_tasks,
        "terminal_provider_tasks": terminal_provider_tasks,
        "provider_terminal_counts": provider_terminal,
        "provider_remaining_counts": provider_remaining,
        "provider_error_or_blocked_counts": provider_errors,
        "next_continuation_state": {
            "continue": remaining_tasks > 0,
            "remaining_provider_tasks": remaining_tasks,
            "remaining_species": int(status["continue_direct_acquisition"].sum()),
        },
        "global_fallback_used": False,
    }
    (output_dir / "three_axis_checkpoint_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command()
def build(
    campaign_dir: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(build_checkpoint(campaign_dir, output_dir), indent=2))


if __name__ == "__main__":
    app()
