"""Build a durable species-by-axis checkpoint from public-web shard outputs.

The checkpoint is deliberately conservative: an axis is complete only when a
non-empty direct candidate exists for one of its source traits. Provider state
is retained separately so completed provider/species attempts are never
repeated. Missing values remain unresolved; no global fallback is introduced.
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
    """Yield species, axis, trait, source from heterogeneous cumulative tables."""
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
                        source = _text(row.get("source")) or _text(row.get("source_type")) or path.stem
                        yield species, axis, trait, source
            continue

        # Older public-web tables are wide. Treat a populated source-backed field
        # as direct only when the row is not explicitly marked as inference.
        for row in frame.to_dict("records"):
            species = _text(row.get(species_column))
            if not species:
                continue
            evidence_type = _text(row.get("evidence_type")).casefold()
            if evidence_type == "inference":
                continue
            for axis, traits in AXIS_TRAITS.items():
                for trait in traits:
                    if trait not in frame.columns:
                        continue
                    value = _text(row.get(trait))
                    if value.casefold() in EMPTY_VALUES:
                        continue
                    source = _text(row.get("source_kind")) or _text(row.get("source")) or path.stem
                    yield species, axis, trait, source


def build_checkpoint(campaign_dir: Path, output_dir: Path) -> dict[str, object]:
    species_checkpoint = _read_csv(campaign_dir / "checkpoint.csv")
    provider_checkpoint = _read_csv(campaign_dir / "provider_checkpoint.csv")
    if species_checkpoint.empty or "species" not in species_checkpoint.columns:
        raise typer.BadParameter(f"missing usable species checkpoint: {campaign_dir / 'checkpoint.csv'}")

    evidence_paths = sorted((campaign_dir / "cumulative").rglob("*.csv"))
    evidence_paths += sorted((campaign_dir / "cumulative").rglob("*.csv.gz"))
    evidence: dict[tuple[str, str], dict[str, set[str]]] = {}
    for species, axis, trait, source in _iter_direct_evidence(evidence_paths):
        item = evidence.setdefault((species, axis), {"traits": set(), "sources": set()})
        item["traits"].add(trait)
        item["sources"].add(source)

    provider_terminal: dict[str, int] = {}
    provider_remaining: dict[str, int] = {}
    if not provider_checkpoint.empty and {"species", "provider", "status"}.issubset(provider_checkpoint):
        terminal = {"hit", "no_hit", "skipped_covered", "legacy_completed", "exhausted"}
        for provider, group in provider_checkpoint.groupby("provider"):
            statuses = group["status"].map(_text)
            provider_terminal[_text(provider)] = int(statuses.isin(terminal).sum())
            provider_remaining[_text(provider)] = int((~statuses.isin(terminal)).sum())

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
        row["continue_direct_acquisition"] = bool(missing) and _text(record.get("status")) != "exhausted"
        row["collection_priority"] = (
            0 if "reproductive_assurance" in missing
            else 1 if "floral_complexity" in missing
            else 2 if "colour_vividness" in missing
            else 3
        )
        rows.append(row)

    status = pd.DataFrame(rows).sort_values(
        ["continue_direct_acquisition", "collection_priority", "accepted_species"],
        ascending=[False, True, True],
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    status.to_csv(output_dir / "three_axis_checkpoint.csv.gz", index=False)
    status.loc[status["continue_direct_acquisition"]].to_csv(
        output_dir / "three_axis_next_queue.csv.gz", index=False
    )

    summary = {
        "contract_version": "three_axis_acquisition_checkpoint_v1",
        "n_species": int(len(status)),
        "n_all_three_axes_direct": int(status["all_three_axes_direct"].sum()),
        "n_continue_direct_acquisition": int(status["continue_direct_acquisition"].sum()),
        "axis_direct_counts": {
            axis: int(status[f"{axis}_direct"].sum()) for axis in AXIS_TRAITS
        },
        "provider_terminal_counts": provider_terminal,
        "provider_remaining_counts": provider_remaining,
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
