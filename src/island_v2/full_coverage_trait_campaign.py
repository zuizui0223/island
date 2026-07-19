"""Aggregate direct three-axis evidence and resumable provider checkpoints.

This module is intentionally small. Collection remains in the existing public-web
shard runner; this module owns only the durable long-form ledger, coverage report,
and continuation state. Missing evidence stays unresolved.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

AXES = {
    "flower_colour": {"flower_primary_color", "flower_color"},
    "floral_structural_complexity": {
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
TERMINAL_PROVIDER_STATUSES = {"hit", "no_hit", "skipped_covered", "legacy_completed", "exhausted"}
EMPTY = {"", "unknown", "unresolved", "na", "nan", "none"}
REQUIRED_EVIDENCE_COLUMNS = {
    "accepted_species",
    "matched_source_name",
    "trait_name",
    "reported_value",
    "raw_candidate_value",
    "source_url",
    "provider",
    "source_title",
    "excerpt",
    "retrieval_date",
    "stable_source_id",
}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _read_many(paths: list[Path]) -> pd.DataFrame:
    frames = []
    for path in paths:
        try:
            frames.append(pd.read_csv(path, dtype=str).fillna(""))
        except (OSError, ValueError, pd.errors.ParserError):
            continue
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _normalise_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=sorted(REQUIRED_EVIDENCE_COLUMNS))
    aliases = {
        "species": "accepted_species",
        "matched_name": "matched_source_name",
        "source_name": "matched_source_name",
        "trait": "trait_name",
        "candidate_value": "raw_candidate_value",
        "provisional_candidate_value": "raw_candidate_value",
        "value": "reported_value",
        "doi": "stable_source_id",
        "source_id": "stable_source_id",
        "supporting_text": "excerpt",
    }
    frame = frame.rename(columns={k: v for k, v in aliases.items() if k in frame.columns and v not in frame.columns})
    for column in REQUIRED_EVIDENCE_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    frame = frame[list(sorted(REQUIRED_EVIDENCE_COLUMNS))].copy()
    for column in frame.columns:
        frame[column] = frame[column].map(_text)
    value = frame["reported_value"].where(frame["reported_value"].ne(""), frame["raw_candidate_value"])
    frame = frame.loc[
        frame["accepted_species"].ne("")
        & frame["trait_name"].ne("")
        & ~value.str.casefold().isin(EMPTY)
        & frame["source_url"].ne("")
        & frame["excerpt"].ne("")
    ].copy()
    return frame.drop_duplicates(
        subset=["accepted_species", "matched_source_name", "trait_name", "reported_value", "raw_candidate_value", "source_url", "excerpt"]
    )


def aggregate(campaign_root: Path, output_dir: Path) -> dict[str, object]:
    checkpoint = _read_many(sorted(campaign_root.rglob("checkpoint.csv")))
    provider = _read_many(sorted(campaign_root.rglob("provider_checkpoint.csv")))
    evidence_paths = sorted(campaign_root.rglob("*.csv")) + sorted(campaign_root.rglob("*.csv.gz"))
    evidence_paths = [p for p in evidence_paths if p.name not in {"checkpoint.csv", "provider_checkpoint.csv"}]
    evidence = _normalise_evidence(_read_many(evidence_paths))

    if checkpoint.empty or "species" not in checkpoint.columns:
        raise typer.BadParameter("no usable shard checkpoints found")
    species = sorted({_text(v) for v in checkpoint["species"] if _text(v)})
    species_set = set(species)
    evidence = evidence.loc[evidence["accepted_species"].isin(species_set)].copy()

    coverage_rows = []
    for accepted_species in species:
        subset = evidence.loc[evidence["accepted_species"].eq(accepted_species)]
        row: dict[str, object] = {"accepted_species": accepted_species}
        missing = []
        for axis, traits in AXES.items():
            covered = bool(subset["trait_name"].isin(traits).any())
            row[f"{axis}_direct"] = covered
            if not covered:
                missing.append(axis)
        row["all_three_axes_direct"] = not missing
        row["unresolved_axes"] = "|".join(missing)
        coverage_rows.append(row)
    coverage = pd.DataFrame(coverage_rows)

    provider_summary: dict[str, dict[str, int]] = {}
    if not provider.empty and {"provider", "status"}.issubset(provider.columns):
        for name, group in provider.groupby("provider"):
            statuses = group["status"].map(_text)
            provider_summary[_text(name)] = {
                "terminal": int(statuses.isin(TERMINAL_PROVIDER_STATUSES).sum()),
                "remaining": int((~statuses.isin(TERMINAL_PROVIDER_STATUSES)).sum()),
                "errors_or_blocked": int(statuses.isin({"retry", "exhausted"}).sum()),
            }

    successful_species = set()
    if not provider.empty and {"species", "status"}.issubset(provider.columns):
        successful_species = set(provider.loc[provider["status"].isin({"hit", "no_hit"}), "species"].map(_text))
    remaining_tasks = sum(item["remaining"] for item in provider_summary.values())
    summary = {
        "contract_version": "full_coverage_three_axis_v1",
        "total_species_processed": len(species),
        "species_successfully_queried": len(successful_species),
        "direct_evidence_records_added": int(len(evidence)),
        "direct_coverage": {
            axis: {
                "species": int(coverage[f"{axis}_direct"].sum()),
                "fraction": float(coverage[f"{axis}_direct"].mean()) if len(coverage) else 0.0,
            }
            for axis in AXES
        },
        "species_with_all_three_axes_direct": int(coverage["all_three_axes_direct"].sum()),
        "unresolved_species_count": int(coverage["unresolved_axes"].ne("").sum()),
        "provider_state": provider_summary,
        "provider_errors_and_blocked_sources": {
            name: item["errors_or_blocked"] for name, item in provider_summary.items()
        },
        "next_continuation_state": {
            "remaining_provider_tasks": remaining_tasks,
            "complete": remaining_tasks == 0,
            "resume_from_existing_checkpoints": remaining_tasks > 0,
        },
        "global_fallback_used": False,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output_dir / "all_master_direct_evidence.csv.gz", index=False)
    coverage.to_csv(output_dir / "three_axis_coverage.csv.gz", index=False)
    (output_dir / "collection_round_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


@app.command()
def build(
    campaign_root: Path = typer.Option(..., exists=True, file_okay=False),
    output_dir: Path = typer.Option(...),
) -> None:
    typer.echo(json.dumps(aggregate(campaign_root, output_dir), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    app()
