"""Aggregate three-axis evidence and resumable provider checkpoints.

Evidence quality is deliberately simple:

- high: literature, flora, and structured biodiversity/trait databases
- medium: species- or synonym-level web descriptions, including Wikipedia and
  horticultural pages, when a supporting excerpt and URL are retained
- low: explicit genus inference only

Family inference and global fallback are prohibited. Missing evidence remains
unresolved.
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
DIRECT_SCOPES = {"species_direct", "synonym_direct"}
PROHIBITED_SCOPES = {"family_inference", "global_fallback"}
TERMINAL_PROVIDER_STATUSES = {"hit", "no_hit", "skipped_covered", "legacy_completed", "exhausted"}
EMPTY = {"", "unknown", "unresolved", "na", "nan", "none"}
HIGH_PROVIDER_MARKERS = {
    "crossref",
    "doi",
    "europe_pmc",
    "pubmed",
    "openalex",
    "semantic_scholar",
    "journal",
    "literature",
    "flora",
    "efloras",
    "gbif",
    "traitbank",
    "austraits",
    "gift",
    "usda",
    "database",
}
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
    "evidence_scope",
    "evidence_quality",
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


def _quality(provider: str, source_title: str, scope: str) -> str:
    if scope == "genus_inference":
        return "low"
    haystack = f"{provider} {source_title}".casefold()
    if any(marker in haystack for marker in HIGH_PROVIDER_MARKERS):
        return "high"
    return "medium"


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
        "source_kind": "provider",
        "source_type": "provider",
        "fill_tier": "evidence_scope",
    }
    frame = frame.rename(columns={k: v for k, v in aliases.items() if k in frame.columns and v not in frame.columns})
    for column in REQUIRED_EVIDENCE_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    for column in frame.columns:
        frame[column] = frame[column].map(_text)

    frame["evidence_scope"] = frame["evidence_scope"].replace(
        {"": "species_direct", "direct": "species_direct", "inference": "genus_inference"}
    )
    prohibited = frame["evidence_scope"].isin(PROHIBITED_SCOPES)
    frame = frame.loc[~prohibited].copy()
    frame["evidence_quality"] = [
        _quality(provider, title, scope)
        for provider, title, scope in zip(
            frame["provider"], frame["source_title"], frame["evidence_scope"], strict=False
        )
    ]

    value = frame["reported_value"].where(frame["reported_value"].ne(""), frame["raw_candidate_value"])
    valid_direct = frame["evidence_scope"].isin(DIRECT_SCOPES) & frame["source_url"].ne("") & frame["excerpt"].ne("")
    valid_genus = frame["evidence_scope"].eq("genus_inference")
    frame = frame.loc[
        frame["accepted_species"].ne("")
        & frame["trait_name"].ne("")
        & ~value.str.casefold().isin(EMPTY)
        & (valid_direct | valid_genus)
    ].copy()
    return frame[list(sorted(REQUIRED_EVIDENCE_COLUMNS))].drop_duplicates(
        subset=[
            "accepted_species",
            "matched_source_name",
            "trait_name",
            "reported_value",
            "raw_candidate_value",
            "source_url",
            "excerpt",
            "evidence_scope",
        ]
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
        unresolved = []
        for axis, traits in AXES.items():
            axis_rows = subset.loc[subset["trait_name"].isin(traits)]
            direct = axis_rows["evidence_scope"].isin(DIRECT_SCOPES).any()
            genus = axis_rows["evidence_scope"].eq("genus_inference").any()
            row[f"{axis}_direct"] = bool(direct)
            row[f"{axis}_genus_inferred"] = bool(genus and not direct)
            row[f"{axis}_resolved"] = bool(direct or genus)
            if not direct and not genus:
                unresolved.append(axis)
        row["all_three_axes_direct"] = all(row[f"{axis}_direct"] for axis in AXES)
        row["all_three_axes_resolved"] = all(row[f"{axis}_resolved"] for axis in AXES)
        row["unresolved_axes"] = "|".join(unresolved)
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
        "contract_version": "full_coverage_three_axis_v2",
        "evidence_hierarchy": {"high": "literature_or_database", "medium": "species_web_description", "low": "genus_inference"},
        "family_inference_used": False,
        "total_species_processed": len(species),
        "species_successfully_queried": len(successful_species),
        "evidence_records": int(len(evidence)),
        "evidence_records_by_quality": evidence["evidence_quality"].value_counts().sort_index().astype(int).to_dict(),
        "direct_evidence_records_added": int(evidence["evidence_scope"].isin(DIRECT_SCOPES).sum()),
        "direct_coverage": {
            axis: {
                "species": int(coverage[f"{axis}_direct"].sum()),
                "fraction": float(coverage[f"{axis}_direct"].mean()) if len(coverage) else 0.0,
            }
            for axis in AXES
        },
        "resolved_coverage_including_genus_low": {
            axis: int(coverage[f"{axis}_resolved"].sum()) for axis in AXES
        },
        "species_with_all_three_axes_direct": int(coverage["all_three_axes_direct"].sum()),
        "species_with_all_three_axes_resolved": int(coverage["all_three_axes_resolved"].sum()),
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
    evidence.to_csv(output_dir / "all_master_evidence.csv.gz", index=False)
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
