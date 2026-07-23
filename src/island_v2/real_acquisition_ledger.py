"""Aggregate public-web acquisition outputs into a three-axis evidence ledger.

Evidence hierarchy:
- high: literature, flora, and structured biodiversity/trait databases
- medium: species/synonym web descriptions, including Wikipedia and horticultural sites
- low: explicit genus inference

Family inference and global fallback are rejected.
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
EMPTY = {"", "unknown", "unresolved", "na", "nan", "none"}
TERMINAL = {"hit", "no_hit", "skipped_covered", "legacy_completed", "exhausted"}
HIGH_MARKERS = {
    "crossref",
    "doi",
    "europe_pmc",
    "pubmed",
    "openalex",
    "semantic_scholar",
    "journal",
    "literature",
    "flora",
    "monograph",
    "efloras",
    "gbif",
    "traitbank",
    "austraits",
    "gift",
    "usda",
    "database",
}
TRAIT_ALIASES = {
    "flower_color": "flower_primary_color",
    "flower_shape": "floral_form",
}
LEDGER_COLUMNS = [
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
]
COVERAGE_COLUMNS = [
    "accepted_species",
    *[
        column
        for axis in AXES
        for column in (
            f"{axis}_direct",
            f"{axis}_genus_inferred",
            f"{axis}_resolved",
        )
    ],
    "all_three_axes_direct",
    "all_three_axes_resolved",
    "unresolved_axes",
]


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _read(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, dtype=str).fillna("")
    except (OSError, ValueError, pd.errors.ParserError):
        return pd.DataFrame()


def _read_many(paths: list[Path]) -> pd.DataFrame:
    frames = [_read(path) for path in paths]
    frames = [frame for frame in frames if not frame.empty]
    return pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame()


def _quality(provider: str, title: str, scope: str) -> str:
    if scope == "genus_inference":
        return "low"
    haystack = f"{provider} {title}".casefold()
    return "high" if any(marker in haystack for marker in HIGH_MARKERS) else "medium"


def _normalise(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(columns=LEDGER_COLUMNS)

    aliases = {
        "species": "accepted_species",
        "matched_name": "matched_source_name",
        "source_name": "matched_source_name",
        "trait": "trait_name",
        "field": "trait_name",
        "candidate_value": "raw_candidate_value",
        "provisional_candidate_value": "raw_candidate_value",
        "value": "reported_value",
        "raw_description": "excerpt",
        "source_excerpt": "excerpt",
        "source_text": "excerpt",
        "source_citation": "source_title",
        "source_type": "provider",
        "source_kind": "provider",
        "doi": "stable_source_id",
        "source_id": "stable_source_id",
        "fill_tier": "evidence_scope",
        "evidence_type": "evidence_scope",
    }
    frame = frame.rename(
        columns={
            old: new
            for old, new in aliases.items()
            if old in frame.columns and new not in frame.columns
        }
    ).copy()
    for column in LEDGER_COLUMNS:
        if column not in frame.columns:
            frame[column] = ""
    for column in frame.columns:
        frame[column] = frame[column].map(_text)

    frame["trait_name"] = frame["trait_name"].replace(TRAIT_ALIASES)
    frame["matched_source_name"] = frame["matched_source_name"].where(
        frame["matched_source_name"].ne(""), frame["accepted_species"]
    )
    frame["evidence_scope"] = frame["evidence_scope"].replace(
        {
            "": "species_direct",
            "direct": "species_direct",
            "flora": "species_direct",
            "reported": "species_direct",
            "inference": "genus_inference",
        }
    )
    frame = frame.loc[~frame["evidence_scope"].isin(PROHIBITED_SCOPES)].copy()
    frame["evidence_quality"] = [
        _quality(provider, title, scope)
        for provider, title, scope in zip(
            frame["provider"],
            frame["source_title"],
            frame["evidence_scope"],
            strict=False,
        )
    ]

    value = frame["reported_value"].where(
        frame["reported_value"].ne(""), frame["raw_candidate_value"]
    )
    valid_direct = (
        frame["evidence_scope"].isin(DIRECT_SCOPES)
        & frame["source_url"].ne("")
        & frame["excerpt"].ne("")
    )
    valid_low = frame["evidence_scope"].eq("genus_inference")
    frame = frame.loc[
        frame["accepted_species"].ne("")
        & frame["trait_name"].ne("")
        & ~value.str.casefold().isin(EMPTY)
        & (valid_direct | valid_low)
    ].copy()
    return frame[LEDGER_COLUMNS].drop_duplicates(
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


def _empty_shard_status(campaign_root: Path) -> bool:
    status_files = sorted(campaign_root.rglob("campaign_status.json"))
    if not status_files:
        return False
    try:
        status = json.loads(status_files[0].read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return False
    return bool(status.get("complete")) and int(status.get("n_species_in_shard", -1)) == 0


def aggregate(campaign_root: Path, output_dir: Path) -> dict[str, object]:
    checkpoint = _read_many(sorted(campaign_root.rglob("checkpoint.csv")))
    provider = _read_many(sorted(campaign_root.rglob("provider_checkpoint.csv")))
    empty_shard = checkpoint.empty and _empty_shard_status(campaign_root)

    if checkpoint.empty and not empty_shard:
        raise typer.BadParameter("no usable checkpoint.csv found")
    if not checkpoint.empty and "species" not in checkpoint.columns:
        raise typer.BadParameter("checkpoint.csv is missing species column")

    preferred = sorted(campaign_root.rglob("v2_reported_candidates.csv"))
    fallback = sorted(campaign_root.rglob("trait_evidence.csv"))
    evidence = _normalise(_read_many(preferred if preferred else fallback))

    if empty_shard:
        all_species: list[str] = []
        processed_species: list[str] = []
    else:
        all_species = sorted({_text(value) for value in checkpoint["species"] if _text(value)})
        attempted_mask = pd.Series(False, index=checkpoint.index)
        if "attempts" in checkpoint.columns:
            attempted_mask |= pd.to_numeric(
                checkpoint["attempts"], errors="coerce"
            ).fillna(0).gt(0)
        if "status" in checkpoint.columns:
            attempted_mask |= ~checkpoint["status"].map(_text).isin({"", "pending"})
        processed_species = sorted(
            {
                _text(value)
                for value in checkpoint.loc[attempted_mask, "species"]
                if _text(value)
            }
        )

    processed_set = set(processed_species)
    evidence = evidence.loc[evidence["accepted_species"].isin(processed_set)].copy()

    rows = []
    for species in processed_species:
        subset = evidence.loc[evidence["accepted_species"].eq(species)]
        row: dict[str, object] = {"accepted_species": species}
        missing = []
        for axis, traits in AXES.items():
            axis_rows = subset.loc[subset["trait_name"].isin(traits)]
            direct = bool(axis_rows["evidence_scope"].isin(DIRECT_SCOPES).any())
            genus = bool(axis_rows["evidence_scope"].eq("genus_inference").any())
            row[f"{axis}_direct"] = direct
            row[f"{axis}_genus_inferred"] = genus and not direct
            row[f"{axis}_resolved"] = direct or genus
            if not direct and not genus:
                missing.append(axis)
        row["all_three_axes_direct"] = all(row[f"{axis}_direct"] for axis in AXES)
        row["all_three_axes_resolved"] = all(row[f"{axis}_resolved"] for axis in AXES)
        row["unresolved_axes"] = "|".join(missing)
        rows.append(row)
    coverage = pd.DataFrame(rows, columns=COVERAGE_COLUMNS)

    provider_state: dict[str, dict[str, int]] = {}
    if not provider.empty and {"provider", "status", "species"}.issubset(provider.columns):
        active = provider.loc[provider["species"].map(_text).isin(processed_set)].copy()
        for name, group in active.groupby("provider"):
            statuses = group["status"].map(_text)
            provider_state[_text(name)] = {
                "hit": int(statuses.eq("hit").sum()),
                "no_hit": int(statuses.eq("no_hit").sum()),
                "terminal": int(statuses.isin(TERMINAL).sum()),
                "errors_or_blocked": int(statuses.isin({"retry", "exhausted"}).sum()),
                "remaining": int((~statuses.isin(TERMINAL)).sum()),
            }

    quality_counts = (
        evidence["evidence_quality"]
        .value_counts()
        .reindex(["high", "medium", "low"], fill_value=0)
        .astype(int)
        .to_dict()
    )
    direct_coverage = {
        axis: int(coverage[f"{axis}_direct"].sum()) if not coverage.empty else 0
        for axis in AXES
    }
    summary = {
        "contract_version": "real_acquisition_ledger_v1",
        "evidence_hierarchy": {
            "high": "literature_or_database",
            "medium": "species_or_synonym_web_description",
            "low": "genus_inference",
        },
        "family_inference_used": False,
        "global_fallback_used": False,
        "empty_shard": empty_shard,
        "total_species_in_shard": len(all_species),
        "total_species_processed": len(processed_species),
        "species_with_any_evidence": int(evidence["accepted_species"].nunique()),
        "evidence_records": int(len(evidence)),
        "evidence_records_by_quality": quality_counts,
        "direct_coverage": direct_coverage,
        "species_with_all_three_axes_direct": int(coverage["all_three_axes_direct"].sum())
        if not coverage.empty
        else 0,
        "species_with_all_three_axes_resolved": int(
            coverage["all_three_axes_resolved"].sum()
        )
        if not coverage.empty
        else 0,
        "unresolved_species_count": int(coverage["unresolved_axes"].ne("").sum())
        if not coverage.empty
        else 0,
        "provider_state": provider_state,
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
