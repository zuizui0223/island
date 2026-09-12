"""Deterministic failure-class triage for Database 2.0 alpha2 taxonomy review rows.

This layer does not resolve names, alter alpha1, or promote any taxon. It partitions
the frozen alpha2 review queue into source-scope problems, kingdom conflicts, and
species-binomial names that are eligible for an independent second-backbone audit.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

from island_v2.taxon_intake import classify_name_rank

app = typer.Typer(add_completion=False, no_args_is_help=True)

ALPHA2_REVIEW_SHA256 = "66669f60d96e9356e3d7e6ad6b803de9a15b60be608443e7694b0eee600828bf"
ALPHA2_REVIEW_ROWS = 3_767
ALPHA2_CANONICAL_RUN = 34_668_274_725
ALPHA2_CANONICAL_ARTIFACT = 10_289_759_727


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            frame.to_csv(zipped, index=False)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def triage_class(row: pd.Series) -> str:
    kingdom = _text(row.get("candidate_kingdom", ""))
    if kingdom and kingdom.casefold() != "plantae":
        return "kingdom_conflict_non_plantae"

    rank_class = _text(row["source_name_rank_class"])
    if rank_class == "genus_or_higher_rank":
        return "source_genus_or_higher_rank"
    if rank_class == "infra_or_hybrid_candidate":
        return "source_infra_or_hybrid"
    if rank_class in {"unresolved_name_rank", "indeterminate_or_placeholder", "blank"}:
        return "source_name_not_species_binomial"

    match_type = _text(row.get("match_type", "")).upper()
    matched_status = _text(row.get("matched_status", "")).upper()
    accepted_rank = _text(row.get("candidate_accepted_rank", "")).upper()
    flags = _text(row.get("processing_flags", ""))

    if match_type == "NONE":
        if "MULTIPLE_MATCHES_SAME_CONFIDENCE" in flags:
            return "colxr_none_ambiguous"
        if "LOW_CONFIDENCE" in flags:
            return "colxr_none_low_confidence"
        return "colxr_none_no_match"
    if match_type == "VARIANT":
        return "colxr_variant_species_name"
    if match_type == "HIGHERRANK":
        return "colxr_higher_rank_only"
    if match_type == "EXACT":
        if matched_status == "PROVISIONALLY_ACCEPTED":
            return "colxr_exact_provisionally_accepted"
        if accepted_rank in {"SUBSPECIES", "VARIETY", "FORM"}:
            return "colxr_exact_target_infraspecific"
        if accepted_rank and accepted_rank != "SPECIES":
            return "colxr_exact_target_non_species"
        if matched_status in {"AMBIGUOUS_SYNONYM", "MISAPPLIED"}:
            return "colxr_exact_ambiguous_or_misapplied"
        return "colxr_exact_other_review"
    return "other_review"


def _next_action(triage: str) -> str:
    if triage == "kingdom_conflict_non_plantae":
        return "audit_source_occurrence_taxon_and_remove_from_species_plant_layer_only_if_confirmed"
    if triage.startswith("source_"):
        return "audit_source_name_scope_before_any_species_level_taxonomy_resolution"
    return "independent_second_backbone_audit_then_manual_review_if_discordant"


def build_triage(review: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    required = {
        "submitted_name",
        "match_type",
        "matched_rank",
        "matched_status",
        "candidate_accepted_rank",
        "candidate_accepted_status",
        "candidate_kingdom",
        "processing_flags",
        "resolution_status",
    }
    missing = sorted(required - set(review.columns))
    if missing:
        raise ValueError(f"alpha2 review queue missing columns: {missing}")
    out = review.copy().fillna("")
    if out["submitted_name"].duplicated().any():
        raise ValueError("alpha2 review queue contains duplicate submitted names")
    out["source_name_rank_class"] = out["submitted_name"].map(classify_name_rank)
    out["triage_class"] = out.apply(triage_class, axis=1)
    kingdom_conflict = (
        out["candidate_kingdom"].astype(str).str.strip().ne("")
        & out["candidate_kingdom"].astype(str).str.casefold().ne("plantae")
    )
    out["second_backbone_eligible"] = (
        out["source_name_rank_class"].eq("binomial_candidate") & ~kingdom_conflict
    ).map({True: "true", False: "false"})
    out["source_scope_status"] = "species_binomial_review_candidate"
    out.loc[out["triage_class"].str.startswith("source_"), "source_scope_status"] = (
        "non_species_name_scope_review"
    )
    out.loc[out["triage_class"].eq("kingdom_conflict_non_plantae"), "source_scope_status"] = (
        "kingdom_conflict_review"
    )
    out["next_action"] = out["triage_class"].map(_next_action)

    triage_counts = out["triage_class"].value_counts().sort_index().to_dict()
    rank_counts = out["source_name_rank_class"].value_counts().sort_index().to_dict()
    summary = {
        "n_input_review_rows": int(len(out)),
        "n_second_backbone_eligible": int(out["second_backbone_eligible"].eq("true").sum()),
        "n_source_scope_review": int(out["source_scope_status"].eq("non_species_name_scope_review").sum()),
        "n_kingdom_conflict_review": int(out["source_scope_status"].eq("kingdom_conflict_review").sum()),
        "triage_class_counts": triage_counts,
        "source_name_rank_counts": rank_counts,
    }
    return out, summary


def write_bundle(review_path: Path, output_dir: Path) -> dict[str, Any]:
    actual_sha = sha256_file(review_path)
    if actual_sha != ALPHA2_REVIEW_SHA256:
        raise ValueError(
            f"alpha2 review SHA mismatch: expected {ALPHA2_REVIEW_SHA256}, got {actual_sha}"
        )
    review = pd.read_csv(review_path, dtype=str).fillna("")
    if len(review) != ALPHA2_REVIEW_ROWS:
        raise ValueError(f"alpha2 review row count changed: {len(review)} != {ALPHA2_REVIEW_ROWS}")

    triage, summary = build_triage(review)
    second = triage.loc[triage["second_backbone_eligible"].eq("true")].copy()
    source_scope = triage.loc[triage["source_scope_status"].eq("non_species_name_scope_review")].copy()
    kingdom = triage.loc[triage["source_scope_status"].eq("kingdom_conflict_review")].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "taxonomy_review_triage.csv.gz": triage,
        "second_backbone_queue.csv.gz": second,
        "source_scope_review_queue.csv.gz": source_scope,
        "kingdom_conflict_queue.csv.gz": kingdom,
    }
    file_receipts: dict[str, dict[str, Any]] = {}
    for name, frame in paths.items():
        path = output_dir / name
        _write_gzip_csv(frame, path)
        file_receipts[name] = {"rows": int(len(frame)), "sha256": sha256_file(path)}

    manifest = {
        "schema_version": 1,
        "database_id": "global_island_plant_database",
        "version": "2.0.0-alpha3-taxonomy-triage",
        "source_layer": "2.0.0-alpha2-taxonomy",
        "source_review_sha256": actual_sha,
        "source_canonical_run": ALPHA2_CANONICAL_RUN,
        "source_canonical_artifact": ALPHA2_CANONICAL_ARTIFACT,
        "summary": summary,
        "files": file_receipts,
        "scientific_boundary": {
            "alpha1_mutated": False,
            "alpha2_taxonomy_mutated": False,
            "taxa_promoted": 0,
            "taxa_removed": 0,
            "second_backbone_results_present": False,
            "triage_only": True,
        },
    }
    (output_dir / "TRIAGE_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


@app.command("build")
def build(
    review_queue: Path = typer.Option(..., "--review-queue", exists=True, dir_okay=False),
    output_dir: Path = typer.Option(..., "--output-dir"),
) -> None:
    manifest = write_bundle(review_queue, output_dir)
    typer.echo(json.dumps(manifest["summary"], sort_keys=True))


if __name__ == "__main__":
    app()
