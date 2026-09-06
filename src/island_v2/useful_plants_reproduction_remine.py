"""Promote manually reviewed reproductive statements from cached pages.

No page is fetched here.  The candidate JSONL contains exact excerpts from
previously completed Ken Fern/PFAF pages, while the review table records a
species-trait decision and the underlying cited work.  Provider mirrors and
multiple species statements copied from the same monograph share one source
lineage.  Apomixis, cultivar-only claims and source conflicts remain excluded.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

app = typer.Typer(add_completion=False, no_args_is_help=True)

ALLOWED_TRAITS = {
    "self_incompatibility": "reproductive_assurance",
    "cleistogamy": "reproductive_assurance",
}
REVIEW_COLUMNS = {
    "accepted_species",
    "trait_name",
    "normalized_value",
    "decision",
    "reason",
    "underlying_citation_ref",
    "underlying_citation",
    "reviewer",
    "reviewed_at_utc",
}
GZIP = {"method": "gzip", "mtime": 0}
TEXT_SUFFIXES = {".csv", ".json", ".jsonl", ".md", ".toml", ".yml", ".yaml"}


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def _sha256(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in TEXT_SUFFIXES:
        payload = payload.replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def _lineage(ref: str) -> str:
    return f"useful-plants-underlying-citation:ref-{ref.casefold()}"


def build_package(
    candidates: list[dict[str, object]], review: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply frozen human decisions to cached exact-excerpt candidates."""

    missing = REVIEW_COLUMNS.difference(review.columns)
    if missing:
        raise ValueError(f"reproduction review table missing columns: {sorted(missing)}")
    reviewed = review.fillna("").copy()
    if reviewed[["accepted_species", "trait_name"]].duplicated().any():
        raise ValueError("review table duplicates a species x trait decision")
    review_map = {
        (_text(row["accepted_species"]), _text(row["trait_name"])): row
        for row in reviewed.to_dict("records")
    }
    candidate_keys = {
        (_text(row.get("accepted_species")), _text(row.get("trait_name"))) for row in candidates
    }
    if missing_reviews := candidate_keys.difference(review_map):
        raise ValueError(f"candidate rows lack review decisions: {sorted(missing_reviews)}")
    if extra_reviews := set(review_map).difference(candidate_keys):
        raise ValueError(f"review rows lack cached candidates: {sorted(extra_reviews)}")

    evidence_rows: list[dict[str, str]] = []
    audit_rows: list[dict[str, str]] = []
    for candidate in candidates:
        species = _text(candidate.get("accepted_species"))
        trait = _text(candidate.get("trait_name"))
        reviewed_row = review_map[(species, trait)]
        decision = _text(reviewed_row["decision"]).casefold()
        reason = _text(reviewed_row["reason"])
        expected_axis = ALLOWED_TRAITS.get(trait, "")
        value = _text(reviewed_row["normalized_value"])
        package_reason = "selected"
        if decision != "accept":
            package_reason = reason or "not_review_accepted"
        elif not expected_axis:
            package_reason = "trait_not_allowed"
        elif value != _text(candidate.get("normalized_value")):
            package_reason = "review_candidate_value_mismatch"
        elif not _text(candidate.get("supporting_excerpt")):
            package_reason = "missing_exact_supporting_excerpt"
        elif not _text(reviewed_row["underlying_citation_ref"]):
            package_reason = "missing_underlying_citation_lineage"
        elif not _text(reviewed_row["reviewer"]) or not _text(reviewed_row["reviewed_at_utc"]):
            package_reason = "incomplete_review_provenance"
        audit_rows.append(
            {
                "accepted_species": species,
                "trait_name": trait,
                "candidate_value": _text(candidate.get("normalized_value")),
                "provider": _text(candidate.get("provider")),
                "source_url": _text(candidate.get("source_url")),
                "source_excerpt": _text(candidate.get("supporting_excerpt")),
                "page_sha256": _text(candidate.get("page_sha256")),
                "excerpt_sha256": _text(candidate.get("excerpt_sha256")),
                "content_fingerprint": _text(candidate.get("content_fingerprint")),
                "review_decision": decision,
                "review_reason": reason,
                "package_decision": ("accept" if package_reason == "selected" else "exclude"),
                "package_reason": package_reason,
                "source_lineage": _lineage(_text(reviewed_row["underlying_citation_ref"])),
                "reviewer": _text(reviewed_row["reviewer"]),
                "reviewed_at_utc": _text(reviewed_row["reviewed_at_utc"]),
            }
        )
        if package_reason != "selected":
            continue
        provider = _text(candidate.get("provider"))
        provider_name = (
            "Useful Tropical/Temperate Plants" if provider == "ken_fern" else "Plants For A Future"
        )
        ref = _text(reviewed_row["underlying_citation_ref"])
        evidence_rows.append(
            {
                "accepted_species": species,
                "axis": expected_axis,
                "trait_name": trait,
                "normalized_value": value,
                "quality": "medium",
                "source_group": "useful_plants_citation_remine",
                "source_provider": provider_name,
                "source_url": _text(candidate.get("source_url")),
                "source_record_id": f"{provider}:" + _text(candidate.get("excerpt_sha256"))[:24],
                "source_citation": _text(reviewed_row["underlying_citation"]),
                "source_excerpt": _text(candidate.get("supporting_excerpt")),
                "evidence_scope": "species_direct",
                "name_match_method": "accepted_name_exact",
                "source_lineage": _lineage(ref),
                "lineage_method": "underlying_citation_ref_provider_mirror_collapse",
                "source_run_id": "local-wave36-cached-reproduction-remine-20260827",
                "source_artifact": "all_cached_reproductive_candidate_excerpts.jsonl.gz",
                "source_file": _text(candidate.get("cache_file")),
                "acceptance_contract": "manual_cached_exact_excerpt_reproduction_medium_v1",
            }
        )
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    audit = pd.DataFrame(audit_rows)
    if not evidence.empty:
        evidence = evidence.sort_values(
            ["accepted_species", "trait_name", "source_lineage", "source_provider"]
        ).reset_index(drop=True)
    return evidence, audit


@app.command("build")
def build(
    candidates_jsonl: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    review_csv: Annotated[Path, typer.Option(exists=True, dir_okay=False)],
    output_dir: Annotated[Path, typer.Option(file_okay=False)],
) -> None:
    candidates = [
        json.loads(line) for line in candidates_jsonl.open(encoding="utf-8") if line.strip()
    ]
    review = pd.read_csv(review_csv, dtype=str).fillna("")
    evidence, audit = build_package(candidates, review)
    output_dir.mkdir(parents=True, exist_ok=True)
    evidence_path = output_dir / "useful_plants_remine_direct_evidence.csv.gz"
    audit_path = output_dir / "useful_plants_remine_review_audit.csv.gz"
    evidence.to_csv(evidence_path, index=False, compression=GZIP)
    audit.to_csv(audit_path, index=False, compression=GZIP)
    accepted = evidence[["accepted_species", "trait_name"]].drop_duplicates()
    summary = {
        "contract": "useful_plants_reproduction_remine_v1",
        "candidate_excerpt_rows": len(candidates),
        "reviewed_species_trait": len(review),
        "accepted_evidence_rows_before_lineage_dedup": len(evidence),
        "accepted_species_trait_after_lineage_dedup": len(accepted),
        "accepted_source_lineages": int(evidence["source_lineage"].nunique()),
        "excluded_species_trait": int(review["decision"].str.casefold().ne("accept").sum()),
        "checks": {
            "provider_mirrors_collapsed_by_underlying_citation": True,
            "apomixis_excluded": True,
            "cultivar_conflict_excluded": True,
            "self_fertile_no_not_mapped": True,
            "autonomous_selfing_not_inferred": True,
            "family_inference_absent": True,
            "global_fallback_absent": True,
        },
        "input_sha256": {
            candidates_jsonl.name: _sha256(candidates_jsonl),
            review_csv.name: _sha256(review_csv),
        },
        "artifact_sha256": {
            evidence_path.name: _sha256(evidence_path),
            audit_path.name: _sha256(audit_path),
        },
    }
    (output_dir / "useful_plants_remine_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
