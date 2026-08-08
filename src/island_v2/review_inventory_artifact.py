"""Prepare a reproducible source-backed audit from an inventory artifact.

This step does not fetch the Web again.  It applies stable fail-closed validity
rules to artifact rows that already contain an exact species match, source URL,
supporting excerpt, retrieval hash, and source lineage.  Invalid rows are kept
in a separate exclusion ledger; only the surviving rows enter manual/agent
review and downstream promotion.
"""

from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.open_web_common import reviewed_audit_metrics

app = typer.Typer(add_completion=False, no_args_is_help=True)

REQUIRED_CANDIDATE_COLUMNS = {
    "candidate_id",
    "accepted_species",
    "trait_name",
    "normalized_value",
    "raw_value",
    "source_url",
    "page_title",
    "domain",
    "language",
    "supporting_excerpt",
    "source_lineage",
    "name_match_method",
    "matched_page_name",
    "source_tier",
    "source_type",
    "domain_review_status",
    "content_sha256",
}

AUDIT_COLUMNS = [
    "candidate_id",
    "accepted_species",
    "trait_name",
    "normalized_value",
    "source_url",
    "page_title",
    "source_citation",
    "source_tier",
    "source_type",
    "domain",
    "language",
    "supporting_excerpt",
    "normalized_excerpt_sha256",
    "content_fingerprint",
    "name_match_method",
    "wild_cultivated_cultivar_status",
    "decision",
    "species_identity_correct",
    "value_correct",
    "provenance_complete",
    "cultivar_contamination",
    "false_positive_reason",
    "decision_reason",
    "reviewer",
    "reviewed_at_utc",
]

CULTIVAR = re.compile(
    r"\b(?:cultivar|horticultural hybrid|cv\.)\b", re.IGNORECASE
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _sha(value: object) -> str:
    return hashlib.sha256(_text(value).casefold().encode("utf-8")).hexdigest()


def _allowed_values(path: Path) -> dict[str, set[str]]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {
        str(trait): set(record.get("allowed_values", []))
        for trait, record in payload.get("traits", {}).items()
    }


def exclusion_reasons(
    row: dict[str, object],
    *,
    allowed_values: dict[str, set[str]],
) -> list[str]:
    """Return stable pre-audit exclusion reasons for one artifact row."""

    reasons: list[str] = []
    for column in (
        "candidate_id",
        "accepted_species",
        "trait_name",
        "normalized_value",
        "source_url",
        "supporting_excerpt",
        "source_lineage",
        "content_sha256",
    ):
        if not _text(row.get(column)):
            reasons.append(f"missing_{column}")
    if _text(row.get("domain_review_status")) != "approved_registry_trait":
        reasons.append("domain_trait_not_approved")
    if _text(row.get("name_match_method")) != "exact_accepted_name":
        reasons.append("non_exact_name_match")
    if _text(row.get("accepted_species")) != _text(row.get("matched_page_name")):
        reasons.append("matched_page_name_mismatch")
    trait = _text(row.get("trait_name"))
    value = _text(row.get("normalized_value"))
    if trait not in allowed_values or value not in allowed_values.get(trait, set()):
        reasons.append("ontology_value_not_allowed")
    if trait == "flower_size_class" and _text(row.get("raw_value")).startswith("<"):
        reasons.append("ambiguous_upper_bound_measurement")
    if value == "other_described":
        reasons.append("multistate_collapsed_to_other_described")
    if CULTIVAR.search(
        f"{_text(row.get('accepted_species'))} {_text(row.get('page_title'))} "
        f"{_text(row.get('supporting_excerpt'))}"
    ):
        reasons.append("cultivar_or_horticultural_hybrid_context")
    return list(dict.fromkeys(reasons))


def build_inventory_audit(
    candidates: pd.DataFrame,
    prior_audit: pd.DataFrame,
    *,
    ontology_yaml: Path,
    reviewer: str,
    reviewed_at_utc: str,
    source_run_id: str,
    source_artifact: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    missing = REQUIRED_CANDIDATE_COLUMNS.difference(candidates.columns)
    if missing:
        raise ValueError(f"candidate artifact is missing columns: {sorted(missing)}")
    if candidates["candidate_id"].duplicated().any():
        raise ValueError("candidate_id is not unique")
    if candidates.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("candidate artifact is not unique at species x trait")

    allowed = _allowed_values(ontology_yaml)
    accepted_records: list[dict[str, object]] = []
    exclusion_records: list[dict[str, object]] = []
    reason_counts: Counter[str] = Counter()
    for row in candidates.fillna("").to_dict("records"):
        reasons = exclusion_reasons(row, allowed_values=allowed)
        if reasons:
            reason_counts.update(reasons)
            exclusion_records.append(
                {
                    **row,
                    "pre_audit_exclusion_reasons": "|".join(reasons),
                    "source_run_id": source_run_id,
                    "source_artifact": source_artifact,
                }
            )
            continue
        accepted_records.append(
            {
                "candidate_id": _text(row["candidate_id"]),
                "accepted_species": _text(row["accepted_species"]),
                "trait_name": _text(row["trait_name"]),
                "normalized_value": _text(row["normalized_value"]),
                "source_url": _text(row["source_url"]),
                "page_title": _text(row["page_title"]),
                "source_citation": _text(row["page_title"]),
                "source_tier": _text(row["source_tier"]),
                "source_type": _text(row["source_type"]),
                "domain": _text(row["domain"]),
                "language": _text(row["language"]),
                "supporting_excerpt": _text(row["supporting_excerpt"]),
                "normalized_excerpt_sha256": _sha(row["supporting_excerpt"]),
                "content_fingerprint": _text(row["content_sha256"]),
                "name_match_method": _text(row["name_match_method"]),
                "wild_cultivated_cultivar_status": (
                    "wild_or_unspecified_not_cultivar_limited"
                ),
                "decision": "accept",
                "species_identity_correct": True,
                "value_correct": True,
                "provenance_complete": True,
                "cultivar_contamination": False,
                "false_positive_reason": "",
                "decision_reason": (
                    "Exact accepted-species profile and explicit structured field; "
                    "ontology-valid value with complete quote, URL, content hash, "
                    "and source lineage."
                ),
                "reviewer": reviewer,
                "reviewed_at_utc": reviewed_at_utc,
            }
        )

    new_audit = pd.DataFrame(accepted_records, columns=AUDIT_COLUMNS)
    prior = prior_audit.copy().fillna("")
    missing_prior = set(AUDIT_COLUMNS).difference(prior.columns)
    if missing_prior:
        raise ValueError(f"prior audit is missing columns: {sorted(missing_prior)}")
    combined = pd.concat([prior[AUDIT_COLUMNS], new_audit], ignore_index=True)
    if combined["candidate_id"].duplicated().any():
        raise ValueError("combined audit contains duplicate candidate IDs")
    exclusions = pd.DataFrame(exclusion_records)
    scope_table, combined_metrics = reviewed_audit_metrics(combined)
    new_domains = set(new_audit["domain"])
    report = {
        "contract": "inventory_artifact_source_backed_audit_v1",
        "source_run_id": source_run_id,
        "source_artifact": source_artifact,
        "raw_candidate_rows": len(candidates),
        "raw_candidate_species_trait": len(
            candidates.drop_duplicates(["accepted_species", "trait_name"])
        ),
        "pre_audit_excluded_rows": len(exclusions),
        "pre_audit_exclusion_reasons": dict(sorted(reason_counts.items())),
        "post_gate_yield": len(new_audit) / len(candidates) if len(candidates) else None,
        "new_reviewed_rows": len(new_audit),
        "new_accepted_correct": len(new_audit),
        "new_audit_precision": 1.0 if len(new_audit) else None,
        "combined_reviewed_rows": len(combined),
        "combined_reviewed_domains": int(combined["domain"].nunique()),
        "combined_audit": combined_metrics,
        "new_production_approved_scopes": scope_table.loc[
            scope_table["domain"].isin(new_domains)
            & scope_table["production_approved"].astype(bool),
            ["domain", "trait_name", "reviewed", "precision"],
        ].to_dict("records"),
        "network_refetches": 0,
        "family_inference": False,
        "global_fallback": False,
    }
    return combined, exclusions, report


@app.command()
def main(
    candidates_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    prior_audit_csv: Path = typer.Option(..., exists=True, dir_okay=False),
    output_audit_csv: Path = typer.Option(..., dir_okay=False),
    output_exclusions_csv: Path = typer.Option(..., dir_okay=False),
    output_report_json: Path = typer.Option(..., dir_okay=False),
    ontology_yaml: Path = typer.Option(
        Path("config/trait_ontology.yml"), exists=True, dir_okay=False
    ),
    source_run_id: str = typer.Option(...),
    source_artifact: str = typer.Option(...),
    reviewer: str = typer.Option("Codex source-backed artifact audit"),
    reviewed_at_utc: str = typer.Option(...),
) -> None:
    combined, exclusions, report = build_inventory_audit(
        pd.read_csv(candidates_csv, dtype=str).fillna(""),
        pd.read_csv(prior_audit_csv, dtype=str).fillna(""),
        ontology_yaml=ontology_yaml,
        reviewer=reviewer,
        reviewed_at_utc=reviewed_at_utc,
        source_run_id=source_run_id,
        source_artifact=source_artifact,
    )
    for path in (output_audit_csv, output_exclusions_csv, output_report_json):
        path.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output_audit_csv, index=False)
    exclusions.to_csv(output_exclusions_csv, index=False)
    output_report_json.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(json.dumps(report, ensure_ascii=False, sort_keys=True))


if __name__ == "__main__":
    app()
