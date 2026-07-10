"""Build a review queue from validation-pilot trait/proxy candidates.

The queue is an adjudication worksheet, not curated data. It preserves the
source candidates, adds a stable candidate id and review priority, and leaves
human decision columns blank. Only rows later marked accepted with a final value
should be promoted into curated trait tables.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

DECISION_ALLOWED_VALUES = "accepted|rejected|needs_source_check"

OUTPUT_COLUMNS = [
    "review_candidate_id",
    "review_priority",
    "review_group",
    "review_action",
    "source_lane",
    "accepted_species",
    "trait_layer",
    "trait_name",
    "candidate_value",
    "candidate_class",
    "candidate_kind",
    "evidence_scope",
    "source",
    "source_type",
    "source_url",
    "source_citation",
    "matched_term",
    "basis_traits",
    "basis_family",
    "raw_description",
    "source_review_status",
    "adjudication_decision",
    "final_value",
    "decision_reason",
    "reviewer",
    "review_date",
    "decision_allowed_values",
    "decision_boundary",
]


@app.callback()
def main() -> None:
    """Create review queues for validation-pilot trait candidates."""


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).split())


def _read_csv(path: Path | None) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, dtype=str).fillna("")


def _candidate_id(row: dict[str, str]) -> str:
    parts = [
        row.get("source_lane", ""),
        row.get("accepted_species", ""),
        row.get("trait_name", ""),
        row.get("candidate_value", ""),
        row.get("source_url", ""),
        row.get("raw_description", ""),
        row.get("basis_traits", ""),
    ]
    digest = hashlib.sha1("\x1f".join(parts).encode("utf-8")).hexdigest()[:16]
    return f"vpilot_trait:{digest}"


def _priority_for(source_lane: str, evidence_scope: str, candidate_class: str) -> tuple[int, str, str]:
    if source_lane == "web_reported_scout" and evidence_scope == "species_direct":
        return (
            10,
            "P1_species_direct_reported",
            "Review first; accept only if the cited source explicitly supports the floral trait value.",
        )
    if source_lane == "web_reported_scout" and evidence_scope == "species_indirect":
        return (
            20,
            "P2_species_indirect_source_check",
            "Hold for stricter source check; reject if the source is genus-level, ambiguous, or not about the accepted species.",
        )
    if candidate_class == "proxy":
        return (
            30,
            "P3_proxy_candidate",
            "Review as proxy only; do not treat as an explicit reported trait value.",
        )
    return (
        40,
        "P4_descriptive_scout_source_check",
        "Review after reported/proxy candidates; require a clear floral-organ context before accepting.",
    )


def _base_row(
    *,
    source_lane: str,
    accepted_species: str,
    trait_layer: str,
    trait_name: str,
    candidate_value: str,
    candidate_class: str,
    candidate_kind: str,
    evidence_scope: str,
    source: str,
    source_type: str,
    source_url: str,
    source_citation: str,
    matched_term: str,
    basis_traits: str,
    basis_family: str,
    raw_description: str,
    source_review_status: str,
) -> dict[str, str]:
    priority, group, action = _priority_for(source_lane, evidence_scope, candidate_class)
    row = {
        "review_candidate_id": "",
        "review_priority": str(priority),
        "review_group": group,
        "review_action": action,
        "source_lane": source_lane,
        "accepted_species": accepted_species,
        "trait_layer": trait_layer,
        "trait_name": trait_name,
        "candidate_value": candidate_value,
        "candidate_class": candidate_class,
        "candidate_kind": candidate_kind,
        "evidence_scope": evidence_scope,
        "source": source,
        "source_type": source_type,
        "source_url": source_url,
        "source_citation": source_citation,
        "matched_term": matched_term,
        "basis_traits": basis_traits,
        "basis_family": basis_family,
        "raw_description": raw_description,
        "source_review_status": source_review_status,
        "adjudication_decision": "",
        "final_value": "",
        "decision_reason": "",
        "reviewer": "",
        "review_date": "",
        "decision_allowed_values": DECISION_ALLOWED_VALUES,
        "decision_boundary": (
            "Queue row only. No value is curated until adjudication_decision=accepted "
            "and final_value is filled by a reviewer."
        ),
    }
    row["review_candidate_id"] = _candidate_id(row)
    return row


def _rows_from_web(table: pd.DataFrame) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in table.to_dict("records"):
        species = _text(record.get("accepted_species"))
        trait_name = _text(record.get("trait_name"))
        value = _text(record.get("provisional_candidate_value"))
        if not species or not trait_name or not value:
            continue
        rows.append(
            _base_row(
                source_lane="web_reported_scout",
                accepted_species=species,
                trait_layer=_text(record.get("trait_layer")),
                trait_name=trait_name,
                candidate_value=value,
                candidate_class=_text(record.get("candidate_class")) or "reported",
                candidate_kind="reported_keyword_match",
                evidence_scope=_text(record.get("evidence_scope")),
                source=_text(record.get("source")),
                source_type=_text(record.get("source_type")),
                source_url=_text(record.get("source_url")),
                source_citation="",
                matched_term=_text(record.get("matched_term")),
                basis_traits="",
                basis_family="",
                raw_description=_text(record.get("raw_description")),
                source_review_status=_text(record.get("review_status")),
            )
        )
    return rows


def _rows_from_proxy(table: pd.DataFrame) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in table.to_dict("records"):
        species = _text(record.get("accepted_species"))
        trait_name = _text(record.get("trait_name"))
        value = _text(record.get("proxy_value"))
        if not species or not trait_name or not value:
            continue
        rows.append(
            _base_row(
                source_lane="trait_proxy_generator",
                accepted_species=species,
                trait_layer="proxy",
                trait_name=trait_name,
                candidate_value=value,
                candidate_class=_text(record.get("candidate_class")) or "proxy",
                candidate_kind="rule_based_proxy",
                evidence_scope=_text(record.get("evidence_scope")),
                source="",
                source_type=_text(record.get("source_type")),
                source_url="",
                source_citation="",
                matched_term="",
                basis_traits=_text(record.get("basis_traits")),
                basis_family=_text(record.get("basis_family")),
                raw_description=_text(record.get("raw_description")),
                source_review_status=_text(record.get("review_status")),
            )
        )
    return rows


def _rows_from_descriptive(table: pd.DataFrame) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in table.to_dict("records"):
        species = _text(record.get("accepted_species"))
        trait_name = _text(record.get("trait_name"))
        value = _text(record.get("provisional_candidate_value"))
        if not species or not trait_name or not value:
            continue
        rows.append(
            _base_row(
                source_lane="m0_descriptive_scout",
                accepted_species=species,
                trait_layer=_text(record.get("trait_layer")),
                trait_name=trait_name,
                candidate_value=value,
                candidate_class="reported",
                candidate_kind=_text(record.get("candidate_kind")) or "descriptive_keyword_match",
                evidence_scope="species_description",
                source=_text(record.get("source")),
                source_type=_text(record.get("description_type")),
                source_url=_text(record.get("source_url")),
                source_citation=_text(record.get("source_citation")),
                matched_term=_text(record.get("matched_term")),
                basis_traits="",
                basis_family="",
                raw_description=_text(record.get("evidence_excerpt")),
                source_review_status=_text(record.get("review_status")),
            )
        )
    return rows


def build_review_queue(
    web_reported: pd.DataFrame,
    trait_proxies: pd.DataFrame,
    descriptive: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Return a prioritized adjudication queue for candidate rows."""
    rows = []
    rows.extend(_rows_from_web(web_reported.fillna("") if not web_reported.empty else web_reported))
    rows.extend(_rows_from_proxy(trait_proxies.fillna("") if not trait_proxies.empty else trait_proxies))
    if descriptive is not None and not descriptive.empty:
        rows.extend(_rows_from_descriptive(descriptive.fillna("")))
    queue = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    if queue.empty:
        return queue
    queue = queue.drop_duplicates("review_candidate_id", keep="first")
    queue["_priority_sort"] = queue["review_priority"].astype(int)
    queue = queue.sort_values(
        ["_priority_sort", "accepted_species", "trait_name", "candidate_value", "review_candidate_id"],
        kind="stable",
    )
    return queue.drop(columns=["_priority_sort"]).reset_index(drop=True)


def review_queue_summary(queue: pd.DataFrame) -> dict[str, Any]:
    """Summarise queued work without treating rows as reviewed evidence."""
    if queue.empty:
        return {
            "n_review_candidates": 0,
            "n_species": 0,
            "n_review_decisions_recorded": 0,
            "n_accepted_ready_for_curation": 0,
            "review_group_counts": {},
            "source_lane_counts": {},
            "decision_allowed_values": DECISION_ALLOWED_VALUES.split("|"),
            "curation_rule": "No rows are curated by queue generation.",
        }
    decisions = queue["adjudication_decision"].astype(str).str.strip()
    accepted = queue.loc[decisions.eq("accepted")]
    return {
        "n_review_candidates": int(len(queue)),
        "n_species": int(queue["accepted_species"].nunique()),
        "n_review_decisions_recorded": int(decisions.ne("").sum()),
        "n_accepted_ready_for_curation": int(
            len(accepted.loc[accepted["final_value"].astype(str).str.strip().ne("")])
        ),
        "review_group_counts": {
            str(k): int(v) for k, v in queue["review_group"].value_counts().sort_index().items()
        },
        "source_lane_counts": {
            str(k): int(v) for k, v in queue["source_lane"].value_counts().sort_index().items()
        },
        "decision_allowed_values": DECISION_ALLOWED_VALUES.split("|"),
        "curation_rule": (
            "Queue generation records no biological decisions. Promote only rows with "
            "adjudication_decision=accepted and nonblank final_value."
        ),
    }


def accepted_for_curation(queue: pd.DataFrame) -> pd.DataFrame:
    """Extract accepted rows only, after human review decisions have been filled."""
    columns = [
        "review_candidate_id",
        "accepted_species",
        "trait_name",
        "final_value",
        "candidate_class",
        "source_lane",
        "source_url",
        "source_citation",
        "decision_reason",
        "reviewer",
        "review_date",
    ]
    if queue.empty:
        return pd.DataFrame(columns=columns)
    accepted = queue.loc[
        queue["adjudication_decision"].astype(str).str.strip().eq("accepted")
        & queue["final_value"].astype(str).str.strip().ne("")
    ].copy()
    return accepted[columns].reset_index(drop=True)


@app.command("build")
def build(
    web_reported_csv: Path = typer.Option(..., exists=True, help="web_reported_candidates.csv"),
    trait_proxy_csv: Path = typer.Option(..., exists=True, help="trait_proxy_candidates.csv"),
    output_dir: Path = typer.Option(..., help="Directory for review queue artifacts."),
    descriptive_csv: Path | None = typer.Option(
        None, exists=True, help="Optional m0_descriptive_candidates.csv"
    ),
) -> None:
    """Write a validation-pilot trait candidate review queue and summary."""
    queue = build_review_queue(
        _read_csv(web_reported_csv),
        _read_csv(trait_proxy_csv),
        _read_csv(descriptive_csv),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    queue.to_csv(output_dir / "trait_candidate_review_queue.csv", index=False)
    summary = review_queue_summary(queue)
    (output_dir / "trait_candidate_review_queue_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    typer.echo(
        f"Built {summary['n_review_candidates']} review candidate(s) "
        f"across {summary['n_species']} species."
    )


if __name__ == "__main__":
    app()
