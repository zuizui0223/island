#!/usr/bin/env python3
"""Materialize machine-wave candidates as a review queue without promoting evidence.

This is an acquisition aid only.  It deliberately does not write to the strict
species-axis ledger.  A row can become strict evidence only after source review,
identity reconciliation, lineage deduplication, and conflict resolution.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

STRICT_REPRODUCTIVE_TRAITS = {
    "self_incompatibility",
    "autonomous_selfing_capacity",
    "mating_system",
    "cleistogamy",
}


def truthy(value: str) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def read_rows(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_rows(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build(input_csv: Path, output_dir: Path, wave_id: str) -> dict[str, object]:
    rows = read_rows(input_csv)
    if not rows:
        raise ValueError("machine candidate file is empty")
    required = {
        "accepted_species",
        "trait_name",
        "candidate_value",
        "source_url",
        "source_citation",
        "source_excerpt",
        "evidence_scope",
        "target_for_task",
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"machine candidate schema missing: {sorted(missing)}")

    eligible: list[dict[str, str]] = []
    non_strict: list[dict[str, str]] = []
    for row in rows:
        row = {key: str(value or "").strip() for key, value in row.items()}
        row["wave_id"] = wave_id
        row["strict_axis_candidate"] = "false"
        row["review_status"] = "unreviewed_machine_candidate"
        row["promotion_allowed"] = "false"
        row["promotion_note"] = (
            "Requires exact-source review, species identity gate, provenance lineage, "
            "and strict-ledger conflict resolution."
        )
        is_target = truthy(row.get("target_for_task", ""))
        species_direct = row.get("evidence_scope", "") == "species_direct"
        reproductive = row.get("trait_name", "") in STRICT_REPRODUCTIVE_TRAITS
        has_source = bool(row.get("source_url")) and bool(row.get("source_excerpt"))
        if is_target and species_direct and reproductive and has_source:
            row["strict_axis_candidate"] = "true"
            eligible.append(row)
        else:
            non_strict.append(row)

    key = lambda row: (
        row.get("trait_name", ""),
        row.get("accepted_species", ""),
        row.get("source_url", ""),
        row.get("candidate_value", ""),
    )
    eligible.sort(key=key)
    non_strict.sort(key=key)
    all_columns = list(rows[0]) + [
        "wave_id",
        "strict_axis_candidate",
        "review_status",
        "promotion_allowed",
        "promotion_note",
    ]
    write_rows(output_dir / "strict_reproductive_review_queue.csv", eligible, all_columns)
    write_rows(output_dir / "non_strict_or_incomplete_candidates.csv", non_strict, all_columns)

    summary: dict[str, object] = {
        "wave_id": wave_id,
        "input_candidates": len(rows),
        "strict_reproductive_review_candidates": len(eligible),
        "strict_candidate_species": len({r["accepted_species"] for r in eligible}),
        "strict_candidate_by_trait": dict(Counter(r["trait_name"] for r in eligible)),
        "strict_candidate_by_value": dict(Counter(r["candidate_value"] for r in eligible)),
        "non_strict_or_incomplete": len(non_strict),
        "promotion_policy": "review_queue_only_no_strict_promotion",
    }
    (output_dir / "review_queue_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--wave-id", required=True)
    args = parser.parse_args()
    print(json.dumps(build(args.input_csv, args.output_dir, args.wave_id), indent=2))


if __name__ == "__main__":
    main()
