"""Convert a reviewed source-scale CSV packet into the strict restart JSON contract.

This is deliberately narrower than a generic CSV importer.  A source packet is
eligible only after source-level methodology review has established a direct
mapping to one strict trait.  Candidate discovery tables, mixed semantic fields,
and rows lacking original source references fail closed.

The adapter does not change coverage.  It only makes a completed source-scale
packet consumable by the existing add-only strict batch integrator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
REVIEW_STATUS = "source_methodology_reviewed_reference_backed_strict_direct"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def false_value(value: object) -> bool:
    return text(value).casefold() in {"false", "0", "no"}


def record_id(row: dict[str, object]) -> str:
    payload = "|".join(
        text(row.get(key))
        for key in (
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_lineage",
            "source_excel_row",
        )
    )
    return "source-packet:" + hashlib.sha256(payload.encode("utf-8")).hexdigest()[:24]


def convert(frame: pd.DataFrame, review_date: str, reviewer: str) -> list[dict[str, object]]:
    frame = frame.fillna("")
    required = {
        "accepted_species",
        "axis",
        "trait_name",
        "normalized_value",
        "quality",
        "source_url",
        "source_lineage",
        "source_reference_raw",
        "source_column",
        "source_raw_value",
        "name_match_method",
        "review_status",
        "promotion_allowed",
        "genus_rule_training_allowed",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"source packet missing columns: {sorted(missing)}")
    if not review_date or not reviewer:
        raise ValueError("review_date and reviewer are required")

    records: list[dict[str, object]] = []
    seen: set[tuple[str, str]] = set()
    for raw in frame.to_dict("records"):
        species = text(raw["accepted_species"])
        trait = text(raw["trait_name"])
        value = text(raw["normalized_value"])
        axis = text(raw["axis"])
        quality = text(raw["quality"]).casefold()
        refs = text(raw["source_reference_raw"])
        if not BINOMIAL.fullmatch(species):
            raise ValueError(f"non-binomial source-packet species: {species}")
        if not trait or not value or not axis:
            raise ValueError(f"missing trait/value/axis for {species}")
        if quality not in {"high", "medium"}:
            raise ValueError(f"unreviewed quality for {species}: {quality}")
        if text(raw["review_status"]) != REVIEW_STATUS:
            raise ValueError(f"source methodology not reviewed for {species}")
        if not false_value(raw["promotion_allowed"]):
            raise ValueError("source packet must remain unpromoted before batch integration")
        if not false_value(raw["genus_rule_training_allowed"]):
            raise ValueError("source packet cannot train genus rules at this stage")
        if not refs or not text(raw["source_url"]) or not text(raw["source_lineage"]):
            raise ValueError(f"missing source provenance for {species}")
        key = (species, trait)
        if key in seen:
            raise ValueError(f"duplicate species-trait in source packet: {key}")
        seen.add(key)
        excerpt = (
            f"{text(raw['source_column'])}={text(raw['source_raw_value'])}; "
            f"references={refs}"
        )
        records.append(
            {
                "record_id": record_id(raw),
                "accepted_species": species,
                "axis": axis,
                "trait_name": trait,
                "normalized_value": value,
                "state_set": [value],
                "quality": quality,
                "source_url": text(raw["source_url"]),
                "source_excerpt": excerpt,
                "source_lineage": text(raw["source_lineage"]),
                "review_date": review_date,
                "reviewer": reviewer,
                "name_match_method": text(raw["name_match_method"]),
                "cultivar_status": "species_database_record",
                "genus_rule_training_allowed": False,
                "review_status": "accepted_direct_statement",
                "source_packet_review_status": REVIEW_STATUS,
                "source_reference_raw": refs,
                "source_article_url": text(raw.get("source_article_url")),
            }
        )
    return records


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--review-date", required=True)
    p.add_argument("--reviewer", required=True)
    args = p.parse_args()
    frame = pd.read_csv(args.input, dtype=str).fillna("")
    records = convert(frame, args.review_date, args.reviewer)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(records, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"records": len(records), "species": len({r['accepted_species'] for r in records})}, indent=2))


if __name__ == "__main__":
    main()
