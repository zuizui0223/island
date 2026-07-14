"""Audit an open bulk plant-trait source against the island v2 species master.

The pilot is intentionally conservative: it performs exact normalized-name matching
and explicit source-trait mapping only. It does not fuzzy-match taxa and does not
convert continuous source values into categorical v2 traits.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable


DEFAULT_TRAIT_MAP = {
    "flower_length": "flower_length_mm",
    "flower_diameter": "flower_diameter_mm",
    "flower_perianth_fusion": "flower_perianth_fusion_fraction",
    "flower_structural_sex_type": "flower_structural_sex_type_raw",
}


def normalize_scientific_name(value: str) -> str:
    """Normalize whitespace only; deliberately avoid silent taxonomic guessing."""
    return " ".join((value or "").strip().split())


def read_master_names(path: Path, name_column: str) -> set[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or name_column not in reader.fieldnames:
            raise ValueError(f"master column {name_column!r} not found")
        return {
            name
            for row in reader
            if (name := normalize_scientific_name(row.get(name_column, "")))
        }


def audit_bulk_source(
    source_path: Path,
    master_names: set[str],
    *,
    taxon_column: str,
    trait_column: str,
    trait_map: dict[str, str] | None = None,
) -> dict[str, object]:
    mapping = trait_map or DEFAULT_TRAIT_MAP
    source_taxa: set[str] = set()
    matched_taxa: set[str] = set()
    mapped_records = Counter()
    unmapped_traits = Counter()
    input_records = 0

    with source_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = {taxon_column, trait_column} - fields
        if missing:
            raise ValueError(f"source columns not found: {sorted(missing)}")

        for row in reader:
            input_records += 1
            taxon = normalize_scientific_name(row.get(taxon_column, ""))
            trait = (row.get(trait_column, "") or "").strip()
            if taxon:
                source_taxa.add(taxon)
                if taxon in master_names:
                    matched_taxa.add(taxon)
            if trait in mapping:
                mapped_records[mapping[trait]] += 1
            elif trait:
                unmapped_traits[trait] += 1

    unmatched_taxa = source_taxa - master_names
    match_rate = len(matched_taxa) / len(source_taxa) if source_taxa else 0.0
    mapped_total = sum(mapped_records.values())

    return {
        "input_records": input_records,
        "unique_source_taxa": len(source_taxa),
        "exact_master_matches": len(matched_taxa),
        "unmatched_source_taxa": len(unmatched_taxa),
        "exact_taxon_match_rate": match_rate,
        "mapped_records": mapped_total,
        "mapped_records_by_target_trait": dict(sorted(mapped_records.items())),
        "unmapped_trait_records": sum(unmapped_traits.values()),
        "top_unmapped_source_traits": dict(unmapped_traits.most_common(25)),
        "policy": {
            "taxon_matching": "exact_after_whitespace_normalization",
            "fuzzy_matching": False,
            "continuous_to_category_derivation": False,
        },
    }


def load_trait_map(path: Path | None) -> dict[str, str]:
    if path is None:
        return dict(DEFAULT_TRAIT_MAP)
    with path.open(encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in payload.items()
    ):
        raise ValueError("trait map must be a JSON object of string-to-string mappings")
    return payload


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-csv", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--source-taxon-column", default="taxon_name")
    parser.add_argument("--source-trait-column", default="trait_name")
    parser.add_argument("--master-name-column", default="scientific_name")
    parser.add_argument("--trait-map-json", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    master_names = read_master_names(args.master_csv, args.master_name_column)
    report = audit_bulk_source(
        args.source_csv,
        master_names,
        taxon_column=args.source_taxon_column,
        trait_column=args.source_trait_column,
        trait_map=load_trait_map(args.trait_map_json),
    )
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
