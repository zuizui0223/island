#!/usr/bin/env python3
"""Rebuild the v2 trait base from source-backed species-level records only.

This deliberately excludes genus/family/global inference.  It preserves raw source
values and emits a species x target-trait gap ledger so acquisition can continue
until every eligible master species has either direct evidence or an explicit
not-found audit record.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

TARGET_TRAITS = (
    "flower_primary_color",
    "floral_form",
    "floral_symmetry",
    "pollination_functional_guild",
    "pollen_vector_mode",
    "sex_system",
    "mating_system",
    "self_incompatibility",
    "autonomous_selfing_capacity",
)

BIEN_MAP = {
    "flower color": "flower_primary_color",
    "flower pollination syndrome": "pollen_vector_mode",
    "whole plant sexual system": "sex_system",
}

AUSTRAITS_MAP = {
    "flower_colour": "flower_primary_color",
    "perianth_colour": "flower_primary_color",
    "flower_perianth_symmetry": "floral_symmetry",
    "flower_structural_sex_type": "sex_system",
}


def norm(value: object) -> str:
    return " ".join(str(value or "").strip().split())


def master_names(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or "accepted_species" not in reader.fieldnames:
            raise ValueError("accepted_species missing from master")
        return sorted({norm(row.get("accepted_species")) for row in reader if norm(row.get("accepted_species"))})


def iter_csv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else path.open
    kwargs = {"mode": "rt", "encoding": "utf-8", "newline": ""} if path.suffix == ".gz" else {"mode": "r", "encoding": "utf-8", "newline": ""}
    with opener(**kwargs) as handle:
        yield from csv.DictReader(handle)


def collect(root: Path, master: set[str]):
    evidence: dict[tuple[str, str, str, str, str], dict[str, str]] = {}
    raw_species_by_source: dict[str, set[str]] = defaultdict(set)

    for path in root.rglob("master_matched_trait_records*.csv.gz"):
        for row in iter_csv(path):
            species = norm(row.get("scrubbed_species_binomial") or row.get("accepted_species"))
            raw_trait = norm(row.get("bien_query_trait") or row.get("trait_name"))
            target = BIEN_MAP.get(raw_trait.lower())
            if not species or species not in master or not target:
                continue
            value = norm(row.get("trait_value"))
            if not value:
                continue
            raw_species_by_source["BIEN"].add(species)
            key = (species, target, value, "BIEN", raw_trait)
            evidence[key] = {
                "accepted_species": species,
                "trait_name": target,
                "candidate_value": value,
                "source_name": "BIEN",
                "source_trait": raw_trait,
                "source_url": norm(row.get("url_source")),
                "source_record_id": norm(row.get("id")),
                "evidence_scope": "species_direct",
                "evidence_status": "source_backed_unreviewed",
            }

    for path in root.rglob("austraits_all.csv"):
        for row in iter_csv(path):
            species = norm(row.get("taxon_name"))
            raw_trait = norm(row.get("trait_name"))
            target = AUSTRAITS_MAP.get(raw_trait)
            if not species or species not in master or not target:
                continue
            value = norm(row.get("trait_value"))
            if not value:
                continue
            raw_species_by_source["AusTraits"].add(species)
            key = (species, target, value, "AusTraits", raw_trait)
            evidence[key] = {
                "accepted_species": species,
                "trait_name": target,
                "candidate_value": value,
                "source_name": "AusTraits",
                "source_trait": raw_trait,
                "source_url": "",
                "source_record_id": norm(row.get("source_record_id")),
                "evidence_scope": "species_direct",
                "evidence_status": "source_backed_unreviewed",
            }
    return list(evidence.values()), raw_species_by_source


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--artifact-root", type=Path, required=True)
    parser.add_argument("--master-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    names = master_names(args.master_csv)
    master = set(names)
    rows, by_source = collect(args.artifact_root, master)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    fields = ["accepted_species", "trait_name", "candidate_value", "source_name", "source_trait", "source_url", "source_record_id", "evidence_scope", "evidence_status"]
    with gzip.open(args.output_dir / "species_direct_evidence.csv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda r: (r["accepted_species"], r["trait_name"], r["source_name"], r["candidate_value"])))

    cells = {(r["accepted_species"], r["trait_name"]) for r in rows}
    covered_species = {s for s, _ in cells}
    by_trait = {trait: len({s for s, t in cells if t == trait}) for trait in TARGET_TRAITS}

    with gzip.open(args.output_dir / "species_trait_gap_ledger.csv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accepted_species", "trait_name", "direct_evidence_present", "acquisition_status"])
        for species in names:
            for trait in TARGET_TRAITS:
                present = (species, trait) in cells
                writer.writerow([species, trait, str(present).lower(), "covered" if present else "pending_direct_search"])

    report = {
        "master_species": len(names),
        "target_traits": len(TARGET_TRAITS),
        "direct_evidence_rows": len(rows),
        "direct_species_trait_cells": len(cells),
        "species_with_any_mapped_target_trait": len(covered_species),
        "species_with_any_mapped_target_trait_rate": len(covered_species) / len(names) if names else 0,
        "species_by_source_for_mapped_targets": {k: len(v) for k, v in sorted(by_source.items())},
        "species_direct_by_target_trait": by_trait,
        "policy": {
            "taxonomic_resolution": "species_direct_only",
            "inference_included": False,
            "unmapped_source_traits_retained_in_raw_source_artifacts": True,
            "missing_cells": "pending_direct_search; never biological absence",
        },
    }
    (args.output_dir / "species_direct_coverage.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
