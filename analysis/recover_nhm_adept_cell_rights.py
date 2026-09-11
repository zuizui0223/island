"""Recover cell-scoped NHM ADEPT rights overrides for Database 1.0.

The immutable Chapter 1 ledger deliberately preserves original eFloras treatment
lineages even when the structured trait value was imported through the licensed
NHM ADEPT dataset. A global lineage override would therefore be too broad.

This audit emits only (target species, axis, exact source lineage) overrides that
can be linked to a reviewed NHM ADEPT row. It also propagates those exact NHM
support lineages to derived Validated-Low target cells through the already-closed
historical provenance audit. It never changes scientific values or grants rights
outside those cell-lineage keys.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


PINNED_EVIDENCE_SHA256 = (
    "3b9204acecd8b66e1f58dd9ae7165a4b3c43375d8fe2fed5ecf30ea03a55e298"
)
NHM_SOURCE_FAMILY = "dataset:nhm_adept_2026"
RESOLVED_QUALITIES = {"high", "medium", "low"}
DIRECT_QUALITIES = {"high", "medium"}
DIRECT_SCOPES = {"species_direct", "synonym_direct"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tokens(value: object) -> list[str]:
    return [token.strip() for token in str(value or "").split("|") if token.strip()]


def recover(
    species_axis_path: Path,
    validated_low_recovery_path: Path,
    nhm_evidence_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    actual_sha = _sha256(nhm_evidence_path)
    if actual_sha != PINNED_EVIDENCE_SHA256:
        raise ValueError(
            "reviewed scale-source evidence SHA256 changed: "
            f"{actual_sha} != {PINNED_EVIDENCE_SHA256}"
        )

    evidence = pd.read_csv(nhm_evidence_path, dtype=str).fillna("")
    required_evidence = {
        "accepted_species",
        "axis",
        "quality",
        "evidence_scope",
        "source_provider",
        "source_lineage",
    }
    missing = required_evidence.difference(evidence.columns)
    if missing:
        raise ValueError(f"NHM evidence missing columns: {sorted(missing)}")

    nhm = evidence.loc[
        evidence["source_provider"].str.startswith("nhm_adept_")
        & evidence["quality"].str.casefold().isin(DIRECT_QUALITIES)
        & evidence["evidence_scope"].str.casefold().isin(DIRECT_SCOPES)
        & evidence["accepted_species"].ne("")
        & evidence["axis"].ne("")
        & evidence["source_lineage"].ne("")
    ].copy()
    if nhm.empty:
        raise ValueError("no reviewed NHM ADEPT evidence rows found")

    nhm = nhm.drop_duplicates(
        ["accepted_species", "axis", "source_lineage", "source_provider"]
    )
    nhm_lineages_by_axis: dict[str, set[str]] = defaultdict(set)
    providers_by_axis_lineage: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in nhm.to_dict("records"):
        axis = str(row["axis"])
        lineage = str(row["source_lineage"])
        provider = str(row["source_provider"])
        nhm_lineages_by_axis[axis].add(lineage)
        providers_by_axis_lineage[(axis, lineage)].add(provider)

    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    required_coverage = {"accepted_species", "axis", "quality", "source_lineages"}
    missing = required_coverage.difference(coverage.columns)
    if missing:
        raise ValueError(f"species-axis ledger missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("species-axis ledger contains duplicate species x axis cells")
    resolved = coverage.loc[
        coverage["quality"].str.casefold().isin(RESOLVED_QUALITIES)
    ].copy()
    coverage_lookup = {
        (str(row["accepted_species"]), str(row["axis"])): set(
            _tokens(row["source_lineages"])
        )
        for row in resolved.to_dict("records")
    }

    rows: list[dict[str, str]] = []

    # Direct cells require exact species + axis + lineage agreement with both
    # the reviewed NHM row and the immutable final Database 1.0 cell.
    for row in nhm.to_dict("records"):
        key = (str(row["accepted_species"]), str(row["axis"]))
        lineage = str(row["source_lineage"])
        if lineage not in coverage_lookup.get(key, set()):
            continue
        rows.append(
            {
                "accepted_species": key[0],
                "axis": key[1],
                "source_lineage": lineage,
                "source_family": NHM_SOURCE_FAMILY,
                "evidence_role": "direct_cell",
                "nhm_source_providers": str(row["source_provider"]),
            }
        )

    # Derived cells use exact support lineages from the closed historical audit.
    # A support lineage is reclassified only if the reviewed NHM packet contains
    # that same lineage on the same response axis.
    recovery = pd.read_csv(validated_low_recovery_path, dtype=str).fillna("")
    required_recovery = {
        "accepted_species",
        "axis",
        "support_source_lineages",
        "all_rules_source_provenance_recovered",
    }
    missing = required_recovery.difference(recovery.columns)
    if missing:
        raise ValueError(
            f"Validated-Low recovery missing columns: {sorted(missing)}"
        )
    for row in recovery.to_dict("records"):
        if str(row["all_rules_source_provenance_recovered"]).casefold() not in {
            "1",
            "true",
            "yes",
            "y",
        }:
            continue
        key = (str(row["accepted_species"]), str(row["axis"]))
        if key not in coverage_lookup:
            continue
        for lineage in _tokens(row["support_source_lineages"]):
            if lineage not in nhm_lineages_by_axis.get(key[1], set()):
                continue
            rows.append(
                {
                    "accepted_species": key[0],
                    "axis": key[1],
                    "source_lineage": lineage,
                    "source_family": NHM_SOURCE_FAMILY,
                    "evidence_role": "validated_low_support",
                    "nhm_source_providers": "|".join(
                        sorted(providers_by_axis_lineage[(key[1], lineage)])
                    ),
                }
            )

    result = pd.DataFrame(
        rows,
        columns=[
            "accepted_species",
            "axis",
            "source_lineage",
            "source_family",
            "evidence_role",
            "nhm_source_providers",
        ],
    )
    if result.empty:
        raise ValueError("NHM ADEPT evidence produced zero final-cell overrides")

    conflicts = (
        result.groupby(["accepted_species", "axis", "source_lineage"])[
            "source_family"
        ]
        .nunique()
        .loc[lambda series: series > 1]
    )
    if not conflicts.empty:
        raise ValueError(f"NHM cell map has {len(conflicts)} conflicting keys")

    grouped: list[dict[str, str]] = []
    for key, frame in result.groupby(
        ["accepted_species", "axis", "source_lineage"], sort=True
    ):
        grouped.append(
            {
                "accepted_species": key[0],
                "axis": key[1],
                "source_lineage": key[2],
                "source_family": NHM_SOURCE_FAMILY,
                "evidence_role": "|".join(sorted(set(frame["evidence_role"]))),
                "nhm_source_providers": "|".join(
                    sorted(
                        {
                            token
                            for value in frame["nhm_source_providers"]
                            for token in _tokens(value)
                        }
                    )
                ),
            }
        )
    mapped = pd.DataFrame(grouped)

    output_dir.mkdir(parents=True, exist_ok=True)
    mapped.to_csv(output_dir / "NHM_ADEPT_CELL_LINEAGE_FAMILY_MAP.csv", index=False)

    role_counts: Counter[str] = Counter()
    for value in mapped["evidence_role"]:
        role_counts.update(_tokens(value))
    provider_counts: Counter[str] = Counter()
    for value in mapped["nhm_source_providers"]:
        provider_counts.update(_tokens(value))

    summary: dict[str, object] = {
        "contract": "chapter1_database_v1_nhm_adept_cell_rights_v1",
        "pinned_evidence_sha256": PINNED_EVIDENCE_SHA256,
        "reviewed_nhm_evidence_rows": int(len(nhm)),
        "cell_lineage_override_rows": int(len(mapped)),
        "distinct_target_cells": int(
            mapped[["accepted_species", "axis"]].drop_duplicates().shape[0]
        ),
        "distinct_nhm_source_lineages": int(mapped["source_lineage"].nunique()),
        "role_counts": dict(role_counts),
        "provider_override_rows": dict(provider_counts),
        "source_family": NHM_SOURCE_FAMILY,
        "scientific_database_modified": False,
        "blanket_efloras_rights_granted": False,
    }
    (output_dir / "NHM_ADEPT_CELL_RIGHTS_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--validated-low-recovery", type=Path, required=True)
    parser.add_argument("--nhm-evidence", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    recover(
        args.species_axis,
        args.validated_low_recovery,
        args.nhm_evidence,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
