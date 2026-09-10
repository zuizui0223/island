"""Evaluate redistribution rights for every resolved Chapter 1 Database 1.0 cell.

This is a release audit only. It never changes scientific trait values or quality
labels. Direct cells are evaluated from their stored source lineages. Cells carrying
``validated-low:*`` derived lineages replace those synthetic lineage tokens with the
upstream source families reconstructed by the Validated-Low provenance audit.

A cell is publicly redistributable only when:

1. every derived rule has upstream source provenance;
2. every required upstream/direct source family is explicitly redistributable and
   carries a non-empty license decision in the source policy; and
3. no unresolved rule-semantic mismatch remains for that cell.

The final Database 1.0 release is ready only when every resolved cell passes.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import yaml

from island_v2.chapter1_database_release import (
    family_for_lineage,
    load_lineage_family_map,
    source_decision,
)


RESOLVED_QUALITIES = {"high", "medium", "low"}


def _tokens(value: object) -> list[str]:
    return [token.strip() for token in str(value or "").split("|") if token.strip()]


def _truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def evaluate(
    species_axis_path: Path,
    validated_low_recovery_path: Path,
    source_policy_path: Path,
    output_dir: Path,
    lineage_family_map_path: Path | None = None,
) -> dict[str, object]:
    policy = yaml.safe_load(source_policy_path.read_text(encoding="utf-8")) or {}
    overrides = load_lineage_family_map(lineage_family_map_path)

    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    required = {"accepted_species", "axis", "quality", "source_lineages"}
    missing = required.difference(coverage.columns)
    if missing:
        raise ValueError(f"species-axis ledger missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("species-axis ledger contains duplicate species x axis cells")
    resolved = coverage.loc[
        coverage["quality"].str.strip().str.casefold().isin(RESOLVED_QUALITIES)
    ].copy()

    recovery = pd.read_csv(validated_low_recovery_path, dtype=str).fillna("")
    required_recovery = {
        "accepted_species",
        "axis",
        "derived_lineages",
        "support_source_families",
        "all_rules_source_provenance_recovered",
        "all_rules_provenance_recovered",
    }
    missing = required_recovery.difference(recovery.columns)
    if missing:
        raise ValueError(f"Validated-Low recovery missing columns: {sorted(missing)}")
    if recovery.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("Validated-Low recovery contains duplicate species x axis cells")
    recovery_lookup = {
        (row["accepted_species"], row["axis"]): row
        for row in recovery.to_dict("records")
    }

    decision_cache: dict[str, tuple[str, str, str]] = {}

    def decision(family: str) -> tuple[str, str, str]:
        if family not in decision_cache:
            decision_cache[family] = source_decision(family, policy)
        return decision_cache[family]

    rows: list[dict[str, object]] = []
    blocked_family_cells: Counter[str] = Counter()
    blocked_family_examples: dict[str, tuple[str, str]] = {}
    status_counts: Counter[str] = Counter()
    derived_cells = 0
    derived_source_complete = 0
    semantic_mismatch_cells = 0

    for row in resolved.to_dict("records"):
        key = (str(row["accepted_species"]), str(row["axis"]))
        stored = _tokens(row["source_lineages"])
        derived = [token for token in stored if token.startswith("validated-low:")]
        direct = [token for token in stored if not token.startswith("validated-low:")]
        direct_families = {family_for_lineage(token, overrides) for token in direct}

        support_families: set[str] = set()
        source_provenance_complete = True
        semantic_exact = True
        if derived:
            derived_cells += 1
            recovered = recovery_lookup.get(key)
            if recovered is None:
                source_provenance_complete = False
                semantic_exact = False
            else:
                recovered_derived = set(_tokens(recovered["derived_lineages"]))
                if not set(derived).issubset(recovered_derived):
                    source_provenance_complete = False
                support_families.update(_tokens(recovered["support_source_families"]))
                source_provenance_complete = (
                    source_provenance_complete
                    and _truthy(recovered["all_rules_source_provenance_recovered"])
                    and bool(support_families)
                )
                semantic_exact = _truthy(recovered["all_rules_provenance_recovered"])
            if source_provenance_complete:
                derived_source_complete += 1
            if not semantic_exact:
                semantic_mismatch_cells += 1

        required_families = direct_families | support_families
        family_decisions = {family: decision(family) for family in required_families}
        blocked_families = sorted(
            family
            for family, (status, license_id, _note) in family_decisions.items()
            if status != "redistributable" or not license_id
        )

        if not source_provenance_complete:
            release_status = "provenance_review_required"
        elif not semantic_exact:
            release_status = "semantic_review_required"
        elif not required_families:
            release_status = "provenance_review_required"
        elif blocked_families:
            release_status = "rights_review_required"
        else:
            release_status = "redistributable"
        status_counts[release_status] += 1

        if release_status == "rights_review_required":
            for family in blocked_families:
                blocked_family_cells[family] += 1
                blocked_family_examples.setdefault(family, key)

        rows.append(
            {
                "accepted_species": key[0],
                "axis": key[1],
                "quality": row["quality"],
                "has_validated_low_lineage": bool(derived),
                "derived_rule_count": len(derived),
                "source_provenance_complete": source_provenance_complete,
                "exact_rule_semantics": semantic_exact,
                "required_source_family_count": len(required_families),
                "required_source_families": "|".join(sorted(required_families)),
                "blocked_source_family_count": len(blocked_families),
                "blocked_source_families": "|".join(blocked_families),
                "release_status": release_status,
            }
        )

    cell_audit = pd.DataFrame(rows)
    release_ready = bool(len(cell_audit)) and bool(
        cell_audit["release_status"].eq("redistributable").all()
    )

    family_rows: list[dict[str, object]] = []
    for family, count in blocked_family_cells.most_common():
        status, license_id, note = decision(family)
        example = blocked_family_examples[family]
        family_rows.append(
            {
                "source_family": family,
                "blocked_resolved_cells": int(count),
                "example_accepted_species": example[0],
                "example_axis": example[1],
                "redistribution_status": status,
                "source_license": license_id,
                "policy_note": note,
            }
        )
    family_blockers = pd.DataFrame(
        family_rows,
        columns=[
            "source_family",
            "blocked_resolved_cells",
            "example_accepted_species",
            "example_axis",
            "redistribution_status",
            "source_license",
            "policy_note",
        ],
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    cell_audit.to_csv(
        output_dir / "CELL_RELEASE_RIGHTS.csv.gz", index=False, compression="gzip"
    )
    cell_audit.loc[cell_audit["has_validated_low_lineage"]].to_csv(
        output_dir / "VALIDATED_LOW_CELL_RIGHTS.csv.gz", index=False, compression="gzip"
    )
    family_blockers.to_csv(output_dir / "CELL_RIGHTS_BLOCKER_FAMILIES.csv", index=False)

    summary: dict[str, object] = {
        "contract": "chapter1_database_cell_level_redistribution_gate_v1",
        "resolved_cells": int(len(cell_audit)),
        "redistributable_cells": int(status_counts["redistributable"]),
        "blocked_cells": int(len(cell_audit) - status_counts["redistributable"]),
        "release_status_counts": dict(status_counts),
        "validated_low_cells": int(derived_cells),
        "validated_low_cells_with_complete_upstream_source_provenance": int(
            derived_source_complete
        ),
        "semantic_mismatch_cells": int(semantic_mismatch_cells),
        "distinct_blocker_source_families": int(len(family_blockers)),
        "synthetic_validated_low_family_used_as_rights_decision": False,
        "release_ready_for_public_zenodo": release_ready,
        "scientific_database_modified": False,
        "rights_granted_by_provenance_recovery": False,
    }
    (output_dir / "FINAL_RELEASE_GATE.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# Chapter 1 Database 1.0 cell-level release gate",
        "",
        f"- resolved cells evaluated: **{len(cell_audit):,}**",
        f"- redistributable now: **{status_counts['redistributable']:,}**",
        f"- blocked/review-required: **{len(cell_audit) - status_counts['redistributable']:,}**",
        f"- Validated-Low cells evaluated through recovered upstream sources: **{derived_cells:,}**",
        f"- Validated-Low cells with complete upstream provenance: **{derived_source_complete:,} / {derived_cells:,}**",
        f"- semantic mismatch cells: **{semantic_mismatch_cells:,}**",
        f"- distinct source families currently blocking cells: **{len(family_blockers):,}**",
        f"- public Zenodo release ready: **{release_ready}**",
        "",
        "The synthetic `derived:validated_low_without_direct_source_lineage` family is not used as a final rights decision. Each derived Low cell is evaluated against all reconstructed upstream source families.",
        "",
        "## Cell status",
        "",
    ]
    for status, count in status_counts.most_common():
        lines.append(f"- `{status}` — {count:,} cells")
    lines += ["", "## Largest blocker source families", ""]
    for _, item in family_blockers.head(30).iterrows():
        lines.append(
            f"- `{item['source_family']}` — {int(item['blocked_resolved_cells']):,} blocked cells"
        )
    (output_dir / "FINAL_RELEASE_GATE.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--validated-low-recovery", type=Path, required=True)
    parser.add_argument("--source-policy", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lineage-family-map", type=Path)
    args = parser.parse_args()
    evaluate(
        args.species_axis,
        args.validated_low_recovery,
        args.source_policy,
        args.output_dir,
        lineage_family_map_path=args.lineage_family_map,
    )


if __name__ == "__main__":
    main()
