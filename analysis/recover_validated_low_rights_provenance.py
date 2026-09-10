"""Recover upstream source lineages behind immutable Database 1.0 Validated-Low cells.

This is a rights/provenance audit only. It never changes a trait value, quality tier,
or the immutable species-axis ledger. Database 1.0 inherited its Validated-Low cells
from the pinned historical rule sidecar. That sidecar records the direct source
lineages supporting each species x axis prediction; this script joins those supports
back to the final low cells and reports the provider/source families that must be
cleared before redistribution.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from island_v2.chapter1_database_release import (
    family_for_lineage,
    load_lineage_family_map,
)

LOW_SIDECAR_SHA256 = "b62d5ae133c7029ce7a67fb57afa230b18160ea6ea27c9fea5beb14a557b60a1"


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _tokens(value: object) -> list[str]:
    return [part.strip() for part in str(value or "").split("|") if part.strip()]


def _json_list(value: object) -> list[str]:
    try:
        decoded = json.loads(str(value or ""))
    except json.JSONDecodeError:
        return []
    if not isinstance(decoded, list):
        return []
    return [str(item).strip() for item in decoded if str(item).strip()]


def recover(
    species_axis_path: Path,
    sidecar_path: Path,
    output_dir: Path,
    lineage_family_map_path: Path | None = None,
) -> dict[str, object]:
    if sha256_file(sidecar_path) != LOW_SIDECAR_SHA256:
        raise ValueError("Pinned historical Validated-Low sidecar checksum mismatch")

    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    required_coverage = {"accepted_species", "axis", "quality", "source_lineages"}
    missing = required_coverage.difference(coverage.columns)
    if missing:
        raise ValueError(f"species-axis ledger missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("species-axis ledger contains duplicate species x axis cells")

    sidecar = pd.read_csv(sidecar_path, dtype=str).fillna("")
    required_sidecar = {
        "accepted_species",
        "axis",
        "support_source_lineages",
        "family_inference",
        "global_fallback",
    }
    missing = required_sidecar.difference(sidecar.columns)
    if missing:
        raise ValueError(f"historical sidecar missing columns: {sorted(missing)}")
    if sidecar.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("historical sidecar contains duplicate species x axis cells")

    forbidden = (
        sidecar["family_inference"].str.casefold().ne("false")
        | sidecar["global_fallback"].str.casefold().ne("false")
    )
    if forbidden.any():
        raise ValueError("historical sidecar contains forbidden family/global fallback rows")

    low = coverage.loc[coverage["quality"].str.casefold().eq("low")].copy()
    low["derived_lineages"] = low["source_lineages"].map(
        lambda value: "|".join(t for t in _tokens(value) if t.startswith("validated-low:"))
    )
    target = low.loc[low["derived_lineages"].ne("")].copy()
    if target.empty:
        raise ValueError("final Database 1.0 contains no derived Validated-Low lineages")

    joined = target.merge(
        sidecar[
            [
                "accepted_species",
                "axis",
                "support_source_lineages",
                "trait_names",
                "predicted_state_sets",
            ]
        ],
        on=["accepted_species", "axis"],
        how="left",
        validate="one_to_one",
    )
    joined["support_source_lineages"] = joined["support_source_lineages"].fillna("")
    joined["sidecar_matched"] = joined["support_source_lineages"].ne("")

    overrides = load_lineage_family_map(lineage_family_map_path)
    family_cell_counts: Counter[str] = Counter()
    family_support_lineages: dict[str, set[str]] = defaultdict(set)
    recovered_rows: list[dict[str, object]] = []
    rule_cells: Counter[str] = Counter()
    rule_support_families: dict[str, set[str]] = defaultdict(set)
    rule_support_lineages: dict[str, set[str]] = defaultdict(set)

    for row in joined.to_dict("records"):
        supports = _json_list(row.get("support_source_lineages"))
        support_families: set[str] = set()
        for lineage in supports:
            family = family_for_lineage(lineage, overrides)
            support_families.add(family)
            family_support_lineages[family].add(lineage)
        family_cell_counts.update(support_families)

        derived = _tokens(row.get("derived_lineages"))
        for lineage in derived:
            rule_cells[lineage] += 1
            rule_support_families[lineage].update(support_families)
            rule_support_lineages[lineage].update(supports)

        recovered_rows.append(
            {
                "accepted_species": row["accepted_species"],
                "axis": row["axis"],
                "derived_lineages": row.get("derived_lineages", ""),
                "sidecar_matched": bool(row.get("sidecar_matched")),
                "support_source_lineage_count": len(set(supports)),
                "support_source_lineages": json.dumps(sorted(set(supports))),
                "support_source_family_count": len(support_families),
                "support_source_families": "|".join(sorted(support_families)),
                "trait_names": row.get("trait_names", ""),
                "predicted_state_sets": row.get("predicted_state_sets", ""),
            }
        )

    recovery = pd.DataFrame(recovered_rows)
    provider_rows = [
        {
            "source_family": family,
            "validated_low_cell_mentions": int(count),
            "distinct_support_source_lineages": len(family_support_lineages[family]),
            "example_support_lineage": sorted(family_support_lineages[family])[0]
            if family_support_lineages[family]
            else "",
        }
        for family, count in family_cell_counts.most_common()
    ]
    provider_summary = pd.DataFrame(provider_rows)

    rule_rows = []
    for lineage, cell_count in rule_cells.most_common():
        families = sorted(rule_support_families[lineage])
        supports = sorted(rule_support_lineages[lineage])
        rule_rows.append(
            {
                "derived_lineage": lineage,
                "target_cells": int(cell_count),
                "support_source_family_count": len(families),
                "support_source_families": "|".join(families),
                "support_source_lineage_count": len(supports),
                "support_source_lineages": json.dumps(supports),
            }
        )
    rule_summary = pd.DataFrame(rule_rows)

    matched = int(recovery["sidecar_matched"].sum())
    zero_support = int(recovery["support_source_lineage_count"].eq(0).sum())
    multi_family = int(recovery["support_source_family_count"].gt(1).sum())
    stable_rule_count = 0
    # A rule is considered family-stable when every target cell carrying it maps to
    # one identical upstream family set. The union above is retained regardless.
    for derived_lineage in rule_cells:
        cell_sets = {
            value
            for value in recovery.loc[
                recovery["derived_lineages"].map(
                    lambda x, token=derived_lineage: token in _tokens(x)
                ),
                "support_source_families",
            ]
        }
        if len(cell_sets) == 1:
            stable_rule_count += 1

    output_dir.mkdir(parents=True, exist_ok=True)
    recovery.to_csv(output_dir / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv", index=False)
    provider_summary.to_csv(output_dir / "VALIDATED_LOW_SUPPORT_FAMILY_SUMMARY.csv", index=False)
    rule_summary.to_csv(output_dir / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv", index=False)

    summary: dict[str, object] = {
        "contract": "chapter1_database_v1_validated_low_upstream_provenance_audit_v1",
        "historical_sidecar_sha256": LOW_SIDECAR_SHA256,
        "final_low_cells": int(len(low)),
        "derived_low_target_cells": int(len(target)),
        "sidecar_matched_cells": matched,
        "sidecar_match_rate": matched / len(target),
        "zero_support_cells": zero_support,
        "multi_family_support_cells": multi_family,
        "distinct_derived_lineages": int(len(rule_summary)),
        "family_stable_derived_lineages": int(stable_rule_count),
        "distinct_support_source_families": int(len(provider_summary)),
        "lineage_family_override_rows_used": int(len(overrides)),
        "scientific_database_modified": False,
        "rights_granted": False,
    }
    (output_dir / "VALIDATED_LOW_PROVENANCE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        "# Validated-Low upstream provenance recovery",
        "",
        f"- final Low cells: **{len(low):,}**",
        f"- cells carrying derived `validated-low:*` lineages: **{len(target):,}**",
        f"- matched to pinned historical sidecar: **{matched:,} / {len(target):,} ({matched / len(target):.1%})**",
        f"- cells with zero recovered direct support lineages: **{zero_support:,}**",
        f"- cells supported by >1 source family: **{multi_family:,}**",
        f"- distinct derived rule lineages: **{len(rule_summary):,}**",
        f"- rule lineages with one stable family set across target cells: **{stable_rule_count:,}**",
        f"- distinct upstream source families: **{len(provider_summary):,}**",
        "",
        "## Largest upstream source families",
        "",
    ]
    for _, row in provider_summary.head(30).iterrows():
        lines.append(
            f"- `{row['source_family']}` — {int(row['validated_low_cell_mentions']):,} Low-cell mentions; "
            f"{int(row['distinct_support_source_lineages']):,} exact support lineages"
        )
    (output_dir / "VALIDATED_LOW_PROVENANCE_SUMMARY.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--sidecar", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lineage-family-map", type=Path)
    args = parser.parse_args()
    recover(
        args.species_axis,
        args.sidecar,
        args.output_dir,
        lineage_family_map_path=args.lineage_family_map,
    )


if __name__ == "__main__":
    main()
