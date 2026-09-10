"""Recover upstream direct-source provenance behind Database 1.0 Validated-Low rules.

This is a rights/provenance audit only. It does not alter the immutable Database 1.0
species-axis ledger. The final source-scale rebuild retained the Wave53 Low layer and
added only reviewed direct evidence to previously empty cells. For each remaining
``validated-low:<genus>:<trait>:<state-set>`` lineage, this script therefore verifies
the corresponding eligible Wave53 ``current_min3`` rule and reconstructs its upstream
provenance from every direct species x trait record used by that genus rule.

All direct species in the genus x trait training cell are treated as provenance
dependencies, including counterexamples, because dominance and leave-one-out gates
were computed from the full direct set. This is conservative for redistribution.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from island_v2.chapter1_database_release import (
    family_for_lineage,
    load_lineage_family_map,
)


def _tokens(value: object) -> list[str]:
    return [part.strip() for part in str(value or "").split("|") if part.strip()]


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().str.casefold().isin({"1", "true", "yes", "y"})


def _canonical_state_set(value: object) -> str:
    text = str(value or "").strip()
    try:
        decoded = json.loads(text)
    except json.JSONDecodeError as exc:
        raise ValueError(f"invalid state-set JSON: {text}") from exc
    if not isinstance(decoded, list) or not decoded:
        raise ValueError(f"state set must be a non-empty JSON list: {text}")
    return json.dumps(sorted({str(item) for item in decoded}), separators=(",", ":"))


def _parse_validated_low(lineage: str) -> tuple[str, str, str]:
    prefix = "validated-low:"
    if not lineage.startswith(prefix):
        raise ValueError(f"not a Validated-Low lineage: {lineage}")
    rest = lineage[len(prefix) :]
    parts = rest.split(":", 2)
    if len(parts) != 3 or not parts[0] or not parts[1]:
        raise ValueError(f"malformed Validated-Low lineage: {lineage}")
    return parts[0], parts[1], _canonical_state_set(parts[2])


def recover(
    species_axis_path: Path,
    wave53_direct_path: Path,
    wave53_rules_path: Path,
    output_dir: Path,
    lineage_family_map_path: Path | None = None,
) -> dict[str, object]:
    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    required_coverage = {"accepted_species", "axis", "quality", "source_lineages"}
    missing = required_coverage.difference(coverage.columns)
    if missing:
        raise ValueError(f"species-axis ledger missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("species-axis ledger contains duplicate species x axis cells")

    low = coverage.loc[coverage["quality"].str.casefold().eq("low")].copy()
    low["derived_lineages"] = low["source_lineages"].map(
        lambda value: "|".join(t for t in _tokens(value) if t.startswith("validated-low:"))
    )
    target = low.loc[low["derived_lineages"].ne("")].copy()
    if target.empty:
        raise ValueError("final Database 1.0 contains no derived Validated-Low lineages")

    derived_lineages = sorted(
        {token for value in target["derived_lineages"] for token in _tokens(value)}
    )
    parsed = []
    for lineage in derived_lineages:
        genus, trait, state_set = _parse_validated_low(lineage)
        parsed.append(
            {
                "derived_lineage": lineage,
                "genus": genus,
                "trait_name": trait,
                "inferred_state_set": state_set,
            }
        )
    derived_rules = pd.DataFrame(parsed)

    rules = pd.read_csv(wave53_rules_path, dtype=str).fillna("")
    required_rules = {
        "setting",
        "genus",
        "trait_name",
        "inferred_state_set",
        "n_direct_species",
        "eligible",
        "diagnostic_only",
    }
    missing = required_rules.difference(rules.columns)
    if missing:
        raise ValueError(f"Wave53 rule audit missing columns: {sorted(missing)}")
    eligible = rules.loc[
        rules["setting"].eq("current_min3")
        & _truthy(rules["eligible"])
        & ~_truthy(rules["diagnostic_only"])
    ].copy()
    eligible["inferred_state_set"] = eligible["inferred_state_set"].map(_canonical_state_set)
    rule_keys = ["genus", "trait_name", "inferred_state_set"]
    if eligible.duplicated(rule_keys).any():
        raise ValueError("Wave53 contains duplicate eligible current_min3 rule keys")

    verified_rules = derived_rules.merge(
        eligible,
        on=rule_keys,
        how="left",
        validate="one_to_one",
        indicator=True,
        suffixes=("", "_wave53"),
    )
    unmatched = verified_rules.loc[verified_rules["_merge"].ne("both")].copy()

    direct = pd.read_csv(wave53_direct_path, dtype=str).fillna("")
    required_direct = {"accepted_species", "trait_name", "source_lineages"}
    missing = required_direct.difference(direct.columns)
    if missing:
        raise ValueError(f"Wave53 direct ledger missing columns: {sorted(missing)}")
    direct["genus"] = direct["accepted_species"].str.split().str[0].fillna("")
    direct = direct.loc[direct["genus"].ne("") & direct["trait_name"].ne("")].copy()

    overrides = load_lineage_family_map(lineage_family_map_path)
    direct_groups = {
        key: group.copy()
        for key, group in direct.groupby(["genus", "trait_name"], sort=False)
    }

    rule_rows: list[dict[str, object]] = []
    rule_support_lineages: dict[str, set[str]] = {}
    rule_support_families: dict[str, set[str]] = {}
    rules_with_direct_count_mismatch = 0
    rules_with_zero_support = 0

    for row in verified_rules.to_dict("records"):
        lineage = str(row["derived_lineage"])
        matched = row["_merge"] == "both"
        group = direct_groups.get((row["genus"], row["trait_name"]), pd.DataFrame())
        support_lineages: set[str] = set()
        direct_species = 0
        if not group.empty:
            direct_species = int(group["accepted_species"].nunique())
            for value in group["source_lineages"]:
                support_lineages.update(_tokens(value))
        support_families = {
            family_for_lineage(source, overrides) for source in support_lineages
        }
        expected_direct = int(row["n_direct_species"]) if matched and str(row.get("n_direct_species", "")) else 0
        count_matches = matched and direct_species == expected_direct
        if matched and not count_matches:
            rules_with_direct_count_mismatch += 1
        if matched and not support_lineages:
            rules_with_zero_support += 1
        rule_support_lineages[lineage] = support_lineages
        rule_support_families[lineage] = support_families
        rule_rows.append(
            {
                "derived_lineage": lineage,
                "genus": row["genus"],
                "trait_name": row["trait_name"],
                "inferred_state_set": row["inferred_state_set"],
                "wave53_rule_matched": matched,
                "expected_n_direct_species": expected_direct,
                "reconstructed_n_direct_species": direct_species,
                "direct_species_count_matches": count_matches,
                "support_source_lineage_count": len(support_lineages),
                "support_source_lineages": json.dumps(sorted(support_lineages)),
                "support_source_family_count": len(support_families),
                "support_source_families": "|".join(sorted(support_families)),
            }
        )
    rule_summary = pd.DataFrame(rule_rows)

    family_cell_counts: Counter[str] = Counter()
    family_support_lineages: dict[str, set[str]] = defaultdict(set)
    cell_rows: list[dict[str, object]] = []
    zero_support_cells = 0
    multi_family_cells = 0
    unmatched_rule_cells = 0

    matched_rule_set = set(rule_summary.loc[rule_summary["wave53_rule_matched"], "derived_lineage"])
    for row in target.to_dict("records"):
        derived = _tokens(row["derived_lineages"])
        supports: set[str] = set()
        families: set[str] = set()
        for lineage in derived:
            supports.update(rule_support_lineages.get(lineage, set()))
            families.update(rule_support_families.get(lineage, set()))
        if any(lineage not in matched_rule_set for lineage in derived):
            unmatched_rule_cells += 1
        if not supports:
            zero_support_cells += 1
        if len(families) > 1:
            multi_family_cells += 1
        family_cell_counts.update(families)
        for source in supports:
            family_support_lineages[family_for_lineage(source, overrides)].add(source)
        cell_rows.append(
            {
                "accepted_species": row["accepted_species"],
                "axis": row["axis"],
                "derived_lineages": row["derived_lineages"],
                "derived_rule_count": len(derived),
                "all_rules_matched_wave53": all(lineage in matched_rule_set for lineage in derived),
                "support_source_lineage_count": len(supports),
                "support_source_lineages": json.dumps(sorted(supports)),
                "support_source_family_count": len(families),
                "support_source_families": "|".join(sorted(families)),
            }
        )
    recovery = pd.DataFrame(cell_rows)

    provider_summary = pd.DataFrame(
        [
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
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    recovery.to_csv(output_dir / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv", index=False)
    provider_summary.to_csv(output_dir / "VALIDATED_LOW_SUPPORT_FAMILY_SUMMARY.csv", index=False)
    rule_summary.to_csv(output_dir / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv", index=False)
    unmatched.to_csv(output_dir / "VALIDATED_LOW_UNMATCHED_RULES.csv", index=False)

    matched_rules = int(rule_summary["wave53_rule_matched"].sum())
    summary: dict[str, object] = {
        "contract": "chapter1_database_v1_validated_low_wave53_rule_provenance_audit_v2",
        "final_low_cells": int(len(low)),
        "derived_low_target_cells": int(len(target)),
        "distinct_derived_lineages": int(len(rule_summary)),
        "wave53_rule_matched_lineages": matched_rules,
        "wave53_rule_match_rate": matched_rules / len(rule_summary),
        "unmatched_rule_cells": unmatched_rule_cells,
        "rules_with_direct_species_count_mismatch": rules_with_direct_count_mismatch,
        "rules_with_zero_support_lineages": rules_with_zero_support,
        "zero_support_cells": zero_support_cells,
        "multi_family_support_cells": multi_family_cells,
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
        f"- distinct derived rule lineages: **{len(rule_summary):,}**",
        f"- matched to eligible Wave53 `current_min3` rules: **{matched_rules:,} / {len(rule_summary):,} ({matched_rules / len(rule_summary):.1%})**",
        f"- rules with direct-species count mismatch: **{rules_with_direct_count_mismatch:,}**",
        f"- rules with zero recovered direct support lineages: **{rules_with_zero_support:,}**",
        f"- target cells with zero support: **{zero_support_cells:,}**",
        f"- target cells supported by >1 source family: **{multi_family_cells:,}**",
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
    parser.add_argument("--wave53-direct", type=Path, required=True)
    parser.add_argument("--wave53-rules", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lineage-family-map", type=Path)
    args = parser.parse_args()
    recover(
        args.species_axis,
        args.wave53_direct,
        args.wave53_rules,
        args.output_dir,
        lineage_family_map_path=args.lineage_family_map,
    )


if __name__ == "__main__":
    main()
