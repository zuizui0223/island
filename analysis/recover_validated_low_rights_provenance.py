"""Recover upstream direct-source provenance behind Database 1.0 Validated-Low rules.

Rights/provenance audit only. The immutable Database 1.0 species-axis ledger is
never changed and no redistribution right is granted here.

The final Database 1.0 retained a historical Validated-Low layer. Provenance is
therefore reconstructed against the exact rule generations that can be verified:

1. Wave53 eligible non-diagnostic ``current_min3`` rules, trained on the union of
   ordinary direct evidence and the Wave53 external-congener direct ledger.
2. The earlier formal integrated baseline (Run 32932103226) for rules that were
   still present in Database 1.0 but were no longer eligible/current in Wave53.

Rules not exactly recovered from either verified generation remain fail-closed.
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
    return (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.casefold()
        .isin({"1", "true", "yes", "y"})
    )


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


def _load_rules(path: Path, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    rules = pd.read_csv(path, dtype=str).fillna("")
    required = {
        "setting",
        "genus",
        "trait_name",
        "inferred_state_set",
        "n_direct_species",
        "eligible",
        "diagnostic_only",
    }
    missing = required.difference(rules.columns)
    if missing:
        raise ValueError(f"{label} rule audit missing columns: {sorted(missing)}")
    current = rules.loc[rules["setting"].eq("current_min3")].copy()
    current["inferred_state_set"] = current["inferred_state_set"].map(_canonical_state_set)
    eligible = current.loc[
        _truthy(current["eligible"]) & ~_truthy(current["diagnostic_only"])
    ].copy()
    keys = ["genus", "trait_name", "inferred_state_set"]
    if eligible.duplicated(keys).any():
        raise ValueError(f"{label} contains duplicate eligible current_min3 rule keys")
    return current, eligible


def _load_direct(
    paths: list[Path],
    label: str,
) -> tuple[pd.DataFrame, dict[tuple[str, str], pd.DataFrame]]:
    frames: list[pd.DataFrame] = []
    required = {"accepted_species", "trait_name", "source_lineages"}
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(
                f"{label} direct ledger {path} missing columns: {sorted(missing)}"
            )
        frame = frame.copy()
        frame["_provenance_input"] = path.name
        frames.append(frame)
    direct = pd.concat(frames, ignore_index=True, sort=False).fillna("")
    direct["genus"] = direct["accepted_species"].str.split().str[0].fillna("")
    direct = direct.loc[
        direct["genus"].ne("") & direct["trait_name"].ne("")
    ].copy()
    groups = {
        key: group.copy()
        for key, group in direct.groupby(["genus", "trait_name"], sort=False)
    }
    return direct, groups


def _diagnose_wave53_unmatched(
    row: dict[str, object],
    wave53_current: pd.DataFrame,
) -> str:
    same_gt = wave53_current.loc[
        wave53_current["genus"].eq(str(row["genus"]))
        & wave53_current["trait_name"].eq(str(row["trait_name"]))
    ]
    if same_gt.empty:
        return "no_genus_trait_in_wave53"
    same_state = same_gt.loc[
        same_gt["inferred_state_set"].eq(str(row["inferred_state_set"]))
    ]
    if not same_state.empty:
        return "current_min3_same_state_ineligible"
    return "genus_trait_present_state_diff"


def recover(
    species_axis_path: Path,
    wave53_direct_path: Path,
    wave53_external_congener_path: Path,
    wave53_rules_path: Path,
    output_dir: Path,
    lineage_family_map_path: Path | None = None,
    run329_direct_path: Path | None = None,
    run329_rules_path: Path | None = None,
) -> dict[str, object]:
    if (run329_direct_path is None) != (run329_rules_path is None):
        raise ValueError("Run329 direct and rule files must be supplied together")

    coverage = pd.read_csv(species_axis_path, dtype=str).fillna("")
    required_coverage = {"accepted_species", "axis", "quality", "source_lineages"}
    missing = required_coverage.difference(coverage.columns)
    if missing:
        raise ValueError(f"species-axis ledger missing columns: {sorted(missing)}")
    if coverage.duplicated(["accepted_species", "axis"]).any():
        raise ValueError("species-axis ledger contains duplicate species x axis cells")

    low = coverage.loc[coverage["quality"].str.casefold().eq("low")].copy()
    low["derived_lineages"] = low["source_lineages"].map(
        lambda value: "|".join(
            t for t in _tokens(value) if t.startswith("validated-low:")
        )
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
    keys = ["genus", "trait_name", "inferred_state_set"]

    wave53_current, wave53_eligible = _load_rules(wave53_rules_path, "Wave53")
    _, wave53_groups = _load_direct(
        [wave53_direct_path, wave53_external_congener_path],
        "Wave53",
    )

    run329_eligible = pd.DataFrame()
    run329_groups: dict[tuple[str, str], pd.DataFrame] = {}
    if run329_rules_path is not None and run329_direct_path is not None:
        _, run329_eligible = _load_rules(run329_rules_path, "Run329")
        _, run329_groups = _load_direct([run329_direct_path], "Run329")

    wave53_lookup = {
        tuple(row[k] for k in keys): row
        for row in wave53_eligible.to_dict("records")
    }
    run329_lookup = {
        tuple(row[k] for k in keys): row
        for row in run329_eligible.to_dict("records")
    }

    overrides = load_lineage_family_map(lineage_family_map_path)
    rule_rows: list[dict[str, object]] = []
    rule_support_lineages: dict[str, set[str]] = {}
    rule_support_families: dict[str, set[str]] = {}
    tier_counts: Counter[str] = Counter()
    mismatch_reasons: Counter[str] = Counter()
    rules_with_direct_count_mismatch = 0
    rules_with_zero_support = 0

    for row in derived_rules.to_dict("records"):
        lineage = str(row["derived_lineage"])
        key = tuple(str(row[k]) for k in keys)
        tier = "unresolved_historical_secondary"
        matched_rule: dict[str, object] | None = None
        groups: dict[tuple[str, str], pd.DataFrame] = {}

        if key in wave53_lookup:
            tier = "wave53_current_min3"
            matched_rule = wave53_lookup[key]
            groups = wave53_groups
        elif key in run329_lookup:
            tier = "run329_formal_current_min3"
            matched_rule = run329_lookup[key]
            groups = run329_groups

        tier_counts[tier] += 1
        diagnosis = ""
        if matched_rule is None:
            diagnosis = _diagnose_wave53_unmatched(row, wave53_current)
            mismatch_reasons[diagnosis] += 1

        group = groups.get((str(row["genus"]), str(row["trait_name"])), pd.DataFrame())
        support_lineages: set[str] = set()
        direct_species = 0
        if matched_rule is not None and not group.empty:
            direct_species = int(group["accepted_species"].nunique())
            for value in group["source_lineages"]:
                support_lineages.update(_tokens(value))
        support_families = {
            family_for_lineage(source, overrides) for source in support_lineages
        }
        expected_direct = (
            int(str(matched_rule["n_direct_species"]))
            if matched_rule is not None and str(matched_rule.get("n_direct_species", ""))
            else 0
        )
        count_matches = matched_rule is not None and direct_species == expected_direct
        if matched_rule is not None and not count_matches:
            rules_with_direct_count_mismatch += 1
        if matched_rule is not None and not support_lineages:
            rules_with_zero_support += 1

        rule_support_lineages[lineage] = support_lineages
        rule_support_families[lineage] = support_families
        rule_rows.append(
            {
                "derived_lineage": lineage,
                "genus": row["genus"],
                "trait_name": row["trait_name"],
                "inferred_state_set": row["inferred_state_set"],
                "provenance_rule_tier": tier,
                "rule_recovered": matched_rule is not None,
                "unmatched_wave53_reason": diagnosis,
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
    unresolved_rule_cells = 0

    recovered_rule_set = set(
        rule_summary.loc[rule_summary["rule_recovered"], "derived_lineage"]
    )
    tier_by_lineage = dict(
        zip(
            rule_summary["derived_lineage"],
            rule_summary["provenance_rule_tier"],
            strict=True,
        )
    )
    for row in target.to_dict("records"):
        derived = _tokens(row["derived_lineages"])
        supports: set[str] = set()
        families: set[str] = set()
        tiers: set[str] = set()
        for lineage in derived:
            supports.update(rule_support_lineages.get(lineage, set()))
            families.update(rule_support_families.get(lineage, set()))
            if lineage in tier_by_lineage:
                tiers.add(tier_by_lineage[lineage])
        all_recovered = all(lineage in recovered_rule_set for lineage in derived)
        if not all_recovered:
            unresolved_rule_cells += 1
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
                "all_rules_provenance_recovered": all_recovered,
                "provenance_rule_tiers": "|".join(sorted(tiers)),
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
                "example_support_lineage": (
                    sorted(family_support_lineages[family])[0]
                    if family_support_lineages[family]
                    else ""
                ),
            }
            for family, count in family_cell_counts.most_common()
        ]
    )

    unmatched = rule_summary.loc[~rule_summary["rule_recovered"]].copy()

    output_dir.mkdir(parents=True, exist_ok=True)
    recovery.to_csv(output_dir / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv", index=False)
    provider_summary.to_csv(
        output_dir / "VALIDATED_LOW_SUPPORT_FAMILY_SUMMARY.csv", index=False
    )
    rule_summary.to_csv(
        output_dir / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv", index=False
    )
    unmatched.to_csv(output_dir / "VALIDATED_LOW_UNMATCHED_RULES.csv", index=False)

    recovered_rules = int(rule_summary["rule_recovered"].sum())
    summary: dict[str, object] = {
        "contract": "chapter1_database_v1_validated_low_tiered_rule_provenance_audit_v3",
        "final_low_cells": int(len(low)),
        "derived_low_target_cells": int(len(target)),
        "distinct_derived_lineages": int(len(rule_summary)),
        "recovered_rule_lineages": recovered_rules,
        "rule_recovery_rate": recovered_rules / len(rule_summary),
        "rule_tier_counts": dict(tier_counts),
        "wave53_unmatched_reason_counts_after_tiered_recovery": dict(mismatch_reasons),
        "unresolved_rule_cells": unresolved_rule_cells,
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
        f"- recovered from verified rule generations: **{recovered_rules:,} / {len(rule_summary):,} ({recovered_rules / len(rule_summary):.1%})**",
        f"- Wave53 current rules: **{tier_counts['wave53_current_min3']:,}**",
        f"- earlier Run329 formal rules: **{tier_counts['run329_formal_current_min3']:,}**",
        f"- unresolved historical/secondary rules: **{tier_counts['unresolved_historical_secondary']:,}**",
        f"- recovered rules with direct-species count mismatch: **{rules_with_direct_count_mismatch:,}**",
        f"- recovered rules with zero direct support lineages: **{rules_with_zero_support:,}**",
        f"- target cells with unresolved rule provenance: **{unresolved_rule_cells:,}**",
        f"- target cells with zero recovered support: **{zero_support_cells:,}**",
        f"- target cells supported by >1 source family: **{multi_family_cells:,}**",
        f"- distinct upstream source families: **{len(provider_summary):,}**",
        "",
        "## Remaining unresolved rules by Wave53 status",
        "",
    ]
    for reason, count in mismatch_reasons.most_common():
        lines.append(f"- `{reason}` — {count:,} rules")
    lines += ["", "## Largest upstream source families", ""]
    for _, row in provider_summary.head(30).iterrows():
        lines.append(
            f"- `{row['source_family']}` — "
            f"{int(row['validated_low_cell_mentions']):,} Low-cell mentions; "
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
    parser.add_argument("--wave53-external-congener", type=Path, required=True)
    parser.add_argument("--wave53-rules", type=Path, required=True)
    parser.add_argument("--run329-direct", type=Path)
    parser.add_argument("--run329-rules", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--lineage-family-map", type=Path)
    args = parser.parse_args()
    recover(
        args.species_axis,
        args.wave53_direct,
        args.wave53_external_congener,
        args.wave53_rules,
        args.output_dir,
        lineage_family_map_path=args.lineage_family_map,
        run329_direct_path=args.run329_direct,
        run329_rules_path=args.run329_rules,
    )


if __name__ == "__main__":
    main()
