"""Complete Database 1.0 Validated-Low provenance from immutable historical rules.

Rights/provenance audit only. This script never modifies the frozen scientific
Database 1.0 ledger and never grants redistribution rights.

The current audit first reconstructs rules from Wave53 and Run329. A small
append-only historical Low layer predates those exact rule states. This module
recovers that layer from the immutable Wave34/35/39/40/41/42/43 artifacts that
actually generated or carried the rules. Exact state-set matches are required
for scientific rule recovery. When the upstream support is recoverable but the
historical state set differs from the stored Database 1.0 lineage, source
provenance is retained while the semantic mismatch remains explicit.
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


def tokens(value: object) -> list[str]:
    return [part.strip() for part in str(value or "").split("|") if part.strip()]


def canonical_state(value: object) -> str:
    text = str(value or "").strip()
    decoded = json.loads(text)
    if not isinstance(decoded, list) or not decoded:
        raise ValueError(f"invalid state set: {text}")
    return json.dumps(sorted({str(item) for item in decoded}), separators=(",", ":"))


def truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in {"1", "true", "yes", "y"}


def parse_json_list(value: object) -> list[str]:
    if isinstance(value, list):
        return [str(item) for item in value]
    text = str(value or "").strip()
    if not text:
        return []
    decoded = json.loads(text)
    if not isinstance(decoded, list):
        raise ValueError(f"expected list JSON: {text}")
    return [str(item) for item in decoded]


def find_one(root: Path, patterns: list[str]) -> Path:
    hits: list[Path] = []
    for pattern in patterns:
        hits.extend(root.rglob(pattern))
    unique = sorted(set(hits))
    if len(unique) != 1:
        raise ValueError(
            f"expected exactly one match under {root} for {patterns}, got {unique}"
        )
    return unique[0]


def load_rule_table(path: Path, state_col: str) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    required = {"genus", "trait_name", state_col, "n_direct_species", "eligible"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"{path} missing {sorted(missing)}")
    if "setting" in frame.columns:
        frame = frame.loc[frame["setting"].eq("current_min3")].copy()
    if "diagnostic_only" in frame.columns:
        frame = frame.loc[~frame["diagnostic_only"].map(truthy)].copy()
    frame = frame.loc[frame["eligible"].map(truthy)].copy()
    frame["_state"] = frame[state_col].map(canonical_state)
    return frame


def conflict_index(paths: list[Path]) -> dict[tuple[str, str], dict[str, object]]:
    frames: list[pd.DataFrame] = []
    required = {"accepted_species", "trait_name", "resolution_status", "source_lineages"}
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path} missing {sorted(missing)}")
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True, sort=False).fillna("")
    combined = combined.loc[combined["resolution_status"].eq("resolved")].copy()
    combined["genus"] = combined["accepted_species"].str.split().str[0].fillna("")
    index: dict[tuple[str, str], dict[str, object]] = {}
    for key, group in combined.groupby(["genus", "trait_name"], sort=False):
        sources: set[str] = set()
        for value in group["source_lineages"]:
            sources.update(tokens(value))
        index[key] = {
            "n_species": int(group["accepted_species"].nunique()),
            "sources": sources,
        }
    return index


def discover_source(
    label: str, root: Path
) -> tuple[pd.DataFrame, dict[tuple[str, str], dict[str, object]], str]:
    if label == "wave34":
        rule_path = find_one(root, ["trait_specific_rule_frontier.csv.gz"])
        rules = load_rule_table(rule_path, "predicted_state_set")
        index: dict[tuple[str, str], dict[str, object]] = {}
        for _, row in rules.iterrows():
            sources = set(parse_json_list(row.get("support_source_lineages", "[]")))
            index[(row["genus"], row["trait_name"])] = {
                "n_species": int(row["n_direct_species"]),
                "sources": sources,
            }
        return rules, index, "wave34_frontier"

    if label == "wave35":
        rule_path = find_one(root, ["wave35_provider_touched_new_rule_audit.csv.gz"])
        conflicts = find_one(root, ["wave35_source_lineage_conflicts.csv.gz"])
        return (
            load_rule_table(rule_path, "inferred_state_set"),
            conflict_index([conflicts]),
            "wave35",
        )

    if label == "wave39":
        rule_path = find_one(root, ["wave39_new_trait_specific_genus_rule_audit.csv.gz"])
        conflicts = find_one(root, ["wave39_touched_source_lineage_conflicts.csv.gz"])
        return (
            load_rule_table(rule_path, "inferred_state_set"),
            conflict_index([conflicts]),
            "wave39",
        )

    if label == "wave40":
        rule_path = find_one(root, ["wave40_new_trait_specific_genus_rule_audit.csv.gz"])
        conflicts = find_one(root, ["wave40-all-evidence-audit/source_lineage_conflicts.csv.gz"])
        return (
            load_rule_table(rule_path, "inferred_state_set"),
            conflict_index([conflicts]),
            "wave40",
        )

    if label in {"wave41", "wave42", "wave43"}:
        rule_path = find_one(root, [f"{label}_new_trait_specific_genus_rule_audit.csv.gz"])
        ordinary = find_one(
            root, [f"{label}-all-evidence-audit/source_lineage_conflicts.csv.gz"]
        )
        external = find_one(
            root,
            [f"{label}-all-evidence-audit/external_congener_source_lineage_conflicts.csv.gz"],
        )
        return (
            load_rule_table(rule_path, "inferred_state_set"),
            conflict_index([ordinary, external]),
            label,
        )

    raise ValueError(f"unsupported historical generation: {label}")


def apply(
    audit_dir: Path,
    source_specs: list[tuple[str, Path]],
    lineage_family_map_path: Path | None = None,
) -> dict[str, object]:
    rule_path = audit_dir / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv"
    cell_path = audit_dir / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv"
    rules = pd.read_csv(rule_path, dtype=str).fillna("")
    cells = pd.read_csv(cell_path, dtype=str).fillna("")
    overrides = load_lineage_family_map(lineage_family_map_path)

    for column in (
        "source_provenance_recovered",
        "exact_semantic_match",
        "semantic_mismatch_receipt_state_set",
    ):
        if column not in rules.columns:
            rules[column] = ""
    already_recovered = rules["rule_recovered"].map(truthy)
    rules.loc[already_recovered, "source_provenance_recovered"] = "true"
    rules.loc[already_recovered, "exact_semantic_match"] = "true"

    sources = []
    for label, root in source_specs:
        rule_table, index, tier = discover_source(label, root)
        exact_lookup = {
            (row["genus"], row["trait_name"], row["_state"]): row
            for _, row in rule_table.iterrows()
        }
        genus_trait_lookup: dict[tuple[str, str], list[pd.Series]] = defaultdict(list)
        for _, row in rule_table.iterrows():
            genus_trait_lookup[(row["genus"], row["trait_name"])].append(row)
        sources.append((label, tier, exact_lookup, genus_trait_lookup, index))

    receipts: list[dict[str, object]] = []
    unresolved = rules.loc[~rules["rule_recovered"].map(truthy)].copy()
    for row_index, row in unresolved.iterrows():
        key = (row["genus"], row["trait_name"], canonical_state(row["inferred_state_set"]))
        recovered = False
        for label, tier, exact_lookup, _, index in sources:
            if key not in exact_lookup:
                continue
            historical_rule = exact_lookup[key]
            support = index.get(
                (row["genus"], row["trait_name"]), {"n_species": 0, "sources": set()}
            )
            if int(historical_rule["n_direct_species"]) != int(support["n_species"]):
                raise ValueError(
                    f"{key} {label} direct count mismatch "
                    f"{historical_rule['n_direct_species']} != {support['n_species']}"
                )
            if not support["sources"]:
                raise ValueError(f"{key} {label} has zero support lineages")
            families = sorted(
                {family_for_lineage(source, overrides) for source in support["sources"]}
            )
            rules.at[row_index, "provenance_rule_tier"] = tier
            rules.at[row_index, "rule_recovered"] = "true"
            rules.at[row_index, "source_provenance_recovered"] = "true"
            rules.at[row_index, "exact_semantic_match"] = "true"
            rules.at[row_index, "expected_n_direct_species"] = str(
                historical_rule["n_direct_species"]
            )
            rules.at[row_index, "reconstructed_n_direct_species"] = str(
                support["n_species"]
            )
            rules.at[row_index, "direct_species_count_matches"] = "true"
            rules.at[row_index, "support_source_lineage_count"] = str(
                len(support["sources"])
            )
            rules.at[row_index, "support_source_lineages"] = json.dumps(
                sorted(support["sources"]), separators=(",", ":")
            )
            rules.at[row_index, "support_source_family_count"] = str(len(families))
            rules.at[row_index, "support_source_families"] = "|".join(families)
            receipts.append(
                {
                    "derived_lineage": row["derived_lineage"],
                    "receipt_source": label,
                    "exact_semantic_match": True,
                }
            )
            recovered = True
            break
        if recovered:
            continue

        candidates = []
        for label, tier, _, genus_trait_lookup, index in sources:
            for historical_rule in genus_trait_lookup.get(
                (row["genus"], row["trait_name"]), []
            ):
                support = index.get(
                    (row["genus"], row["trait_name"]),
                    {"n_species": 0, "sources": set()},
                )
                if (
                    int(historical_rule["n_direct_species"]) == int(support["n_species"])
                    and support["sources"]
                ):
                    candidates.append((label, tier, historical_rule, support))
        if len(candidates) == 1:
            label, tier, historical_rule, support = candidates[0]
            families = sorted(
                {family_for_lineage(source, overrides) for source in support["sources"]}
            )
            rules.at[row_index, "provenance_rule_tier"] = f"{tier}_semantic_mismatch"
            rules.at[row_index, "source_provenance_recovered"] = "true"
            rules.at[row_index, "exact_semantic_match"] = "false"
            rules.at[row_index, "semantic_mismatch_receipt_state_set"] = historical_rule[
                "_state"
            ]
            rules.at[row_index, "expected_n_direct_species"] = str(
                historical_rule["n_direct_species"]
            )
            rules.at[row_index, "reconstructed_n_direct_species"] = str(
                support["n_species"]
            )
            rules.at[row_index, "direct_species_count_matches"] = "true"
            rules.at[row_index, "support_source_lineage_count"] = str(
                len(support["sources"])
            )
            rules.at[row_index, "support_source_lineages"] = json.dumps(
                sorted(support["sources"]), separators=(",", ":")
            )
            rules.at[row_index, "support_source_family_count"] = str(len(families))
            rules.at[row_index, "support_source_families"] = "|".join(families)
            receipts.append(
                {
                    "derived_lineage": row["derived_lineage"],
                    "receipt_source": label,
                    "exact_semantic_match": False,
                }
            )

    support_by: dict[str, set[str]] = {}
    family_by: dict[str, set[str]] = {}
    exact_by: dict[str, bool] = {}
    source_by: dict[str, bool] = {}
    tier_by: dict[str, str] = {}
    for _, row in rules.iterrows():
        lineage = row["derived_lineage"]
        supports = set(parse_json_list(row["support_source_lineages"]))
        support_by[lineage] = supports
        family_by[lineage] = {
            family_for_lineage(source, overrides) for source in supports
        }
        exact_by[lineage] = truthy(row["rule_recovered"])
        source_by[lineage] = truthy(row["source_provenance_recovered"]) or truthy(
            row["rule_recovered"]
        )
        tier_by[lineage] = row["provenance_rule_tier"]

    aggregates = []
    family_cell_counts: Counter[str] = Counter()
    family_support_lineages: dict[str, set[str]] = defaultdict(set)
    for _, row in cells.iterrows():
        derived = tokens(row["derived_lineages"])
        supports: set[str] = set()
        families: set[str] = set()
        tiers: set[str] = set()
        for lineage in derived:
            supports.update(support_by.get(lineage, set()))
            families.update(family_by.get(lineage, set()))
            if tier_by.get(lineage):
                tiers.add(tier_by[lineage])
        family_cell_counts.update(families)
        for source in supports:
            family_support_lineages[family_for_lineage(source, overrides)].add(source)
        aggregates.append(
            (
                all(exact_by.get(lineage, False) for lineage in derived),
                all(source_by.get(lineage, False) for lineage in derived),
                "|".join(sorted(tiers)),
                len(supports),
                json.dumps(sorted(supports), separators=(",", ":")),
                len(families),
                "|".join(sorted(families)),
            )
        )

    columns = [
        "all_rules_provenance_recovered",
        "all_rules_source_provenance_recovered",
        "provenance_rule_tiers",
        "support_source_lineage_count",
        "support_source_lineages",
        "support_source_family_count",
        "support_source_families",
    ]
    for column_index, column in enumerate(columns):
        cells[column] = [value[column_index] for value in aggregates]

    provider = pd.DataFrame(
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
    unmatched = rules.loc[~rules["rule_recovered"].map(truthy)].copy()
    receipts_frame = pd.DataFrame(receipts)

    rules.to_csv(rule_path, index=False)
    cells.to_csv(cell_path, index=False)
    provider.to_csv(audit_dir / "VALIDATED_LOW_SUPPORT_FAMILY_SUMMARY.csv", index=False)
    unmatched.to_csv(audit_dir / "VALIDATED_LOW_UNMATCHED_RULES.csv", index=False)
    receipts_frame.to_csv(audit_dir / "VALIDATED_LOW_HISTORICAL_RECEIPTS.csv", index=False)

    exact = int(rules["rule_recovered"].map(truthy).sum())
    source_recovered = int(
        (
            rules["source_provenance_recovered"].map(truthy)
            | rules["rule_recovered"].map(truthy)
        ).sum()
    )
    semantic_mismatches = int(
        (
            (
                rules["source_provenance_recovered"].map(truthy)
                | rules["rule_recovered"].map(truthy)
            )
            & ~rules["rule_recovered"].map(truthy)
        ).sum()
    )
    unresolved_source_cells = int(
        (~cells["all_rules_source_provenance_recovered"].map(truthy)).sum()
    )
    zero_support_cells = int(
        (cells["support_source_lineage_count"].astype(int) == 0).sum()
    )
    summary: dict[str, object] = {
        "contract": "chapter1_database_v1_validated_low_historical_rule_provenance_audit_v4",
        "distinct_derived_lineages": int(len(rules)),
        "exact_rule_semantics_recovered": exact,
        "exact_rule_semantics_recovery_rate": exact / len(rules),
        "source_provenance_recovered": source_recovered,
        "source_provenance_recovery_rate": source_recovered / len(rules),
        "semantic_mismatch_rule_lineages": semantic_mismatches,
        "target_cells_with_unresolved_source_provenance": unresolved_source_cells,
        "target_cells_with_zero_support": zero_support_cells,
        "scientific_database_modified": False,
        "rights_granted": False,
    }
    (audit_dir / "VALIDATED_LOW_HISTORICAL_PROVENANCE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "# Validated-Low historical provenance completion",
        "",
        f"- exact rule semantics recovered: **{exact:,} / {len(rules):,} ({exact / len(rules):.2%})**",
        f"- upstream source provenance recovered: **{source_recovered:,} / {len(rules):,} ({source_recovered / len(rules):.2%})**",
        f"- semantic mismatch rules retained fail-closed: **{semantic_mismatches:,}**",
        f"- target cells with unresolved source provenance: **{unresolved_source_cells:,}**",
        f"- target cells with zero support lineage: **{zero_support_cells:,}**",
        "",
        "The historical completion does not modify Database 1.0 or grant redistribution rights.",
    ]
    (audit_dir / "VALIDATED_LOW_HISTORICAL_PROVENANCE_SUMMARY.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--lineage-family-map", type=Path)
    for label in ("wave34", "wave35", "wave39", "wave40", "wave41", "wave42", "wave43"):
        parser.add_argument(f"--{label}-dir", type=Path, required=True)
    args = parser.parse_args()
    source_specs = [
        (label, getattr(args, f"{label}_dir"))
        for label in ("wave34", "wave35", "wave39", "wave40", "wave41", "wave42", "wave43")
    ]
    print(
        json.dumps(
            apply(args.audit_dir, source_specs, args.lineage_family_map), sort_keys=True
        )
    )


if __name__ == "__main__":
    main()
