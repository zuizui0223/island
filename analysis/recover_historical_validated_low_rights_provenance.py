"""Recover historical Database 1.0 Validated-Low rule provenance.

Rights/provenance audit only. The immutable scientific species-axis ledger is never
changed and no redistribution right is granted. This extends the Wave53/Run329
first pass by examining only the still-unresolved rule keys against pinned historical
artifacts. Reads of large direct ledgers are restricted to the handful of required
genus × trait keys so the audit remains CI-sized.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from island_v2.chapter1_database_release import family_for_lineage, load_lineage_family_map

KEYS = ["genus", "trait_name", "inferred_state_set"]
_TRUE = {"1", "true", "yes", "y"}


def _truthy(value: object) -> bool:
    return str(value or "").strip().casefold() in _TRUE


def _tokens(value: object) -> list[str]:
    return [p.strip() for p in str(value or "").split("|") if p.strip()]


def _json_list(value: object) -> list[str]:
    out = json.loads(str(value or "[]"))
    if not isinstance(out, list):
        raise ValueError(f"expected JSON list: {value}")
    return [str(v) for v in out]


def _canonical_state_set(value: object) -> str:
    vals = _json_list(value)
    if not vals:
        raise ValueError("state set must be non-empty")
    return json.dumps(sorted(set(vals)), separators=(",", ":"))


def _load_rules(path: Path, state_col: str = "inferred_state_set") -> pd.DataFrame:
    rules = pd.read_csv(path, dtype=str).fillna("")
    required = {"genus", "trait_name", state_col, "n_direct_species", "eligible"}
    missing = required.difference(rules.columns)
    if missing:
        raise ValueError(f"rule audit {path} missing columns: {sorted(missing)}")
    if "setting" in rules.columns:
        rules = rules.loc[rules["setting"].eq("current_min3")].copy()
    rules = rules.loc[rules["eligible"].map(_truthy)].copy()
    if "diagnostic_only" in rules.columns:
        rules = rules.loc[~rules["diagnostic_only"].map(_truthy)].copy()
    rules["inferred_state_set"] = rules[state_col].map(_canonical_state_set)
    return rules


def _targeted_support_index(
    paths: list[Path],
    genus_traits: set[tuple[str, str]],
) -> dict[tuple[str, str], tuple[set[str], set[str]]]:
    species: dict[tuple[str, str], set[str]] = defaultdict(set)
    lineages: dict[tuple[str, str], set[str]] = defaultdict(set)
    if not genus_traits:
        return {}
    target_genera = {g for g, _ in genus_traits}
    target_traits = {t for _, t in genus_traits}
    for path in paths:
        header = pd.read_csv(path, nrows=0)
        required = {"accepted_species", "trait_name", "source_lineages"}
        missing = required.difference(header.columns)
        if missing:
            raise ValueError(f"direct ledger {path} missing columns: {sorted(missing)}")
        for chunk in pd.read_csv(
            path,
            dtype=str,
            usecols=["accepted_species", "trait_name", "source_lineages"],
            chunksize=40000,
        ):
            chunk = chunk.fillna("")
            chunk["genus"] = chunk["accepted_species"].str.split().str[0].fillna("")
            chunk = chunk.loc[
                chunk["genus"].isin(target_genera)
                & chunk["trait_name"].isin(target_traits)
            ]
            for row in chunk.itertuples(index=False):
                key = (row.genus, row.trait_name)
                if key not in genus_traits:
                    continue
                species[key].add(row.accepted_species)
                lineages[key].update(_tokens(row.source_lineages))
    return {key: (species[key], lineages[key]) for key in genus_traits}


def _targeted_coverage_support(
    path: Path,
    genus_traits: set[tuple[str, str]],
) -> dict[tuple[str, str], tuple[set[str], set[str]]]:
    """Narrow bridge for historical High/Medium single-trait axis cells."""
    species: dict[tuple[str, str], set[str]] = defaultdict(set)
    lineages: dict[tuple[str, str], set[str]] = defaultdict(set)
    if not genus_traits:
        return {}
    genera = {g for g, _ in genus_traits}
    usecols = ["accepted_species", "trait_names", "source_lineages", "quality"]
    for chunk in pd.read_csv(path, dtype=str, usecols=usecols, chunksize=40000):
        chunk = chunk.fillna("")
        chunk["genus"] = chunk["accepted_species"].str.split().str[0].fillna("")
        chunk = chunk.loc[chunk["genus"].isin(genera)]
        for row in chunk.itertuples(index=False):
            if (
                row.quality.casefold() not in {"high", "medium"}
                or "validated-low:" in row.source_lineages
            ):
                continue
            raw = row.trait_names.strip()
            if not raw:
                continue
            try:
                names = _json_list(raw) if raw.startswith("[") else _tokens(raw)
            except json.JSONDecodeError:
                names = _tokens(raw)
            if len(names) != 1:
                continue
            key = (row.genus, names[0])
            if key not in genus_traits:
                continue
            species[key].add(row.accepted_species)
            lineages[key].update(_tokens(row.source_lineages))
    return {key: (species[key], lineages[key]) for key in genus_traits}


def _merge_support_indexes(*indexes):
    species: dict[tuple[str, str], set[str]] = defaultdict(set)
    lineages: dict[tuple[str, str], set[str]] = defaultdict(set)
    for index in indexes:
        for key, (sp, src) in index.items():
            species[key].update(sp)
            lineages[key].update(src)
    return {
        key: (species[key], lineages[key])
        for key in set(species).union(lineages)
    }


def _frontier_receipts(path: Path, target_keys: set[tuple[str, str, str]]):
    rules = _load_rules(path, state_col="predicted_state_set")
    needed = {"support_species", "support_source_lineages"}
    if missing := needed.difference(rules.columns):
        raise ValueError(f"Wave34 frontier missing columns: {sorted(missing)}")
    out = {}
    for row in rules.to_dict("records"):
        key = tuple(str(row[k]) for k in KEYS)
        if key not in target_keys:
            continue
        sp = set(_json_list(row["support_species"]))
        src = set(_json_list(row["support_source_lineages"]))
        expected = int(row["n_direct_species"])
        if len(sp) == expected and src:
            out[key] = ("wave34_frontier", expected, len(sp), src)
    return out


def _rule_receipts(
    path: Path,
    direct_paths: list[Path],
    tier: str,
    target_keys: set[tuple[str, str, str]],
    coverage_path: Path | None = None,
):
    rules = _load_rules(path)
    rules = rules.loc[
        rules.apply(
            lambda row: (
                row["genus"],
                row["trait_name"],
                row["inferred_state_set"],
            ) in target_keys,
            axis=1,
        )
    ]
    genus_traits = {(row.genus, row.trait_name) for row in rules.itertuples(index=False)}
    support = _targeted_support_index(direct_paths, genus_traits)
    if coverage_path is not None:
        support = _merge_support_indexes(
            support,
            _targeted_coverage_support(coverage_path, genus_traits),
        )
    out = {}
    for row in rules.to_dict("records"):
        key = tuple(str(row[k]) for k in KEYS)
        sp, src = support.get((row["genus"], row["trait_name"]), (set(), set()))
        expected = int(row["n_direct_species"])
        if len(sp) == expected and src:
            out[key] = (tier, expected, len(sp), src)
    return out


def _secondary_mismatch(row: dict[str, str], sidecar: pd.DataFrame):
    candidates = sidecar.loc[sidecar["accepted_species"].eq(row["genus"])]
    for item in candidates.to_dict("records"):
        try:
            traits = _json_list(item.get("trait_names", "[]"))
            nested = _json_list(item.get("predicted_state_sets", "[]"))
        except (ValueError, json.JSONDecodeError):
            continue
        if traits != [row["trait_name"]] or len(nested) != 1:
            continue
        historical_state = _canonical_state_set(nested[0])
        if historical_state == row["inferred_state_set"]:
            continue
        sources = set(_json_list(item.get("support_source_lineages", "[]")))
        expected = int(str(item.get("min_direct_species_support", "0") or "0"))
        if sources and expected > 0:
            return (
                "wave33_secondary_state_mismatch",
                expected,
                expected,
                sources,
                historical_state,
            )
    return None


def extend(
    rights_dir: Path,
    lineage_family_map_path: Path | None,
    run329_direct: Path,
    wave34_frontier: Path,
    wave33_secondary_low: Path,
    wave35_rules: Path,
    wave35_coverage: Path,
    wave35_direct: Path,
    wave39_rules: Path,
    wave39_direct: Path,
    wave40_rules: Path,
    wave40_direct: Path,
    wave41_rules: Path,
    wave41_direct: Path,
    wave41_external: Path,
    wave42_rules: Path,
    wave42_direct: Path,
    wave42_external: Path,
    wave43_rules: Path,
    wave43_direct: Path,
    wave43_external: Path,
) -> dict[str, object]:
    rules_path = rights_dir / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv"
    cells_path = rights_dir / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv"
    rules = pd.read_csv(rules_path, dtype=str).fillna("")
    cells = pd.read_csv(cells_path, dtype=str).fillna("")
    rules["inferred_state_set"] = rules["inferred_state_set"].map(_canonical_state_set)
    rules["source_provenance_recovered"] = (
        rules["rule_recovered"].map(_truthy).astype(object)
    )
    rules["semantic_status"] = rules["rule_recovered"].map(
        lambda value: "exact_rule" if _truthy(value) else ""
    )
    rules["historical_state_set"] = ""

    unresolved_idx = rules.index[~rules["rule_recovered"].map(_truthy)].tolist()
    target_keys = {tuple(rules.loc[i, k] for k in KEYS) for i in unresolved_idx}
    receipts = {}

    def add(batch):
        receipts.update({key: value for key, value in batch.items() if key not in receipts})

    add(_frontier_receipts(wave34_frontier, target_keys))
    remaining = target_keys - set(receipts)
    add(
        _rule_receipts(
            wave35_rules,
            [run329_direct, wave35_direct],
            "wave35",
            remaining,
            wave35_coverage,
        )
    )
    remaining = target_keys - set(receipts)
    add(_rule_receipts(wave39_rules, [run329_direct, wave39_direct], "wave39", remaining))
    remaining = target_keys - set(receipts)
    add(_rule_receipts(wave40_rules, [wave40_direct], "wave40", remaining))
    remaining = target_keys - set(receipts)
    add(
        _rule_receipts(
            wave41_rules,
            [wave41_direct, wave41_external],
            "wave41",
            remaining,
        )
    )
    remaining = target_keys - set(receipts)
    add(
        _rule_receipts(
            wave42_rules,
            [wave42_direct, wave42_external],
            "wave42",
            remaining,
        )
    )
    remaining = target_keys - set(receipts)
    add(
        _rule_receipts(
            wave43_rules,
            [wave43_direct, wave43_external],
            "wave43",
            remaining,
        )
    )

    overrides = load_lineage_family_map(lineage_family_map_path)
    sidecar = pd.read_csv(wave33_secondary_low, dtype=str).fillna("")
    exact_added = 0
    mismatch_added = 0
    for idx in unresolved_idx:
        row = rules.loc[idx].to_dict()
        key = tuple(row[k] for k in KEYS)
        receipt = receipts.get(key)
        historical_state = ""
        semantic_status = "exact_rule"
        exact = receipt is not None
        if receipt is None:
            mismatch = _secondary_mismatch(row, sidecar)
            if mismatch is None:
                continue
            tier, expected, reconstructed, sources, historical_state = mismatch
            semantic_status = "historical_state_set_mismatch"
            mismatch_added += 1
        else:
            tier, expected, reconstructed, sources = receipt
            exact_added += 1
        families = {family_for_lineage(source, overrides) for source in sources}
        rules.at[idx, "provenance_rule_tier"] = tier
        rules.at[idx, "rule_recovered"] = str(exact)
        rules.at[idx, "source_provenance_recovered"] = True
        rules.at[idx, "semantic_status"] = semantic_status
        rules.at[idx, "historical_state_set"] = historical_state
        rules.at[idx, "expected_n_direct_species"] = str(expected)
        rules.at[idx, "reconstructed_n_direct_species"] = str(reconstructed)
        rules.at[idx, "direct_species_count_matches"] = "True"
        rules.at[idx, "support_source_lineage_count"] = str(len(sources))
        rules.at[idx, "support_source_lineages"] = json.dumps(
            sorted(sources), separators=(",", ":")
        )
        rules.at[idx, "support_source_family_count"] = str(len(families))
        rules.at[idx, "support_source_families"] = "|".join(sorted(families))

    exact_mask = rules["rule_recovered"].map(_truthy)
    source_mask = rules["source_provenance_recovered"].map(_truthy)
    support_by_rule = {
        row["derived_lineage"]: set(_json_list(row["support_source_lineages"]))
        for row in rules.to_dict("records")
        if _truthy(row["source_provenance_recovered"])
    }
    families_by_rule = {
        lineage: {family_for_lineage(source, overrides) for source in sources}
        for lineage, sources in support_by_rule.items()
    }
    exact_rules = set(rules.loc[exact_mask, "derived_lineage"])
    source_rules = set(rules.loc[source_mask, "derived_lineage"])
    tier_by_rule = dict(
        zip(rules["derived_lineage"], rules["provenance_rule_tier"], strict=True)
    )

    family_cell_counts: Counter[str] = Counter()
    family_support_lineages: dict[str, set[str]] = defaultdict(set)
    rows = []
    unresolved_exact = 0
    unresolved_source = 0
    zero_support = 0
    multi_family = 0
    for row in cells.to_dict("records"):
        derived = _tokens(row["derived_lineages"])
        sources: set[str] = set()
        families: set[str] = set()
        tiers: set[str] = set()
        for lineage in derived:
            sources.update(support_by_rule.get(lineage, set()))
            families.update(families_by_rule.get(lineage, set()))
            tier = tier_by_rule.get(lineage, "")
            if tier:
                tiers.add(tier)
        exact_ok = all(lineage in exact_rules for lineage in derived)
        source_ok = all(lineage in source_rules for lineage in derived)
        unresolved_exact += int(not exact_ok)
        unresolved_source += int(not source_ok)
        zero_support += int(not sources)
        multi_family += int(len(families) > 1)
        family_cell_counts.update(families)
        for source in sources:
            family_support_lineages[family_for_lineage(source, overrides)].add(source)
        rows.append(
            (
                str(exact_ok),
                str(source_ok),
                "|".join(sorted(tiers)),
                str(len(sources)),
                json.dumps(sorted(sources), separators=(",", ":")),
                str(len(families)),
                "|".join(sorted(families)),
            )
        )

    columns = list(zip(*rows, strict=True))
    names = [
        "all_rules_exact_rule_recovered",
        "all_rules_source_provenance_recovered",
        "provenance_rule_tiers",
        "support_source_lineage_count",
        "support_source_lineages",
        "support_source_family_count",
        "support_source_families",
    ]
    for name, values in zip(names, columns, strict=True):
        cells[name] = values
    cells["all_rules_provenance_recovered"] = cells["all_rules_exact_rule_recovered"]

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
    unmatched = rules.loc[~exact_mask].copy()
    mismatch = rules.loc[
        rules["semantic_status"].eq("historical_state_set_mismatch")
    ].copy()

    rules.to_csv(rules_path, index=False, lineterminator="\n")
    cells.to_csv(cells_path, index=False, lineterminator="\n")
    provider.to_csv(
        rights_dir / "VALIDATED_LOW_SUPPORT_FAMILY_SUMMARY.csv",
        index=False,
        lineterminator="\n",
    )
    unmatched.to_csv(
        rights_dir / "VALIDATED_LOW_UNMATCHED_RULES.csv",
        index=False,
        lineterminator="\n",
    )
    mismatch.to_csv(
        rights_dir / "VALIDATED_LOW_SEMANTIC_MISMATCHES.csv",
        index=False,
        lineterminator="\n",
    )

    summary = {
        "contract": "chapter1_database_v1_validated_low_historical_rule_provenance_audit_v4",
        "distinct_derived_lineages": len(rules),
        "exact_rule_generations_recovered": int(exact_mask.sum()),
        "source_provenance_recovered_rule_lineages": int(source_mask.sum()),
        "new_historical_exact_rule_recoveries": exact_added,
        "historical_state_set_mismatch_rules": mismatch_added,
        "unresolved_exact_rule_cells": unresolved_exact,
        "unresolved_source_provenance_cells": unresolved_source,
        "zero_support_cells": zero_support,
        "multi_family_support_cells": multi_family,
        "distinct_support_source_families": len(provider),
        "derived_low_target_cells": len(cells),
        "rule_tier_counts": dict(Counter(rules["provenance_rule_tier"])),
    }
    (rights_dir / "VALIDATED_LOW_PROVENANCE_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    lines = [
        "# Validated-Low upstream provenance recovery",
        "",
        f"- exact rule semantics: **{summary['exact_rule_generations_recovered']:,}/{len(rules):,}**",
        f"- upstream source provenance: **{summary['source_provenance_recovered_rule_lineages']:,}/{len(rules):,}**",
        f"- unresolved exact-rule cells: **{unresolved_exact:,}**",
        f"- unresolved source-provenance cells: **{unresolved_source:,}**",
        f"- semantic mismatch rules: **{mismatch_added:,}**",
    ]
    (rights_dir / "VALIDATED_LOW_PROVENANCE_SUMMARY.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    for arg in [
        "rights-dir",
        "run329-direct",
        "wave34-frontier",
        "wave33-secondary-low",
        "wave35-rules",
        "wave35-coverage",
        "wave35-direct",
        "wave39-rules",
        "wave39-direct",
        "wave40-rules",
        "wave40-direct",
        "wave41-rules",
        "wave41-direct",
        "wave41-external",
        "wave42-rules",
        "wave42-direct",
        "wave42-external",
        "wave43-rules",
        "wave43-direct",
        "wave43-external",
    ]:
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--lineage-family-map", type=Path)
    args = parser.parse_args()
    kwargs = vars(args)
    kwargs["lineage_family_map_path"] = kwargs.pop("lineage_family_map")
    print(json.dumps(extend(**kwargs), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
