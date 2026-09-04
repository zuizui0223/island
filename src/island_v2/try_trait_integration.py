"""Integrate prepared private TRY evidence into the current strict v2 ledger.

This script deliberately starts from ``try_common_direct_evidence.csv.gz``
created by :mod:`island_v2.try_traits`; it never needs the raw TRY request.
Only genus x trait cells touched by TRY are recomputed, including their
Validated-Low rules, while unaffected direct/low cells are preserved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import pandas as pd

from island_v2.all_evidence_trait_audit import (
    apply_genus_rules,
    build_rule_audit,
    coverage_snapshot,
    dedupe_direct_lineages,
    direct_evidence_from_integrated,
    load_ontology,
    resolve_direct_cells,
    species_axis_coverage,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

CONTRACT = "try_v6_private_species_direct_integration_v1"


def text(value: Any) -> str:
    return "" if value is None or pd.isna(value) else " ".join(str(value).split())


def keyed(frame: pd.DataFrame, columns: list[str]) -> set[tuple[str, ...]]:
    return set(frame[columns].fillna("").itertuples(index=False, name=None))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_gz(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression={"method": "gzip", "mtime": 0})


def load_try_common(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, dtype=str).fillna("")
    missing = set(EVIDENCE_COLUMNS).difference(frame.columns)
    if missing:
        raise ValueError(f"TRY common evidence missing columns: {sorted(missing)}")
    frame = frame.reindex(columns=EVIDENCE_COLUMNS).fillna("")
    if not frame["source_group"].eq("try").all():
        raise ValueError("TRY common evidence must use source_group=try")
    if not frame["evidence_scope"].isin({"species_direct", "synonym_direct"}).all():
        raise ValueError("TRY common evidence contains non-direct rows")
    if frame["source_lineage"].eq("").any():
        raise ValueError("TRY common evidence must preserve original-source lineage")
    return frame


def build_master(
    current_coverage: pd.DataFrame,
    master_genus_map: pd.DataFrame,
    current_low: pd.DataFrame,
) -> pd.DataFrame:
    genus = master_genus_map[["accepted_species", "genus"]].fillna("").copy()
    genus = genus.loc[genus.genus.ne("")].drop_duplicates("accepted_species", keep="first")
    low_genus = pd.DataFrame(columns=["accepted_species", "genus"])
    if {"accepted_species", "genus"}.issubset(current_low.columns):
        low_genus = (
            current_low.loc[current_low["genus"].fillna("").ne(""), ["accepted_species", "genus"]]
            .drop_duplicates("accepted_species", keep="first")
        )
    master = current_coverage[["accepted_species"]].drop_duplicates().copy()
    master = master.merge(
        genus.rename(columns={"genus": "formal_genus"}),
        on="accepted_species",
        how="left",
    )
    master = master.merge(
        low_genus.rename(columns={"genus": "low_genus"}),
        on="accepted_species",
        how="left",
    )
    master["genus"] = master["formal_genus"].fillna("")
    master["genus"] = master["genus"].where(master["genus"].ne(""), master["low_genus"].fillna(""))
    master["genus"] = master["genus"].where(
        master["genus"].ne(""), master["accepted_species"].str.split().str[0]
    )
    return master[["accepted_species", "genus"]]


def integrate(
    *,
    try_common: pd.DataFrame,
    formal_direct_evidence: pd.DataFrame,
    additional_common: list[pd.DataFrame],
    current_direct: pd.DataFrame,
    current_low: pd.DataFrame,
    current_coverage: pd.DataFrame,
    master_genus_map: pd.DataFrame,
    ontology: dict[str, set[str]],
) -> dict[str, Any]:
    if try_common.empty:
        raise ValueError("TRY common evidence is empty")
    master = build_master(current_coverage, master_genus_map, current_low)
    master_species = set(master["accepted_species"])
    input_try_rows = len(try_common)
    input_try_species = int(try_common["accepted_species"].nunique())
    in_scope = try_common["accepted_species"].isin(master_species)
    out_of_scope_try = try_common.loc[~in_scope].copy()
    try_common = try_common.loc[in_scope].copy()
    if try_common.empty:
        raise ValueError("TRY common evidence has no rows in the current strict coverage universe")
    raw = pd.concat(
        [formal_direct_evidence, *additional_common, try_common],
        ignore_index=True,
    ).fillna("")
    genus_map = dict(master.itertuples(index=False, name=None))

    lineages, lineage_duplicates = dedupe_direct_lineages(raw, ontology)
    lineages["genus"] = lineages["accepted_species"].map(genus_map).fillna(
        lineages["accepted_species"].str.split().str[0]
    )
    current_direct = current_direct.copy().fillna("")
    current_low = current_low.copy().fillna("")
    current_direct["genus"] = current_direct["accepted_species"].map(genus_map).fillna(
        current_direct["accepted_species"].str.split().str[0]
    )
    if "genus" not in current_low.columns:
        current_low["genus"] = current_low["accepted_species"].map(genus_map).fillna(
            current_low["accepted_species"].str.split().str[0]
        )

    affected = {
        (genus_map.get(sp, sp.split()[0]), trait)
        for sp, trait in try_common[["accepted_species", "trait_name"]].itertuples(
            index=False, name=None
        )
    }
    lineage_mask = [
        (genus, trait) in affected
        for genus, trait in zip(lineages["genus"], lineages["trait_name"], strict=True)
    ]
    affected_lineages = lineages.loc[lineage_mask].copy()
    recomputed_direct, conflict_audit = resolve_direct_cells(affected_lineages)
    recomputed_direct["genus"] = recomputed_direct["accepted_species"].map(genus_map).fillna(
        recomputed_direct["accepted_species"].str.split().str[0]
    )

    existing_direct_mask = [
        (genus, trait) in affected
        for genus, trait in zip(current_direct["genus"], current_direct["trait_name"], strict=True)
    ]
    updated_direct = pd.concat(
        [current_direct.loc[[not x for x in existing_direct_mask]], recomputed_direct],
        ignore_index=True,
    ).fillna("")
    updated_direct = updated_direct.drop_duplicates(["accepted_species", "trait_name"], keep="last")

    affected_direct_mask = [
        (genus, trait) in affected
        for genus, trait in zip(updated_direct["genus"], updated_direct["trait_name"], strict=True)
    ]
    affected_low_mask = [
        (genus, trait) in affected
        for genus, trait in zip(current_low["genus"], current_low["trait_name"], strict=True)
    ]
    rules = build_rule_audit(
        updated_direct.loc[affected_direct_mask],
        affected_lineages,
        current_low.loc[affected_low_mask],
    )
    replacement_low = apply_genus_rules(master, updated_direct, rules, "current_min3")
    updated_low = pd.concat(
        [current_low.loc[[not x for x in affected_low_mask]], replacement_low],
        ignore_index=True,
    ).fillna("")
    updated_low = updated_low.drop_duplicates(["accepted_species", "trait_name"], keep="last")
    direct_keys = pd.MultiIndex.from_frame(updated_direct[["accepted_species", "trait_name"]])
    low_keys = pd.MultiIndex.from_frame(updated_low[["accepted_species", "trait_name"]])
    updated_low = updated_low.loc[~low_keys.isin(direct_keys)].reset_index(drop=True)
    updated_coverage = species_axis_coverage(master, updated_direct, updated_low)

    before_direct = keyed(current_direct, ["accepted_species", "trait_name"])
    after_direct = keyed(updated_direct, ["accepted_species", "trait_name"])
    before_low = keyed(current_low, ["accepted_species", "trait_name"])
    after_low = keyed(updated_low, ["accepted_species", "trait_name"])
    before_axis = keyed(
        current_coverage.loc[current_coverage["quality"].fillna("").ne("")],
        ["accepted_species", "axis"],
    )
    after_axis = keyed(
        updated_coverage.loc[updated_coverage["quality"].fillna("").ne("")],
        ["accepted_species", "axis"],
    )
    summary = {
        "contract": CONTRACT,
        "try": {
            "input_common_rows": input_try_rows,
            "input_species": input_try_species,
            "common_rows": len(try_common),
            "species": int(try_common.accepted_species.nunique()),
            "species_trait": int(
                try_common[["accepted_species", "trait_name"]]
                .drop_duplicates()
                .shape[0]
            ),
            "out_of_scope_rows_dropped": int(len(out_of_scope_try)),
            "out_of_scope_species_dropped": int(out_of_scope_try["accepted_species"].nunique()),
            "by_trait": try_common.groupby("trait_name").size().astype(int).to_dict(),
        },
        "direct": {
            "before_species_trait": len(before_direct),
            "after_species_trait": len(after_direct),
            "added_species_trait": len(after_direct - before_direct),
            "removed_species_trait": len(before_direct - after_direct),
            "conflicted_cells_excluded": (
                int(conflict_audit["resolution_status"].eq("excluded").sum())
                if len(conflict_audit)
                else 0
            ),
        },
        "validated_low": {
            "before_species_trait": len(before_low),
            "after_species_trait": len(after_low),
            "added_species_trait": len(after_low - before_low),
            "invalidated_species_trait": len(before_low - after_low),
        },
        "coverage": {
            "before": coverage_snapshot(current_coverage),
            "after": coverage_snapshot(updated_coverage),
            "gross_gain": len(after_axis - before_axis),
            "loss": len(before_axis - after_axis),
            "net_change": len(after_axis) - len(before_axis),
        },
        "lineage": {
            "rows_after_dedup": len(lineages),
            "redistribution_duplicates": len(lineage_duplicates),
            "try_lineages": int(try_common.source_lineage.nunique()),
        },
        "affected_genus_trait": len(affected),
        "safeguards": {
            "raw_try_redistributed": False,
            "family_inference": False,
            "corolla_fusion_mixed_into_floral_form": False,
            "original_source_lineage_deduplicated": True,
            "non_colour_within_source_conflicts_preexcluded": True,
        },
    }
    return {
        "lineages": lineages,
        "lineage_duplicates": lineage_duplicates,
        "direct": updated_direct,
        "low": updated_low,
        "coverage": updated_coverage,
        "rules": rules,
        "conflicts": conflict_audit,
        "summary": summary,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--try-common", type=Path, required=True)
    parser.add_argument("--formal-lineage", type=Path, required=True)
    parser.add_argument("--additional-common", type=Path, action="append", default=[])
    parser.add_argument("--current-direct", type=Path, required=True)
    parser.add_argument("--current-low", type=Path, required=True)
    parser.add_argument("--current-coverage", type=Path, required=True)
    parser.add_argument("--master-genus-map", type=Path, required=True)
    parser.add_argument("--ontology", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    try_common = load_try_common(args.try_common)
    formal = direct_evidence_from_integrated(args.formal_lineage)
    extra = [
        pd.read_csv(path, dtype=str)
        .fillna("")
        .reindex(columns=EVIDENCE_COLUMNS)
        .fillna("")
        for path in args.additional_common
    ]
    result = integrate(
        try_common=try_common,
        formal_direct_evidence=formal,
        additional_common=extra,
        current_direct=pd.read_csv(args.current_direct, dtype=str).fillna(""),
        current_low=pd.read_csv(args.current_low, dtype=str).fillna(""),
        current_coverage=pd.read_csv(args.current_coverage, dtype=str).fillna(""),
        master_genus_map=pd.read_csv(args.master_genus_map, dtype=str).fillna(""),
        ontology=load_ontology(args.ontology),
    )
    write_gz(result["lineages"], args.output / "direct_source_lineages_after_try.csv.gz")
    write_gz(result["lineage_duplicates"], args.output / "try_lineage_duplicate_audit.csv.gz")
    write_gz(result["direct"], args.output / "strict_species_direct_after_try.csv.gz")
    write_gz(result["low"], args.output / "validated_low_after_try.csv.gz")
    write_gz(result["coverage"], args.output / "strict_species_axis_coverage_after_try.csv.gz")
    write_gz(result["rules"], args.output / "affected_genus_rules_after_try.csv.gz")
    write_gz(result["conflicts"], args.output / "try_direct_conflict_audit.csv.gz")
    summary = result["summary"]
    summary["input_sha256"] = {
        "try_common": sha256(args.try_common),
        "formal_lineage": sha256(args.formal_lineage),
        "current_direct": sha256(args.current_direct),
        "current_low": sha256(args.current_low),
        "current_coverage": sha256(args.current_coverage),
        "master_genus_map": sha256(args.master_genus_map),
        "ontology": sha256(args.ontology),
    }
    (args.output / "try_integration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
