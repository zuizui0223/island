"""Rank Meyer et al. source blocks against the current strict reproductive gaps.

Diagnostic only. This script does not classify or promote any row. It reads the
pinned public Input_Data_MS.xlsx workbook, inventories each original `ref` block,
and measures exact-binomial overlap with the current fixed 106,295-species
reproductive coverage. Raw mating-system labels, quantitative fields, notes, and
citations are retained so source-specific semantics can be reviewed before any
block becomes a strict packet.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
EXPECTED_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
SOURCE_REPO = "elenacicada/interesting_flowers"
SOURCE_COMMIT = "3d57315fe77097fe582025b884917550139801e8"
SOURCE_FILE = "Input_Data_MS.xlsx"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _sample(values: pd.Series, n: int = 8) -> str:
    uniq: list[str] = []
    seen: set[str] = set()
    for value in values.map(text):
        if not value or value in seen:
            continue
        seen.add(value)
        uniq.append(value)
        if len(uniq) >= n:
            break
    return " || ".join(uniq)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-xlsx", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    observed = sha256(args.input_xlsx)
    if observed != EXPECTED_SHA256:
        raise ValueError(f"pinned Meyer workbook hash mismatch: {observed}")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].copy()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    universe = set(rep["accepted_species"])
    unresolved = set(rep.loc[rep["quality"].eq(""), "accepted_species"])

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    direct_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=False))

    book = pd.ExcelFile(args.input_xlsx, engine="openpyxl")
    sheet_rows: list[dict[str, object]] = []
    block_rows: list[dict[str, object]] = []
    gap_rows: list[dict[str, object]] = []

    for sheet in book.sheet_names:
        frame = pd.read_excel(args.input_xlsx, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        sheet_rows.append({
            "sheet": sheet,
            "rows": len(frame),
            "columns": "|".join(map(str, frame.columns)),
            "has_ref": "ref" in lookup,
            "has_genus_species": "genus_species" in lookup,
            "has_mating_system": "mating_system" in lookup,
            "has_quant": "quant" in lookup,
            "has_notes": "notes" in lookup,
            "has_citation": "citation" in lookup,
        })
        if "ref" not in lookup or "genus_species" not in lookup:
            continue

        ref_col = lookup["ref"]
        species_col = lookup["genus_species"]
        mating_col = lookup.get("mating_system")
        quant_col = lookup.get("quant")
        notes_col = lookup.get("notes")
        citation_col = lookup.get("citation")
        work = frame.copy()
        work["_source_row"] = work.index + 2
        work["_ref"] = work[ref_col].map(text)
        work["_species"] = work[species_col].map(text).str.replace("_", " ", regex=False)
        work["_exact_binomial"] = work["_species"].map(lambda x: bool(BINOMIAL.fullmatch(x)))
        work["_in_universe"] = work["_species"].isin(universe)
        work["_current_gap"] = work["_species"].isin(unresolved)

        for ref, group in work.groupby("_ref", sort=True, dropna=False):
            if not ref:
                continue
            exact = group.loc[group["_exact_binomial"]]
            fixed = exact.loc[exact["_in_universe"]]
            gap = fixed.loc[fixed["_current_gap"]]
            block_rows.append({
                "sheet": sheet,
                "ref": ref,
                "rows": len(group),
                "raw_taxon_labels": group["_species"].nunique(),
                "exact_binomial_rows": int(group["_exact_binomial"].sum()),
                "exact_fixed_universe_species": int(fixed["_species"].nunique()),
                "exact_current_gap_species": int(gap["_species"].nunique()),
                "current_gap_rows": len(gap),
                "mating_system_nonempty_rows": int(group[mating_col].map(text).ne("").sum()) if mating_col is not None else 0,
                "mating_system_examples": _sample(group[mating_col], 12) if mating_col is not None else "",
                "current_gap_mating_system_examples": _sample(gap[mating_col], 12) if mating_col is not None and len(gap) else "",
                "quant_nonempty_rows": int(group[quant_col].map(text).ne("").sum()) if quant_col is not None else 0,
                "quant_numeric_rows": int(pd.to_numeric(group[quant_col], errors="coerce").notna().sum()) if quant_col is not None else 0,
                "quant_examples": _sample(group[quant_col]) if quant_col is not None else "",
                "notes_examples": _sample(group[notes_col], 5) if notes_col is not None else "",
                "citation_examples": _sample(group[citation_col], 5) if citation_col is not None else "",
                "all_columns": "|".join(map(str, frame.columns)),
                "formal_gain": 0,
                "promotion_allowed": False,
            })
            if len(gap):
                for row in gap.to_dict("records"):
                    species = text(row["_species"])
                    gap_rows.append({
                        "sheet": sheet,
                        "ref": ref,
                        "source_row": int(row["_source_row"]),
                        "accepted_species_candidate": species,
                        "mating_system_raw": text(row.get(mating_col, "")) if mating_col is not None else "",
                        "quant": text(row.get(quant_col, "")) if quant_col is not None else "",
                        "notes": text(row.get(notes_col, "")) if notes_col is not None else "",
                        "citation": text(row.get(citation_col, "")) if citation_col is not None else "",
                        "self_incompatibility_already_direct": (species, "self_incompatibility") in direct_pairs,
                        "mating_system_already_direct": (species, "mating_system") in direct_pairs,
                        "status": "current_axis_gap_source_row_only_not_evidence",
                    })

    sheets = pd.DataFrame(sheet_rows)
    blocks = pd.DataFrame(block_rows)
    gaps = pd.DataFrame(gap_rows)
    if len(blocks):
        blocks = blocks.sort_values(
            ["exact_current_gap_species", "exact_fixed_universe_species", "rows", "ref"],
            ascending=[False, False, False, True],
            kind="stable",
        ).reset_index(drop=True)
    sheets.to_csv(args.output / "meyer_sheet_inventory.csv", index=False)
    blocks.to_csv(args.output / "meyer_source_block_roi.csv", index=False)
    gaps.to_csv(args.output / "meyer_current_gap_source_rows.csv", index=False)

    summary = {
        "contract": "meyer_source_block_current_gap_audit_v2",
        "source_repo": SOURCE_REPO,
        "source_commit": SOURCE_COMMIT,
        "source_file": SOURCE_FILE,
        "source_sha256": observed,
        "current_reproductive_unresolved": len(unresolved),
        "sheets": len(sheets),
        "source_blocks": len(blocks),
        "blocks_with_current_gap": int(blocks["exact_current_gap_species"].gt(0).sum()) if len(blocks) else 0,
        "current_gap_source_rows": len(gaps),
        "unique_current_gap_species_seen_any_block": int(gaps["accepted_species_candidate"].nunique()) if len(gaps) else 0,
        "top_blocks": blocks.head(20)[[
            "sheet", "ref", "rows", "exact_fixed_universe_species", "exact_current_gap_species",
            "current_gap_rows", "mating_system_examples", "current_gap_mating_system_examples",
            "quant_numeric_rows", "quant_examples", "notes_examples", "citation_examples"
        ]].to_dict("records") if len(blocks) else [],
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "source-specific methodology review and species-level aggregation before any strict packet",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
