"""Rank all multi-species review blocks in the Meyer mating-system workbook.

This is a source-selection audit only. It reads every worksheet containing a
`ref` and `genus_species` column, groups complete rows by source reference, and
measures exact overlap with the current unresolved reproductive-assurance
universe. No row is promoted and no trait value is inferred here.
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
SOURCE_REPO = "elenacicada/interesting_flowers"
SOURCE_COMMIT = "3d57315fe77097fe582025b884917550139801e8"
SOURCE_BLOB = "6c1f7de0f3e5a6c62d1fca2647d7adf33a19c380"


def text(v: object) -> str:
    if v is None or pd.isna(v):
        return ""
    return " ".join(str(v).strip().split())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-xlsx", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected 106295 unique reproductive species, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)

    book = pd.ExcelFile(args.input_xlsx)
    blocks: list[pd.DataFrame] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.input_xlsx, sheet_name=sheet, dtype=object).fillna("")
        lookup = {text(c).casefold(): c for c in frame.columns}
        if "ref" not in lookup or "genus_species" not in lookup:
            continue
        frame = frame.copy()
        frame["_source_sheet"] = sheet
        frame["_source_row"] = frame.index + 2
        frame["_ref"] = frame[lookup["ref"]].map(text)
        frame["_species"] = frame[lookup["genus_species"]].map(text).str.replace("_", " ", regex=False)
        frame = frame.loc[frame["_ref"].ne("") & frame["_species"].ne("")].copy()
        if frame.empty:
            continue
        # Normalize optional common columns without changing source values.
        for canonical in ("mating_system", "quant", "notes", "citation"):
            actual = next((c for k, c in lookup.items() if k == canonical), None)
            frame[f"_{canonical}"] = frame[actual].map(text) if actual is not None else ""
        blocks.append(frame[["_source_sheet", "_source_row", "_ref", "_species", "_mating_system", "_quant", "_notes", "_citation"]])

    if not blocks:
        raise ValueError("no ref/genus_species review blocks found")
    rows = pd.concat(blocks, ignore_index=True)
    rows["exact_binomial"] = rows["_species"].map(lambda x: bool(BINOMIAL.fullmatch(x)))
    rows["in_fixed_universe_exact"] = rows["_species"].isin(universe)
    rows["current_reproductive_quality"] = rows["_species"].map(quality).fillna("")
    rows["currently_unresolved_reproductive"] = rows["in_fixed_universe_exact"] & rows["current_reproductive_quality"].eq("")
    rows["quant_numeric"] = pd.to_numeric(rows["_quant"], errors="coerce")
    rows["quant_in_unit_interval"] = rows["quant_numeric"].between(0, 1, inclusive="both")
    rows.to_csv(args.output / "meyer_review_all_rows.csv.gz", index=False, compression={"method": "gzip", "mtime": 0})

    ranked: list[dict[str, object]] = []
    for ref, g in rows.groupby("_ref", sort=False):
        exact = g.loc[g["exact_binomial"]]
        fixed = g.loc[g["in_fixed_universe_exact"]]
        unresolved = g.loc[g["currently_unresolved_reproductive"]]
        qnum = g["quant_numeric"].dropna()
        mating_values = sorted(set(g["_mating_system"]) - {""})
        sheets = sorted(set(g["_source_sheet"]))
        ranked.append({
            "source_reference": ref,
            "source_sheets": "|".join(sheets),
            "source_rows": int(len(g)),
            "raw_unique_taxa": int(g["_species"].nunique()),
            "exact_binomial_species": int(exact["_species"].nunique()),
            "exact_fixed_universe_species": int(fixed["_species"].nunique()),
            "exact_current_unresolved_reproductive_species": int(unresolved["_species"].nunique()),
            "rows_with_mating_system": int(g["_mating_system"].ne("").sum()),
            "mating_system_raw_values": "|".join(mating_values[:30]),
            "rows_with_numeric_quant": int(qnum.size),
            "quant_in_unit_interval_rows": int(g["quant_in_unit_interval"].sum()),
            "quant_min": float(qnum.min()) if len(qnum) else "",
            "quant_max": float(qnum.max()) if len(qnum) else "",
            "rows_with_citation": int(g["_citation"].ne("").sum()),
            "selection_priority": int(unresolved["_species"].nunique()),
            "promotion_allowed": False,
            "status": "source_scale_screen_only",
        })
    rank = pd.DataFrame(ranked).sort_values(
        ["selection_priority", "exact_fixed_universe_species", "source_rows"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    rank.insert(0, "rank", range(1, len(rank) + 1))
    rank.to_csv(args.output / "meyer_review_source_roi.csv", index=False)

    top = rank.iloc[0].to_dict() if len(rank) else {}
    summary = {
        "contract": "meyer_review_source_block_roi_v1",
        "source_repo": SOURCE_REPO,
        "source_commit": SOURCE_COMMIT,
        "source_blob": SOURCE_BLOB,
        "source_sha256": sha256(args.input_xlsx),
        "review_blocks": int(len(rank)),
        "source_rows_scanned": int(len(rows)),
        "current_reproductive_unresolved": int(rep["quality"].eq("").sum()),
        "top_source_reference": top.get("source_reference", ""),
        "top_exact_unresolved_reproductive_species": int(top.get("exact_current_unresolved_reproductive_species", 0) or 0),
        "formal_gain": 0,
        "promotion_allowed": False,
    }
    (args.output / "meyer_review_source_roi_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(rank.head(20).to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
