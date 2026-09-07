"""Inventory and extract the complete Ackerman et al. orchid pollination workbook.

This is a source-scale acquisition stage, not a coverage promotion step.  The
entire workbook is inspected before any species-axis accounting.  Exact fixed-
universe matches are retained with sheet/row/column provenance.  Only explicit
self-(in)compatibility or autonomous self-pollination statements are normalized;
all other breeding-system labels remain unreviewed raw source values.
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
HEADER_HINT = re.compile(
    r"species|taxon|scientific|breeding|mating|self|compat|autog|agamo|pollinat|repro",
    re.I,
)
BREEDING_HINT = re.compile(
    r"breeding|mating|self|compat|autog|agamo|pollinator.?depend|repro",
    re.I,
)
VALUE_HINT = re.compile(
    r"self.?compat|self.?incompat|autonomous self|automatic self|autopollin|auto.?pollin",
    re.I,
)
SOURCE_DOI = "10.5281/zenodo.14601785"
SOURCE_URL = "https://zenodo.org/records/14601785/files/Pollination%20List%20Thru%202024.xlsx?download=1"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def find_header(raw: pd.DataFrame) -> int:
    best = (float("-inf"), 0)
    for idx in range(min(25, len(raw))):
        vals = [text(v) for v in raw.iloc[idx].tolist()]
        nonempty = sum(bool(v) for v in vals)
        hints = sum(bool(HEADER_HINT.search(v)) for v in vals if v)
        unique = len(set(v for v in vals if v))
        score = hints * 20 + nonempty + unique * 0.1
        if score > best[0]:
            best = (score, idx)
    return best[1]


def frame_from_sheet(path: Path, sheet: str) -> tuple[pd.DataFrame, int]:
    raw = pd.read_excel(path, sheet_name=sheet, header=None, dtype=object)
    header = find_header(raw)
    cols: list[str] = []
    seen: dict[str, int] = {}
    for i, value in enumerate(raw.iloc[header].tolist()):
        name = text(value) or f"unnamed_{i}"
        seen[name] = seen.get(name, 0) + 1
        cols.append(name if seen[name] == 1 else f"{name}__{seen[name]}")
    frame = raw.iloc[header + 1 :].copy()
    frame.columns = cols
    frame = frame.dropna(how="all").fillna("")
    frame["_source_excel_row"] = frame.index + 1
    return frame, header + 1


def column_inventory(sheet: str, frame: pd.DataFrame, header_excel_row: int) -> list[dict]:
    rows = []
    for col in frame.columns:
        if col == "_source_excel_row":
            continue
        vals = frame[col].map(text)
        nonempty = vals[vals.ne("")]
        samples = list(dict.fromkeys(nonempty.astype(str)))[:8]
        rows.append(
            {
                "sheet": sheet,
                "header_excel_row": header_excel_row,
                "column": col,
                "nonempty_rows": int(len(nonempty)),
                "breeding_name_hint": bool(BREEDING_HINT.search(col)),
                "breeding_value_hint_rows": int(nonempty.map(lambda x: bool(VALUE_HINT.search(x))).sum()),
                "sample_values": " | ".join(samples),
            }
        )
    return rows


def detect_species(frame: pd.DataFrame) -> tuple[pd.Series, str]:
    candidates: list[tuple[float, str]] = []
    for col in frame.columns:
        if col == "_source_excel_row":
            continue
        vals = frame[col].map(text)
        if vals.empty:
            continue
        binomial_rate = vals.map(lambda x: bool(BINOMIAL.fullmatch(x))).mean()
        name_bonus = 1.0 if re.search(r"species|scientific|taxon", col, re.I) else 0.0
        candidates.append((binomial_rate * 10 + name_bonus, col))
    if candidates:
        score, col = max(candidates)
        vals = frame[col].map(text)
        if score >= 1.5 and vals.map(lambda x: bool(BINOMIAL.fullmatch(x))).any():
            return vals, col

    genus_col = next((c for c in frame if re.fullmatch(r"genus", c.strip(), re.I)), None)
    species_col = next(
        (c for c in frame if re.fullmatch(r"species|specific epithet|epithet", c.strip(), re.I)),
        None,
    )
    if genus_col and species_col:
        combined = (frame[genus_col].map(text) + " " + frame[species_col].map(text)).str.strip()
        return combined, f"{genus_col}+{species_col}"
    return pd.Series([""] * len(frame), index=frame.index), ""


def map_explicit(raw: str) -> tuple[str, str, str]:
    folded = raw.casefold().replace("‐", "-").replace("–", "-")
    # Test incompatibility before compatibility because the latter is a substring.
    if re.search(r"\bself[- ]?incompatib(?:le|ility)\b", folded):
        return "self_incompatibility", "SI", "explicit_self_incompatible"
    if re.search(r"\bself[- ]?compatib(?:le|ility)\b", folded):
        return "self_incompatibility", "SC", "explicit_self_compatible"
    if re.search(r"\b(?:autonomous|automatic)\s+self[- ]?pollinat", folded) or re.search(
        r"\bauto[- ]?pollinat", folded
    ):
        return "autonomous_selfing_capacity", "autonomous", "explicit_autonomous_self_pollination"
    return "", "", "unmapped_raw_breeding_value"


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--workbook", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    reproductive = coverage.loc[coverage["axis"].eq(AXIS)].copy()
    if len(reproductive) != 106_295 or reproductive["accepted_species"].nunique() != 106_295:
        raise ValueError("expected one reproductive row for each of 106,295 fixed species")
    fixed = set(reproductive["accepted_species"])
    unresolved = set(reproductive.loc[reproductive["quality"].eq(""), "accepted_species"])

    xls = pd.ExcelFile(args.workbook)
    inventory: list[dict] = []
    raw_rows: list[dict] = []
    mapped_rows: list[dict] = []
    sheet_summary: list[dict] = []

    for sheet in xls.sheet_names:
        frame, header_excel_row = frame_from_sheet(args.workbook, sheet)
        inventory.extend(column_inventory(sheet, frame, header_excel_row))
        species, species_method = detect_species(frame)
        breeding_cols = [
            c
            for c in frame.columns
            if c != "_source_excel_row"
            and (
                BREEDING_HINT.search(c)
                or frame[c].map(text).map(lambda x: bool(VALUE_HINT.search(x))).sum() >= 3
            )
        ]
        exact_species_rows = 0
        for idx, row in frame.iterrows():
            source_species = text(species.loc[idx])
            if source_species not in fixed:
                continue
            exact_species_rows += 1
            for col in breeding_cols:
                raw = text(row[col])
                if not raw:
                    continue
                trait, value, mapping = map_explicit(raw)
                base = {
                    "source_species_name": source_species,
                    "accepted_species": source_species,
                    "axis": AXIS,
                    "sheet": sheet,
                    "source_excel_row": int(row["_source_excel_row"]),
                    "source_column": col,
                    "raw_value": raw,
                    "name_match_method": "exact_fixed_universe_name",
                    "species_detection_method": species_method,
                    "source_lineage": f"zenodo:{SOURCE_DOI}:row:{sheet}:{int(row['_source_excel_row'])}",
                    "source_url": SOURCE_URL,
                    "unresolved_reproductive_before_source": str(source_species in unresolved).lower(),
                    "promotion_allowed": "false",
                    "genus_rule_training_allowed": "false",
                }
                raw_rows.append({**base, "mapping_status": mapping})
                if trait and value:
                    mapped_rows.append(
                        {
                            **base,
                            "trait_name": trait,
                            "normalized_value": value,
                            "quality": "unreviewed",
                            "mapping_status": mapping,
                            "review_status": "source_scale_mapping_needs_record_review",
                        }
                    )
        sheet_summary.append(
            {
                "sheet": sheet,
                "rows_after_header": int(len(frame)),
                "header_excel_row": header_excel_row,
                "species_detection_method": species_method,
                "breeding_columns": " | ".join(breeding_cols),
                "exact_fixed_species_rows": exact_species_rows,
            }
        )

    inv = pd.DataFrame(inventory)
    raw = pd.DataFrame(raw_rows)
    mapped = pd.DataFrame(mapped_rows)
    sheets = pd.DataFrame(sheet_summary)
    inv.to_csv(args.output / "orchid_column_inventory.csv", index=False)
    raw.to_csv(args.output / "orchid_reproductive_raw_rows.csv.gz", index=False, compression="gzip")
    mapped.to_csv(args.output / "orchid_explicit_mapping_review_queue.csv.gz", index=False, compression="gzip")
    sheets.to_csv(args.output / "orchid_sheet_inventory.csv", index=False)

    unresolved_mapped_species = (
        mapped.loc[mapped["unresolved_reproductive_before_source"].eq("true"), "accepted_species"].nunique()
        if not mapped.empty
        else 0
    )
    summary = {
        "contract": "ackerman_orchid_pollination_source_scale_v2",
        "source_doi": SOURCE_DOI,
        "source_url": SOURCE_URL,
        "workbook_sha256": sha256(args.workbook),
        "source_scale_complete": True,
        "sheets": len(xls.sheet_names),
        "sheet_rows_scanned": int(sheets["rows_after_header"].sum()) if not sheets.empty else 0,
        "exact_fixed_species_rows": int(sheets["exact_fixed_species_rows"].sum()) if not sheets.empty else 0,
        "raw_reproductive_cells_retained": int(len(raw)),
        "explicit_mapped_review_rows": int(len(mapped)),
        "explicit_mapped_species": int(mapped["accepted_species"].nunique()) if not mapped.empty else 0,
        "exact_unresolved_reproductive_species_overlap": int(unresolved_mapped_species),
        "formal_gain": 0,
        "promotion_allowed": False,
        "integration_deferred_to_source_batch": True,
    }
    summary["mapped_by_trait_value"] = {
        f"{trait}:{value}": int(n)
        for (trait, value), n in (
            mapped.groupby(["trait_name", "normalized_value"]).size().items()
            if not mapped.empty
            else []
        )
    }
    (args.output / "orchid_source_scale_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
