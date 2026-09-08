"""Audit the complete Goodwillie et al. (2010) floral-display/outcrossing source.

This is a source-scale diagnostic adapter.  It inventories the public workbook,
identifies species and outcrossing-rate fields, exact-matches the fixed species
universe, and reports current reproductive-assurance overlap.  It deliberately
keeps formal promotion disabled until the source workbook/README semantics are
verified against the paper's species-mean outcrossing-rate aggregation contract.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
OUTCROSS_RE = re.compile(r"outcross|(?:^|_)tm(?:$|_)|mating", re.I)
SPECIES_KEYS = (
    "accepted_species", "scientific_name", "scientificname", "species_name",
    "speciesname", "taxon_name", "taxon", "species",
)


def norm(value: object) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(value).strip().casefold()).strip("_")


def clean_species(value: object) -> str:
    return " ".join(str(value).replace("_", " ").strip().split())


def find_species(frame: pd.DataFrame) -> tuple[pd.Series, str]:
    cols = {norm(c): c for c in frame.columns}
    for key in SPECIES_KEYS:
        if key in cols:
            s = frame[cols[key]].map(clean_species)
            # A column named merely Species may contain only epithets; require
            # a meaningful fraction of binomials before accepting it alone.
            frac = s.map(lambda x: bool(BINOMIAL.fullmatch(x))).mean() if len(s) else 0
            if frac >= 0.5:
                return s, str(cols[key])
    if "genus" in cols and "species" in cols:
        s = frame[cols["genus"]].map(clean_species) + " " + frame[cols["species"]].map(clean_species)
        return s.str.strip(), f"{cols['genus']}+{cols['species']}"
    return pd.Series([""] * len(frame), index=frame.index, dtype=object), ""


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input-xls", type=Path, required=True)
    p.add_argument("--readme", type=Path)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage.axis.eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep.accepted_species.duplicated().any():
        raise ValueError(f"expected fixed 106295-species reproductive coverage, got {len(rep)}")
    quality = dict(zip(rep.accepted_species, rep.quality, strict=True))
    universe = set(quality)
    unresolved = {s for s, q in quality.items() if not q}

    book = pd.ExcelFile(args.input_xls, engine="xlrd")
    inventories: list[dict[str, object]] = []
    row_frames: list[pd.DataFrame] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.input_xls, sheet_name=sheet, dtype=object, engine="xlrd").fillna("")
        species, species_source = find_species(frame)
        out_cols = [str(c) for c in frame.columns if OUTCROSS_RE.search(norm(c))]
        inventories.append({
            "sheet": sheet,
            "rows": len(frame),
            "columns": "|".join(map(str, frame.columns)),
            "species_source": species_source,
            "exact_binomial_rows": int(species.map(lambda x: bool(BINOMIAL.fullmatch(x))).sum()),
            "outcrossing_columns": "|".join(out_cols),
        })
        if species_source and out_cols:
            work = frame.copy()
            work["source_sheet"] = sheet
            work["source_row_in_sheet"] = work.index + 2
            work["accepted_species_candidate"] = species
            work["species_source_column"] = species_source
            work["in_fixed_universe_exact"] = work.accepted_species_candidate.isin(universe)
            work["currently_unresolved_reproductive"] = work.accepted_species_candidate.isin(unresolved)
            for c in out_cols:
                work[f"numeric__{norm(c)}"] = pd.to_numeric(work[c], errors="coerce")
            row_frames.append(work)

    inv = pd.DataFrame(inventories)
    inv.to_csv(args.output / "goodwillie2010_workbook_inventory.csv", index=False)
    if row_frames:
        rows = pd.concat(row_frames, ignore_index=True, sort=False)
    else:
        rows = pd.DataFrame()
    rows.to_csv(
        args.output / "goodwillie2010_reproductive_rows.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    if len(rows):
        exact = rows.loc[rows.in_fixed_universe_exact]
        gap = rows.loc[rows.currently_unresolved_reproductive]
        exact_species = sorted(set(exact.accepted_species_candidate))
        gap_species = sorted(set(gap.accepted_species_candidate))
    else:
        exact_species = []
        gap_species = []

    pd.DataFrame({"accepted_species": gap_species}).to_csv(
        args.output / "goodwillie2010_unresolved_overlap.csv", index=False
    )
    readme_text = ""
    if args.readme and args.readme.exists():
        readme_text = args.readme.read_text(encoding="utf-8", errors="replace")
        (args.output / "goodwillie2010_readme.txt").write_text(readme_text, encoding="utf-8")

    summary = {
        "contract": "goodwillie2010_source_scale_diagnostic_v1",
        "source_doi": "10.5061/dryad.bn5v3cm5",
        "article_doi": "10.1111/j.1469-8137.2009.03043.x",
        "workbook_sheets": len(book.sheet_names),
        "inventory": inventories,
        "rows_in_reproductive_tables": int(len(rows)),
        "exact_fixed_universe_species": len(exact_species),
        "exact_current_unresolved_reproductive_species": len(gap_species),
        "current_unresolved_species": gap_species,
        "readme_present": bool(readme_text),
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "verify workbook field semantics and whether stored values are the paper-defined species-mean outcrossing rates before applying the established tm mating-system thresholds",
    }
    (args.output / "goodwillie2010_source_scale_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
