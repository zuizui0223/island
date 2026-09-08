"""Audit the Whitehead et al. Meyer source block for strict mating-system reuse.

Diagnostic only.  The pinned Meyer input retains Whitehead population-level
numeric values in the raw ``mating_system`` column while the companion ``quant``
field identifies them as ``TM``.  This script inventories the complete block,
measures current-gap overlap, preserves all rows for the target species, searches
all workbook sheets for Whitehead metadata, and computes descriptive species-level
tm summaries.  It does not promote a mating-system class until source methodology
confirms the appropriate species aggregation contract.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
REF = "Whitehead et al"
EXPECTED_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
BINOMIAL = re.compile(r"^[A-Z][A-Za-z.-]+ [a-z][A-Za-z.-]+$")
TARGET = "Xanthorrhoea johnsonii"
EXPECTED_TM = [0.901, 0.948, 0.956, 1.0]


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
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError("current reproductive coverage is not the fixed 106295-species ledger")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)
    unresolved = {s for s, q in quality.items() if not q}

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    target_direct = direct.loc[
        direct["accepted_species"].eq(TARGET)
        & direct["trait_name"].isin(["mating_system", "self_incompatibility"])
    ].copy()
    if len(target_direct):
        raise ValueError(f"{TARGET} already has direct reproductive evidence in the current ledger")

    book = pd.ExcelFile(args.input_xlsx, engine="openpyxl")
    source_hits: list[pd.DataFrame] = []
    metadata_rows: list[dict[str, object]] = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.input_xlsx, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        if {"ref", "genus_species", "mating_system"}.issubset(lookup):
            subset = frame.loc[frame[lookup["ref"]].map(text).eq(REF)].copy()
            if len(subset):
                subset["_sheet"] = sheet
                subset["_sheet_row"] = subset.index + 2
                source_hits.append(subset)
        for idx, row in frame.iterrows():
            values = [text(value) for value in row.tolist()]
            joined = " | ".join(v for v in values if v)
            if "whitehead" in joined.casefold():
                metadata_rows.append({
                    "sheet": sheet,
                    "sheet_row": int(idx + 2),
                    "row_text": joined,
                })

    if len(source_hits) != 1:
        raise ValueError(f"expected one Whitehead block, found {len(source_hits)}")
    source = source_hits[0]
    lookup = {str(c).strip().casefold(): c for c in source.columns if not str(c).startswith("_")}
    species_col = lookup["genus_species"]
    value_col = lookup["mating_system"]
    quant_col = lookup.get("quant")
    notes_col = lookup.get("notes")
    citation_col = lookup.get("citation")

    source["accepted_species_candidate"] = source[species_col].map(text).str.replace("_", " ", regex=False)
    source["raw_tm"] = pd.to_numeric(source[value_col], errors="coerce")
    source["quant_label"] = source[quant_col].map(text) if quant_col is not None else ""
    source["exact_binomial"] = source["accepted_species_candidate"].map(lambda x: bool(BINOMIAL.fullmatch(x)))
    source["in_fixed_universe_exact"] = source["accepted_species_candidate"].isin(universe)
    source["currently_unresolved_reproductive"] = source["accepted_species_candidate"].isin(unresolved)

    target = source.loc[source["accepted_species_candidate"].eq(TARGET)].copy()
    if len(target) != 4:
        raise ValueError(f"expected four {TARGET} Whitehead rows, found {len(target)}")
    observed_tm = sorted(float(x) for x in target["raw_tm"].dropna())
    if observed_tm != EXPECTED_TM:
        raise ValueError(f"{TARGET} tm values changed: {observed_tm}")
    if not target["quant_label"].str.casefold().str.replace(" ", "", regex=False).eq("quant=tm").all():
        raise ValueError(f"{TARGET} rows are not all labelled quant=TM: {target['quant_label'].tolist()}")
    if TARGET not in unresolved:
        raise ValueError(f"{TARGET} is no longer a current reproductive gap")

    target_export = pd.DataFrame({
        "accepted_species": target["accepted_species_candidate"],
        "sheet": target["_sheet"],
        "sheet_row": target["_sheet_row"],
        "raw_tm": target["raw_tm"],
        "quant_label": target["quant_label"],
        "notes": target[notes_col].map(text) if notes_col is not None else "",
        "citation": target[citation_col].map(text) if citation_col is not None else "",
    })
    target_export.to_csv(args.output / "whitehead_target_rows.csv", index=False)
    pd.DataFrame(metadata_rows).to_csv(args.output / "whitehead_workbook_metadata_hits.csv", index=False)
    source.to_csv(
        args.output / "whitehead_full_source_audit.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    numeric = source.loc[source["raw_tm"].notna()].copy()
    exact_fixed = numeric.loc[source["exact_binomial"] & source["in_fixed_universe_exact"]]
    current_gap = exact_fixed.loc[exact_fixed["currently_unresolved_reproductive"]]
    summary = {
        "contract": "whitehead_tm_source_scale_audit_v1",
        "source_ref_label": REF,
        "source_rows": int(len(source)),
        "numeric_tm_rows": int(len(numeric)),
        "raw_unique_taxon_labels": int(source["accepted_species_candidate"].nunique()),
        "exact_fixed_universe_species_with_numeric_tm": int(exact_fixed["accepted_species_candidate"].nunique()),
        "exact_current_gap_species_with_numeric_tm": int(current_gap["accepted_species_candidate"].nunique()),
        "current_gap_species": sorted(current_gap["accepted_species_candidate"].unique().tolist()),
        "target_species": TARGET,
        "target_tm_values": observed_tm,
        "target_tm_mean": float(target["raw_tm"].mean()),
        "target_tm_min": float(target["raw_tm"].min()),
        "target_tm_max": float(target["raw_tm"].max()),
        "target_provisional_repository_threshold_class": "predominantly_outcrossing" if float(target["raw_tm"].mean()) > 0.8 else "not_predominantly_outcrossing",
        "quant_labels": sorted(set(source["quant_label"].map(text)) - {""}),
        "metadata_hits": metadata_rows,
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "identify the Whitehead source publication/data contract and verify source-defined species aggregation before mapping mean tm to mating_system",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
