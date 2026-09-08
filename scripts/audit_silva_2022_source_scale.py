"""Diagnostic audit for the Meyer workbook block labelled Silva et al (2022).

This script is intentionally non-promoting. It inventories the complete source
block, exact current reproductive gaps, and every workbook row containing the
Silva reference token so that the original publication/data source can be
identified before any SI/SC value is accepted.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"
REF = "Silva et al (2022)"
EXPECTED_SHA256 = "d6f30a8ad728a3b1b3c8ac5a5143fb935ffc64a38ce06708f6f2df4331737647"
TARGETS = ["Nectandra megapotamica", "Dinebra decipiens"]


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
        raise ValueError(f"Meyer workbook hash mismatch: {observed}")

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected 106295 reproductive rows, got {len(rep)}")
    universe = set(rep["accepted_species"])
    unresolved = set(rep.loc[rep["quality"].eq(""), "accepted_species"])

    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    direct_pairs = set(zip(direct["accepted_species"], direct["trait_name"], strict=False))

    book = pd.ExcelFile(args.input_xlsx, engine="openpyxl")
    block_parts = []
    metadata_hits = []
    for sheet in book.sheet_names:
        frame = pd.read_excel(args.input_xlsx, sheet_name=sheet, dtype=object, engine="openpyxl").fillna("")
        lookup = {str(c).strip().casefold(): c for c in frame.columns}
        # Search every cell/row for Silva or either target. This catches metadata
        # outside the main data sheet without assuming a schema.
        for idx, row in frame.iterrows():
            values = [text(v) for v in row.tolist()]
            joined = " | ".join(v for v in values if v)
            folded = joined.casefold()
            if "silva" in folded or any(t.casefold() in folded for t in TARGETS):
                metadata_hits.append({"sheet": sheet, "sheet_row": int(idx + 2), "row_text": joined})
        if "ref" not in lookup or "genus_species" not in lookup:
            continue
        ref_col = lookup["ref"]
        subset = frame.loc[frame[ref_col].map(text).eq(REF)].copy()
        if subset.empty:
            continue
        subset.insert(0, "source_sheet", sheet)
        subset.insert(1, "source_row", subset.index + 2)
        block_parts.append(subset)

    if len(block_parts) != 1:
        raise ValueError(f"expected one Silva block, found {len(block_parts)}")
    block = block_parts[0].copy()
    lookup = {str(c).strip().casefold(): c for c in block.columns}
    species_col = lookup["genus_species"]
    mating_col = lookup.get("mating_system")
    if mating_col is None:
        raise ValueError("Silva block lacks mating_system column")

    block["accepted_species_candidate"] = block[species_col].map(text).str.replace("_", " ", regex=False)
    block["mating_system_raw"] = block[mating_col].map(text)
    block["in_fixed_universe_exact"] = block["accepted_species_candidate"].isin(universe)
    block["current_reproductive_gap"] = block["accepted_species_candidate"].isin(unresolved)
    block["self_incompatibility_already_direct"] = block["accepted_species_candidate"].map(
        lambda s: (s, "self_incompatibility") in direct_pairs
    )
    block["mating_system_already_direct"] = block["accepted_species_candidate"].map(
        lambda s: (s, "mating_system") in direct_pairs
    )
    block.to_csv(
        args.output / "silva_2022_full_block.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    target = block.loc[block["accepted_species_candidate"].isin(TARGETS)].copy()
    target.to_csv(args.output / "silva_2022_target_rows.csv", index=False)
    pd.DataFrame(metadata_hits).to_csv(args.output / "silva_2022_workbook_metadata_hits.csv", index=False)

    gap = block.loc[block["current_reproductive_gap"]]
    state_counts = block["mating_system_raw"].replace("", "<blank>").value_counts().to_dict()
    gap_state_counts = gap["mating_system_raw"].replace("", "<blank>").value_counts().to_dict()
    target_states = {
        s: sorted(set(target.loc[target["accepted_species_candidate"].eq(s), "mating_system_raw"]) - {""})
        for s in TARGETS
    }
    summary = {
        "contract": "silva_2022_meyer_source_block_diagnostic_v1",
        "meyer_workbook_sha256": observed,
        "source_ref_label": REF,
        "source_rows": int(len(block)),
        "raw_unique_taxon_labels": int(block["accepted_species_candidate"].nunique()),
        "state_counts": {str(k): int(v) for k, v in state_counts.items()},
        "exact_fixed_universe_species": int(block.loc[block["in_fixed_universe_exact"], "accepted_species_candidate"].nunique()),
        "exact_current_gap_species": int(gap["accepted_species_candidate"].nunique()),
        "current_gap_state_counts": {str(k): int(v) for k, v in gap_state_counts.items()},
        "target_states": target_states,
        "target_rows": int(len(target)),
        "workbook_metadata_hit_rows": int(len(metadata_hits)),
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "identify original Silva et al. 2022 publication/data and verify that its species-level SI/SC states are direct compatibility evidence rather than a proxy or secondary recoding",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
