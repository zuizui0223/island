"""Audit the committed Zenodo 18208949 orchid pending queue at source scale.

This is a selection/review audit only.  It measures how many reference-backed
SI/SC rows from the already acquired queue still hit empty reproductive cells
at the current cumulative checkpoint.  It never treats mixed pollination,
selfing evidence, or autonomous selfing/agamospermy as strict mating traits.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

AXIS = "reproductive_assurance"


def text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return " ".join(str(value).strip().split())


def present(value: object) -> bool:
    return text(value).casefold() not in {"", "0", "na", "n/a", "none", "?", "nan", "false"}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--pending", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--zenodo-metadata", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    pending = pd.read_csv(args.pending, dtype=str).fillna("")
    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[coverage["axis"].eq(AXIS), ["accepted_species", "quality"]].drop_duplicates()
    if len(rep) != 106_295 or rep["accepted_species"].duplicated().any():
        raise ValueError(f"expected 106295 reproductive species, got {len(rep)}")
    quality = dict(zip(rep["accepted_species"], rep["quality"], strict=True))
    universe = set(quality)
    unresolved = {s for s, q in quality.items() if not q}

    required = {"accepted_species", "SI", "SC", "references", "source_record", "source_sha256", "review_status"}
    missing = required.difference(pending.columns)
    if missing:
        raise ValueError(f"pending queue missing columns: {sorted(missing)}")

    rows: list[dict[str, object]] = []
    for raw in pending.to_dict("records"):
        species = text(raw["accepted_species"])
        si = present(raw["SI"])
        sc = present(raw["SC"])
        refs = text(raw["references"])
        if not (si or sc):
            continue
        value = "mixed_or_variable" if si and sc else ("SI" if si else "SC")
        rows.append({
            "accepted_species": species,
            "axis": AXIS,
            "trait_name": "self_incompatibility",
            "provisional_normalized_value": value,
            "si_raw": text(raw["SI"]),
            "sc_raw": text(raw["SC"]),
            "references": refs,
            "source_excel_row": text(raw.get("source_excel_row")),
            "source_record": text(raw["source_record"]),
            "source_sha256": text(raw["source_sha256"]),
            "prior_review_status": text(raw["review_status"]),
            "in_fixed_universe_exact": species in universe,
            "currently_unresolved_reproductive": species in unresolved,
            "reference_present": bool(refs),
            "strict_state_present": True,
            "promotion_allowed": False,
            "review_status": "needs_source_version_and_methodology_equivalence_review",
        })
    audit = pd.DataFrame(rows)
    if audit.empty:
        audit = pd.DataFrame(columns=[
            "accepted_species", "axis", "trait_name", "provisional_normalized_value",
            "si_raw", "sc_raw", "references", "source_excel_row", "source_record",
            "source_sha256", "prior_review_status", "in_fixed_universe_exact",
            "currently_unresolved_reproductive", "reference_present", "strict_state_present",
            "promotion_allowed", "review_status",
        ])
    audit.to_csv(args.output / "orchid_18208949_si_sc_audit.csv.gz", index=False, compression={"method":"gzip","mtime":0})

    eligible = audit.loc[
        audit["in_fixed_universe_exact"].astype(bool)
        & audit["currently_unresolved_reproductive"].astype(bool)
        & audit["reference_present"].astype(bool)
    ].copy()
    eligible.to_csv(args.output / "orchid_18208949_current_gap_review_queue.csv", index=False)

    meta = json.loads(args.zenodo_metadata.read_text(encoding="utf-8"))
    md = meta.get("metadata", {}) if isinstance(meta, dict) else {}
    title = text(md.get("title"))
    doi = text(md.get("doi")) or text(meta.get("doi"))
    files = meta.get("files", []) if isinstance(meta, dict) else []
    file_inventory = [
        {
            "key": text(f.get("key")),
            "size": f.get("size"),
            "checksum": text(f.get("checksum")),
        }
        for f in files if isinstance(f, dict)
    ]
    (args.output / "zenodo_18208949_file_inventory.json").write_text(
        json.dumps(file_inventory, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    summary = {
        "contract": "orchid_18208949_pending_source_roi_v1",
        "zenodo_record": 18208949,
        "zenodo_doi": doi,
        "zenodo_title": title,
        "pending_rows_total": int(len(pending)),
        "pending_unique_species": int(pending["accepted_species"].nunique()),
        "si_sc_rows": int(len(audit)),
        "si_sc_species": int(audit["accepted_species"].nunique()) if len(audit) else 0,
        "exact_fixed_universe_si_sc_species": int(audit.loc[audit["in_fixed_universe_exact"].astype(bool), "accepted_species"].nunique()) if len(audit) else 0,
        "current_unresolved_reference_backed_si_sc_rows": int(len(eligible)),
        "current_unresolved_reference_backed_si_sc_species": int(eligible["accepted_species"].nunique()) if len(eligible) else 0,
        "current_unresolved_state_counts": eligible["provisional_normalized_value"].value_counts().sort_index().astype(int).to_dict() if len(eligible) else {},
        "source_record_values": sorted(set(pending["source_record"].map(text)) - {""}),
        "source_sha256_values": sorted(set(pending["source_sha256"].map(text)) - {""}),
        "zenodo_files": len(file_inventory),
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "compare Zenodo 18208949 source identity/methods to reviewed Ackerman 14601785 contract before promotion",
    }
    (args.output / "orchid_18208949_roi_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
