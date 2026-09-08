from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

TARGET = "Wilkesia hobdyi"


def text(v: object) -> str:
    if v is None or pd.isna(v):
        return ""
    return " ".join(str(v).strip().split())


def read_table(path: Path) -> pd.DataFrame | None:
    try:
        if path.suffix == ".csv" or path.name.endswith(".csv.gz"):
            return pd.read_csv(path, dtype=str).fillna("")
    except Exception:
        return None
    return None


def match_rows(frame: pd.DataFrame) -> pd.DataFrame:
    name_cols = [
        c for c in frame.columns
        if c.casefold() in {
            "accepted_species", "species", "scientific_name", "submitted_species",
            "source_scientific_name", "matched_page_name", "taxon", "name"
        }
    ]
    if not name_cols:
        return frame.iloc[0:0]
    mask = pd.Series(False, index=frame.index)
    for c in name_cols:
        mask |= frame[c].map(text).eq(TARGET)
    return frame.loc[mask].copy()


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--coverage", type=Path, required=True)
    p.add_argument("--direct-ledger", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    coverage = pd.read_csv(args.coverage, dtype=str).fillna("")
    rep = coverage.loc[
        coverage["accepted_species"].eq(TARGET)
        & coverage["axis"].eq("reproductive_assurance")
    ]
    if len(rep) != 1:
        raise ValueError("target reproductive coverage row missing or duplicated")
    direct = pd.read_csv(args.direct_ledger, dtype=str).fillna("")
    current_target_direct = direct.loc[direct["accepted_species"].eq(TARGET)].copy()

    hits = []
    file_summaries = []
    for path in sorted(args.root.rglob("*")):
        if not path.is_file() or not (path.suffix == ".csv" or path.name.endswith(".csv.gz")):
            continue
        frame = read_table(path)
        if frame is None:
            continue
        rows = match_rows(frame)
        if rows.empty:
            continue
        file_summaries.append({
            "file": str(path),
            "rows": int(len(frame)),
            "target_rows": int(len(rows)),
            "columns": list(frame.columns),
        })
        for idx, row in rows.iterrows():
            record = {str(k): text(v) for k, v in row.to_dict().items()}
            hits.append({"file": str(path), "row_index": int(idx), "record": record})

    # Compact table preserving fields useful for strict-evidence decisions.
    useful = [
        "accepted_species", "axis", "trait_name", "normalized_value", "quality",
        "evidence_quality", "evidence_scope", "review_status", "evidence_status",
        "decision", "accepted_correct", "species_identity_correct", "value_correct",
        "provenance_complete", "cultivar_contamination", "name_match_method",
        "source_provider", "source_url", "source_record_id", "source_citation",
        "source_excerpt", "source_lineage", "lineage_method", "analysis_role",
        "classification", "resolution_status", "state_set", "source_lineages",
        "source_groups", "n_independent_lineages",
    ]
    compact = []
    for hit in hits:
        rec = hit["record"]
        out = {"file": hit["file"], "row_index": hit["row_index"]}
        for c in useful:
            if c in rec and rec[c] != "":
                out[c] = rec[c]
        compact.append(out)

    summary = {
        "contract": "wilkesia_hobdyi_residual_provenance_audit_v1",
        "target": TARGET,
        "current_reproductive_quality": text(rep.iloc[0]["quality"]),
        "current_reproductive_unresolved": text(rep.iloc[0]["quality"]) == "",
        "current_direct_traits": sorted(current_target_direct["trait_name"].unique().tolist()),
        "files_with_target": len(file_summaries),
        "target_rows_found": len(hits),
        "files": file_summaries,
        "formal_gain": 0,
        "promotion_allowed": False,
        "next_gate": "identify an original species-direct SI row plus completed review/provenance contract; aggregated resolved rows alone are insufficient",
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (args.output / "target_rows.json").write_text(json.dumps(compact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    pd.DataFrame(compact).fillna("").to_csv(args.output / "target_rows.csv", index=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(compact, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
