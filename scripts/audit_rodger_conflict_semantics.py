from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

TARGETS = {"Chasmanthe aethiopica", "Ferraria crispa"}
AUTO_PREFIX = "auto."


def text(v: object) -> str:
    if v is None or pd.isna(v):
        return ""
    return " ".join(str(v).strip().split())


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--snapshot", type=Path, required=True)
    p.add_argument("--evidence", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    src = pd.read_csv(args.snapshot, dtype=str).fillna("")
    ev = pd.read_csv(args.evidence, dtype=str).fillna("")
    src["resolved_name"] = src["genus.species"].str.replace("_", " ", regex=False)
    sub = src.loc[src["resolved_name"].isin(TARGETS)].copy()
    if sub.empty:
        raise SystemExit("target species absent from Rodger snapshot")

    sub.to_csv(args.output / "target_source_rows.csv", index=False)
    evsub = ev.loc[ev["accepted_species"].isin(TARGETS)].copy()
    evsub.to_csv(args.output / "target_evidence_rows.csv", index=False)

    auto_cols = [c for c in src.columns if c.startswith(AUTO_PREFIX)]
    non_auto_cols = [c for c in src.columns if c not in auto_cols and c != "source_row"]

    species_summaries = []
    diff_rows = []
    for species, g in sub.groupby("resolved_name", sort=True):
        varying_non_auto = []
        for c in non_auto_cols:
            vals = sorted({text(x) for x in g[c] if text(x)})
            if len(vals) > 1:
                varying_non_auto.append({"column": c, "values": vals})
        varying_auto = []
        for c in auto_cols:
            vals = sorted({text(x) for x in g[c] if text(x)})
            if len(vals) > 1:
                varying_auto.append({"column": c, "values": vals})
        states = []
        for _, row in g.iterrows():
            nums = []
            for c in auto_cols:
                raw = text(row[c])
                if not raw:
                    continue
                try:
                    nums.append(float(raw))
                except ValueError:
                    pass
            state = "autonomous" if nums and any(v > 0 for v in nums) else ("absent" if nums else "missing")
            states.append({
                "source_row": text(row.get("source_row")),
                "pc.unique.number": text(row.get("pc.unique.number")),
                "source.dataset": text(row.get("source.dataset")),
                "study": text(row.get("study")),
                "taxon": text(row.get("taxon")),
                "state": state,
            })
        species_summaries.append({
            "species": species,
            "n_rows": int(len(g)),
            "source_datasets": sorted(set(g["source.dataset"].map(text))) if "source.dataset" in g else [],
            "studies": sorted(set(g["study"].map(text))) if "study" in g else [],
            "pc_unique_numbers": sorted(set(g["pc.unique.number"].map(text))) if "pc.unique.number" in g else [],
            "taxa": sorted(set(g["taxon"].map(text))) if "taxon" in g else [],
            "varying_non_auto_fields": varying_non_auto,
            "varying_auto_fields": varying_auto,
            "row_states": states,
        })
        for item in varying_non_auto:
            diff_rows.append({"species": species, **item})

    summary = {
        "contract": "rodger_within_lineage_conflict_semantics_v1",
        "snapshot_rows": int(len(src)),
        "snapshot_columns": list(src.columns),
        "auto_columns": auto_cols,
        "targets": species_summaries,
        "decision_rule": (
            "Do not resolve absent/autonomous disagreement as biological variability unless the source "
            "contains an explicit non-measurement field distinguishing population/site/treatment or otherwise "
            "documents within-species variation. Different row IDs alone are insufficient."
        ),
        "formal_gain": 0,
        "promotion_allowed": False,
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    pd.DataFrame(diff_rows).to_csv(args.output / "varying_non_auto_fields.csv", index=False)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
