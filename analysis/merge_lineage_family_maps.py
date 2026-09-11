"""Merge audit-only source-lineage family maps with conflict detection.

The merged map changes only release-rights grouping. It never edits the immutable
Chapter 1 scientific database and never grants redistribution permission.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def merge_maps(
    paths: tuple[Path, ...],
    output_path: Path,
    summary_path: Path | None = None,
) -> dict[str, object]:
    if not paths:
        raise ValueError("at least one lineage-family map is required")

    frames: list[pd.DataFrame] = []
    input_rows: dict[str, int] = {}
    for path in paths:
        frame = pd.read_csv(path, dtype=str).fillna("")
        required = {"source_lineage", "source_family"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(
                f"lineage-family map {path} missing columns: {sorted(missing)}"
            )
        frame = frame[["source_lineage", "source_family"]].copy()
        if frame["source_lineage"].eq("").any() or frame["source_family"].eq("").any():
            raise ValueError(f"lineage-family map {path} contains blank lineage/family")
        frames.append(frame)
        input_rows[str(path)] = int(len(frame))

    combined = pd.concat(frames, ignore_index=True)
    conflicts = combined.groupby("source_lineage")["source_family"].nunique()
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        details = combined.loc[combined["source_lineage"].isin(conflicts.index)]
        pairs = details.groupby("source_lineage")["source_family"].agg(
            lambda values: "|".join(sorted(set(values)))
        )
        raise ValueError(f"conflicting lineage-family assignments: {pairs.to_dict()}")

    merged = (
        combined[["source_lineage", "source_family"]]
        .drop_duplicates()
        .sort_values(["source_lineage", "source_family"])
        .reset_index(drop=True)
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output_path, index=False)

    summary = {
        "contract": "chapter1_lineage_family_map_merge_v1",
        "input_maps": [str(path) for path in paths],
        "input_rows": input_rows,
        "merged_rows": int(len(merged)),
        "duplicate_same_family_rows_removed": int(len(combined) - len(merged)),
        "conflicting_lineages": 0,
        "scientific_database_modified": False,
        "rights_granted": False,
    }
    if summary_path is not None:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--map", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path)
    args = parser.parse_args()
    merge_maps(tuple(args.map), args.output, args.summary)


if __name__ == "__main__":
    main()
