"""Recover provider-level provenance for Database 1.0 FlorML source lineages.

This is an audit helper only. It never grants redistribution rights and never
changes the immutable Database 1.0 ledger. It joins exact ``florml-content:*``
lineages from the rights audit back to committed WFO acquisition evidence so
rights review can operate on flora/provider units rather than content hashes.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

DEFAULT_EVIDENCE = [
    Path("data/v2/staging/traits/direct_llm_pilot/20260810_wfo_high_yield_source_acquisition/wfo_high_yield_incremental_evidence.csv.gz"),
    Path("data/v2/staging/traits/direct_llm_pilot/20260810_wfo_fna_foc_africa_source_acquisition/wfo_combined_high_yield_evidence.csv.gz"),
]


def _provider_column(frame: pd.DataFrame) -> str:
    for name in ("source_provider", "provider"):
        if name in frame.columns:
            return name
    raise ValueError("evidence file has no provider/source_provider column")


def _provider_family(provider: str) -> str:
    normalized = re.sub(r"[^a-z0-9]+", "_", provider.strip().lower()).strip("_")
    if not normalized:
        raise ValueError("cannot create source family from blank provider")
    return f"provider_treatment:{normalized}"


def recover(lineage_details: Path, output_dir: Path, evidence_paths: list[Path]) -> pd.DataFrame:
    details = pd.read_csv(lineage_details, dtype=str).fillna("")
    target = details.loc[details["source_lineage"].str.startswith("florml-content:")].copy()
    if target.empty:
        raise ValueError("rights audit contains no florml-content lineages")

    maps: list[pd.DataFrame] = []
    for path in evidence_paths:
        if not path.exists():
            continue
        frame = pd.read_csv(path, dtype=str).fillna("")
        if "source_lineage" not in frame.columns:
            continue
        provider_col = _provider_column(frame)
        columns = ["source_lineage", provider_col]
        if "source_url" in frame.columns:
            columns.append("source_url")
        keep = frame.loc[frame["source_lineage"].str.startswith("florml-content:"), columns].copy()
        keep = keep.rename(columns={provider_col: "provider"})
        if "source_url" not in keep.columns:
            keep["source_url"] = ""
        maps.append(keep)
    if not maps:
        raise ValueError("no FlorML evidence provenance files were found")

    mapping = pd.concat(maps, ignore_index=True).drop_duplicates()
    conflicts = mapping.groupby("source_lineage")["provider"].nunique()
    conflicts = conflicts[conflicts > 1]
    if not conflicts.empty:
        raise ValueError(f"provider conflict for {len(conflicts)} FlorML lineages")
    mapping = mapping.sort_values(["source_lineage", "provider"]).drop_duplicates("source_lineage")

    recovered = target.merge(mapping, on="source_lineage", how="left", validate="one_to_one")
    recovered["provider"] = recovered["provider"].fillna("")
    recovered["source_url"] = recovered["source_url"].fillna("")
    recovered["provider_recovered"] = recovered["provider"].ne("")
    recovered["recovered_source_family"] = ""
    mask = recovered["provider_recovered"]
    recovered.loc[mask, "recovered_source_family"] = recovered.loc[mask, "provider"].map(
        _provider_family
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    recovered.to_csv(output_dir / "FLORML_PROVENANCE_RECOVERY.csv", index=False)

    override_map = recovered.loc[
        recovered["provider_recovered"],
        ["source_lineage", "recovered_source_family", "provider", "source_url"],
    ].rename(columns={"recovered_source_family": "source_family"})
    if override_map["source_lineage"].duplicated().any():
        raise ValueError("FlorML lineage-family override map contains duplicate lineages")
    override_map.to_csv(output_dir / "FLORML_LINEAGE_FAMILY_MAP.csv", index=False)

    summary = (
        recovered.groupby("provider", dropna=False)
        .agg(
            exact_lineages=("source_lineage", "nunique"),
            lineage_mentions=(
                "lineage_mentions",
                lambda x: pd.to_numeric(x, errors="coerce").fillna(0).sum(),
            ),
        )
        .reset_index()
        .sort_values("lineage_mentions", ascending=False)
    )
    summary.to_csv(output_dir / "FLORML_PROVIDER_SUMMARY.csv", index=False)

    recovered_n = int(recovered["provider_recovered"].sum())
    total_n = int(len(recovered))
    lines = [
        "# FlorML provenance recovery",
        "",
        f"- exact FlorML lineages in Database 1.0 rights audit: **{total_n:,}**",
        f"- provider recovered: **{recovered_n:,}** ({recovered_n / total_n:.1%})",
        f"- still unresolved: **{total_n - recovered_n:,}**",
        f"- audit-family overrides emitted: **{len(override_map):,}**",
        "",
        "## Provider units",
        "",
    ]
    for _, row in summary.iterrows():
        provider = str(row["provider"]) or "UNRESOLVED"
        lines.append(
            f"- `{provider}` — {int(row['exact_lineages']):,} exact lineages; "
            f"{int(row['lineage_mentions']):,} lineage mentions"
        )
    (output_dir / "FLORML_PROVENANCE_SUMMARY.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return recovered


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineage-details", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, action="append", default=[])
    args = parser.parse_args()
    evidence = args.evidence or DEFAULT_EVIDENCE
    recover(args.lineage_details, args.output_dir, evidence)


if __name__ == "__main__":
    main()
