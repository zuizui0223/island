"""Build non-duplicating evidence layers for PR138 syndrome analyses.

The maximal coverage layer uses every *analysis-eligible* evidence source available in
our locked trait artifact without double-counting a species x trait combination:

1. accepted species-direct/synonym-direct evidence has priority;
2. validated genus-consensus evidence fills only species x trait gaps lacking direct evidence.

Raw unreviewed candidates and diagnostic-only low evidence are not silently promoted to
biological observations. Separate direct-only and genus-consensus-only ledgers are emitted
so evidence-layer dependence is always visible.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer

app = typer.Typer(add_completion=False, no_args_is_help=True)

KEY = ["accepted_species", "trait_name"]


def _truthy(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.lower().isin({"true", "1", "yes", "y"})


def _normalize_direct(direct: pd.DataFrame) -> pd.DataFrame:
    required = {"accepted_species", "trait_name", "resolution_status", "normalized_value"}
    missing = required - set(direct.columns)
    if missing:
        raise typer.BadParameter(f"direct ledger missing columns: {sorted(missing)}")
    work = direct.loc[
        direct["resolution_status"].fillna("").astype(str).str.lower().eq("resolved")
    ].copy()
    work["accepted_species"] = work["accepted_species"].fillna("").astype(str).str.strip()
    work["trait_name"] = work["trait_name"].fillna("").astype(str).str.strip()
    work["normalized_value"] = work["normalized_value"].fillna("").astype(str).str.strip()
    work = work.loc[
        work["accepted_species"].ne("")
        & work["trait_name"].ne("")
        & work["normalized_value"].ne("")
    ].copy()
    work = work.drop_duplicates(KEY, keep="first")
    work["state_set"] = work.get("state_set", work["normalized_value"]).fillna("").astype(str)
    work["evidence_scope"] = "species_direct"
    work["evidence_quality"] = work.get("selected_quality", work.get("quality", "medium"))
    work["evidence_quality"] = work["evidence_quality"].fillna("medium").astype(str)
    work["evidence_origin"] = "direct_species_trait_ledger"
    return work[[*KEY, "normalized_value", "state_set", "evidence_scope", "evidence_quality", "evidence_origin"]]


def _normalize_validated_low(low: pd.DataFrame) -> pd.DataFrame:
    required = {
        "accepted_species",
        "trait_name",
        "normalized_value",
        "eligible",
        "diagnostic_only",
    }
    missing = required - set(low.columns)
    if missing:
        raise typer.BadParameter(f"validated-low ledger missing columns: {sorted(missing)}")
    work = low.loc[_truthy(low["eligible"]) & ~_truthy(low["diagnostic_only"])].copy()
    work["accepted_species"] = work["accepted_species"].fillna("").astype(str).str.strip()
    work["trait_name"] = work["trait_name"].fillna("").astype(str).str.strip()
    work["normalized_value"] = work["normalized_value"].fillna("").astype(str).str.strip()
    work = work.loc[
        work["accepted_species"].ne("")
        & work["trait_name"].ne("")
        & work["normalized_value"].ne("")
    ].copy()
    work = work.drop_duplicates(KEY, keep="first")
    work["state_set"] = work.get("state_set", work["normalized_value"]).fillna("").astype(str)
    work["evidence_scope"] = work.get("evidence_scope", "genus_consensus").fillna("genus_consensus").astype(str)
    work["evidence_quality"] = work.get("evidence_quality", "low").fillna("low").astype(str)
    work["evidence_origin"] = "validated_genus_consensus"
    return work[[*KEY, "normalized_value", "state_set", "evidence_scope", "evidence_quality", "evidence_origin"]]


def build_evidence_ladder(
    direct: pd.DataFrame,
    validated_low: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    direct_norm = _normalize_direct(direct)
    low_norm = _normalize_validated_low(validated_low)

    direct_keys = pd.MultiIndex.from_frame(direct_norm[KEY])
    low_keys = pd.MultiIndex.from_frame(low_norm[KEY])
    gap_mask = ~low_keys.isin(direct_keys)
    low_gap = low_norm.loc[gap_mask].copy()

    expanded = pd.concat([direct_norm, low_gap], ignore_index=True)
    expanded["evidence_mode"] = "all_analysis_eligible"
    direct_only = direct_norm.copy()
    direct_only["evidence_mode"] = "direct_only"
    genus_only = low_norm.copy()
    genus_only["evidence_mode"] = "validated_genus_consensus_only"

    summary_rows: list[dict[str, Any]] = []
    for mode, frame in (
        ("all_analysis_eligible", expanded),
        ("direct_only", direct_only),
        ("validated_genus_consensus_only", genus_only),
    ):
        for trait_name, group in frame.groupby("trait_name", sort=True):
            summary_rows.append(
                {
                    "evidence_mode": mode,
                    "trait_name": str(trait_name),
                    "n_species": int(group["accepted_species"].nunique()),
                    "n_rows": int(len(group)),
                }
            )
    summary = pd.DataFrame(summary_rows)
    return expanded, direct_only, genus_only, summary


@app.command("run")
def run(
    direct_ledger_csv: Path = typer.Option(..., exists=True),
    validated_low_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    expanded, direct_only, genus_only, summary = build_evidence_ladder(
        pd.read_csv(direct_ledger_csv),
        pd.read_csv(validated_low_csv),
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    expanded.to_csv(output_dir / "all_analysis_eligible_trait_ledger.csv.gz", index=False)
    direct_only.to_csv(output_dir / "direct_only_trait_ledger.csv.gz", index=False)
    genus_only.to_csv(output_dir / "validated_genus_consensus_only_trait_ledger.csv.gz", index=False)
    summary.to_csv(output_dir / "evidence_ladder_coverage.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_evidence_ladder_v1",
        "direct_priority": True,
        "validated_genus_consensus_fills_only_direct_gaps": True,
        "raw_unreviewed_candidates_promoted": False,
        "n_expanded_species_trait_rows": int(len(expanded)),
        "n_direct_rows": int(len(direct_only)),
        "n_genus_consensus_rows": int(len(genus_only)),
    }
    (output_dir / "evidence_ladder_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
