"""Within-region H1 tests for WCVP floristic-status partitions.

The existing WCVP partition diagnostic only tested the North--Tropical difference.
This module instead runs the current seven-response H1 within each of the four
geographic regions for each partition, including the combined
regionally_incompatible_or_introduced group.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import (
    build_broad_counts,
    run_probability_analysis,
)
from island_v2.chapter1_wcvp_partition_diagnostic import classify_wcvp_partitions

app = typer.Typer(add_completion=False, no_args_is_help=True)

PARTITIONS = {
    "source_native": {"source_native"},
    "source_introduced": {"source_introduced"},
    "wcvp_compatible_unresolved": {"wcvp_compatible_unresolved"},
    "wcvp_incompatible_unresolved": {"wcvp_incompatible_unresolved"},
    "wcvp_unclassifiable_unresolved": {"wcvp_unclassifiable_unresolved"},
    "regional_native_compatible": {"source_native", "wcvp_compatible_unresolved"},
    "regionally_incompatible_or_introduced": {
        "source_introduced",
        "wcvp_incompatible_unresolved",
    },
}


def run_partition_within_h1(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    wcvp_ranges: pd.DataFrame,
    island_tdwg: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    classified = classify_wcvp_partitions(status_flora, wcvp_ranges, island_tdwg)
    slope_parts: list[pd.DataFrame] = []
    omnibus_parts: list[pd.DataFrame] = []
    coverage_rows: list[dict[str, Any]] = []

    for label, members in PARTITIONS.items():
        part = classified.loc[classified["wcvp_partition"].isin(members)].copy()
        coverage_rows.append(
            {
                "partition": label,
                "n_flora_rows": int(len(part)),
                "n_islands": int(part["island_id"].nunique()),
                "n_species": int(part["accepted_species"].nunique()),
            }
        )
        if part.empty:
            continue
        cfg = dict(config)
        cfg["strata"] = ["all_observed"]
        counts = build_broad_counts(part, state_audit, cfg)
        if counts.empty:
            continue
        counts["stratum"] = label
        cfg["strata"] = [label]
        within_slopes, _, within, _ = run_probability_analysis(
            counts,
            covariates,
            cfg,
        )
        if not within_slopes.empty:
            within_slopes.insert(0, "partition", label)
            slope_parts.append(within_slopes)
        if not within.empty:
            within.insert(0, "partition", label)
            omnibus_parts.append(within)

    slopes = (
        pd.concat(slope_parts, ignore_index=True)
        if slope_parts
        else pd.DataFrame()
    )
    omnibus = (
        pd.concat(omnibus_parts, ignore_index=True)
        if omnibus_parts
        else pd.DataFrame()
    )
    return slopes, omnibus, pd.DataFrame(coverage_rows)


@app.command()
def main(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    wcvp_ranges_csv: Path = typer.Option(..., exists=True),
    island_tdwg_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_dir: Path = typer.Option(...),
) -> None:
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    slopes, omnibus, coverage = run_partition_within_h1(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        pd.read_csv(covariates_csv),
        pd.read_csv(wcvp_ranges_csv),
        pd.read_csv(island_tdwg_csv),
        cfg,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    for frame in (slopes, omnibus, coverage):
        if not frame.empty:
            frame.insert(0, "evidence_scope", evidence_scope)
    slopes.to_csv(output_dir / "within_slopes.csv", index=False)
    omnibus.to_csv(output_dir / "within_omnibus.csv", index=False)
    coverage.to_csv(output_dir / "coverage.csv", index=False)
    typer.echo(omnibus.to_csv(index=False))


if __name__ == "__main__":
    app()
