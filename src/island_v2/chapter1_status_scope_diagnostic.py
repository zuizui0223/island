"""Diagnose whether the broad all-observed branch is a sample-size or status-composition result.

The key comparison holds the island set fixed. It asks whether the North--Tropical
beta-binomial vector difference remains when all observed species versus only native
species are counted on the same status-supported islands. Unresolved and introduced
rows are then opened only as diagnostics; they are never relabelled as native.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import (
    build_broad_counts,
    run_probability_analysis,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _one_contrast(
    counts: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> dict[str, Any]:
    cfg = dict(config)
    cfg["strata"] = ["all_observed"]
    cfg["contexts"] = ["northern_midlatitude", "tropical"]
    cfg["between_contexts"] = [["northern_midlatitude", "tropical"]]
    _, slopes, _, omnibus = run_probability_analysis(counts, covariates, cfg)
    if omnibus.empty:
        return {"status": "not_testable"}
    row = omnibus.iloc[0].to_dict()
    if not slopes.empty:
        component = {
            str(r["outcome"]): {
                "estimate": float(r["slope_difference_b_minus_a"]),
                "se": float(r["cluster_robust_se"]),
                "p_value": float(r["p_value"]),
            }
            for _, r in slopes.iterrows()
        }
    else:
        component = {}
    row["component_results_json"] = json.dumps(component, sort_keys=True)
    return row


def run_status_scope_diagnostic(
    status_flora: pd.DataFrame,
    state_audit: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    counts = build_broad_counts(status_flora, state_audit, config)
    context = covariates[["island_id", config["context_column"]]].drop_duplicates("island_id")
    native = counts.loc[counts["stratum"].eq("all_native")].merge(
        context, on="island_id", how="left", validate="many_to_one"
    )
    focal_contexts = {"northern_midlatitude", "tropical"}
    supported_ids = set(
        native.loc[native[config["context_column"]].isin(focal_contexts), "island_id"].astype(str)
    )

    cases: list[tuple[str, pd.DataFrame]] = []
    full_observed = counts.loc[counts["stratum"].eq("all_observed")].copy()
    cases.append(("all_observed_full", full_observed))
    cases.append((
        "all_observed_same_islands_as_native",
        full_observed.loc[full_observed["island_id"].astype(str).isin(supported_ids)].copy(),
    ))

    for source_stratum, label in (
        ("all_native", "native_only_same_support_frame"),
        ("native_nonendemic", "native_nonendemic_same_support_frame"),
    ):
        part = counts.loc[
            counts["stratum"].eq(source_stratum)
            & counts["island_id"].astype(str).isin(supported_ids)
        ].copy()
        part["stratum"] = "all_observed"
        cases.append((label, part))

    flora = status_flora.copy()
    flora["origin_status"] = flora["origin_status"].fillna("unresolved").astype(str)
    flora = flora.loc[flora["island_id"].astype(str).isin(supported_ids)].copy()
    status_filters = (
        ("unresolved_only_same_islands", flora["origin_status"].eq("unresolved")),
        ("introduced_only_same_islands", flora["origin_status"].eq("introduced")),
        ("non_native_status_same_islands", ~flora["origin_status"].eq("native")),
    )
    for label, mask in status_filters:
        part_flora = flora.loc[mask].copy()
        part_counts = build_broad_counts(
            part_flora,
            state_audit,
            {**config, "strata": ["all_observed"]},
        )
        cases.append((label, part_counts))

    rows: list[dict[str, Any]] = []
    for label, part in cases:
        result = _one_contrast(part, covariates, config)
        rows.append(
            {
                "scope": label,
                "n_status_supported_island_ids": len(supported_ids),
                **result,
            }
        )
    return pd.DataFrame(rows)


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    state_audit_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    evidence_scope: str = typer.Option(...),
    output_csv: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    result = run_status_scope_diagnostic(
        pd.read_csv(status_flora_csv),
        pd.read_csv(state_audit_csv),
        pd.read_csv(covariates_csv),
        config,
    )
    result.insert(0, "evidence_scope", evidence_scope)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_csv, index=False)
    typer.echo(result.to_csv(index=False))


if __name__ == "__main__":
    app()
