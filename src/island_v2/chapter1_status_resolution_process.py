"""Audit geographic structure in Chapter 1 floristic-status resolution.

This is an observation-process diagnostic, not a biological response model. It asks
whether the probability that an observed island-species record has resolved origin
status varies with the same geography/context covariates used by the biological
analysis. No trait outcome enters this model.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_all_data_probability import run_probability_analysis

app = typer.Typer(add_completion=False, no_args_is_help=True)


def build_resolution_counts(status_flora: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "accepted_species", "origin_status"}
    missing = required - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    work = status_flora[list(required)].copy()
    work["island_id"] = work["island_id"].astype(str)
    work["accepted_species"] = work["accepted_species"].astype(str)
    work["origin_status"] = work["origin_status"].fillna("unresolved").astype(str)
    work = work.drop_duplicates(["island_id", "accepted_species"])
    work["resolved"] = work["origin_status"].ne("unresolved").astype(int)
    out = (
        work.groupby("island_id", as_index=False)
        .agg(successes=("resolved", "sum"), trials=("resolved", "size"))
    )
    out["outcome"] = "status_resolved"
    out["stratum"] = "all_observed"
    out["resolution_fraction"] = out["successes"] / out["trials"]
    return out


def _audit_config(config: dict[str, Any]) -> dict[str, Any]:
    out = dict(config)
    out["model_outcomes"] = ["status_resolved"]
    out["minimum_outcomes_per_vector"] = 1
    out["strata"] = ["all_observed"]
    return out


def run_status_resolution_audit(
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    counts = build_resolution_counts(status_flora)
    cfg = _audit_config(config)
    within_slopes, between_slopes, within, between = run_probability_analysis(
        counts, covariates, cfg
    )
    context = str(config["context_column"])
    merged = counts.merge(
        covariates[["island_id", context]].drop_duplicates("island_id"),
        on="island_id", how="left", validate="one_to_one",
    )
    raw = (
        merged.groupby(context, dropna=False, as_index=False)
        .agg(
            n_islands=("island_id", "nunique"),
            n_status_resolved=("successes", "sum"),
            n_observed_species=("trials", "sum"),
        )
    )
    raw["resolved_fraction"] = raw["n_status_resolved"] / raw["n_observed_species"]
    return raw, within_slopes, between_slopes, within, between


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    raw, within_slopes, between_slopes, within, between = run_status_resolution_audit(
        pd.read_csv(status_flora_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    raw.to_csv(output_dir / "status_resolution_by_context.csv", index=False)
    within_slopes.to_csv(output_dir / "status_resolution_within_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "status_resolution_between_slopes.csv", index=False)
    within.to_csv(output_dir / "status_resolution_within_omnibus.csv", index=False)
    between.to_csv(output_dir / "status_resolution_between_omnibus.csv", index=False)
    typer.echo(raw.to_csv(index=False))
    typer.echo(between.to_csv(index=False))


if __name__ == "__main__":
    app()
