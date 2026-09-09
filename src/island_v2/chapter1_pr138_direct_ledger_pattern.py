"""Direct-ledger frontend for the PR138 observed assemblage syndrome test.

The locked direct species-trait ledger contains more syndrome-relevant axes than the
three-trait Chapter 1 state-audit table. This frontend builds the same fail-closed
island x stratum broad counts from that ledger, then delegates all inference to
``chapter1_pr138_biogeographic_pattern.run_observed_pattern``.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_biogeographic_pattern import (
    _classify_signature,
    run_observed_pattern,
)
from island_v2.flora_status_support import stratum_mask

app = typer.Typer(add_completion=False, no_args_is_help=True)


def build_direct_ledger_broad_counts(
    status_flora: pd.DataFrame,
    trait_ledger: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    required_flora = {
        "island_id",
        "accepted_species",
        "origin_status",
        "endemic_status",
        "floristic_status",
    }
    required_ledger = {
        "accepted_species",
        "trait_name",
        "resolution_status",
        "normalized_value",
    }
    missing = required_flora - set(status_flora.columns)
    if missing:
        raise typer.BadParameter(f"status flora missing columns: {sorted(missing)}")
    missing = required_ledger - set(trait_ledger.columns)
    if missing:
        raise typer.BadParameter(f"direct trait ledger missing columns: {sorted(missing)}")

    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)

    ledger = trait_ledger.loc[
        trait_ledger["resolution_status"].fillna("").astype(str).str.lower().eq("resolved"),
        ["accepted_species", "trait_name", "normalized_value"],
    ].copy()
    ledger["accepted_species"] = ledger["accepted_species"].astype(str)
    ledger["trait_name"] = ledger["trait_name"].astype(str)
    ledger["normalized_value"] = ledger["normalized_value"].fillna("").astype(str)

    # Fail closed if a supposedly direct/resolved ledger contains more than one distinct
    # normalized signature for the same species x trait.
    collapsed = (
        ledger.groupby(["accepted_species", "trait_name"], as_index=False)
        .agg(
            n_signatures=("normalized_value", "nunique"),
            normalized_value=("normalized_value", "first"),
        )
    )
    collapsed = collapsed.loc[collapsed["n_signatures"].eq(1)].copy()

    parts: list[pd.DataFrame] = []
    for outcome, spec in config["broad_outcomes"].items():
        trait_name = str(spec["trait_name"])
        positive = {str(x) for x in spec["positive_states"]}
        negative = {str(x) for x in spec["negative_states"]}
        states = collapsed.loc[
            collapsed["trait_name"].eq(trait_name),
            ["accepted_species", "normalized_value"],
        ].copy()
        states["state"] = [
            _classify_signature(value, positive, negative)
            for value in states["normalized_value"]
        ]
        states = states.dropna(subset=["state"])[["accepted_species", "state"]]
        if states.empty:
            continue
        joined = flora.merge(states, on="accepted_species", how="inner", validate="many_to_one")
        for stratum in [str(x) for x in config["strata"]]:
            subset = joined.loc[
                stratum_mask(joined, stratum),
                ["island_id", "accepted_species", "state"],
            ].drop_duplicates(["island_id", "accepted_species"])
            if subset.empty:
                continue
            summary = (
                subset.groupby("island_id", as_index=False)
                .agg(successes=("state", "sum"), trials=("state", "size"))
            )
            summary["share"] = summary["successes"] / summary["trials"]
            summary["outcome"] = str(outcome)
            summary["stratum"] = stratum
            parts.append(summary)

    if not parts:
        return pd.DataFrame(
            columns=["island_id", "successes", "trials", "share", "outcome", "stratum"]
        )
    return pd.concat(parts, ignore_index=True)


@app.command("run")
def run(
    status_flora_csv: Path = typer.Option(..., exists=True),
    trait_ledger_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_pr136_biogeographic_pattern.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    counts = build_direct_ledger_broad_counts(
        pd.read_csv(status_flora_csv),
        pd.read_csv(trait_ledger_csv),
        config,
    )
    within_slopes, between_slopes, within, between = run_observed_pattern(
        counts,
        pd.read_csv(covariates_csv),
        config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    counts.to_csv(output_dir / "observed_broad_syndrome_counts.csv.gz", index=False)
    within_slopes.to_csv(output_dir / "observed_within_outcome_slopes.csv", index=False)
    between_slopes.to_csv(output_dir / "observed_between_outcome_slope_differences.csv", index=False)
    within.to_csv(output_dir / "observed_within_region_omnibus.csv", index=False)
    between.to_csv(output_dir / "observed_between_region_omnibus.csv", index=False)
    support = (
        counts.groupby(["stratum", "outcome"], as_index=False)
        .agg(n_islands=("island_id", "nunique"), n_direct_trials=("trials", "sum"))
    )
    support.to_csv(output_dir / "observed_broad_syndrome_support.csv", index=False)
    manifest = {
        "contract": str(config["contract"]),
        "input": "locked direct species-trait ledger",
        "primary_response": "observed status-stratified broad trait composition",
        "broad_outcomes": list(config["broad_outcomes"]),
        "pollinator_predictors_in_primary_model": False,
        "secondary_layer": "genus-fixed residual decomposition",
        "discussion_layer": "pollination-syndrome concordance only",
    }
    (output_dir / "observed_biogeographic_pattern_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
