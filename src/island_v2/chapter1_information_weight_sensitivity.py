"""Information-weight robustness for Chapter 1 when/where inference.

The canonical grouped-binomial model gives more likelihood information to islands
with more directly trait-resolved species. This module preserves each island's
observed trait share while varying only its effective trial count:

- canonical species-count weighting;
- capped at 100, 50, or 20 effective trials;
- equal-island weighting (one effective trial per island × atomic response).

These are pseudo-likelihood robustness checks, not ecological abundance models.
They test whether highly trait-resolved islands are necessary for the WHEN/WHERE
response-vector result.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_when_where_omnibus import run_when_where_omnibus

app = typer.Typer(add_completion=False, no_args_is_help=True)


def reweight_composition(composition: pd.DataFrame, mode: str) -> pd.DataFrame:
    work = composition.copy()
    trials = pd.to_numeric(work["trials"], errors="coerce")
    successes = pd.to_numeric(work["successes"], errors="coerce")
    valid = trials.gt(0) & successes.ge(0) & successes.le(trials)
    work = work.loc[valid].copy()
    trials = trials.loc[valid].astype(float)
    successes = successes.loc[valid].astype(float)
    share = successes / trials

    if mode == "canonical":
        effective_trials = trials
    elif mode.startswith("cap_"):
        cap = float(mode.split("_", 1)[1])
        effective_trials = np.minimum(trials, cap)
    elif mode == "equal_island":
        effective_trials = pd.Series(1.0, index=work.index)
    else:
        raise typer.BadParameter(f"unknown weighting mode: {mode}")

    work["original_trials"] = trials.to_numpy()
    work["original_successes"] = successes.to_numpy()
    work["trials"] = np.asarray(effective_trials, dtype=float)
    work["successes"] = share.to_numpy() * work["trials"].to_numpy(float)
    work["information_weight_mode"] = mode
    return work


def _headline_summary(
    within: pd.DataFrame,
    between: pd.DataFrame,
    *,
    mode: str,
    strata: tuple[str, ...],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for stratum in strata:
        w = within.loc[
            within["stratum"].eq(stratum)
            & within["support_tier"].eq("confirmatory")
        ]
        north = w.loc[w["context"].eq("northern_midlatitude")]
        tropical = w.loc[w["context"].eq("tropical")]
        b = between.loc[
            between["stratum"].eq(stratum)
            & between["support_tier"].eq("confirmatory")
            & (
                (
                    between["context_a"].eq("northern_midlatitude")
                    & between["context_b"].eq("tropical")
                )
                | (
                    between["context_a"].eq("tropical")
                    & between["context_b"].eq("northern_midlatitude")
                )
            )
        ]

        north_testable = not north.empty
        tropical_testable = not tropical.empty
        between_testable = not b.empty
        north_supported = bool(north_testable and north.iloc[0]["where_supported"])
        tropical_supported = bool(
            tropical_testable and tropical.iloc[0]["where_supported"]
        )
        between_supported = bool(
            between_testable and b.iloc[0]["between_where_supported"]
        )
        testable = north_testable and tropical_testable and between_testable
        rows.append(
            {
                "information_weight_mode": mode,
                "stratum": stratum,
                "north_testable": north_testable,
                "north_supported": north_supported,
                "north_q": float(north.iloc[0]["q_within_stratum_tier"])
                if north_testable
                else np.nan,
                "tropical_testable": tropical_testable,
                "tropical_supported": tropical_supported,
                "tropical_q": float(tropical.iloc[0]["q_within_stratum_tier"])
                if tropical_testable
                else np.nan,
                "between_testable": between_testable,
                "between_supported": between_supported,
                "between_q": float(b.iloc[0]["q_within_stratum_tier"])
                if between_testable
                else np.nan,
                "headline_testable": bool(testable),
                "headline_replication": bool(
                    testable
                    and north_supported
                    and tropical_supported
                    and between_supported
                ),
            }
        )
    return rows


def run_information_weight_sensitivity(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    omnibus = dict(config["omnibus"])
    modes = [str(x) for x in config["weighting_modes"]]
    headline_strata = tuple(str(x) for x in config["headline_strata"])
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    for mode in modes:
        weighted = reweight_composition(composition, mode)
        within, between, _ = run_when_where_omnibus(weighted, covariates, omnibus)
        if not within.empty:
            within.insert(0, "information_weight_mode", mode)
            within_parts.append(within)
        if not between.empty:
            between.insert(0, "information_weight_mode", mode)
            between_parts.append(between)
        summary_rows.extend(
            _headline_summary(within, between, mode=mode, strata=headline_strata)
        )

    within_all = (
        pd.concat(within_parts, ignore_index=True) if within_parts else pd.DataFrame()
    )
    between_all = (
        pd.concat(between_parts, ignore_index=True) if between_parts else pd.DataFrame()
    )
    return within_all, between_all, pd.DataFrame(summary_rows)


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_information_weight_sensitivity.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within, between, summary = run_information_weight_sensitivity(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "information_weight_within.csv", index=False)
    between.to_csv(output_dir / "information_weight_between.csv", index=False)
    summary.to_csv(output_dir / "information_weight_headline_summary.csv", index=False)

    testable = summary.loc[summary["headline_testable"].fillna(False)]
    manifest = {
        "contract": "chapter1_information_weight_sensitivity_v1",
        "weighting_modes": [str(x) for x in config["weighting_modes"]],
        "canonical_species_count_weighting_is_only_one_scenario": True,
        "trait_shares_preserved_within_island": True,
        "headline_testable": int(len(testable)),
        "headline_replications": int(
            testable["headline_replication"].fillna(False).sum()
        ),
    }
    (output_dir / "chapter1_information_weight_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
