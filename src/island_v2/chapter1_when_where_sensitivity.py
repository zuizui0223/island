"""Functional-form robustness checks for the Chapter 1 when/where result.

All scenarios retain the full island universe. The focal mainland-distance gradient is
never thresholded or truncated because the gradient can encode both dispersal
limitation and accessibility to a mainland/source species pool.

The checks ask whether the northern-midlatitude/tropical WHERE and BETWEEN-WHERE
conclusions persist under alternative monotone representations of the same geographic
gradient: log1p(distance), sqrt(distance), and raw distance.

Pollinator variables never enter these scenarios.
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


def _scenario_covariates(
    covariates: pd.DataFrame,
    scenario: dict[str, Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, str]:
    """Return the unchanged island universe with an alternative distance transform."""
    work = covariates.copy()
    raw_column = str(config["canonical_distance_column"])
    raw = pd.to_numeric(work[raw_column], errors="coerce")
    form = str(scenario.get("distance_form", "log1p"))

    if form == "log1p":
        gradient_column = str(config["canonical_gradient_column"])
    elif form == "sqrt":
        gradient_column = "sqrt_distance_to_continent_km"
        work[gradient_column] = np.sqrt(np.clip(raw, 0, None))
    elif form == "raw":
        gradient_column = raw_column
    else:
        raise typer.BadParameter(f"unknown distance_form: {form}")

    # The scientific contract forbids dropping islands by mainland-distance thresholds.
    if int(work["island_id"].nunique()) != int(covariates["island_id"].nunique()):
        raise RuntimeError("distance sensitivity changed the island universe")
    return work, gradient_column


def _supported_flag(
    table: pd.DataFrame,
    *,
    stratum: str,
    context: str | None = None,
    context_a: str | None = None,
    context_b: str | None = None,
    flag: str,
) -> tuple[bool, bool, int, float]:
    subset = table.loc[
        table["stratum"].eq(stratum)
        & table["support_tier"].eq("confirmatory")
    ].copy()
    if context is not None:
        subset = subset.loc[subset["context"].eq(context)]
    if context_a is not None and context_b is not None:
        direct = subset.loc[
            subset["context_a"].eq(context_a) & subset["context_b"].eq(context_b)
        ]
        reverse = subset.loc[
            subset["context_a"].eq(context_b) & subset["context_b"].eq(context_a)
        ]
        subset = pd.concat([direct, reverse], ignore_index=True)
    if subset.empty:
        return False, False, 0, float("nan")
    row = subset.iloc[0]
    return (
        True,
        bool(row.get(flag, False)),
        int(row.get("n_responses", 0)),
        float(row.get("q_within_stratum_tier", np.nan)),
    )


def run_when_where_sensitivity(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    within_parts: list[pd.DataFrame] = []
    between_parts: list[pd.DataFrame] = []
    persistence_parts: list[pd.DataFrame] = []
    summary_rows: list[dict[str, Any]] = []

    base_config = {
        "baseline_covariates": [str(x) for x in config["baseline_covariates"]],
        "cluster_column": str(config["cluster_column"]),
        "context_column": str(config["canonical_context_column"]),
        "contexts": [str(x) for x in config["contexts"]],
        "strata": [str(x) for x in config["strata"]],
        "support_tiers": {str(k): int(v) for k, v in config["support_tiers"].items()},
    }
    full_n = int(covariates["island_id"].nunique())

    for scenario in config["scenarios"]:
        name = str(scenario["name"])
        scenario_cov, gradient_column = _scenario_covariates(covariates, scenario, config)
        fit_config = {
            **base_config,
            "isolation_column": gradient_column,
        }
        within, between, persistence = run_when_where_omnibus(
            composition, scenario_cov, fit_config
        )
        for table in (within, between, persistence):
            if not table.empty:
                table.insert(0, "scenario", name)
        within_parts.append(within)
        between_parts.append(between)
        persistence_parts.append(persistence)

        for stratum in ("all_native", "native_nonendemic"):
            north_testable, north, north_n, north_q = _supported_flag(
                within,
                stratum=stratum,
                context="northern_midlatitude",
                flag="where_supported",
            )
            tropical_testable, tropical, tropical_n, tropical_q = _supported_flag(
                within,
                stratum=stratum,
                context="tropical",
                flag="where_supported",
            )
            contrast_testable, contrast, common_n, contrast_q = _supported_flag(
                between,
                stratum=stratum,
                context_a="northern_midlatitude",
                context_b="tropical",
                flag="between_where_supported",
            )
            headline_testable = bool(
                north_testable and tropical_testable and contrast_testable
            )
            n_covariate_islands = int(scenario_cov["island_id"].nunique())
            summary_rows.append(
                {
                    "scenario": name,
                    "stratum": stratum,
                    "distance_form": scenario.get("distance_form", "log1p"),
                    "n_covariate_islands": n_covariate_islands,
                    "all_islands_retained": n_covariate_islands == full_n,
                    "north_where_testable": north_testable,
                    "north_where_supported": north,
                    "north_n_responses": north_n,
                    "north_q": north_q,
                    "tropical_where_testable": tropical_testable,
                    "tropical_where_supported": tropical,
                    "tropical_n_responses": tropical_n,
                    "tropical_q": tropical_q,
                    "north_vs_tropical_testable": contrast_testable,
                    "north_vs_tropical_supported": contrast,
                    "common_n_responses": common_n,
                    "north_vs_tropical_q": contrast_q,
                    "headline_testable": headline_testable,
                    "headline_replication": bool(
                        headline_testable and north and tropical and contrast
                    ),
                }
            )

    def combine(parts: list[pd.DataFrame]) -> pd.DataFrame:
        good = [x for x in parts if not x.empty]
        return pd.concat(good, ignore_index=True) if good else pd.DataFrame()

    return (
        combine(within_parts),
        combine(between_parts),
        combine(persistence_parts),
        pd.DataFrame(summary_rows),
    )


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
    config_path: Path = typer.Option(
        Path("config/chapter1_when_where_sensitivity.yml"), exists=True
    ),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    within, between, persistence, summary = run_when_where_sensitivity(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    within.to_csv(output_dir / "when_where_sensitivity_within.csv", index=False)
    between.to_csv(output_dir / "when_where_sensitivity_between.csv", index=False)
    persistence.to_csv(output_dir / "when_where_sensitivity_persistence.csv", index=False)
    summary.to_csv(output_dir / "when_where_sensitivity_summary.csv", index=False)

    manuscript = summary.loc[summary["stratum"].isin(["all_native", "native_nonendemic"])]
    testable = manuscript.loc[manuscript["headline_testable"].fillna(False)]
    manifest = {
        "contract": "chapter1_when_where_sensitivity_v3_all_islands",
        "scenarios": [str(x["name"]) for x in config["scenarios"]],
        "distance_thresholds_used": False,
        "all_islands_required": True,
        "distance_gradient_interpretation": (
            "composite geographic axis: dispersal limitation plus mainland/source-pool accessibility"
        ),
        "n_summary_rows": int(len(summary)),
        "headline_testable_opportunities": int(len(testable)),
        "headline_replications": int(testable["headline_replication"].fillna(False).sum()),
        "headline_untestable_opportunities": int(len(manuscript) - len(testable)),
        "pollinator_predictors": False,
    }
    (output_dir / "chapter1_when_where_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
