"""Predeclared robustness checks for the Chapter 1 when/where result.

These checks do not redefine the canonical result. They ask whether the supported
northern-midlatitude/tropical WHERE and BETWEEN-WHERE conclusions persist when:

- latitude boundaries are moved inward and buffered;
- the isolation functional form changes from log1p(distance) to sqrt(distance);
- zero-distance islands are removed; or
- only islands at least 50 km from a major continent are retained.

A sensitivity that cannot meet the declared response-support threshold is recorded
as not testable, not as a biological contradiction. Pollinator variables never enter
these scenarios.
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


def _core_latitude_context(
    latitude: pd.Series,
    *,
    tropical_abs_max: float,
    north_mid_min: float,
    north_mid_max: float,
    north_high_min: float,
    south_mid_min_abs: float,
    south_mid_max_abs: float,
) -> pd.Series:
    lat = pd.to_numeric(latitude, errors="coerce")
    out = pd.Series("unresolved", index=lat.index, dtype="object")
    out.loc[lat.abs() < tropical_abs_max] = "tropical"
    out.loc[(lat >= north_mid_min) & (lat < north_mid_max)] = "northern_midlatitude"
    out.loc[lat >= north_high_min] = "northern_high_latitude"
    out.loc[
        (lat <= -south_mid_min_abs) & (lat > -south_mid_max_abs)
    ] = "southern_extratropical"
    return out


def _scenario_covariates(
    covariates: pd.DataFrame,
    scenario: dict[str, Any],
    config: dict[str, Any],
) -> tuple[pd.DataFrame, str, str]:
    work = covariates.copy()
    distance = pd.to_numeric(work["distance_to_continent_km"], errors="coerce")

    filter_name = str(scenario.get("distance_filter", "all"))
    if filter_name == "positive_only":
        work = work.loc[distance > 0].copy()
    elif filter_name == "oceanic_50km":
        work = work.loc[distance >= 50].copy()
    elif filter_name != "all":
        raise typer.BadParameter(f"unknown distance_filter: {filter_name}")

    context_rule = str(scenario.get("context_rule", "canonical"))
    if context_rule == "canonical":
        context_column = str(config["canonical_context_column"])
    elif context_rule == "core_latitudes":
        context_column = "sensitivity_context"
        core = config["core_latitude_rule"]
        work[context_column] = _core_latitude_context(
            work["island_latitude"],
            tropical_abs_max=float(core["tropical_abs_max"]),
            north_mid_min=float(core["north_mid_min"]),
            north_mid_max=float(core["north_mid_max"]),
            north_high_min=float(core["north_high_min"]),
            south_mid_min_abs=float(core["south_mid_min_abs"]),
            south_mid_max_abs=float(core["south_mid_max_abs"]),
        )
    else:
        raise typer.BadParameter(f"unknown context_rule: {context_rule}")

    isolation_form = str(scenario.get("isolation_form", "log1p"))
    if isolation_form == "log1p":
        isolation_column = str(config["canonical_isolation_column"])
    elif isolation_form == "sqrt":
        isolation_column = "sqrt_distance_to_continent_km"
        work[isolation_column] = np.sqrt(
            np.clip(pd.to_numeric(work["distance_to_continent_km"], errors="coerce"), 0, None)
        )
    else:
        raise typer.BadParameter(f"unknown isolation_form: {isolation_form}")

    return work, context_column, isolation_column


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
        "contexts": [str(x) for x in config["contexts"]],
        "strata": [str(x) for x in config["strata"]],
        "support_tiers": {str(k): int(v) for k, v in config["support_tiers"].items()},
    }

    for scenario in config["scenarios"]:
        name = str(scenario["name"])
        scenario_cov, context_column, isolation_column = _scenario_covariates(
            covariates, scenario, config
        )
        fit_config = {
            **base_config,
            "isolation_column": isolation_column,
            "context_column": context_column,
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
            summary_rows.append(
                {
                    "scenario": name,
                    "stratum": stratum,
                    "context_rule": scenario.get("context_rule", "canonical"),
                    "isolation_form": scenario.get("isolation_form", "log1p"),
                    "distance_filter": scenario.get("distance_filter", "all"),
                    "n_covariate_islands": int(scenario_cov["island_id"].nunique()),
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
        "contract": "chapter1_when_where_sensitivity_v2",
        "scenarios": [str(x["name"]) for x in config["scenarios"]],
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
