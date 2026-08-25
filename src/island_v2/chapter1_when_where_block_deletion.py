"""Leave-one-spatial-block leverage sensitivity for Chapter 1 when/where."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer
import yaml

from island_v2.chapter1_when_where_omnibus import run_when_where_omnibus

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _extract_headline(
    within: pd.DataFrame,
    between: pd.DataFrame,
    *,
    stratum: str,
    north: str,
    tropical: str,
) -> dict[str, bool]:
    def within_flag(context: str) -> tuple[bool, bool]:
        x = within.loc[
            within["stratum"].eq(stratum)
            & within["support_tier"].eq("confirmatory")
            & within["context"].eq(context)
        ]
        return (not x.empty, bool(not x.empty and x.iloc[0]["where_supported"]))

    north_testable, north_supported = within_flag(north)
    tropical_testable, tropical_supported = within_flag(tropical)
    x = between.loc[
        between["stratum"].eq(stratum)
        & between["support_tier"].eq("confirmatory")
        & (
            (between["context_a"].eq(north) & between["context_b"].eq(tropical))
            | (between["context_a"].eq(tropical) & between["context_b"].eq(north))
        )
    ]
    between_testable = not x.empty
    between_supported = bool(not x.empty and x.iloc[0]["between_where_supported"])
    headline_testable = bool(north_testable and tropical_testable and between_testable)
    return {
        "north_testable": north_testable,
        "north_supported": north_supported,
        "tropical_testable": tropical_testable,
        "tropical_supported": tropical_supported,
        "between_testable": between_testable,
        "between_supported": between_supported,
        "headline_testable": headline_testable,
        "headline_replication": bool(
            headline_testable and north_supported and tropical_supported and between_supported
        ),
    }


def run_block_deletion(
    composition: pd.DataFrame,
    covariates: pd.DataFrame,
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    north = str(config.get("north_context", "northern_midlatitude"))
    tropical = str(config.get("tropical_context", "tropical"))
    context_column = str(config["context_column"])
    cluster_column = str(config["cluster_column"])
    strata = [str(x) for x in config.get("strata", ["all_native", "native_nonendemic"])]
    fit_config = {
        "baseline_covariates": [str(x) for x in config["baseline_covariates"]],
        "isolation_column": str(config["isolation_column"]),
        "context_column": context_column,
        "cluster_column": cluster_column,
        "contexts": [north, tropical],
        "strata": strata,
        "support_tiers": {"confirmatory": int(config.get("confirmatory_min_islands", 50))},
    }

    relevant_ids = set(composition.loc[composition["stratum"].isin(strata), "island_id"].astype(str))
    cov = covariates.copy()
    cov["island_id"] = cov["island_id"].astype(str)
    relevant = cov.loc[
        cov["island_id"].isin(relevant_ids)
        & cov[context_column].isin([north, tropical])
        & cov[cluster_column].notna()
    ].copy()
    blocks = sorted(relevant[cluster_column].astype(str).unique())

    rows: list[dict[str, object]] = []
    for block in blocks:
        reduced = cov.loc[cov[cluster_column].astype(str).ne(block)].copy()
        within, between, _ = run_when_where_omnibus(composition, reduced, fit_config)
        block_rows = relevant.loc[relevant[cluster_column].astype(str).eq(block)]
        block_contexts = "|".join(sorted(block_rows[context_column].astype(str).unique()))
        for stratum in strata:
            flags = _extract_headline(
                within,
                between,
                stratum=stratum,
                north=north,
                tropical=tropical,
            )
            rows.append(
                {
                    "deleted_block": block,
                    "deleted_block_contexts": block_contexts,
                    "deleted_block_islands": int(block_rows["island_id"].nunique()),
                    "stratum": stratum,
                    **flags,
                }
            )

    detail = pd.DataFrame(rows)
    summary_rows: list[dict[str, object]] = []
    for stratum, group in detail.groupby("stratum", sort=True):
        testable = group.loc[group["headline_testable"]]
        summary_rows.append(
            {
                "stratum": stratum,
                "n_deleted_blocks": int(group["deleted_block"].nunique()),
                "n_headline_testable": int(len(testable)),
                "n_headline_replications": int(testable["headline_replication"].sum()),
                "headline_replication_fraction": (
                    float(testable["headline_replication"].mean()) if len(testable) else float("nan")
                ),
                "n_headline_untestable": int(len(group) - len(testable)),
                "north_supported_fraction": float(group["north_supported"].mean()),
                "tropical_supported_fraction": float(group["tropical_supported"].mean()),
                "between_supported_fraction": float(group["between_supported"].mean()),
            }
        )
    return detail, pd.DataFrame(summary_rows)


@app.command("run")
def run(
    composition_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    detail, summary = run_block_deletion(
        pd.read_csv(composition_csv), pd.read_csv(covariates_csv), config
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    detail.to_csv(output_dir / "when_where_leave_one_block_detail.csv", index=False)
    summary.to_csv(output_dir / "when_where_leave_one_block_summary.csv", index=False)
    manifest = {
        "contract": "chapter1_when_where_leave_one_block_v1",
        "n_blocks": int(detail["deleted_block"].nunique()) if not detail.empty else 0,
        "strata": sorted(detail["stratum"].unique().tolist()) if not detail.empty else [],
        "pollinator_predictors": False,
    }
    (output_dir / "chapter1_when_where_leave_one_block_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
