"""Leave-one-spatial-block robustness for PR138 syndrome conclusions.

The canonical island-level syndrome scores are held fixed. Each represented spatial
block is removed in turn and the full confirmatory FDR family is refit. This tests
whether one archipelago/spatial cluster is necessary for either the northern classic
projection or the direct northern-midlatitude versus tropical syndrome-vector contrast.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr138_syndrome_analysis import (
    _between_contexts,
    _bh,
    _prepare,
    _within_context,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _fit_confirmatory_family(
    data: pd.DataFrame,
    *,
    stratum: str,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    contexts = [str(x) for x in pattern_config["contexts"]]
    threshold = int(pattern_config["support_tiers"]["confirmatory"])
    slopes_parts: list[pd.DataFrame] = []
    within_rows: list[dict[str, Any]] = []
    between_rows: list[dict[str, Any]] = []

    for context in contexts:
        slopes, result = _within_context(
            data,
            stratum=stratum,
            context_value=context,
            support_tier="confirmatory",
            threshold=threshold,
            pattern_config=pattern_config,
            syndrome_config=syndrome_config,
        )
        if not slopes.empty:
            slopes_parts.append(slopes)
        within_rows.append(result)

    for context_a, context_b in itertools.combinations(contexts, 2):
        between_rows.append(
            _between_contexts(
                data,
                stratum=stratum,
                context_a=context_a,
                context_b=context_b,
                support_tier="confirmatory",
                threshold=threshold,
                pattern_config=pattern_config,
            )
        )

    slopes = pd.concat(slopes_parts, ignore_index=True) if slopes_parts else pd.DataFrame()
    within = pd.DataFrame(within_rows)
    between = pd.DataFrame(between_rows)
    if "northern_classic_projection_one_sided_p" in within.columns:
        within["q_classic_projection"] = _bh(
            within["northern_classic_projection_one_sided_p"]
        )
    if "p_value" in between.columns:
        between["q_vector_difference"] = _bh(between["p_value"])
    return slopes, within, between


def _extract_headline(
    slopes: pd.DataFrame,
    within: pd.DataFrame,
    between: pd.DataFrame,
    *,
    north: str,
    tropical: str,
) -> dict[str, Any]:
    north_row = within.loc[within["context"].eq(north)]
    tropical_row = within.loc[within["context"].eq(tropical)]
    contrast = between.loc[
        between["context_a"].eq(north) & between["context_b"].eq(tropical)
    ]
    if contrast.empty:
        contrast = between.loc[
            between["context_a"].eq(tropical) & between["context_b"].eq(north)
        ]

    north_testable = bool(not north_row.empty and north_row.iloc[0]["status"] == "fit")
    between_testable = bool(not contrast.empty and contrast.iloc[0]["status"] == "fit")
    north_q = (
        float(north_row.iloc[0].get("q_classic_projection", np.nan))
        if north_testable
        else np.nan
    )
    between_q = (
        float(contrast.iloc[0].get("q_vector_difference", np.nan))
        if between_testable
        else np.nan
    )
    north_supported = bool(north_testable and np.isfinite(north_q) and north_q <= 0.05)
    between_supported = bool(
        between_testable and np.isfinite(between_q) and between_q <= 0.05
    )
    tropical_supported = bool(
        not tropical_row.empty
        and tropical_row.iloc[0]["status"] == "fit"
        and float(tropical_row.iloc[0].get("q_classic_projection", np.nan)) <= 0.05
    )

    north_slopes = slopes.loc[slopes["context"].eq(north)]

    def slope(name: str) -> float:
        hit = north_slopes.loc[north_slopes["syndrome"].eq(name), "distance_slope"]
        return float(hit.iloc[0]) if not hit.empty else np.nan

    large_bee = slope("large_bee_like")
    generalized = slope("generalized_accessible")
    return {
        "north_testable": north_testable,
        "north_classic_projection": (
            float(north_row.iloc[0].get("northern_classic_projection", np.nan))
            if north_testable
            else np.nan
        ),
        "north_classic_q": north_q,
        "north_classic_supported": north_supported,
        "tropical_classic_supported": tropical_supported,
        "north_tropical_testable": between_testable,
        "north_tropical_vector_q": between_q,
        "north_tropical_vector_difference_supported": between_supported,
        "north_large_bee_slope": large_bee,
        "north_generalized_slope": generalized,
        "north_attraction_direction_supported": bool(
            np.isfinite(large_bee)
            and np.isfinite(generalized)
            and large_bee < 0
            and generalized > 0
        ),
        "headline_testable": bool(north_testable and between_testable),
        "headline_replication": bool(
            north_testable and between_testable and north_supported and between_supported
        ),
    }


def run_syndrome_block_deletion(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    syndrome_config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    north = "northern_midlatitude"
    tropical = "tropical"
    cluster = str(pattern_config["cluster_column"])
    context = str(pattern_config["context_column"])
    strata = [
        s
        for s in ("all_native", "native_nonendemic")
        if s in {str(x) for x in pattern_config["strata"]}
    ]

    prepared = _prepare(island_scores, covariates, pattern_config)
    relevant = prepared.loc[
        prepared["stratum"].isin(strata)
        & prepared[context].isin([north, tropical])
        & prepared[cluster].ne("")
    ].copy()
    blocks = sorted(relevant[cluster].astype(str).unique())

    rows: list[dict[str, Any]] = []
    for block in blocks:
        reduced = prepared.loc[prepared[cluster].astype(str).ne(block)].copy()
        block_rows = relevant.loc[relevant[cluster].astype(str).eq(block)]
        block_contexts = "|".join(sorted(block_rows[context].astype(str).unique()))
        for stratum in strata:
            slopes, within, between = _fit_confirmatory_family(
                reduced,
                stratum=stratum,
                pattern_config=pattern_config,
                syndrome_config=syndrome_config,
            )
            rows.append(
                {
                    "deleted_block": block,
                    "deleted_block_contexts": block_contexts,
                    "deleted_block_islands": int(block_rows["island_id"].nunique()),
                    "stratum": stratum,
                    **_extract_headline(
                        slopes,
                        within,
                        between,
                        north=north,
                        tropical=tropical,
                    ),
                }
            )

    detail = pd.DataFrame(rows)
    summary_rows: list[dict[str, Any]] = []
    if not detail.empty:
        for stratum, group in detail.groupby("stratum", sort=True):
            testable = group.loc[group["headline_testable"]]
            summary_rows.append(
                {
                    "stratum": stratum,
                    "n_deleted_blocks": int(group["deleted_block"].nunique()),
                    "n_headline_testable": int(len(testable)),
                    "n_headline_replications": int(testable["headline_replication"].sum()),
                    "headline_replication_fraction": (
                        float(testable["headline_replication"].mean())
                        if len(testable)
                        else np.nan
                    ),
                    "north_classic_supported_fraction": float(
                        group["north_classic_supported"].mean()
                    ),
                    "north_attraction_direction_fraction": float(
                        group["north_attraction_direction_supported"].mean()
                    ),
                    "north_tropical_difference_supported_fraction": float(
                        group["north_tropical_vector_difference_supported"].mean()
                    ),
                    "n_untestable": int(len(group) - len(testable)),
                }
            )
    return detail, pd.DataFrame(summary_rows)


@app.command("run")
def run(
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    syndrome_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    syndrome_config = yaml.safe_load(syndrome_config_path.read_text(encoding="utf-8"))
    detail, summary = run_syndrome_block_deletion(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pattern_config,
        syndrome_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    detail.to_csv(output_dir / "syndrome_leave_one_block_detail.csv", index=False)
    summary.to_csv(output_dir / "syndrome_leave_one_block_summary.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_syndrome_leave_one_block_v1",
        "n_blocks": int(detail["deleted_block"].nunique()) if not detail.empty else 0,
        "strata": sorted(detail["stratum"].unique().tolist()) if not detail.empty else [],
        "canonical_island_scores_held_fixed": True,
        "full_confirmatory_fdr_family_refit_each_deletion": True,
    }
    (output_dir / "syndrome_leave_one_block_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
