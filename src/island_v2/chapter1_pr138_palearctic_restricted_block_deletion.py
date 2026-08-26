"""Exhaustive Palearctic spatial-block deletion for exact reproductive restrictions.

This is a spatial-leverage sensitivity for the source-adjusted PR138 attraction response.
Every Palearctic spatial block is deleted in turn, with the block list derived from the
frozen geography table only. No syndrome value, fitted effect, or p-value selects which
blocks are examined.

The response and covariate model are identical to the unweighted exact-restriction fits
used by ``chapter1_pr138_palearctic_restricted_ipw``. The purpose is descriptive/inferential
robustness: determine whether the positive SI/outcrossing-restricted Palearctic response
is uniquely dependent on the preidentified Aegean/eastern-Mediterranean block or whether
its sign is broadly retained across leave-one-block-out fits.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_pr138_outcrossing_selection_stress import (
    PREIDENTIFIED_AEGEAN_BLOCK,
    _bh,
    build_restricted_attraction_scores,
)
from island_v2.chapter1_pr138_palearctic_restricted_ipw import (
    _fit_response,
    _palearctic_universe,
)

FULL_SENTINEL = "__full_palearctic__"


def _unweighted_scores(scores: pd.DataFrame) -> pd.DataFrame:
    islands = scores[["island_id"]].drop_duplicates().copy()
    islands["selection_mode"] = "unweighted"
    islands["analysis_weight"] = 1.0
    return islands


def run_palearctic_restricted_block_deletion(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    attraction = build_restricted_attraction_scores(adjusted_scores)
    cluster = str(pattern_config["cluster_column"])
    full_universe = _palearctic_universe(
        covariates,
        realm_assignment,
        pattern_config,
        excluded_block=None,
    )
    blocks = sorted(
        x
        for x in full_universe[cluster].dropna().astype(str).unique()
        if x
    )
    restrictions = sorted(attraction["restriction"].dropna().astype(str).unique())
    source_modes = sorted(attraction["source_mode"].dropna().astype(str).unique())
    strata = [x for x in ["all_native", "native_nonendemic"] if x in set(attraction["stratum"])]

    rows: list[dict[str, Any]] = []
    specs: list[tuple[str, str | None]] = [(FULL_SENTINEL, None)] + [(b, b) for b in blocks]
    for label, deleted_block in specs:
        universe = _palearctic_universe(
            covariates,
            realm_assignment,
            pattern_config,
            excluded_block=deleted_block,
        )
        valid_ids = set(universe["island_id"].astype(str))
        for restriction in restrictions:
            for source_mode in source_modes:
                scenario = attraction.loc[
                    attraction["restriction"].eq(restriction)
                    & attraction["source_mode"].eq(source_mode)
                    & attraction["island_id"].astype(str).isin(valid_ids)
                ].copy()
                for stratum in strata:
                    scores = scenario.loc[
                        scenario["stratum"].eq(stratum)
                        & scenario["syndrome_score"].notna()
                    ].copy()
                    response = _fit_response(
                        scores,
                        _unweighted_scores(scores),
                        universe,
                        pattern_config,
                        "unweighted",
                    )
                    rows.append(
                        {
                            "deleted_block": label,
                            "is_full_palearctic": deleted_block is None,
                            "restriction": restriction,
                            "source_mode": source_mode,
                            "stratum": stratum,
                            **response,
                        }
                    )

    result = pd.DataFrame(rows)
    result["distance_q_across_source_modes"] = np.nan
    fit_mask = result["status"].eq("fit")
    if fit_mask.any():
        family = ["deleted_block", "restriction", "stratum"]
        result.loc[fit_mask, "distance_q_across_source_modes"] = (
            result.loc[fit_mask]
            .groupby(family, group_keys=False)["distance_p"]
            .apply(_bh)
        )

    full = result.loc[result["deleted_block"].eq(FULL_SENTINEL)].copy()
    full = full.rename(
        columns={
            "distance_estimate": "full_distance_estimate",
            "distance_q_across_source_modes": "full_distance_q",
            "n_unique_islands": "full_n_unique_islands",
        }
    )[
        [
            "restriction",
            "source_mode",
            "stratum",
            "full_distance_estimate",
            "full_distance_q",
            "full_n_unique_islands",
        ]
    ]
    influence = result.loc[~result["deleted_block"].eq(FULL_SENTINEL)].merge(
        full,
        on=["restriction", "source_mode", "stratum"],
        how="left",
        validate="many_to_one",
    )
    influence["deletion_delta_estimate"] = (
        pd.to_numeric(influence["distance_estimate"], errors="coerce")
        - pd.to_numeric(influence["full_distance_estimate"], errors="coerce")
    )
    influence["absolute_deletion_delta"] = influence["deletion_delta_estimate"].abs()
    influence["influence_rank_desc"] = influence.groupby(
        ["restriction", "source_mode", "stratum"], group_keys=False
    )["absolute_deletion_delta"].rank(method="min", ascending=False)
    influence["is_preidentified_aegean"] = influence["deleted_block"].eq(
        PREIDENTIFIED_AEGEAN_BLOCK
    )

    fitted_deletions = influence.loc[influence["status"].eq("fit")].copy()
    if fitted_deletions.empty:
        summary = pd.DataFrame()
    else:
        fitted_deletions["positive_direction"] = (
            pd.to_numeric(fitted_deletions["distance_estimate"], errors="coerce") > 0
        )
        fitted_deletions["fdr_supported"] = (
            pd.to_numeric(
                fitted_deletions["distance_q_across_source_modes"], errors="coerce"
            )
            < float(pattern_config.get("alpha", 0.05))
        )
        summary = (
            fitted_deletions.groupby(["restriction", "stratum"], as_index=False)
            .agg(
                n_fitted_block_source_scenarios=("distance_estimate", "size"),
                n_positive=("positive_direction", "sum"),
                n_fdr_supported=("fdr_supported", "sum"),
                estimate_min=("distance_estimate", "min"),
                estimate_max=("distance_estimate", "max"),
                q_min=("distance_q_across_source_modes", "min"),
                q_max=("distance_q_across_source_modes", "max"),
                minimum_islands=("n_unique_islands", "min"),
            )
        )

    block_influence = (
        influence.loc[influence["status"].eq("fit")]
        .groupby("deleted_block", as_index=False)
        .agg(
            mean_absolute_delta=("absolute_deletion_delta", "mean"),
            max_absolute_delta=("absolute_deletion_delta", "max"),
            median_absolute_delta=("absolute_deletion_delta", "median"),
            n_scenarios=("absolute_deletion_delta", "count"),
            n_positive=("distance_estimate", lambda x: int((pd.to_numeric(x, errors="coerce") > 0).sum())),
        )
        .sort_values(["mean_absolute_delta", "max_absolute_delta"], ascending=False)
        .reset_index(drop=True)
    )
    if not block_influence.empty:
        block_influence["overall_influence_rank"] = np.arange(1, len(block_influence) + 1)
        block_influence["is_preidentified_aegean"] = block_influence["deleted_block"].eq(
            PREIDENTIFIED_AEGEAN_BLOCK
        )

    manifest = {
        "contract": "chapter1_pr138_palearctic_restricted_block_deletion_v1",
        "target_population": "Palearctic islands",
        "response": "source-adjusted attraction_shift within exact reproductive restriction",
        "block_selection": "all Palearctic spatial blocks from frozen covariates",
        "outcome_values_used_to_select_blocks": False,
        "effect_estimates_used_to_select_blocks": False,
        "p_values_used_to_select_blocks": False,
        "n_palearctic_blocks": int(len(blocks)),
        "preidentified_aegean_block": PREIDENTIFIED_AEGEAN_BLOCK,
        "selection_weighting_in_this_module": False,
        "claim_ceiling": "Spatial leverage sensitivity only; block-specific causation is not identified.",
    }
    return {
        "results": result,
        "influence": influence,
        "summary": summary,
        "block_influence": block_influence,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adjusted-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pattern = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    outputs = run_palearctic_restricted_block_deletion(
        pd.read_csv(args.adjusted_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        pattern,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs["results"].to_csv(
        args.output_dir / "palearctic_restricted_block_deletion_results.csv", index=False
    )
    outputs["influence"].to_csv(
        args.output_dir / "palearctic_restricted_block_deletion_influence.csv", index=False
    )
    outputs["summary"].to_csv(
        args.output_dir / "palearctic_restricted_block_deletion_summary.csv", index=False
    )
    outputs["block_influence"].to_csv(
        args.output_dir / "palearctic_restricted_block_influence_rank.csv", index=False
    )
    (args.output_dir / "palearctic_restricted_block_deletion_manifest.json").write_text(
        json.dumps(outputs["manifest"], indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(outputs["manifest"], indent=2))


if __name__ == "__main__":
    main()
