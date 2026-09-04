"""Exhaustive Palearctic spatial-block deletion for exact reproductive restrictions.

Every Palearctic spatial block is deleted in turn. For each deletion this module repeats
both the source-adjusted exact-restriction response and the existing within-Palearctic
outcome-blind observation-selection weighting. The block list comes only from the frozen
geography table: no syndrome value, effect estimate, or p-value selects a deletion.

The purpose is to distinguish a spatially distributed restricted floral response from a
single-block result, and to ask whether any apparent block dependence is specific to the
preidentified Aegean/eastern-Mediterranean block or common under selection weighting.
This is a robustness sensitivity, not causal identification of pollinator loss.
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
    build_within_palearctic_weights,
)

FULL_SENTINEL = "__full_palearctic__"


def run_palearctic_restricted_block_deletion(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    selection_config: dict[str, Any],
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
        x for x in full_universe[cluster].dropna().astype(str).unique() if x
    )
    restrictions = sorted(attraction["restriction"].dropna().astype(str).unique())
    source_modes = sorted(attraction["source_mode"].dropna().astype(str).unique())
    strata = [
        x
        for x in ["all_native", "native_nonendemic"]
        if x in set(attraction["stratum"])
    ]
    selection_modes = [str(x) for x in selection_config["weight_modes"]]

    rows: list[dict[str, Any]] = []
    fit_rows: list[dict[str, Any]] = []
    diagnostic_parts: list[pd.DataFrame] = []
    specs: list[tuple[str, str | None]] = [(FULL_SENTINEL, None)] + [
        (block, block) for block in blocks
    ]
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
                    observed_ids = set(scores["island_id"].astype(str))
                    weights, fit, diagnostics = build_within_palearctic_weights(
                        observed_ids,
                        universe,
                        pattern_config,
                        selection_config,
                    )
                    fit_rows.append(
                        {
                            "deleted_block": label,
                            "is_full_palearctic": deleted_block is None,
                            "restriction": restriction,
                            "source_mode": source_mode,
                            "stratum": stratum,
                            **fit,
                        }
                    )
                    diagnostics.insert(0, "stratum", stratum)
                    diagnostics.insert(0, "source_mode", source_mode)
                    diagnostics.insert(0, "restriction", restriction)
                    diagnostics.insert(0, "deleted_block", label)
                    diagnostic_parts.append(diagnostics)
                    for selection_mode in selection_modes:
                        response = _fit_response(
                            scores,
                            weights,
                            universe,
                            pattern_config,
                            selection_mode,
                        )
                        rows.append(
                            {
                                "deleted_block": label,
                                "is_full_palearctic": deleted_block is None,
                                "restriction": restriction,
                                "source_mode": source_mode,
                                "stratum": stratum,
                                "selection_mode": selection_mode,
                                **response,
                            }
                        )

    result = pd.DataFrame(rows)
    result["distance_q_across_source_modes"] = np.nan
    fit_mask = result["status"].eq("fit")
    if fit_mask.any():
        family = ["deleted_block", "restriction", "stratum", "selection_mode"]
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
            "selection_mode",
            "full_distance_estimate",
            "full_distance_q",
            "full_n_unique_islands",
        ]
    ]
    influence = result.loc[~result["deleted_block"].eq(FULL_SENTINEL)].merge(
        full,
        on=["restriction", "source_mode", "stratum", "selection_mode"],
        how="left",
        validate="many_to_one",
    )
    influence["deletion_delta_estimate"] = (
        pd.to_numeric(influence["distance_estimate"], errors="coerce")
        - pd.to_numeric(influence["full_distance_estimate"], errors="coerce")
    )
    influence["absolute_deletion_delta"] = influence["deletion_delta_estimate"].abs()
    influence["influence_rank_desc"] = influence.groupby(
        ["restriction", "source_mode", "stratum", "selection_mode"],
        group_keys=False,
    )["absolute_deletion_delta"].rank(method="min", ascending=False)
    influence["is_preidentified_aegean"] = influence["deleted_block"].eq(
        PREIDENTIFIED_AEGEAN_BLOCK
    )

    fitted_deletions = influence.loc[influence["status"].eq("fit")].copy()
    if fitted_deletions.empty:
        summary = pd.DataFrame()
        block_support = pd.DataFrame()
        block_influence = pd.DataFrame()
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
            fitted_deletions.groupby(
                ["restriction", "stratum", "selection_mode"], as_index=False
            )
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
        block_support = (
            fitted_deletions.groupby(
                ["restriction", "stratum", "selection_mode", "deleted_block"],
                as_index=False,
            )
            .agg(
                n_source_modes=("source_mode", "nunique"),
                n_supported_source_modes=("fdr_supported", "sum"),
                all_positive=("positive_direction", "all"),
                estimate_min=("distance_estimate", "min"),
                q_max=("distance_q_across_source_modes", "max"),
            )
        )
        block_support["all_source_modes_supported"] = block_support[
            "n_supported_source_modes"
        ].eq(block_support["n_source_modes"])
        block_influence = (
            fitted_deletions.groupby("deleted_block", as_index=False)
            .agg(
                mean_absolute_delta=("absolute_deletion_delta", "mean"),
                max_absolute_delta=("absolute_deletion_delta", "max"),
                median_absolute_delta=("absolute_deletion_delta", "median"),
                n_scenarios=("absolute_deletion_delta", "count"),
                n_positive=("positive_direction", "sum"),
                n_fdr_supported=("fdr_supported", "sum"),
            )
            .sort_values(["mean_absolute_delta", "max_absolute_delta"], ascending=False)
            .reset_index(drop=True)
        )
        block_influence["overall_influence_rank"] = np.arange(
            1, len(block_influence) + 1
        )
        block_influence["is_preidentified_aegean"] = block_influence[
            "deleted_block"
        ].eq(PREIDENTIFIED_AEGEAN_BLOCK)

    manifest = {
        "contract": "chapter1_pr138_palearctic_restricted_block_deletion_v2",
        "target_population": "Palearctic islands",
        "response": "source-adjusted attraction_shift within exact reproductive restriction",
        "block_selection": "all Palearctic spatial blocks from frozen covariates",
        "outcome_values_used_to_select_blocks": False,
        "effect_estimates_used_to_select_blocks": False,
        "p_values_used_to_select_blocks": False,
        "n_palearctic_blocks": int(len(blocks)),
        "preidentified_aegean_block": PREIDENTIFIED_AEGEAN_BLOCK,
        "selection_modes": selection_modes,
        "selection_weights_recomputed_after_each_block_deletion": True,
        "claim_ceiling": (
            "Spatial leverage plus measured observation-selection sensitivity only; "
            "block-specific causation and random observation are not identified."
        ),
    }
    return {
        "results": result,
        "selection_fits": pd.DataFrame(fit_rows),
        "diagnostics": (
            pd.concat(diagnostic_parts, ignore_index=True)
            if diagnostic_parts
            else pd.DataFrame()
        ),
        "influence": influence,
        "summary": summary,
        "block_support": block_support,
        "block_influence": block_influence,
        "manifest": manifest,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--adjusted-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--selection-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pattern = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    selection = yaml.safe_load(args.selection_config_path.read_text(encoding="utf-8"))
    outputs = run_palearctic_restricted_block_deletion(
        pd.read_csv(args.adjusted_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        pattern,
        selection,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs["results"].to_csv(
        args.output_dir / "palearctic_restricted_block_deletion_results.csv", index=False
    )
    outputs["selection_fits"].to_csv(
        args.output_dir / "palearctic_restricted_block_selection_fits.csv", index=False
    )
    outputs["diagnostics"].to_csv(
        args.output_dir / "palearctic_restricted_block_ipw_diagnostics.csv", index=False
    )
    outputs["influence"].to_csv(
        args.output_dir / "palearctic_restricted_block_deletion_influence.csv", index=False
    )
    outputs["summary"].to_csv(
        args.output_dir / "palearctic_restricted_block_deletion_summary.csv", index=False
    )
    outputs["block_support"].to_csv(
        args.output_dir / "palearctic_restricted_block_support.csv", index=False
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
