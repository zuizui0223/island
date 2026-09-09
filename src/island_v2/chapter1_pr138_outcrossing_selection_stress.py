"""Joint observation-selection and Aegean-leverage stress for reproductive restrictions.

This module starts from the source-adjusted SI/outcrossing-restriction artifact. It does
not recompute or tune the floral response. For each reproductive restriction and source
mode it:

1. forms the predeclared attraction shift from the two source-adjusted syndrome axes;
2. estimates outcome-blind score-availability weights with the existing PR138 selection
   model over the complete island universe;
3. refits the distance response with those weights;
4. repeats the complete propensity + response fit after removing the preidentified
   Aegean/eastern-Mediterranean spatial block ``lat12_lon20``.

The goal is to ask whether the strongest mechanistic result -- persistence of the
Palearctic floral response among SI / predominantly-outcrossing taxa -- survives source
pool adjustment, observation-selection weighting, and the known spatial leverage point
simultaneously. This is still not causal identification of pollinator loss.
"""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_global_branch_selection import build_selection_weights
from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design


PREIDENTIFIED_AEGEAN_BLOCK = "lat12_lon20"
AXIS_NAME = "restricted_attraction"


def _bh(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    ok = p.notna()
    if not ok.any():
        return out
    x = p.loc[ok].to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    n = len(ranked)
    adjusted = np.minimum.accumulate((ranked * n / np.arange(1, n + 1))[::-1])[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def _standardize(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def build_restricted_attraction_scores(adjusted: pd.DataFrame) -> pd.DataFrame:
    required = {
        "restriction",
        "source_mode",
        "island_id",
        "stratum",
        "syndrome",
        "syndrome_score",
    }
    missing = required - set(adjusted.columns)
    if missing:
        raise ValueError(f"restricted adjusted scores missing columns: {sorted(missing)}")
    keep = adjusted.loc[
        adjusted["syndrome"].isin(["large_bee_like", "generalized_accessible"]),
        [
            "restriction",
            "source_mode",
            "island_id",
            "stratum",
            "syndrome",
            "syndrome_score",
        ],
    ].copy()
    pivot = (
        keep.drop_duplicates(
            ["restriction", "source_mode", "island_id", "stratum", "syndrome"]
        )
        .pivot(
            index=["restriction", "source_mode", "island_id", "stratum"],
            columns="syndrome",
            values="syndrome_score",
        )
        .reset_index()
    )
    for column in ["large_bee_like", "generalized_accessible"]:
        if column not in pivot.columns:
            pivot[column] = np.nan
    pivot["syndrome_score"] = (
        -pd.to_numeric(pivot["large_bee_like"], errors="coerce")
        + pd.to_numeric(pivot["generalized_accessible"], errors="coerce")
    ) / 2.0
    pivot["syndrome"] = AXIS_NAME
    return pivot


def _branching_config(pattern_config: dict[str, Any], realm_assignment: pd.DataFrame) -> dict[str, Any]:
    broad_contexts = [str(x) for x in pattern_config["contexts"]]
    realms = sorted(
        realm_assignment["biogeographic_realm"]
        .fillna("")
        .astype(str)
        .loc[lambda x: x.ne("")]
        .unique()
    )
    return {
        "axis_sets": {AXIS_NAME: {"axes": [AXIS_NAME]}},
        "context_layers": {
            "analysis_regime": {
                "column": str(pattern_config["context_column"]),
                "contexts": broad_contexts,
            },
            "biogeographic_realm": {
                "column": "biogeographic_realm",
                "contexts": realms,
            },
        },
    }


def _fit_context(
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    weights: pd.DataFrame,
    *,
    stratum: str,
    context_column: str,
    context: str,
    pattern_config: dict[str, Any],
    selection_mode: str,
    threshold: int,
) -> dict[str, Any]:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    work = scores.loc[scores["stratum"].eq(stratum)].copy()
    w = weights.loc[
        weights["stratum"].eq(stratum)
        & weights["syndrome"].eq(AXIS_NAME)
        & weights["selection_mode"].eq(selection_mode),
        ["island_id", "analysis_weight"],
    ].copy()
    work = work.merge(w, on="island_id", how="inner", validate="one_to_one")
    needed_cov = ["island_id", geography, context_column, cluster, *baseline]
    work = work.merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["syndrome_score", "analysis_weight", geography, *baseline]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[context_column] = work[context_column].fillna("").astype(str)
    work[cluster] = work[cluster].fillna("").astype(str)
    complete = work.loc[work[context_column].eq(context)].dropna(
        subset=["syndrome_score", "analysis_weight", geography, cluster, *baseline]
    ).copy()
    n_islands = int(complete["island_id"].nunique())
    base: dict[str, Any] = {
        "stratum": stratum,
        "context": context,
        "selection_mode": selection_mode,
        "threshold": threshold,
        "n_unique_islands": n_islands,
    }
    if n_islands < threshold:
        return {**base, "status": "not_testable"}
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    design = np.column_stack(
        [
            np.ones(len(complete), dtype=float),
            _standardize(complete[geography]),
            *[_standardize(complete[x]) for x in baseline],
        ]
    )
    coef, _, fit = _fit_weighted_clustered_design(
        complete["syndrome_score"].to_numpy(float),
        complete["analysis_weight"].to_numpy(float),
        design,
        names,
        complete[cluster].astype(str).to_numpy(),
    )
    if coef.empty:
        return {**base, "status": str(fit.get("status", "fit_failed"))}
    distance = coef.set_index("predictor").loc[f"z_{geography}"]
    weights_arr = complete["analysis_weight"].to_numpy(float)
    sw = float(np.sum(weights_arr))
    ess = sw**2 / float(np.sum(np.square(weights_arr)))
    return {
        **base,
        "status": "fit",
        "n_clusters": int(fit["n_clusters"]),
        "distance_estimate": float(distance["estimate"]),
        "distance_se": float(distance["cluster_robust_se"]),
        "distance_p": float(distance["p_value"]),
        "effective_sample_size": ess,
        "effective_sample_fraction": ess / len(complete),
        "weight_max": float(np.max(weights_arr)),
    }


def run_stress(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
    selection_config: dict[str, Any],
) -> dict[str, pd.DataFrame | dict[str, Any]]:
    attraction = build_restricted_attraction_scores(adjusted_scores)
    restrictions = sorted(attraction["restriction"].dropna().astype(str).unique())
    source_modes = sorted(attraction["source_mode"].dropna().astype(str).unique())
    strata = [
        x
        for x in ["all_native", "native_nonendemic"]
        if x in set(attraction["stratum"].astype(str))
    ]
    threshold = int(pattern_config["support_tiers"]["confirmatory"])
    branching = _branching_config(pattern_config, realm_assignment)
    sel_config = deepcopy(selection_config)
    sel_config["primary_axis_set"] = AXIS_NAME

    result_rows: list[dict[str, Any]] = []
    diagnostics_parts: list[pd.DataFrame] = []
    fit_parts: list[pd.DataFrame] = []
    universe_specs = {
        "full": None,
        "drop_aegean_lat12_lon20": PREIDENTIFIED_AEGEAN_BLOCK,
    }
    for universe_mode, excluded_block in universe_specs.items():
        if excluded_block is None:
            universe_cov = covariates.copy()
        else:
            universe_cov = covariates.loc[
                covariates[str(pattern_config["cluster_column"])].astype(str).ne(excluded_block)
            ].copy()
        valid_ids = set(universe_cov["island_id"].astype(str))
        for restriction in restrictions:
            for source_mode in source_modes:
                scenario = attraction.loc[
                    attraction["restriction"].eq(restriction)
                    & attraction["source_mode"].eq(source_mode)
                    & attraction["island_id"].astype(str).isin(valid_ids)
                ].copy()
                branch_scores = scenario[
                    ["island_id", "stratum", "syndrome", "syndrome_score"]
                ].copy()
                for layer_name, layer_spec in branching["context_layers"].items():
                    weights, fits, diagnostics = build_selection_weights(
                        branch_scores,
                        universe_cov,
                        realm_assignment,
                        pattern_config,
                        branching,
                        sel_config,
                        layer_name=layer_name,
                    )
                    if not fits.empty:
                        fits.insert(0, "source_mode", source_mode)
                        fits.insert(0, "restriction", restriction)
                        fits.insert(0, "universe_mode", universe_mode)
                        fit_parts.append(fits)
                    if not diagnostics.empty:
                        diagnostics.insert(0, "source_mode", source_mode)
                        diagnostics.insert(0, "restriction", restriction)
                        diagnostics.insert(0, "universe_mode", universe_mode)
                        diagnostics_parts.append(diagnostics)

                    if layer_name == "analysis_regime":
                        fit_cov = universe_cov
                    else:
                        fit_cov = universe_cov.merge(
                            realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates(
                                "island_id"
                            ),
                            on="island_id",
                            how="left",
                            validate="one_to_one",
                        )
                    context_column = str(layer_spec["column"])
                    present_contexts = set(fit_cov[context_column].fillna("").astype(str))
                    contexts = [
                        str(x) for x in layer_spec["contexts"] if str(x) in present_contexts
                    ]
                    for stratum in strata:
                        for selection_mode in sel_config["weight_modes"]:
                            for context in contexts:
                                row = _fit_context(
                                    scenario,
                                    fit_cov,
                                    weights,
                                    stratum=stratum,
                                    context_column=context_column,
                                    context=context,
                                    pattern_config=pattern_config,
                                    selection_mode=str(selection_mode),
                                    threshold=threshold,
                                )
                                row.update(
                                    {
                                        "universe_mode": universe_mode,
                                        "excluded_block": excluded_block or "",
                                        "restriction": restriction,
                                        "source_mode": source_mode,
                                        "context_layer": layer_name,
                                    }
                                )
                                result_rows.append(row)

    results = pd.DataFrame(result_rows)
    results["distance_q"] = np.nan
    fit_mask = results["status"].eq("fit")
    if fit_mask.any():
        family = [
            "universe_mode",
            "restriction",
            "source_mode",
            "selection_mode",
            "stratum",
            "context_layer",
        ]
        results.loc[fit_mask, "distance_q"] = (
            results.loc[fit_mask]
            .groupby(family, group_keys=False)["distance_p"]
            .apply(_bh)
        )
    manifest = {
        "contract": "chapter1_pr138_outcrossing_selection_stress_v1",
        "response": "source-adjusted (-large_bee_like + generalized_accessible)/2 within exact reproductive restrictions",
        "selection_weighting": "existing outcome-blind PR138 stabilized-IPW model, refit separately for each restriction/source/stratum/layer",
        "universe_modes": list(universe_specs),
        "preidentified_excluded_block": PREIDENTIFIED_AEGEAN_BLOCK,
        "score_values_used_in_propensity_model": False,
        "source_expectation_recomputed_within_reproductive_subset_upstream": True,
        "causal_pollinator_claimed": False,
    }
    return {
        "results": results,
        "selection_fits": pd.concat(fit_parts, ignore_index=True) if fit_parts else pd.DataFrame(),
        "selection_diagnostics": (
            pd.concat(diagnostics_parts, ignore_index=True) if diagnostics_parts else pd.DataFrame()
        ),
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
    outputs = run_stress(
        pd.read_csv(args.adjusted_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8")),
        yaml.safe_load(args.selection_config_path.read_text(encoding="utf-8")),
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in {
        "results": "outcrossing_selection_stress_results.csv",
        "selection_fits": "outcrossing_selection_stress_selection_fits.csv",
        "selection_diagnostics": "outcrossing_selection_stress_diagnostics.csv",
    }.items():
        frame = outputs[key]
        assert isinstance(frame, pd.DataFrame)
        frame.to_csv(args.output_dir / filename, index=False)
    manifest = outputs["manifest"]
    assert isinstance(manifest, dict)
    (args.output_dir / "outcrossing_selection_stress_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
