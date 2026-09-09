"""Test whether the isolation-associated attraction shift depends on reproductive assurance.

This is a mechanistic-concordance sensitivity, not causal mediation. It fits the already
source-adjusted plant-side responses as

    attraction_shift ~ distance + selfing_core + distance:selfing_core + covariates

within each declared biogeographic context. No threshold on selfing is introduced: the
interaction uses the continuous strict reproductive-assurance score.
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

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design


def _standardize(series: pd.Series) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(values))
    sd = float(np.std(values, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (values - mean) / sd


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


def build_interaction_table(adjusted_scores: pd.DataFrame) -> pd.DataFrame:
    required = {"source_mode", "island_id", "stratum", "syndrome", "syndrome_score"}
    missing = required - set(adjusted_scores.columns)
    if missing:
        raise ValueError(f"source-adjusted scores missing columns: {sorted(missing)}")
    keep = adjusted_scores.loc[
        adjusted_scores["syndrome"].isin(
            ["large_bee_like", "generalized_accessible", "selfing_core"]
        ),
        ["source_mode", "island_id", "stratum", "syndrome", "syndrome_score"],
    ].copy()
    pivot = (
        keep.drop_duplicates(["source_mode", "island_id", "stratum", "syndrome"])
        .pivot(
            index=["source_mode", "island_id", "stratum"],
            columns="syndrome",
            values="syndrome_score",
        )
        .reset_index()
    )
    for column in ["large_bee_like", "generalized_accessible", "selfing_core"]:
        if column not in pivot.columns:
            pivot[column] = np.nan
    pivot["attraction_shift"] = (
        -pd.to_numeric(pivot["large_bee_like"], errors="coerce")
        + pd.to_numeric(pivot["generalized_accessible"], errors="coerce")
    ) / 2.0
    return pivot


def _fit_interaction(
    data: pd.DataFrame,
    *,
    geography: str,
    baseline: list[str],
    cluster: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    needed = ["attraction_shift", geography, "selfing_core", *baseline, cluster]
    complete = data.dropna(subset=needed).copy()
    if complete.empty:
        return pd.DataFrame(), {"status": "no_complete_rows", "n_unique_islands": 0}

    z_distance = _standardize(complete[geography])
    z_selfing = _standardize(complete["selfing_core"])
    names = ["intercept", f"z_{geography}", "z_selfing_core", "z_distance:z_selfing_core"]
    columns = [
        np.ones(len(complete), dtype=float),
        z_distance,
        z_selfing,
        z_distance * z_selfing,
    ]
    for predictor in baseline:
        names.append(f"z_{predictor}")
        columns.append(_standardize(complete[predictor]))

    coef, _, fit = _fit_weighted_clustered_design(
        pd.to_numeric(complete["attraction_shift"], errors="coerce").to_numpy(float),
        np.ones(len(complete), dtype=float),
        np.column_stack(columns),
        names,
        complete[cluster].astype(str).to_numpy(),
    )
    return coef, {
        **fit,
        "n_unique_islands": int(complete["island_id"].nunique()),
    }


def run_interaction_analysis(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
    *,
    context_column: str,
    contexts: list[str],
    context_layer: str,
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    strata = [str(x) for x in pattern_config["strata"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}

    table = build_interaction_table(adjusted_scores)
    needed_cov = ["island_id", geography, context_column, cluster, *baseline]
    missing = set(needed_cov) - set(covariates.columns)
    if missing:
        raise ValueError(f"covariates missing columns: {sorted(missing)}")
    data = table.merge(
        covariates[needed_cov].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    for column in ["attraction_shift", "selfing_core", geography, *baseline]:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[context_column] = data[context_column].fillna("").astype(str)
    data[cluster] = data[cluster].fillna("").astype(str)

    rows: list[dict[str, Any]] = []
    for source_mode in data["source_mode"].dropna().astype(str).drop_duplicates():
        source_data = data.loc[data["source_mode"].eq(source_mode)].copy()
        for stratum in strata:
            for context in contexts:
                base = source_data.loc[
                    source_data["stratum"].eq(stratum)
                    & source_data[context_column].eq(context)
                ].copy()
                complete_n = int(
                    base.dropna(
                        subset=["attraction_shift", geography, "selfing_core", *baseline, cluster]
                    )["island_id"].nunique()
                )
                for tier, threshold in tiers.items():
                    if complete_n < threshold:
                        rows.append(
                            {
                                "context_layer": context_layer,
                                "source_mode": source_mode,
                                "stratum": stratum,
                                "context": context,
                                "support_tier": tier,
                                "threshold": threshold,
                                "status": "not_testable",
                                "n_unique_islands": complete_n,
                            }
                        )
                        continue
                    coef, fit = _fit_interaction(
                        base,
                        geography=geography,
                        baseline=baseline,
                        cluster=cluster,
                    )
                    if coef.empty:
                        rows.append(
                            {
                                "context_layer": context_layer,
                                "source_mode": source_mode,
                                "stratum": stratum,
                                "context": context,
                                "support_tier": tier,
                                "threshold": threshold,
                                "status": str(fit.get("status", "fit_failed")),
                                "n_unique_islands": int(fit.get("n_unique_islands", complete_n)),
                            }
                        )
                        continue
                    indexed = coef.set_index("predictor")
                    distance = indexed.loc[f"z_{geography}"]
                    selfing = indexed.loc["z_selfing_core"]
                    interaction = indexed.loc["z_distance:z_selfing_core"]
                    rows.append(
                        {
                            "context_layer": context_layer,
                            "source_mode": source_mode,
                            "stratum": stratum,
                            "context": context,
                            "support_tier": tier,
                            "threshold": threshold,
                            "status": "fit",
                            "n_unique_islands": int(fit["n_unique_islands"]),
                            "n_clusters": int(fit["n_clusters"]),
                            "distance_estimate": float(distance["estimate"]),
                            "distance_se": float(distance["cluster_robust_se"]),
                            "distance_p": float(distance["p_value"]),
                            "selfing_core_estimate": float(selfing["estimate"]),
                            "selfing_core_se": float(selfing["cluster_robust_se"]),
                            "selfing_core_p": float(selfing["p_value"]),
                            "distance_x_selfing_estimate": float(interaction["estimate"]),
                            "distance_x_selfing_se": float(interaction["cluster_robust_se"]),
                            "distance_x_selfing_p": float(interaction["p_value"]),
                        }
                    )

    result = pd.DataFrame(rows)
    if result.empty:
        return result
    fit_mask = result["status"].eq("fit")
    result["distance_q"] = np.nan
    result["distance_x_selfing_q"] = np.nan
    if fit_mask.any():
        grouping = ["context_layer", "source_mode", "stratum", "support_tier"]
        result.loc[fit_mask, "distance_q"] = (
            result.loc[fit_mask].groupby(grouping, group_keys=False)["distance_p"].apply(_bh)
        )
        result.loc[fit_mask, "distance_x_selfing_q"] = (
            result.loc[fit_mask]
            .groupby(grouping, group_keys=False)["distance_x_selfing_p"]
            .apply(_bh)
        )
    return result


def run_all_contexts(
    adjusted_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    realm_assignment: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> pd.DataFrame:
    regime_config = deepcopy(pattern_config)
    regime = run_interaction_analysis(
        adjusted_scores,
        covariates,
        regime_config,
        context_column=str(regime_config["context_column"]),
        contexts=[str(x) for x in regime_config["contexts"]],
        context_layer="analysis_regime",
    )

    realm_cov = covariates.merge(
        realm_assignment[["island_id", "biogeographic_realm"]].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    realm_contexts = sorted(
        realm_cov["biogeographic_realm"].dropna().astype(str).loc[lambda x: x.ne("")].unique()
    )
    realm = run_interaction_analysis(
        adjusted_scores,
        realm_cov,
        pattern_config,
        context_column="biogeographic_realm",
        contexts=realm_contexts,
        context_layer="biogeographic_realm",
    )
    return pd.concat([regime, realm], ignore_index=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-adjusted-scores-csv", type=Path, required=True)
    parser.add_argument("--covariates-csv", type=Path, required=True)
    parser.add_argument("--realm-assignment-csv", type=Path, required=True)
    parser.add_argument("--pattern-config-path", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    config = yaml.safe_load(args.pattern_config_path.read_text(encoding="utf-8"))
    result = run_all_contexts(
        pd.read_csv(args.source_adjusted_scores_csv),
        pd.read_csv(args.covariates_csv),
        pd.read_csv(args.realm_assignment_csv),
        config,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / "source_adjusted_selfing_interaction.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_source_adjusted_selfing_interaction_v1",
        "model": "attraction_shift ~ distance + selfing_core + distance:selfing_core + baseline_covariates",
        "attraction_shift": "(-large_bee_like + generalized_accessible) / 2",
        "selfing_core_is_continuous": True,
        "selfing_threshold_introduced": False,
        "causal_mediation_claimed": False,
        "interpretation": {
            "positive_interaction": "distance-associated attraction/accessibility shift strengthens as reproductive assurance increases",
            "null_interaction": "distance-associated attraction/accessibility shift is not detectably contingent on measured reproductive assurance",
            "negative_interaction": "distance-associated attraction/accessibility shift weakens as reproductive assurance increases",
        },
    }
    (args.output_dir / "selfing_interaction_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
