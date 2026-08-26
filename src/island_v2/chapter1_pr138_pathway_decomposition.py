"""Conditional decomposition of PR138 attraction and selfing syndrome responses.

This analysis distinguishes two correlated responses without claiming causal mediation:

1. reproductive selfing core (SC / mating system / autonomous selfing) versus distance;
2. floral attraction/accessibility shift versus distance;
3. the same attraction shift versus distance after conditioning on selfing core.

If the distance association of attraction shift persists after selfing-core adjustment,
the floral/pollination-syndrome pattern is not reducible to the measured reproductive
selfing core. Independent pollinator data are still required for causal attribution.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import typer
import yaml

from island_v2.chapter1_pr136_biogeographic_residual import _fit_weighted_clustered_design

app = typer.Typer(add_completion=False, no_args_is_help=True)


def _standardize(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


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


def build_pathway_table(island_scores: pd.DataFrame) -> pd.DataFrame:
    required = {"island_id", "stratum", "syndrome", "syndrome_score"}
    missing = required - set(island_scores.columns)
    if missing:
        raise typer.BadParameter(f"island syndrome scores missing columns: {sorted(missing)}")
    pivot = (
        island_scores.loc[
            island_scores["syndrome"].isin(
                ["large_bee_like", "generalized_accessible", "selfing_core"]
            ),
            ["island_id", "stratum", "syndrome", "syndrome_score"],
        ]
        .drop_duplicates(["island_id", "stratum", "syndrome"])
        .pivot(index=["island_id", "stratum"], columns="syndrome", values="syndrome_score")
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


def _fit_one(
    work: pd.DataFrame,
    *,
    outcome: str,
    predictors: list[str],
    cluster_column: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    complete = work.dropna(subset=[outcome, *predictors, cluster_column]).copy()
    if complete.empty:
        return pd.DataFrame(), {"status": "no_complete_rows"}
    names = ["intercept", *[f"z_{x}" for x in predictors]]
    columns = [np.ones(len(complete), dtype=float)]
    for predictor in predictors:
        columns.append(_standardize(complete[predictor]))
    coef, _, fit = _fit_weighted_clustered_design(
        pd.to_numeric(complete[outcome], errors="coerce").to_numpy(float),
        np.ones(len(complete), dtype=float),
        np.column_stack(columns),
        names,
        complete[cluster_column].astype(str).to_numpy(),
    )
    return coef, {
        **fit,
        "n_unique_islands": int(complete["island_id"].nunique()),
    }


def run_pathway_decomposition(
    island_scores: pd.DataFrame,
    covariates: pd.DataFrame,
    pattern_config: dict[str, Any],
) -> pd.DataFrame:
    geography = str(pattern_config["geography_column"])
    context_column = str(pattern_config["context_column"])
    cluster = str(pattern_config["cluster_column"])
    baseline = [str(x) for x in pattern_config["baseline_covariates"]]
    contexts = [str(x) for x in pattern_config["contexts"]]
    strata = [str(x) for x in pattern_config["strata"]]
    tiers = {str(k): int(v) for k, v in pattern_config["support_tiers"].items()}

    table = build_pathway_table(island_scores)
    needed_cov = ["island_id", geography, context_column, cluster, *baseline]
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
    model_specs = {
        "selfing_core_distance": ("selfing_core", [geography, *baseline]),
        "attraction_unconditional": ("attraction_shift", [geography, *baseline]),
        "attraction_conditional_on_selfing_core": (
            "attraction_shift",
            [geography, "selfing_core", *baseline],
        ),
    }
    for stratum in strata:
        for context in contexts:
            base = data.loc[
                data["stratum"].eq(stratum) & data[context_column].eq(context)
            ].copy()
            for tier, threshold in tiers.items():
                for model, (outcome, predictors) in model_specs.items():
                    complete = base.dropna(subset=[outcome, *predictors, cluster]).copy()
                    if complete["island_id"].nunique() < threshold:
                        rows.append(
                            {
                                "stratum": stratum,
                                "context": context,
                                "support_tier": tier,
                                "threshold": threshold,
                                "model": model,
                                "status": "not_testable",
                                "n_unique_islands": int(complete["island_id"].nunique()),
                            }
                        )
                        continue
                    coef, fit = _fit_one(
                        complete,
                        outcome=outcome,
                        predictors=predictors,
                        cluster_column=cluster,
                    )
                    if coef.empty:
                        rows.append(
                            {
                                "stratum": stratum,
                                "context": context,
                                "support_tier": tier,
                                "threshold": threshold,
                                "model": model,
                                "status": str(fit.get("status", "fit_failed")),
                                "n_unique_islands": int(fit.get("n_unique_islands", 0)),
                            }
                        )
                        continue
                    indexed = coef.set_index("predictor")
                    geo = indexed.loc[f"z_{geography}"]
                    row: dict[str, Any] = {
                        "stratum": stratum,
                        "context": context,
                        "support_tier": tier,
                        "threshold": threshold,
                        "model": model,
                        "status": "fit",
                        "n_unique_islands": int(fit["n_unique_islands"]),
                        "n_clusters": int(fit["n_clusters"]),
                        "distance_estimate": float(geo["estimate"]),
                        "distance_se": float(geo["cluster_robust_se"]),
                        "distance_p": float(geo["p_value"]),
                    }
                    if "selfing_core" in predictors:
                        selfing = indexed.loc["z_selfing_core"]
                        row.update(
                            {
                                "selfing_core_estimate": float(selfing["estimate"]),
                                "selfing_core_se": float(selfing["cluster_robust_se"]),
                                "selfing_core_p": float(selfing["p_value"]),
                            }
                        )
                    rows.append(row)
    result = pd.DataFrame(rows)
    if not result.empty and "distance_p" in result.columns:
        fit_mask = result["status"].eq("fit")
        result["distance_q"] = np.nan
        result.loc[fit_mask, "distance_q"] = (
            result.loc[fit_mask]
            .groupby(["stratum", "support_tier", "model"], group_keys=False)["distance_p"]
            .apply(_bh)
        )
    return result


@app.command("run")
def run(
    island_scores_csv: Path = typer.Option(..., exists=True),
    covariates_csv: Path = typer.Option(..., exists=True),
    pattern_config_path: Path = typer.Option(..., exists=True),
    output_dir: Path = typer.Option(...),
) -> None:
    pattern_config = yaml.safe_load(pattern_config_path.read_text(encoding="utf-8"))
    result = run_pathway_decomposition(
        pd.read_csv(island_scores_csv),
        pd.read_csv(covariates_csv),
        pattern_config,
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_dir / "attraction_selfing_pathway_decomposition.csv", index=False)
    manifest = {
        "contract": "chapter1_pr138_pathway_decomposition_v1",
        "attraction_shift": "(-large_bee_like + generalized_accessible) / 2",
        "selfing_core_excludes_flower_size": True,
        "causal_mediation_claimed": False,
    }
    (output_dir / "pathway_decomposition_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    typer.echo(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    app()
