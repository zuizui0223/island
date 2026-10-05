"""Post-hoc hinge diagnostic for regional pollination-syndrome nonlinearity.

This diagnostic asks whether the corrected architecture|colour isolation associations
change after the upper edge of four-region common isolation support. It is not part of
the frozen Chapter 1 submission inference.

For each evidence scope, region, and frozen colour-architecture combination, fit:

  linear:
    architecture ~ log(1 + distance) + selfing_core + area + climate

  hinge:
    architecture ~ log(1 + distance)
                 + max(0, log(1 + distance) - log(1 + hinge_km))
                 + selfing_core + area + climate

The coefficient on log(1 + distance) is the slope below the hinge. The hinge
coefficient is the change in slope above the hinge; their sum is the post-hinge slope.
Cluster-robust covariance follows the frozen raw colour audit implementation.

The hinge is fixed from the four-region 5th--95th percentile common-support upper
edge and is never tuned to syndrome results.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design
from island_v2.chapter1_v13_raw_colour_audit import _selfing_table
from island_v2.chapter1_v13_raw_colour_coupling_audit import COUPLING_SPECS


FOCUS = {
    "north_mid_yellow_large_bee_form": (
        "northern_midlatitude",
        "yellow_orange__large_bee_form_given_colour",
    ),
    "north_high_blue_butterfly_form": (
        "northern_high_latitude",
        "blue_purple__butterfly_form_given_colour",
    ),
    "north_high_blue_large_bee_deep": (
        "northern_high_latitude",
        "blue_purple__large_bee_deep_tube_given_colour",
    ),
    "tropical_yellow_butterfly_deep": (
        "tropical",
        "yellow_orange__butterfly_deep_tube_given_colour",
    ),
    "south_yellow_bird_form": (
        "southern_extratropical",
        "yellow_orange__bird_form_given_colour",
    ),
    "south_yellow_bird_deep": (
        "southern_extratropical",
        "yellow_orange__bird_deep_tube_given_colour",
    ),
}


def _p_two_sided(z: float) -> float:
    return math.erfc(abs(float(z)) / math.sqrt(2.0)) if math.isfinite(z) else float("nan")


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
    adjusted = np.minimum.accumulate(
        (ranked * n / np.arange(1, n + 1))[::-1]
    )[::-1]
    restored = np.empty(n, dtype=float)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out.loc[ok] = restored
    return out


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _load_hinge_km(path: Path) -> float:
    windows = pd.read_csv(path)
    row = windows.loc[windows["window"].astype(str).eq("common_05_95")]
    if len(row) != 1:
        raise RuntimeError("expected exactly one common_05_95 window")
    hinge = float(row.iloc[0]["upper_km"])
    if not np.isfinite(hinge) or hinge <= 0:
        raise RuntimeError("invalid hinge distance")
    return hinge


def _assemble_data(
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    geography = str(cfg["geography_column"])
    context = str(cfg["context_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]
    needed = [
        "island_id",
        "distance_to_continent_km",
        geography,
        context,
        cluster,
        *baseline,
    ]
    missing = set(needed) - set(covariates.columns)
    if missing:
        raise RuntimeError(f"covariates missing columns: {sorted(missing)}")
    work = counts.loc[
        counts["stratum"].astype(str).eq("all_observed")
    ].copy()
    work["island_id"] = work["island_id"].astype(str)
    work = work.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    ).merge(
        _selfing_table(scores),
        on=["island_id", "stratum"],
        how="left",
        validate="many_to_one",
    )
    for column in [
        "successes",
        "trials",
        "distance_to_continent_km",
        geography,
        "selfing_core",
        *baseline,
    ]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[context] = work[context].fillna("").astype(str)
    work[cluster] = work[cluster].fillna("").astype(str)
    return work


def _fit_one(
    part: pd.DataFrame,
    *,
    cfg: dict,
    hinge_km: float,
) -> dict[str, object]:
    geography = str(cfg["geography_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]
    required = [
        "successes",
        "trials",
        "distance_to_continent_km",
        geography,
        "selfing_core",
        cluster,
        *baseline,
    ]
    complete = part.dropna(subset=required).copy()
    complete = complete.loc[
        complete["trials"].gt(0)
        & complete["successes"].ge(0)
        & complete["successes"].le(complete["trials"])
        & complete[cluster].ne("")
    ].copy()

    n_islands = int(complete["island_id"].nunique())
    n_clusters = int(complete[cluster].nunique())
    minimum = int(cfg["support_tiers"]["confirmatory"])
    if n_islands < minimum:
        return {
            "status": "not_testable",
            "n_islands": n_islands,
            "n_clusters": n_clusters,
        }

    x = pd.to_numeric(complete[geography], errors="coerce").to_numpy(float)
    hinge_log = float(np.log1p(hinge_km))
    tail = np.maximum(0.0, x - hinge_log)
    tail_mask = pd.to_numeric(
        complete["distance_to_continent_km"], errors="coerce"
    ).to_numpy(float) > hinge_km

    control_names = ["z_selfing_core", *[f"z_{x}" for x in baseline]]
    controls = [_z(complete["selfing_core"])]
    controls.extend(_z(complete[x]) for x in baseline)

    linear_names = ["intercept", "log1p_distance", *control_names]
    linear_design = np.column_stack(
        [np.ones(len(complete), dtype=float), x, *controls]
    )
    hinge_names = ["intercept", "log1p_distance", "tail_hinge", *control_names]
    hinge_design = np.column_stack(
        [np.ones(len(complete), dtype=float), x, tail, *controls]
    )

    if np.linalg.matrix_rank(linear_design) != linear_design.shape[1]:
        return {
            "status": "linear_rank_deficient",
            "n_islands": n_islands,
            "n_clusters": n_clusters,
        }
    if np.linalg.matrix_rank(hinge_design) != hinge_design.shape[1]:
        return {
            "status": "hinge_rank_deficient",
            "n_islands": n_islands,
            "n_clusters": n_clusters,
        }

    y = complete["successes"].to_numpy(float)
    n = complete["trials"].to_numpy(float)
    clusters = complete[cluster].astype(str).to_numpy()

    linear_coef, linear_fit, _ = _fit_grouped_binomial_design(
        y, n, linear_design, linear_names, clusters, max_iter=500
    )
    hinge_coef, hinge_fit, hinge_cov = _fit_grouped_binomial_design(
        y, n, hinge_design, hinge_names, clusters, max_iter=500
    )

    if linear_fit["status"] != "fit" or hinge_fit["status"] != "fit":
        return {
            "status": "fit_failure",
            "linear_status": linear_fit["status"],
            "hinge_status": hinge_fit["status"],
            "n_islands": n_islands,
            "n_clusters": n_clusters,
            "n_tail_islands": int(np.sum(tail_mask)),
            "n_tail_clusters": int(
                complete.loc[tail_mask, cluster].nunique()
            ),
        }

    hi = hinge_coef.set_index("predictor")
    pre = hi.loc["log1p_distance"]
    change = hi.loc["tail_hinge"]
    i_pre = hinge_names.index("log1p_distance")
    i_change = hinge_names.index("tail_hinge")
    post_estimate = float(
        pre["estimate_log_odds"] + change["estimate_log_odds"]
    )
    post_var = float(
        hinge_cov[i_pre, i_pre]
        + hinge_cov[i_change, i_change]
        + 2.0 * hinge_cov[i_pre, i_change]
    )
    post_se = float(np.sqrt(max(post_var, 0.0)))
    post_z = post_estimate / post_se if post_se > 0 else float("nan")

    lin = linear_coef.set_index("predictor").loc["log1p_distance"]

    return {
        "status": "fit",
        "n_islands": n_islands,
        "n_clusters": n_clusters,
        "n_tail_islands": int(np.sum(tail_mask)),
        "n_tail_clusters": int(complete.loc[tail_mask, cluster].nunique()),
        "tail_island_fraction": float(np.mean(tail_mask)),
        "hinge_km": float(hinge_km),
        "linear_estimate": float(lin["estimate_log_odds"]),
        "linear_se": float(lin["cluster_robust_se"]),
        "linear_p": float(lin["p_value"]),
        "pre_hinge_estimate": float(pre["estimate_log_odds"]),
        "pre_hinge_se": float(pre["cluster_robust_se"]),
        "pre_hinge_p": float(pre["p_value"]),
        "hinge_change_estimate": float(change["estimate_log_odds"]),
        "hinge_change_se": float(change["cluster_robust_se"]),
        "hinge_change_p": float(change["p_value"]),
        "post_hinge_estimate": post_estimate,
        "post_hinge_se": post_se,
        "post_hinge_p": _p_two_sided(post_z),
        "linear_aic": float(linear_fit["aic"]),
        "hinge_aic": float(hinge_fit["aic"]),
        "delta_aic_hinge_minus_linear": float(
            hinge_fit["aic"] - linear_fit["aic"]
        ),
    }


def _run_scope(
    scope: str,
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    hinge_km: float,
) -> pd.DataFrame:
    data = _assemble_data(counts, scores, covariates, cfg)
    context_col = str(cfg["context_column"])
    rows = []
    for context in [str(x) for x in cfg["contexts"]]:
        for combination in COUPLING_SPECS:
            part = data.loc[
                data[context_col].eq(context)
                & data["combination"].astype(str).eq(combination)
            ].copy()
            rows.append(
                {
                    "evidence_scope": scope,
                    "context": context,
                    "combination": combination,
                    "architecture_label": COUPLING_SPECS[combination]["architecture_label"],
                    "trait_name": COUPLING_SPECS[combination]["trait_name"],
                    **_fit_one(part, cfg=cfg, hinge_km=hinge_km),
                }
            )
    result = pd.DataFrame(rows)
    result["hinge_change_q"] = np.nan
    result["post_hinge_q"] = np.nan
    fit = result["status"].eq("fit")
    for _, index in result.loc[fit].groupby(
        ["evidence_scope", "context"]
    ).groups.items():
        result.loc[index, "hinge_change_q"] = _bh(
            result.loc[index, "hinge_change_p"]
        )
        result.loc[index, "post_hinge_q"] = _bh(
            result.loc[index, "post_hinge_p"]
        )
    return result


def _focus(table: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, (context, combination) in FOCUS.items():
        part = table.loc[
            table["context"].eq(context)
            & table["combination"].eq(combination)
        ]
        for row in part.to_dict("records"):
            rows.append({"signal": label, **row})
    return pd.DataFrame(rows)


def _summary(table: pd.DataFrame) -> pd.DataFrame:
    fit = table.loc[table["status"].eq("fit")].copy()
    rows = []
    for (scope, context), part in fit.groupby(
        ["evidence_scope", "context"], sort=False
    ):
        rows.append(
            {
                "evidence_scope": scope,
                "context": context,
                "n_combinations_fit": int(len(part)),
                "n_hinge_change_p_lt_0_05": int(part["hinge_change_p"].lt(0.05).sum()),
                "n_hinge_change_q_lt_0_05": int(part["hinge_change_q"].lt(0.05).sum()),
                "n_post_hinge_p_lt_0_05": int(part["post_hinge_p"].lt(0.05).sum()),
                "n_post_hinge_q_lt_0_05": int(part["post_hinge_q"].lt(0.05).sum()),
                "n_hinge_aic_better_by_2": int(
                    part["delta_aic_hinge_minus_linear"].lt(-2).sum()
                ),
                "median_tail_fraction": float(part["tail_island_fraction"].median()),
                "median_tail_clusters": float(part["n_tail_clusters"].median()),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--all-counts", type=Path, required=True)
    parser.add_argument("--direct-counts", type=Path, required=True)
    parser.add_argument("--all-scores", type=Path, required=True)
    parser.add_argument("--direct-scores", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    hinge_km = _load_hinge_km(args.windows)

    parts = []
    for scope, count_path, score_path in [
        ("all", args.all_counts, args.all_scores),
        ("direct", args.direct_counts, args.direct_scores),
    ]:
        parts.append(
            _run_scope(
                scope,
                pd.read_csv(count_path, dtype={"island_id": str}),
                pd.read_csv(score_path, dtype={"island_id": str}),
                cov,
                cfg,
                hinge_km,
            )
        )
    result = pd.concat(parts, ignore_index=True)
    focus = _focus(result)
    summary = _summary(result)

    args.output.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output / "syndrome_hinge_models.csv", index=False)
    focus.to_csv(args.output / "syndrome_hinge_focus_signals.csv", index=False)
    summary.to_csv(args.output / "syndrome_hinge_summary.csv", index=False)
    manifest = {
        "contract": "chapter1_pollination_syndrome_hinge_diagnostic_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "hinge_km": hinge_km,
        "hinge_source": (
            "upper bound of four-region common 5th-95th percentile support"
        ),
        "model": (
            "continuous linear spline in log1p mainland distance, conditional on "
            "selfing_core, island area, climate PC1-PC4"
        ),
        "inference": (
            "cluster-robust normal-reference Wald diagnostics; BH within "
            "evidence-scope x region across frozen coupling combinations"
        ),
        "claim_boundary": (
            "A hinge pattern diagnoses distance-regime dependence and possible "
            "nonlinearity; it does not identify realized pollinators or establish "
            "historical causal replacement."
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print("=== HINGE SUMMARY ===")
    print(summary.to_string(index=False))
    print("\n=== FOCUS SIGNALS ===")
    show = [
        "signal",
        "evidence_scope",
        "context",
        "combination",
        "n_islands",
        "n_clusters",
        "n_tail_islands",
        "n_tail_clusters",
        "linear_estimate",
        "linear_p",
        "pre_hinge_estimate",
        "pre_hinge_p",
        "hinge_change_estimate",
        "hinge_change_p",
        "hinge_change_q",
        "post_hinge_estimate",
        "post_hinge_p",
        "post_hinge_q",
        "delta_aic_hinge_minus_linear",
    ]
    print(focus[show].to_string(index=False))


if __name__ == "__main__":
    main()
