"""Post-hoc hinge diagnostic for regional pollination-syndrome nonlinearity.

This analysis is outside the frozen Chapter 1 submission inference. It uses the
already-frozen architecture|colour island counts and corrected geography to ask
whether focal syndrome-concordance signals change slope beyond the upper edge of
the four-region common 5--95% isolation support.

The hinge knot is not tuned to any trait result. It is the common-support upper
bound recorded in results/h1_regional_common_support_20261005/common_support_windows.csv.
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
from island_v2.chapter1_v13_raw_colour_coupling_audit import (
    fit_colour_conditioned_architecture_models,
)

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

KEYS = ["stratum", "context", "combination", "support_tier", "model"]


def _z(series: pd.Series) -> np.ndarray:
    x = pd.to_numeric(series, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _normal_two_sided_p(z_value: float) -> float:
    if not math.isfinite(z_value):
        return float("nan")
    return float(math.erfc(abs(z_value) / math.sqrt(2.0)))


def _load_windows(path: Path) -> tuple[float, float, float]:
    table = pd.read_csv(path)
    common = table.loc[table["window"].eq("common_05_95")]
    iqr = table.loc[table["window"].eq("common_iqr")]
    if len(common) != 1 or len(iqr) != 1:
        raise RuntimeError("expected one common_05_95 row and one common_iqr row")
    return (
        float(common.iloc[0]["lower_km"]),
        float(iqr.iloc[0]["upper_km"]),
        float(common.iloc[0]["upper_km"]),
    )


def _verify_full_replay(
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    reference: pd.DataFrame,
    scope: str,
) -> None:
    local = dict(cfg)
    local["strata"] = ["all_observed"]
    replay = fit_colour_conditioned_architecture_models(counts, scores, covariates, local)
    merged = replay.merge(
        reference,
        on=KEYS,
        suffixes=("_replay", "_reference"),
        validate="one_to_one",
    )
    if len(merged) != len(reference):
        raise RuntimeError(f"{scope}: replay row mismatch")
    for column in ["distance_estimate", "distance_se", "distance_p", "distance_q"]:
        a = pd.to_numeric(merged[f"{column}_replay"], errors="coerce").to_numpy(float)
        b = pd.to_numeric(merged[f"{column}_reference"], errors="coerce").to_numpy(float)
        if not np.allclose(a, b, atol=1e-6, rtol=1e-5, equal_nan=True):
            delta = float(np.nanmax(np.abs(a - b)))
            raise RuntimeError(f"{scope}: replay mismatch {column}: {delta}")


def _prepare_focus(
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    *,
    context: str,
    combination: str,
    lower_km: float,
) -> pd.DataFrame:
    geography = str(cfg["geography_column"])
    context_col = str(cfg["context_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]

    needed_cov = [
        "island_id",
        "distance_to_continent_km",
        geography,
        context_col,
        cluster,
        *baseline,
    ]
    work = (
        counts.loc[
            counts["stratum"].astype(str).eq("all_observed")
            & counts["combination"].astype(str).eq(combination)
        ]
        .copy()
        .merge(
            covariates[needed_cov].drop_duplicates("island_id"),
            on="island_id",
            how="left",
            validate="many_to_one",
        )
        .merge(
            _selfing_table(scores),
            on=["island_id", "stratum"],
            how="left",
            validate="many_to_one",
        )
    )
    for column in [
        "successes",
        "trials",
        "selfing_core",
        "distance_to_continent_km",
        geography,
        *baseline,
    ]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[context_col] = work[context_col].fillna("").astype(str)
    work[cluster] = work[cluster].fillna("").astype(str)
    work = work.loc[work[context_col].eq(context)].dropna(
        subset=[
            "successes",
            "trials",
            "selfing_core",
            "distance_to_continent_km",
            geography,
            *baseline,
            cluster,
        ]
    )
    work = work.loc[
        work["trials"].gt(0)
        & work["distance_to_continent_km"].ge(lower_km)
        & work[cluster].ne("")
    ].copy()
    return work


def _fit_models(
    work: pd.DataFrame,
    cfg: dict,
    *,
    knot_km: float,
) -> dict:
    geography = str(cfg["geography_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]

    log_knot = float(np.log1p(knot_km))
    x = pd.to_numeric(work[geography], errors="coerce").to_numpy(float)
    x_centered = x - log_knot
    hinge = np.maximum(x_centered, 0.0)

    common_names = ["intercept", "log_distance_centered"]
    common_cols = [np.ones(len(work), dtype=float), x_centered]
    for predictor in ["selfing_core", *baseline]:
        common_names.append(f"z_{predictor}")
        common_cols.append(_z(work[predictor]))

    y = work["successes"].to_numpy(float)
    n = work["trials"].to_numpy(float)
    clusters = work[cluster].astype(str).to_numpy()

    linear_X = np.column_stack(common_cols)
    linear_coef, linear_fit, linear_cov = _fit_grouped_binomial_design(
        y, n, linear_X, common_names, clusters
    )

    hinge_names = [*common_names, "hinge_above_common_upper"]
    hinge_X = np.column_stack([*common_cols, hinge])
    hinge_coef, hinge_fit, hinge_cov = _fit_grouped_binomial_design(
        y, n, hinge_X, hinge_names, clusters
    )

    li = linear_coef.set_index("predictor")
    hi = hinge_coef.set_index("predictor")
    below_beta = float(hi.loc["log_distance_centered", "estimate_log_odds"])
    below_se = float(hi.loc["log_distance_centered", "cluster_robust_se"])
    hinge_beta = float(hi.loc["hinge_above_common_upper", "estimate_log_odds"])
    hinge_se = float(hi.loc["hinge_above_common_upper", "cluster_robust_se"])

    b_idx = hinge_names.index("log_distance_centered")
    h_idx = hinge_names.index("hinge_above_common_upper")
    above_beta = below_beta + hinge_beta
    above_var = (
        hinge_cov[b_idx, b_idx]
        + hinge_cov[h_idx, h_idx]
        + 2.0 * hinge_cov[b_idx, h_idx]
    )
    above_se = float(np.sqrt(max(float(above_var), 0.0)))
    above_z = above_beta / above_se if above_se > 0 else float("nan")

    tail = work["distance_to_continent_km"].gt(knot_km)
    below = ~tail

    linear_beta = float(li.loc["log_distance_centered", "estimate_log_odds"])
    linear_se = float(li.loc["log_distance_centered", "cluster_robust_se"])

    return {
        "n_islands": int(work["island_id"].nunique()),
        "n_clusters": int(work[cluster].nunique()),
        "n_below_or_at_knot": int(work.loc[below, "island_id"].nunique()),
        "n_above_knot": int(work.loc[tail, "island_id"].nunique()),
        "n_blocks_below_or_at_knot": int(work.loc[below, cluster].nunique()),
        "n_blocks_above_knot": int(work.loc[tail, cluster].nunique()),
        "linear_status": linear_fit["status"],
        "linear_slope_per_log1p_km": linear_beta,
        "linear_se": linear_se,
        "linear_p": float(li.loc["log_distance_centered", "p_value"]),
        "linear_aic": float(linear_fit["aic"]),
        "hinge_status": hinge_fit["status"],
        "below_slope_per_log1p_km": below_beta,
        "below_se": below_se,
        "below_p": float(hi.loc["log_distance_centered", "p_value"]),
        "hinge_change": hinge_beta,
        "hinge_change_se": hinge_se,
        "hinge_change_p": float(hi.loc["hinge_above_common_upper", "p_value"]),
        "above_slope_per_log1p_km": above_beta,
        "above_se": above_se,
        "above_p": _normal_two_sided_p(above_z),
        "hinge_aic": float(hinge_fit["aic"]),
        "delta_aic_hinge_minus_linear": float(hinge_fit["aic"] - linear_fit["aic"]),
    }


def _fit_three_segment(
    work: pd.DataFrame,
    cfg: dict,
    *,
    first_knot_km: float,
    second_knot_km: float,
) -> dict:
    geography = str(cfg["geography_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]

    log_first = float(np.log1p(first_knot_km))
    log_second = float(np.log1p(second_knot_km))
    x = pd.to_numeric(work[geography], errors="coerce").to_numpy(float)
    x_centered = x - log_first
    hinge_first = np.maximum(x - log_first, 0.0)
    hinge_second = np.maximum(x - log_second, 0.0)

    names = ["intercept", "log_distance_centered"]
    cols = [np.ones(len(work), dtype=float), x_centered]
    for predictor in ["selfing_core", *baseline]:
        names.append(f"z_{predictor}")
        cols.append(_z(work[predictor]))
    names.extend(["hinge_above_iqr", "hinge_above_common_upper"])
    cols.extend([hinge_first, hinge_second])

    coef, fit, covariance = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(cols),
        names,
        work[cluster].astype(str).to_numpy(),
    )
    indexed = coef.set_index("predictor")

    b = float(indexed.loc["log_distance_centered", "estimate_log_odds"])
    h1 = float(indexed.loc["hinge_above_iqr", "estimate_log_odds"])
    h2 = float(indexed.loc["hinge_above_common_upper", "estimate_log_odds"])
    i_b = names.index("log_distance_centered")
    i_h1 = names.index("hinge_above_iqr")
    i_h2 = names.index("hinge_above_common_upper")

    middle = b + h1
    middle_var = (
        covariance[i_b, i_b]
        + covariance[i_h1, i_h1]
        + 2.0 * covariance[i_b, i_h1]
    )
    middle_se = float(np.sqrt(max(float(middle_var), 0.0)))
    middle_p = _normal_two_sided_p(middle / middle_se if middle_se > 0 else float("nan"))

    far = b + h1 + h2
    idx = [i_b, i_h1, i_h2]
    far_var = float(covariance[np.ix_(idx, idx)].sum())
    far_se = float(np.sqrt(max(far_var, 0.0)))
    far_p = _normal_two_sided_p(far / far_se if far_se > 0 else float("nan"))

    d = work["distance_to_continent_km"]
    seg1 = d.le(first_knot_km)
    seg2 = d.gt(first_knot_km) & d.le(second_knot_km)
    seg3 = d.gt(second_knot_km)

    return {
        "three_segment_status": fit["status"],
        "three_segment_aic": float(fit["aic"]),
        "n_segment_le_iqr": int(work.loc[seg1, "island_id"].nunique()),
        "n_segment_iqr_to_common_upper": int(work.loc[seg2, "island_id"].nunique()),
        "n_segment_above_common_upper": int(work.loc[seg3, "island_id"].nunique()),
        "n_blocks_segment_le_iqr": int(work.loc[seg1, cluster].nunique()),
        "n_blocks_segment_iqr_to_common_upper": int(work.loc[seg2, cluster].nunique()),
        "n_blocks_segment_above_common_upper": int(work.loc[seg3, cluster].nunique()),
        "segment1_slope": b,
        "segment1_p": float(indexed.loc["log_distance_centered", "p_value"]),
        "change_at_iqr": h1,
        "change_at_iqr_p": float(indexed.loc["hinge_above_iqr", "p_value"]),
        "segment2_slope": middle,
        "segment2_p": middle_p,
        "change_at_common_upper": h2,
        "change_at_common_upper_p": float(
            indexed.loc["hinge_above_common_upper", "p_value"]
        ),
        "segment3_slope": far,
        "segment3_p": far_p,
    }


def _reference_focus(reference: pd.DataFrame, context: str, combination: str) -> dict:
    row = reference.loc[
        reference["stratum"].astype(str).eq("all_observed")
        & reference["context"].astype(str).eq(context)
        & reference["combination"].astype(str).eq(combination)
        & reference["support_tier"].astype(str).eq("confirmatory")
        & reference["model"].astype(str).eq("conditional_selfing")
    ]
    if len(row) != 1:
        raise RuntimeError(f"reference focus row missing/duplicated: {context} {combination}")
    row = row.iloc[0]
    return {
        "reference_standardized_distance_estimate": float(row["distance_estimate"]),
        "reference_distance_p": float(row["distance_p"]),
        "reference_distance_q": float(row["distance_q"]),
        "reference_n_islands": int(row["n_unique_islands"]),
        "reference_n_clusters": int(row["n_clusters"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--all-counts", type=Path, required=True)
    parser.add_argument("--direct-counts", type=Path, required=True)
    parser.add_argument("--all-scores", type=Path, required=True)
    parser.add_argument("--direct-scores", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--reference-all", type=Path, required=True)
    parser.add_argument("--reference-direct", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    lower_km, iqr_upper_km, knot_km = _load_windows(args.windows)

    rows = []
    for scope, count_path, score_path, reference_path in [
        ("all", args.all_counts, args.all_scores, args.reference_all),
        ("direct", args.direct_counts, args.direct_scores, args.reference_direct),
    ]:
        counts = pd.read_csv(count_path, dtype={"island_id": str})
        scores = pd.read_csv(score_path, dtype={"island_id": str})
        reference = pd.read_csv(reference_path)
        _verify_full_replay(counts, scores, cov, cfg, reference, scope)

        for signal, (context, combination) in FOCUS.items():
            work = _prepare_focus(
                counts,
                scores,
                cov,
                cfg,
                context=context,
                combination=combination,
                lower_km=lower_km,
            )
            if len(work) < 50:
                rows.append(
                    {
                        "signal": signal,
                        "evidence_scope": scope,
                        "context": context,
                        "combination": combination,
                        "status": "fewer_than_50_complete_islands",
                        "n_islands": int(work["island_id"].nunique()),
                    }
                )
                continue
            result = _fit_models(work, cfg, knot_km=knot_km)
            three = _fit_three_segment(
                work,
                cfg,
                first_knot_km=iqr_upper_km,
                second_knot_km=knot_km,
            )
            ref = _reference_focus(reference, context, combination)
            rows.append(
                {
                    "signal": signal,
                    "evidence_scope": scope,
                    "context": context,
                    "combination": combination,
                    "status": "fit",
                    "common_lower_km": lower_km,
                    "first_knot_iqr_upper_km": iqr_upper_km,
                    "hinge_knot_km": knot_km,
                    **ref,
                    **result,
                    **three,
                }
            )

    out = pd.DataFrame(rows)
    args.output.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output / "syndrome_hinge_models.csv", index=False)

    summary_rows = []
    for signal, part in out.loc[out["status"].eq("fit")].groupby("signal", sort=False):
        for row in part.itertuples(index=False):
            summary_rows.append(
                {
                    "signal": signal,
                    "evidence_scope": row.evidence_scope,
                    "context": row.context,
                    "n_islands": row.n_islands,
                    "n_above_knot": row.n_above_knot,
                    "n_blocks_above_knot": row.n_blocks_above_knot,
                    "below_slope": row.below_slope_per_log1p_km,
                    "below_p": row.below_p,
                    "hinge_change": row.hinge_change,
                    "hinge_change_p": row.hinge_change_p,
                    "above_slope": row.above_slope_per_log1p_km,
                    "above_p": row.above_p,
                    "delta_aic_hinge_minus_linear": row.delta_aic_hinge_minus_linear,
                    "delta_aic_three_segment_minus_linear": row.three_segment_aic - row.linear_aic,
                    "segment1_slope": row.segment1_slope,
                    "segment1_p": row.segment1_p,
                    "change_at_iqr": row.change_at_iqr,
                    "change_at_iqr_p": row.change_at_iqr_p,
                    "segment2_slope": row.segment2_slope,
                    "segment2_p": row.segment2_p,
                    "change_at_common_upper": row.change_at_common_upper,
                    "change_at_common_upper_p": row.change_at_common_upper_p,
                    "segment3_slope": row.segment3_slope,
                    "segment3_p": row.segment3_p,
                }
            )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(args.output / "syndrome_hinge_summary.csv", index=False)

    manifest = {
        "contract": "chapter1_pollination_syndrome_hinge_diagnostic_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "lower_bound_km": lower_km,
        "first_knot_iqr_upper_km": iqr_upper_km,
        "hinge_knot_km": knot_km,
        "knot_definition": (
            "support-derived knots: upper bound of four-region common IQR and "
            "upper bound of four-region common 5-95% isolation support"
        ),
        "full_corrected_replay_verified": True,
        "model": (
            "grouped-binomial logit with cluster-robust covariance; selfing_core + "
            "area + climate adjusted; raw log(1+distance) slope plus hinge above knot"
        ),
        "interpretation": (
            "hinge_change tests a slope change beyond the shared regional distance range. "
            "This is exploratory evidence about exposure nonlinearity, not a causal pollinator test."
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print("=== POLLINATION-SYNDROME HINGE SUMMARY ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
