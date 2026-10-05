"""Post-hoc regional common-support diagnostic for final Chapter 1 H1.

This diagnostic does not modify the frozen submission estimand. It asks whether the
four regional H1 slopes depend strongly on different mainland-distance support.

Windows are derived from the PRIMARY all-analysis H1 union, not tuned to trait results:
- common_05_95 = intersection of the four regional 5th--95th percentile ranges;
- common_iqr   = intersection of the four regional interquartile ranges.

The same frozen beta-binomial traitwise model is then refit within each window.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_all_data_probability import _prepare, build_broad_counts
from island_v2.chapter1_h1_traitwise import fit_trait


def _prepare_scope(
    flora: pd.DataFrame,
    audit: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    local = dict(cfg, strata=["all_observed"])
    counts = build_broad_counts(flora, audit, local)
    prepared = _prepare(counts, covariates, local)
    raw = covariates[["island_id", "distance_to_continent_km"]].drop_duplicates("island_id")
    prepared = prepared.merge(raw, on="island_id", how="left", validate="many_to_one")
    prepared["distance_to_continent_km"] = pd.to_numeric(
        prepared["distance_to_continent_km"], errors="coerce"
    )
    return prepared.dropna(subset=["distance_to_continent_km"]).copy()


def _regional_quantiles(prepared: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    context_col = cfg["context_column"]
    rows = []
    for context in cfg["contexts"]:
        part = (
            prepared.loc[
                prepared[context_col].eq(context),
                ["island_id", "distance_to_continent_km"],
            ]
            .drop_duplicates("island_id")
        )
        values = part["distance_to_continent_km"].to_numpy(float)
        rows.append(
            {
                "context": context,
                "n_islands": int(len(part)),
                "q05_km": float(np.quantile(values, 0.05)),
                "q25_km": float(np.quantile(values, 0.25)),
                "median_km": float(np.quantile(values, 0.50)),
                "q75_km": float(np.quantile(values, 0.75)),
                "q95_km": float(np.quantile(values, 0.95)),
            }
        )
    return pd.DataFrame(rows)


def _window_definitions(quantiles: pd.DataFrame) -> list[dict]:
    common_05_95 = {
        "window": "common_05_95",
        "lower_km": float(quantiles["q05_km"].max()),
        "upper_km": float(quantiles["q95_km"].min()),
        "definition": "intersection of regional 5th-95th percentile distance ranges",
    }
    common_iqr = {
        "window": "common_iqr",
        "lower_km": float(quantiles["q25_km"].max()),
        "upper_km": float(quantiles["q75_km"].min()),
        "definition": "intersection of regional interquartile distance ranges",
    }
    if common_05_95["lower_km"] >= common_05_95["upper_km"]:
        raise RuntimeError("no common 5-95 support across regions")
    if common_iqr["lower_km"] >= common_iqr["upper_km"]:
        raise RuntimeError("no common IQR support across regions")
    return [
        {
            "window": "full",
            "lower_km": None,
            "upper_km": None,
            "definition": "full regional support",
        },
        common_05_95,
        common_iqr,
    ]


def _mask_window(frame: pd.DataFrame, lower: float | None, upper: float | None) -> pd.Series:
    mask = pd.Series(True, index=frame.index)
    if lower is not None:
        mask &= frame["distance_to_continent_km"].ge(float(lower))
    if upper is not None:
        mask &= frame["distance_to_continent_km"].le(float(upper))
    return mask


def _fit_scope(
    prepared: pd.DataFrame,
    cfg: dict,
    scope: str,
    windows: list[dict],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    context_col = cfg["context_column"]
    cluster_col = cfg["cluster_column"]
    geo_col = cfg["geography_column"]
    fit_rows = []
    support_rows = []

    full_union = {
        context: set(
            prepared.loc[prepared[context_col].eq(context), "island_id"].astype(str).unique()
        )
        for context in cfg["contexts"]
    }

    for window in windows:
        mask = _mask_window(prepared, window["lower_km"], window["upper_km"])
        sliced = prepared.loc[mask].copy()
        for context in cfg["contexts"]:
            union = (
                sliced.loc[sliced[context_col].eq(context)]
                [["island_id", cluster_col, "distance_to_continent_km"]]
                .drop_duplicates("island_id")
            )
            support_rows.append(
                {
                    "evidence_scope": scope,
                    "window": window["window"],
                    "context": context,
                    "lower_km": window["lower_km"],
                    "upper_km": window["upper_km"],
                    "n_union_islands": int(len(union)),
                    "n_spatial_blocks": int(union[cluster_col].nunique()),
                    "retained_union_fraction": (
                        float(len(union) / len(full_union[context]))
                        if full_union[context]
                        else np.nan
                    ),
                    "distance_median_km": (
                        float(union["distance_to_continent_km"].median())
                        if len(union)
                        else np.nan
                    ),
                }
            )
            for outcome in cfg["model_outcomes"]:
                part = sliced.loc[
                    sliced[context_col].eq(context) & sliced["outcome"].eq(outcome)
                ].copy()
                result = fit_trait(part, outcome, cfg)
                log_sd = (
                    float(np.std(part[geo_col].to_numpy(float), ddof=0))
                    if len(part) > 1
                    else np.nan
                )
                estimate = result.get("estimate", np.nan)
                se = result.get("se", np.nan)
                raw_est = (
                    float(estimate / log_sd)
                    if np.isfinite(estimate) and np.isfinite(log_sd) and log_sd > 0
                    else np.nan
                )
                raw_se = (
                    float(se / log_sd)
                    if np.isfinite(se) and np.isfinite(log_sd) and log_sd > 0
                    else np.nan
                )
                fit_rows.append(
                    {
                        "evidence_scope": scope,
                        "window": window["window"],
                        "context": context,
                        "outcome": outcome,
                        "lower_km": window["lower_km"],
                        "upper_km": window["upper_km"],
                        "log_distance_sd": log_sd,
                        "estimate_per_window_sd": estimate,
                        "se_per_window_sd": se,
                        "estimate_per_log1p_km": raw_est,
                        "se_per_log1p_km": raw_se,
                        "p_two_sided": result.get("p_two_sided", np.nan),
                        "n_islands": int(result.get("n_islands", len(part))),
                        "n_spatial_blocks": int(part[cluster_col].nunique()),
                        "status": result.get("status", "unknown"),
                        "optimizer_success": result.get("optimizer_success", np.nan),
                    }
                )

    fits = pd.DataFrame(fit_rows)
    support = pd.DataFrame(support_rows)

    full = fits.loc[fits["window"].eq("full"), [
        "evidence_scope", "context", "outcome",
        "estimate_per_log1p_km", "p_two_sided", "status"
    ]].rename(
        columns={
            "estimate_per_log1p_km": "full_estimate_per_log1p_km",
            "p_two_sided": "full_p_two_sided",
            "status": "full_status",
        }
    )
    fits = fits.merge(
        full,
        on=["evidence_scope", "context", "outcome"],
        how="left",
        validate="many_to_one",
    )
    fits["same_sign_as_full"] = np.where(
        fits["estimate_per_log1p_km"].notna() & fits["full_estimate_per_log1p_km"].notna(),
        np.sign(fits["estimate_per_log1p_km"]) == np.sign(fits["full_estimate_per_log1p_km"]),
        np.nan,
    )
    fits["raw_slope_change_from_full"] = (
        fits["estimate_per_log1p_km"] - fits["full_estimate_per_log1p_km"]
    )
    return fits, support


def _summarize(fits: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for (scope, window, context), part in fits.groupby(
        ["evidence_scope", "window", "context"], sort=False
    ):
        ok = part.loc[part["status"].eq("ok")].copy()
        restricted = ok.loc[ok["window"].ne("full")]
        rows.append(
            {
                "evidence_scope": scope,
                "window": window,
                "context": context,
                "n_traits_ok": int(len(ok)),
                "n_traits_nominal_p_lt_0_05": int(pd.to_numeric(ok["p_two_sided"], errors="coerce").lt(0.05).sum()),
                "n_same_sign_as_full": (
                    int(pd.Series(restricted["same_sign_as_full"]).fillna(False).sum())
                    if window != "full"
                    else int(len(ok))
                ),
                "median_abs_raw_slope_change": (
                    float(restricted["raw_slope_change_from_full"].abs().median())
                    if len(restricted)
                    else 0.0
                ),
            }
        )
    region_summary = pd.DataFrame(rows)

    vector_rows = []
    contexts = list(dict.fromkeys(fits["context"].astype(str)))
    for (scope, window), part in fits.loc[fits["status"].eq("ok")].groupby(
        ["evidence_scope", "window"], sort=False
    ):
        pivot = part.pivot(index="outcome", columns="context", values="estimate_per_log1p_km")
        for i, a in enumerate(contexts):
            for b in contexts[i + 1 :]:
                if a not in pivot.columns or b not in pivot.columns:
                    continue
                pair = pivot[[a, b]].dropna()
                if pair.empty:
                    continue
                diff = pair[a] - pair[b]
                denom = float(np.linalg.norm(pair[a]) * np.linalg.norm(pair[b]))
                cosine = (
                    float(np.dot(pair[a], pair[b]) / denom)
                    if denom > 0
                    else np.nan
                )
                vector_rows.append(
                    {
                        "evidence_scope": scope,
                        "window": window,
                        "context_a": a,
                        "context_b": b,
                        "n_shared_traits": int(len(pair)),
                        "rms_raw_slope_difference": float(np.sqrt(np.mean(np.square(diff)))),
                        "cosine_similarity": cosine,
                    }
                )
    return region_summary, pd.DataFrame(vector_rows)


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--flora", type=Path, required=True)
    parser.add_argument("--all-audit", type=Path, required=True)
    parser.add_argument("--direct-audit", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8-sig"))
    flora = pd.read_csv(args.flora, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    all_audit = pd.read_csv(args.all_audit)
    direct_audit = pd.read_csv(args.direct_audit)

    all_prepared = _prepare_scope(flora, all_audit, cov, cfg)
    quantiles = _regional_quantiles(all_prepared, cfg)
    windows = _window_definitions(quantiles)

    fit_parts = []
    support_parts = []
    for scope, audit in [("all", all_audit), ("direct", direct_audit)]:
        prepared = _prepare_scope(flora, audit, cov, cfg)
        fits, support = _fit_scope(prepared, cfg, scope, windows)
        fit_parts.append(fits)
        support_parts.append(support)

    fits = pd.concat(fit_parts, ignore_index=True)
    support = pd.concat(support_parts, ignore_index=True)
    region_summary, vector = _summarize(fits)

    args.output.mkdir(parents=True, exist_ok=True)
    quantiles.to_csv(args.output / "regional_distance_quantiles.csv", index=False)
    pd.DataFrame(windows).to_csv(args.output / "common_support_windows.csv", index=False)
    support.to_csv(args.output / "regional_support_retention.csv", index=False)
    fits.to_csv(args.output / "trait_common_support_refits.csv", index=False)
    region_summary.to_csv(args.output / "regional_common_support_summary.csv", index=False)
    vector.to_csv(args.output / "regional_vector_distance.csv", index=False)

    manifest = {
        "contract": "chapter1_h1_regional_common_support_diagnostic_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "window_source": "primary all-analysis H1 union",
        "windows": windows,
        "claim_boundary": (
            "Restricted-support refits diagnose exposure-support dependence only. "
            "They do not replace the frozen H1 estimand or prove a pollinator mechanism."
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print("=== REGIONAL QUANTILES ===")
    print(quantiles.to_string(index=False))
    print("\n=== WINDOWS ===")
    print(pd.DataFrame(windows).to_string(index=False))
    print("\n=== SUPPORT RETENTION ===")
    print(support.to_string(index=False))
    print("\n=== REGIONAL SUMMARY ===")
    print(region_summary.to_string(index=False))
    print("\n=== REGION VECTOR DISTANCE ===")
    print(vector.to_string(index=False))


if __name__ == "__main__":
    main()
