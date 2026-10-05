"""Post-hoc direct tests of region-by-isolation heterogeneity in syndrome concordance.

The preceding diagnostics showed that the four Chapter 1 regions occupy different
mainland-distance ranges. This script asks the stricter question: for the SAME frozen
architecture|colour response, do isolation slopes differ among regions when distance
support is made comparable?

Two designs are diagnostic only:

1. common-support linear model:
   restrict all regions to the four-region 5th--95th percentile overlap and estimate
   one isolation slope per region while allowing region-specific intercepts and
   region-specific nuisance slopes for selfing_core, area, and climate.

2. full-support hinge model:
   retain distances above the common lower bound, use the common-support upper bound
   as a fixed hinge, and estimate region-specific below-hinge slopes and region-specific
   slope changes above the hinge.

No knot or response is tuned to the resulting coefficients. These analyses do not
replace the frozen Chapter 1 submission inference and do not identify realized pollinators.
"""
from __future__ import annotations

import argparse
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_context_analysis import (
    _chi_square_sf_integer_df,
    _fit_grouped_binomial_design,
)
from island_v2.chapter1_v13_raw_colour_audit import _selfing_table
from island_v2.chapter1_v13_raw_colour_coupling_audit import COUPLING_SPECS


CONTEXTS = [
    "northern_midlatitude",
    "northern_high_latitude",
    "tropical",
    "southern_extratropical",
]

FOCUS = {
    "north_mid_vs_tropical_yellow_butterfly_deep": (
        "yellow_orange__butterfly_deep_tube_given_colour",
        "northern_midlatitude",
        "tropical",
    ),
    "north_mid_vs_tropical_yellow_large_bee_form": (
        "yellow_orange__large_bee_form_given_colour",
        "northern_midlatitude",
        "tropical",
    ),
    "north_mid_vs_north_high_blue_large_bee_deep": (
        "blue_purple__large_bee_deep_tube_given_colour",
        "northern_midlatitude",
        "northern_high_latitude",
    ),
    "north_mid_vs_north_high_blue_butterfly_form": (
        "blue_purple__butterfly_form_given_colour",
        "northern_midlatitude",
        "northern_high_latitude",
    ),
}


def _normal_two_sided_p(z: float) -> float:
    if not math.isfinite(z):
        return float("nan")
    return float(math.erfc(abs(float(z)) / math.sqrt(2.0)))


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


def _z_within_context(
    values: pd.Series,
    context: pd.Series,
    contexts: list[str],
) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    c = context.astype(str).to_numpy()
    out = np.zeros(len(x), dtype=float)
    for label in contexts:
        mask = c == label
        vals = x[mask]
        if len(vals) < 2:
            raise ValueError(f"too few rows to standardize control in {label}")
        mean = float(np.mean(vals))
        sd = float(np.std(vals, ddof=0))
        if not np.isfinite(sd) or sd <= 0:
            raise ValueError(f"constant control in {label}")
        out[mask] = (vals - mean) / sd
    return out


def _load_window(path: Path) -> tuple[float, float]:
    windows = pd.read_csv(path)
    row = windows.loc[windows["window"].astype(str).eq("common_05_95")]
    if len(row) != 1:
        raise RuntimeError("expected one common_05_95 row")
    return float(row.iloc[0]["lower_km"]), float(row.iloc[0]["upper_km"])


def _assemble(
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
    work = counts.loc[counts["stratum"].astype(str).eq("all_observed")].copy()
    work["island_id"] = work["island_id"].astype(str)
    work = (
        work.merge(
            covariates[needed].drop_duplicates("island_id"),
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
        "distance_to_continent_km",
        geography,
        "selfing_core",
        *baseline,
    ]:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[context] = work[context].fillna("").astype(str)
    work[cluster] = work[cluster].fillna("").astype(str)
    return work


def _prepare_combination(
    data: pd.DataFrame,
    cfg: dict,
    combination: str,
    *,
    lower_km: float,
    upper_km: float | None,
) -> pd.DataFrame:
    geography = str(cfg["geography_column"])
    context = str(cfg["context_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]
    required = [
        "successes",
        "trials",
        "distance_to_continent_km",
        geography,
        "selfing_core",
        context,
        cluster,
        *baseline,
    ]
    work = data.loc[
        data["combination"].astype(str).eq(combination)
        & data[context].isin(CONTEXTS)
    ].dropna(subset=required).copy()
    work = work.loc[
        work["trials"].gt(0)
        & work["successes"].ge(0)
        & work["successes"].le(work["trials"])
        & work["distance_to_continent_km"].ge(lower_km)
        & work[cluster].ne("")
    ].copy()
    if upper_km is not None:
        work = work.loc[work["distance_to_continent_km"].le(upper_km)].copy()
    return work


def _design(
    work: pd.DataFrame,
    cfg: dict,
    *,
    hinge_km: float | None,
) -> tuple[np.ndarray, list[str]]:
    geography = str(cfg["geography_column"])
    context_col = str(cfg["context_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]
    c = work[context_col].astype(str)
    x = pd.to_numeric(work[geography], errors="coerce").to_numpy(float)
    controls = {
        predictor: _z_within_context(work[predictor], c, CONTEXTS)
        for predictor in ["selfing_core", *baseline]
    }

    columns: list[np.ndarray] = []
    names: list[str] = []
    for context in CONTEXTS:
        dummy = c.eq(context).astype(float).to_numpy()
        columns.append(dummy)
        names.append(f"intercept[{context}]")
        for predictor in ["selfing_core", *baseline]:
            columns.append(dummy * controls[predictor])
            names.append(f"z_{predictor}[{context}]")
        columns.append(dummy * x)
        names.append(f"slope[{context}]")
        if hinge_km is not None:
            hinge = np.maximum(0.0, x - float(np.log1p(hinge_km)))
            columns.append(dummy * hinge)
            names.append(f"hinge[{context}]")

    X = np.column_stack(columns)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError("rank-deficient region-specific design")
    return X, names


def _support(
    work: pd.DataFrame,
    cfg: dict,
    *,
    hinge_km: float | None,
) -> dict[str, object]:
    context = str(cfg["context_column"])
    cluster = str(cfg["cluster_column"])
    result: dict[str, object] = {}
    for label in CONTEXTS:
        part = work.loc[work[context].eq(label)]
        result[f"n_{label}"] = int(part["island_id"].nunique())
        result[f"blocks_{label}"] = int(part[cluster].nunique())
        if hinge_km is not None:
            tail = part.loc[part["distance_to_continent_km"].gt(hinge_km)]
            result[f"tail_n_{label}"] = int(tail["island_id"].nunique())
            result[f"tail_blocks_{label}"] = int(tail[cluster].nunique())
    return result


def _fit(
    work: pd.DataFrame,
    cfg: dict,
    *,
    hinge_km: float | None,
) -> tuple[pd.DataFrame, dict[str, object], np.ndarray, list[str]]:
    cluster = str(cfg["cluster_column"])
    X, names = _design(work, cfg, hinge_km=hinge_km)
    args = (
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        X,
        names,
        work[cluster].astype(str).to_numpy(),
    )
    coef, fit, covariance = _fit_grouped_binomial_design(
        *args,
        max_iter=500,
    )
    if fit["status"] != "fit":
        coef, fit, covariance = _fit_grouped_binomial_design(
            *args,
            max_iter=5000,
            tolerance=1e-10,
        )
        fit = {**fit, "retry_used": True}
    else:
        fit = {**fit, "retry_used": False}
    return coef, fit, covariance, names


def _contrast(
    coef: pd.DataFrame,
    covariance: np.ndarray,
    names: list[str],
    terms: dict[str, float],
) -> dict[str, float]:
    beta = coef.set_index("predictor")["estimate_log_odds"].to_dict()
    weights = np.zeros(len(names), dtype=float)
    estimate = 0.0
    for name, weight in terms.items():
        if name not in beta:
            return {
                "estimate": float("nan"),
                "se": float("nan"),
                "p": float("nan"),
            }
        estimate += float(weight) * float(beta[name])
        weights[names.index(name)] = float(weight)
    variance = float(weights @ covariance @ weights)
    se = float(np.sqrt(max(variance, 0.0)))
    z = estimate / se if se > 0 else float("nan")
    return {"estimate": float(estimate), "se": se, "p": _normal_two_sided_p(z)}


def _joint_heterogeneity(
    coef: pd.DataFrame,
    covariance: np.ndarray,
    names: list[str],
    term_prefix: str,
    *,
    add_hinge: bool = False,
) -> dict[str, float]:
    ref = CONTEXTS[0]
    contrasts = []
    beta = coef.set_index("predictor")["estimate_log_odds"].to_dict()
    estimates = []
    for context in CONTEXTS[1:]:
        w = np.zeros(len(names), dtype=float)
        if add_hinge:
            ref_terms = [f"slope[{ref}]", f"hinge[{ref}]"]
            other_terms = [f"slope[{context}]", f"hinge[{context}]"]
            for name in other_terms:
                w[names.index(name)] += 1.0
            for name in ref_terms:
                w[names.index(name)] -= 1.0
            estimate = sum(float(beta[x]) for x in other_terms) - sum(
                float(beta[x]) for x in ref_terms
            )
        else:
            other = f"{term_prefix}[{context}]"
            base = f"{term_prefix}[{ref}]"
            w[names.index(other)] = 1.0
            w[names.index(base)] = -1.0
            estimate = float(beta[other]) - float(beta[base])
        contrasts.append(w)
        estimates.append(estimate)
    C = np.vstack(contrasts)
    values = np.asarray(estimates, dtype=float)
    V = C @ covariance @ C.T
    rank = int(np.linalg.matrix_rank(V))
    statistic = (
        float(values @ np.linalg.pinv(V) @ values)
        if rank > 0
        else float("nan")
    )
    return {
        "chisq": statistic,
        "df": rank,
        "p": _chi_square_sf_integer_df(statistic, rank),
    }


def _fit_common(
    data: pd.DataFrame,
    cfg: dict,
    combination: str,
    *,
    lower_km: float,
    upper_km: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    work = _prepare_combination(
        data,
        cfg,
        combination,
        lower_km=lower_km,
        upper_km=upper_km,
    )
    context_col = str(cfg["context_column"])
    threshold = int(cfg["support_tiers"]["confirmatory"])
    counts = work.groupby(context_col)["island_id"].nunique()
    if any(int(counts.get(context, 0)) < threshold for context in CONTEXTS):
        return {
            "status": "not_testable",
            "reason": "at_least_one_region_below_confirmatory_island_threshold",
            **_support(work, cfg, hinge_km=None),
        }, []

    coef, fit, cov, names = _fit(work, cfg, hinge_km=None)
    if fit["status"] != "fit":
        return {"status": fit["status"], **_support(work, cfg, hinge_km=None)}, []

    joint = _joint_heterogeneity(coef, cov, names, "slope")
    row: dict[str, object] = {
        "status": "fit",
        "common_lower_km": lower_km,
        "common_upper_km": upper_km,
        "joint_slope_heterogeneity_chisq": joint["chisq"],
        "joint_slope_heterogeneity_df": joint["df"],
        "joint_slope_heterogeneity_p": joint["p"],
        **_support(work, cfg, hinge_km=None),
    }
    pairwise = []
    for a, b in combinations(CONTEXTS, 2):
        result = _contrast(
            coef,
            cov,
            names,
            {f"slope[{a}]": 1.0, f"slope[{b}]": -1.0},
        )
        pairwise.append(
            {
                "analysis": "common_support_linear",
                "context_a": a,
                "context_b": b,
                "contrast": "slope_a_minus_slope_b",
                **result,
            }
        )
    for context in CONTEXTS:
        result = _contrast(
            coef,
            cov,
            names,
            {f"slope[{context}]": 1.0},
        )
        row[f"slope_{context}"] = result["estimate"]
        row[f"slope_p_{context}"] = result["p"]
    return row, pairwise


def _fit_hinge(
    data: pd.DataFrame,
    cfg: dict,
    combination: str,
    *,
    lower_km: float,
    hinge_km: float,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    work = _prepare_combination(
        data,
        cfg,
        combination,
        lower_km=lower_km,
        upper_km=None,
    )
    context_col = str(cfg["context_column"])
    threshold = int(cfg["support_tiers"]["confirmatory"])
    counts = work.groupby(context_col)["island_id"].nunique()
    if any(int(counts.get(context, 0)) < threshold for context in CONTEXTS):
        return {
            "status": "not_testable",
            "reason": "at_least_one_region_below_confirmatory_island_threshold",
            **_support(work, cfg, hinge_km=hinge_km),
        }, []

    coef, fit, cov, names = _fit(work, cfg, hinge_km=hinge_km)
    if fit["status"] != "fit":
        return {
            "status": fit["status"],
            **_support(work, cfg, hinge_km=hinge_km),
        }, []

    pre = _joint_heterogeneity(coef, cov, names, "slope")
    change = _joint_heterogeneity(coef, cov, names, "hinge")
    post = _joint_heterogeneity(coef, cov, names, "slope", add_hinge=True)
    row: dict[str, object] = {
        "status": "fit",
        "lower_km": lower_km,
        "hinge_km": hinge_km,
        "joint_pre_hinge_heterogeneity_chisq": pre["chisq"],
        "joint_pre_hinge_heterogeneity_df": pre["df"],
        "joint_pre_hinge_heterogeneity_p": pre["p"],
        "joint_hinge_change_heterogeneity_chisq": change["chisq"],
        "joint_hinge_change_heterogeneity_df": change["df"],
        "joint_hinge_change_heterogeneity_p": change["p"],
        "joint_post_hinge_heterogeneity_chisq": post["chisq"],
        "joint_post_hinge_heterogeneity_df": post["df"],
        "joint_post_hinge_heterogeneity_p": post["p"],
        "aic": float(fit["aic"]),
        **_support(work, cfg, hinge_km=hinge_km),
    }
    pairwise = []
    for a, b in combinations(CONTEXTS, 2):
        for label, terms in [
            (
                "pre_hinge_slope_a_minus_b",
                {f"slope[{a}]": 1.0, f"slope[{b}]": -1.0},
            ),
            (
                "hinge_change_a_minus_b",
                {f"hinge[{a}]": 1.0, f"hinge[{b}]": -1.0},
            ),
            (
                "post_hinge_slope_a_minus_b",
                {
                    f"slope[{a}]": 1.0,
                    f"hinge[{a}]": 1.0,
                    f"slope[{b}]": -1.0,
                    f"hinge[{b}]": -1.0,
                },
            ),
        ]:
            result = _contrast(coef, cov, names, terms)
            pairwise.append(
                {
                    "analysis": "full_support_hinge",
                    "context_a": a,
                    "context_b": b,
                    "contrast": label,
                    **result,
                }
            )
    for context in CONTEXTS:
        pre_one = _contrast(
            coef, cov, names, {f"slope[{context}]": 1.0}
        )
        change_one = _contrast(
            coef, cov, names, {f"hinge[{context}]": 1.0}
        )
        post_one = _contrast(
            coef,
            cov,
            names,
            {f"slope[{context}]": 1.0, f"hinge[{context}]": 1.0},
        )
        row[f"pre_slope_{context}"] = pre_one["estimate"]
        row[f"pre_p_{context}"] = pre_one["p"]
        row[f"hinge_change_{context}"] = change_one["estimate"]
        row[f"hinge_change_p_{context}"] = change_one["p"]
        row[f"post_slope_{context}"] = post_one["estimate"]
        row[f"post_p_{context}"] = post_one["p"]
    return row, pairwise


def _focus_pairwise(pairwise: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, (combination, a, b) in FOCUS.items():
        part = pairwise.loc[
            pairwise["combination"].eq(combination)
            & (
                (
                    pairwise["context_a"].eq(a)
                    & pairwise["context_b"].eq(b)
                )
                |
                (
                    pairwise["context_a"].eq(b)
                    & pairwise["context_b"].eq(a)
                )
            )
        ].copy()
        if part.empty:
            continue
        part.insert(0, "focus", label)
        rows.append(part)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


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
    lower_km, upper_km = _load_window(args.windows)

    common_rows = []
    hinge_rows = []
    pairwise_rows = []
    for scope, counts_path, scores_path in [
        ("all", args.all_counts, args.all_scores),
        ("direct", args.direct_counts, args.direct_scores),
    ]:
        data = _assemble(
            pd.read_csv(counts_path, dtype={"island_id": str}),
            pd.read_csv(scores_path, dtype={"island_id": str}),
            cov,
            cfg,
        )
        for combination in COUPLING_SPECS:
            common, common_pairs = _fit_common(
                data,
                cfg,
                combination,
                lower_km=lower_km,
                upper_km=upper_km,
            )
            common_rows.append(
                {
                    "evidence_scope": scope,
                    "combination": combination,
                    **common,
                }
            )
            for row in common_pairs:
                pairwise_rows.append(
                    {
                        "evidence_scope": scope,
                        "combination": combination,
                        **row,
                    }
                )

            hinge, hinge_pairs = _fit_hinge(
                data,
                cfg,
                combination,
                lower_km=lower_km,
                hinge_km=upper_km,
            )
            hinge_rows.append(
                {
                    "evidence_scope": scope,
                    "combination": combination,
                    **hinge,
                }
            )
            for row in hinge_pairs:
                pairwise_rows.append(
                    {
                        "evidence_scope": scope,
                        "combination": combination,
                        **row,
                    }
                )

    common = pd.DataFrame(common_rows)
    hinge = pd.DataFrame(hinge_rows)
    pairwise = pd.DataFrame(pairwise_rows)

    common["joint_slope_heterogeneity_q"] = np.nan
    hinge["joint_pre_hinge_heterogeneity_q"] = np.nan
    hinge["joint_hinge_change_heterogeneity_q"] = np.nan
    hinge["joint_post_hinge_heterogeneity_q"] = np.nan
    for scope in ["all", "direct"]:
        idx = common.index[
            common["evidence_scope"].eq(scope) & common["status"].eq("fit")
        ]
        common.loc[idx, "joint_slope_heterogeneity_q"] = _bh(
            common.loc[idx, "joint_slope_heterogeneity_p"]
        )
        idx = hinge.index[
            hinge["evidence_scope"].eq(scope) & hinge["status"].eq("fit")
        ]
        for pcol, qcol in [
            (
                "joint_pre_hinge_heterogeneity_p",
                "joint_pre_hinge_heterogeneity_q",
            ),
            (
                "joint_hinge_change_heterogeneity_p",
                "joint_hinge_change_heterogeneity_q",
            ),
            (
                "joint_post_hinge_heterogeneity_p",
                "joint_post_hinge_heterogeneity_q",
            ),
        ]:
            hinge.loc[idx, qcol] = _bh(hinge.loc[idx, pcol])

    focus = _focus_pairwise(pairwise)

    args.output.mkdir(parents=True, exist_ok=True)
    common.to_csv(
        args.output / "common_support_region_interaction.csv", index=False
    )
    hinge.to_csv(
        args.output / "full_support_hinge_region_interaction.csv", index=False
    )
    pairwise.to_csv(
        args.output / "pairwise_region_contrasts.csv", index=False
    )
    focus.to_csv(
        args.output / "focus_region_contrasts.csv", index=False
    )
    manifest = {
        "contract": "chapter1_pollination_syndrome_region_scale_interaction_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "contexts": CONTEXTS,
        "common_support_km": [lower_km, upper_km],
        "hinge_km": upper_km,
        "nuisance_control": (
            "region-specific intercepts and region-specific standardized slopes for "
            "selfing_core, island area, and climate PC1-PC4"
        ),
        "primary_diagnostic_questions": [
            "do raw architecture|colour isolation slopes differ among regions over shared distance support?",
            "do region-specific slope changes beyond shared support differ?",
        ],
        "claim_boundary": (
            "Direct region-by-isolation heterogeneity is a diagnostic of scale and "
            "context dependence. It does not identify realized pollinators or causal "
            "historical replacement."
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print("=== COMMON SUPPORT REGION INTERACTION ===")
    show_common = [
        "evidence_scope",
        "combination",
        "status",
        "joint_slope_heterogeneity_p",
        "joint_slope_heterogeneity_q",
    ]
    print(common[show_common].to_string(index=False))
    print("\n=== HINGE REGION INTERACTION ===")
    show_hinge = [
        "evidence_scope",
        "combination",
        "status",
        "joint_pre_hinge_heterogeneity_p",
        "joint_pre_hinge_heterogeneity_q",
        "joint_hinge_change_heterogeneity_p",
        "joint_hinge_change_heterogeneity_q",
        "joint_post_hinge_heterogeneity_p",
        "joint_post_hinge_heterogeneity_q",
    ]
    print(hinge[show_hinge].to_string(index=False))
    print("\n=== FOCUS PAIRWISE ===")
    print(focus.to_string(index=False))


if __name__ == "__main__":
    main()
