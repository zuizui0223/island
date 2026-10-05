"""Post-hoc piecewise distance diagnostic for pollination-syndrome concordance.

The hinge location is NOT tuned to syndrome results. It is the upper boundary of the
four-region 5th--95th percentile common-support window derived independently from the
final H1 island union (783.0170568 km).

For selected, previously highlighted architecture|colour combinations, fit:

  logit(p) = b0 + b1 * (log1p(distance)-cutoff)
                 + b2 * max(0, log1p(distance)-cutoff)
                 + selfing_core + area + climate

Below-cutoff slope = b1.
Above-cutoff slope = b1 + b2.
b2 tests a change in slope beyond shared regional support.

This is a post-hoc diagnostic and does not replace the frozen submission estimand.
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
from island_v2.chapter1_v13_raw_colour_audit import _bh, _selfing_table

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


def _z(values: pd.Series) -> np.ndarray:
    x = pd.to_numeric(values, errors="coerce").to_numpy(float)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError("constant or invalid predictor")
    return (x - mean) / sd


def _p_two_sided(z: float) -> float:
    return float(math.erfc(abs(float(z)) / math.sqrt(2.0))) if math.isfinite(z) else np.nan


def _fit_focus(
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    *,
    scope: str,
    cutoff_km: float,
) -> pd.DataFrame:
    geography = str(cfg["geography_column"])
    context_col = str(cfg["context_column"])
    cluster_col = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]
    cutoff_x = float(np.log1p(cutoff_km))

    needed_cov = [
        "island_id",
        "distance_to_continent_km",
        geography,
        context_col,
        cluster_col,
        *baseline,
    ]
    data = (
        counts.merge(
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
    numeric = [
        "successes",
        "trials",
        "distance_to_continent_km",
        geography,
        "selfing_core",
        *baseline,
    ]
    for col in numeric:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    data[context_col] = data[context_col].fillna("").astype(str)
    data[cluster_col] = data[cluster_col].fillna("").astype(str)

    rows: list[dict] = []
    for label, (context, combination) in FOCUS.items():
        part = data.loc[
            data["stratum"].astype(str).eq("all_observed")
            & data[context_col].eq(context)
            & data["combination"].astype(str).eq(combination)
        ].dropna(
            subset=[
                "successes",
                "trials",
                "distance_to_continent_km",
                geography,
                "selfing_core",
                *baseline,
                cluster_col,
            ]
        ).copy()
        part = part.loc[
            part["trials"].gt(0)
            & part["successes"].ge(0)
            & part["successes"].le(part["trials"])
        ].copy()

        n = int(part["island_id"].nunique())
        n_below = int(part.loc[part["distance_to_continent_km"].le(cutoff_km), "island_id"].nunique())
        n_above = int(part.loc[part["distance_to_continent_km"].gt(cutoff_km), "island_id"].nunique())
        g_below = int(part.loc[part["distance_to_continent_km"].le(cutoff_km), cluster_col].nunique())
        g_above = int(part.loc[part["distance_to_continent_km"].gt(cutoff_km), cluster_col].nunique())
        if n < 50 or n_below < 30 or n_above < 30 or g_below < 5 or g_above < 5:
            rows.append(
                {
                    "signal": label,
                    "evidence_scope": scope,
                    "context": context,
                    "combination": combination,
                    "status": "insufficient_two_sided_support",
                    "n_islands": n,
                    "n_below_cutoff": n_below,
                    "n_above_cutoff": n_above,
                    "n_clusters_below": g_below,
                    "n_clusters_above": g_above,
                }
            )
            continue

        x = part[geography].to_numpy(float)
        centered = x - cutoff_x
        hinge = np.maximum(0.0, centered)
        names = [
            "intercept",
            "distance_centered",
            "remote_hinge",
            "z_selfing_core",
            *[f"z_{b}" for b in baseline],
        ]
        design = np.column_stack(
            [
                np.ones(len(part), dtype=float),
                centered,
                hinge,
                _z(part["selfing_core"]),
                *[_z(part[b]) for b in baseline],
            ]
        )
        coef, fit, cov = _fit_grouped_binomial_design(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            design,
            names,
            part[cluster_col].astype(str).to_numpy(),
        )
        idx = {name: i for i, name in enumerate(names)}
        beta = coef.set_index("predictor")["estimate_log_odds"]
        b1 = float(beta["distance_centered"])
        b2 = float(beta["remote_hinge"])
        se1 = float(np.sqrt(max(cov[idx["distance_centered"], idx["distance_centered"]], 0.0)))
        se2 = float(np.sqrt(max(cov[idx["remote_hinge"], idx["remote_hinge"]], 0.0)))
        remote = b1 + b2
        remote_var = (
            cov[idx["distance_centered"], idx["distance_centered"]]
            + cov[idx["remote_hinge"], idx["remote_hinge"]]
            + 2.0 * cov[idx["distance_centered"], idx["remote_hinge"]]
        )
        remote_se = float(np.sqrt(max(remote_var, 0.0)))
        z1 = b1 / se1 if se1 > 0 else np.nan
        z2 = b2 / se2 if se2 > 0 else np.nan
        zr = remote / remote_se if remote_se > 0 else np.nan

        # Linear comparator with no hinge, same controls.
        linear_names = [
            "intercept",
            "distance_centered",
            "z_selfing_core",
            *[f"z_{b}" for b in baseline],
        ]
        linear_design = np.column_stack(
            [
                np.ones(len(part), dtype=float),
                centered,
                _z(part["selfing_core"]),
                *[_z(part[b]) for b in baseline],
            ]
        )
        _, linear_fit, _ = _fit_grouped_binomial_design(
            part["successes"].to_numpy(float),
            part["trials"].to_numpy(float),
            linear_design,
            linear_names,
            part[cluster_col].astype(str).to_numpy(),
        )

        rows.append(
            {
                "signal": label,
                "evidence_scope": scope,
                "context": context,
                "combination": combination,
                "status": str(fit["status"]),
                "cutoff_km": cutoff_km,
                "n_islands": n,
                "n_below_cutoff": n_below,
                "n_above_cutoff": n_above,
                "n_clusters": int(fit["n_clusters"]),
                "n_clusters_below": g_below,
                "n_clusters_above": g_above,
                "below_slope_per_log1p_km": b1,
                "below_slope_se": se1,
                "below_slope_p": _p_two_sided(z1),
                "slope_change_above_cutoff": b2,
                "slope_change_se": se2,
                "slope_change_p": _p_two_sided(z2),
                "above_slope_per_log1p_km": remote,
                "above_slope_se": remote_se,
                "above_slope_p": _p_two_sided(zr),
                "hinge_aic": float(fit["aic"]),
                "linear_aic": float(linear_fit["aic"]),
                "linear_minus_hinge_aic": float(linear_fit["aic"] - fit["aic"]),
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
    parser.add_argument("--cutoff-km", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    parts = []
    for scope, count_path, score_path in [
        ("all", args.all_counts, args.all_scores),
        ("direct", args.direct_counts, args.direct_scores),
    ]:
        counts = pd.read_csv(count_path, dtype={"island_id": str})
        scores = pd.read_csv(score_path, dtype={"island_id": str})
        parts.append(
            _fit_focus(
                counts,
                scores,
                cov,
                cfg,
                scope=scope,
                cutoff_km=float(args.cutoff_km),
            )
        )
    result = pd.concat(parts, ignore_index=True)
    fit_mask = result["status"].eq("fit")
    result["slope_change_q"] = np.nan
    for scope, index in result.loc[fit_mask].groupby("evidence_scope").groups.items():
        result.loc[index, "slope_change_q"] = _bh(
            result.loc[index, "slope_change_p"]
        )

    args.output.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output / "syndrome_distance_hinge.csv", index=False)
    manifest = {
        "contract": "chapter1_pollination_syndrome_distance_hinge_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "cutoff_km": float(args.cutoff_km),
        "cutoff_origin": (
            "upper boundary of four-region 5th-95th percentile common-support "
            "window from final H1 primary island union"
        ),
        "model": (
            "continuous piecewise log1p(distance) slope + selfing_core + area + climate"
        ),
        "claim_boundary": (
            "Tests whether selected syndrome-associated architecture couplings change "
            "slope beyond shared regional isolation support. Does not identify realized "
            "pollinators or causal historical replacement."
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    print("=== POLLINATION-SYNDROME DISTANCE HINGE ===")
    print(result.to_string(index=False))


if __name__ == "__main__":
    main()
