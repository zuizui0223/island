"""Leave-one-tail-spatial-block robustness for the remote-island deep-tube response.

Post-hoc diagnostic only. Uses the same region-specific hinge model as
audit_pollination_syndrome_nonlinearity.py and asks whether the post-hinge
yellow/orange × butterfly-associated deep-tube slope depends on any single
spatial block containing >783 km observations.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from audit_pollination_syndrome_nonlinearity import (
    _fit_models,
    _load_windows,
    _prepare_focus,
)

COMBINATION = "yellow_orange__butterfly_deep_tube_given_colour"
CONTEXTS = ("northern_midlatitude", "tropical")


def _run_one(
    *,
    scope: str,
    context: str,
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    lower_km: float,
    knot_km: float,
) -> tuple[dict, pd.DataFrame]:
    cluster = str(cfg["cluster_column"])
    work = _prepare_focus(
        counts,
        scores,
        covariates,
        cfg,
        context=context,
        combination=COMBINATION,
        lower_km=lower_km,
    )
    base = _fit_models(work, cfg, knot_km=knot_km)
    tail = work.loc[
        pd.to_numeric(work["distance_to_continent_km"], errors="coerce").gt(knot_km)
    ].copy()
    tail_blocks = sorted(tail[cluster].astype(str).unique())

    rows = []
    for block in tail_blocks:
        reduced = work.loc[work[cluster].astype(str).ne(block)].copy()
        fit = _fit_models(reduced, cfg, knot_km=knot_km)
        rows.append(
            {
                "evidence_scope": scope,
                "context": context,
                "excluded_tail_block": block,
                "n_rows_removed": int(work[cluster].astype(str).eq(block).sum()),
                "n_tail_rows_removed": int(tail[cluster].astype(str).eq(block).sum()),
                "n_islands_after": int(fit["n_islands"]),
                "n_clusters_after": int(fit["n_clusters"]),
                "n_above_knot_after": int(fit["n_above_knot"]),
                "n_blocks_above_knot_after": int(fit["n_blocks_above_knot"]),
                "post_hinge_slope": float(fit["above_slope_per_log1p_km"]),
                "post_hinge_se": float(fit["above_se"]),
                "post_hinge_p": float(fit["above_p"]),
                "hinge_change": float(fit["hinge_change"]),
                "hinge_change_p": float(fit["hinge_change_p"]),
                "delta_aic_hinge_minus_linear": float(
                    fit["delta_aic_hinge_minus_linear"]
                ),
                "same_post_hinge_sign_as_full": bool(
                    np.sign(fit["above_slope_per_log1p_km"])
                    == np.sign(base["above_slope_per_log1p_km"])
                ),
            }
        )

    summary = {
        "evidence_scope": scope,
        "context": context,
        "combination": COMBINATION,
        "hinge_km": knot_km,
        "baseline_post_hinge_slope": float(base["above_slope_per_log1p_km"]),
        "baseline_post_hinge_se": float(base["above_se"]),
        "baseline_post_hinge_p": float(base["above_p"]),
        "baseline_hinge_change": float(base["hinge_change"]),
        "baseline_hinge_change_p": float(base["hinge_change_p"]),
        "baseline_n_tail_islands": int(base["n_above_knot"]),
        "baseline_n_tail_blocks": int(base["n_blocks_above_knot"]),
        "n_tail_blocks_leaveout": int(len(tail_blocks)),
    }
    table = pd.DataFrame(rows)
    if not table.empty:
        summary.update(
            {
                "n_same_post_hinge_sign": int(
                    table["same_post_hinge_sign_as_full"].sum()
                ),
                "min_post_hinge_slope": float(table["post_hinge_slope"].min()),
                "max_post_hinge_slope": float(table["post_hinge_slope"].max()),
                "median_post_hinge_slope": float(table["post_hinge_slope"].median()),
                "max_post_hinge_p": float(table["post_hinge_p"].max()),
                "n_post_hinge_p_lt_0_05": int(table["post_hinge_p"].lt(0.05).sum()),
                "min_delta_aic_hinge_minus_linear": float(
                    table["delta_aic_hinge_minus_linear"].min()
                ),
                "max_delta_aic_hinge_minus_linear": float(
                    table["delta_aic_hinge_minus_linear"].max()
                ),
            }
        )
    return summary, table


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
    lower_km, _, knot_km = _load_windows(args.windows)

    summaries = []
    leaveouts = []
    for scope, count_path, score_path in (
        ("all", args.all_counts, args.all_scores),
        ("direct", args.direct_counts, args.direct_scores),
    ):
        counts = pd.read_csv(count_path, dtype={"island_id": str})
        scores = pd.read_csv(score_path, dtype={"island_id": str})
        for context in CONTEXTS:
            summary, table = _run_one(
                scope=scope,
                context=context,
                counts=counts,
                scores=scores,
                covariates=cov,
                cfg=cfg,
                lower_km=lower_km,
                knot_km=knot_km,
            )
            summaries.append(summary)
            leaveouts.append(table)

    summary_df = pd.DataFrame(summaries)
    leaveout_df = pd.concat(leaveouts, ignore_index=True)

    args.output.mkdir(parents=True, exist_ok=True)
    summary_df.to_csv(args.output / "far_tail_block_leaveout_summary.csv", index=False)
    leaveout_df.to_csv(args.output / "far_tail_block_leaveout_rows.csv", index=False)

    manifest = {
        "contract": "chapter1_pollination_syndrome_far_tail_block_leaveout_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "combination": COMBINATION,
        "contexts": list(CONTEXTS),
        "hinge_km": knot_km,
        "leaveout_unit": "spatial blocks containing at least one >hinge observation",
        "model": (
            "same region-specific hinge grouped-binomial model as the parent "
            "nonlinearity diagnostic; selfing_core + area + climate adjusted"
        ),
        "claim_boundary": (
            "leave-one-block stability addresses spatial leverage only; it does not "
            "identify realized pollinators or establish a causal mechanism"
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print("=== FAR-TAIL BLOCK LEAVEOUT SUMMARY ===")
    print(summary_df.to_string(index=False))


if __name__ == "__main__":
    main()
