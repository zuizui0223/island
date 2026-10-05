"""Post-hoc common-isolation-support diagnostic for pollination-syndrome concordance.

Uses frozen raw architecture|colour counts and corrected Chapter 1 geography. The
diagnostic asks whether region-specific colour-architecture coupling results persist
when all four regions are compared over the same mainland-distance support.

This is diagnostic only and must not replace the frozen submission estimand.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_v13_raw_colour_coupling_audit import (
    fit_colour_conditioned_architecture_models,
)


KEYS = ["stratum", "context", "combination", "support_tier", "model"]
FOCUS = {
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
    "north_mid_yellow_large_bee_form": (
        "northern_midlatitude",
        "yellow_orange__large_bee_form_given_colour",
    ),
}


def _load_window(path: Path) -> tuple[float, float]:
    table = pd.read_csv(path)
    row = table.loc[table["window"].eq("common_05_95")]
    if len(row) != 1:
        raise RuntimeError("expected one common_05_95 window")
    return float(row.iloc[0]["lower_km"]), float(row.iloc[0]["upper_km"])


def _fit(
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
) -> pd.DataFrame:
    local = dict(cfg)
    local["strata"] = ["all_observed"]
    return fit_colour_conditioned_architecture_models(counts, scores, covariates, local)


def _verify_replay(replayed: pd.DataFrame, reference: pd.DataFrame, scope: str) -> None:
    merged = replayed.merge(
        reference,
        on=KEYS,
        suffixes=("_replay", "_reference"),
        validate="one_to_one",
    )
    if len(merged) != len(reference):
        raise RuntimeError(f"{scope}: corrected full replay row mismatch")
    for column in ["distance_estimate", "distance_se", "distance_p", "distance_q"]:
        a = pd.to_numeric(merged[f"{column}_replay"], errors="coerce").to_numpy(float)
        b = pd.to_numeric(merged[f"{column}_reference"], errors="coerce").to_numpy(float)
        if not np.allclose(a, b, atol=1e-6, rtol=1e-5, equal_nan=True):
            delta = np.nanmax(np.abs(a - b))
            raise RuntimeError(f"{scope}: full replay mismatch for {column}: {delta}")


def _focus(table: pd.DataFrame) -> pd.DataFrame:
    return table.loc[
        table["support_tier"].eq("confirmatory")
        & table["model"].eq("conditional_selfing")
    ].copy()


def _comparison(full: pd.DataFrame, common: pd.DataFrame, scope: str) -> pd.DataFrame:
    full = _focus(full)
    common = _focus(common)
    keep = [
        "stratum",
        "context",
        "combination",
        "architecture_label",
        "trait_name",
        "status",
        "n_unique_islands",
        "n_clusters",
        "distance_estimate",
        "distance_se",
        "distance_p",
        "distance_q",
    ]
    merged = full[keep].merge(
        common[keep],
        on=["stratum", "context", "combination", "architecture_label", "trait_name"],
        how="outer",
        suffixes=("_full", "_common"),
        validate="one_to_one",
    )
    merged.insert(0, "evidence_scope", scope)
    merged["same_sign"] = np.where(
        merged["distance_estimate_full"].notna() & merged["distance_estimate_common"].notna(),
        np.sign(merged["distance_estimate_full"]) == np.sign(merged["distance_estimate_common"]),
        np.nan,
    )
    merged["common_minus_full"] = (
        merged["distance_estimate_common"] - merged["distance_estimate_full"]
    )
    merged["full_fdr_supported"] = merged["distance_q_full"].lt(0.05)
    merged["common_fdr_supported"] = merged["distance_q_common"].lt(0.05)
    return merged


def _summary(comparison: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (scope, context), part in comparison.groupby(["evidence_scope", "context"]):
        full_fit = part["status_full"].eq("fit")
        common_fit = part["status_common"].eq("fit")
        both = full_fit & common_fit
        rows.append(
            {
                "evidence_scope": scope,
                "context": context,
                "n_confirmatory_combinations_full": int(full_fit.sum()),
                "n_confirmatory_combinations_common": int(common_fit.sum()),
                "n_fdr_supported_full": int(part["full_fdr_supported"].fillna(False).sum()),
                "n_fdr_supported_common": int(part["common_fdr_supported"].fillna(False).sum()),
                "n_same_sign_among_both_fit": int(
                    part.loc[both, "same_sign"].fillna(False).sum()
                ),
                "n_both_fit": int(both.sum()),
                "median_abs_estimate_change": float(
                    part.loc[both, "common_minus_full"].abs().median()
                )
                if both.any()
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _focus_signals(comparison: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, (context, combination) in FOCUS.items():
        part = comparison.loc[
            comparison["context"].eq(context)
            & comparison["combination"].eq(combination)
        ]
        for row in part.itertuples(index=False):
            rows.append(
                {
                    "signal": label,
                    "evidence_scope": row.evidence_scope,
                    "context": context,
                    "combination": combination,
                    "full_estimate": row.distance_estimate_full,
                    "full_p": row.distance_p_full,
                    "full_q": row.distance_q_full,
                    "common_estimate": row.distance_estimate_common,
                    "common_p": row.distance_p_common,
                    "common_q": row.distance_q_common,
                    "common_n_islands": row.n_unique_islands_common,
                    "common_n_clusters": row.n_clusters_common,
                    "same_sign": row.same_sign,
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
    parser.add_argument("--reference-all", type=Path, required=True)
    parser.add_argument("--reference-direct", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    lower, upper = _load_window(args.windows)
    common_cov = cov.loc[
        pd.to_numeric(cov["distance_to_continent_km"], errors="coerce").between(
            lower, upper, inclusive="both"
        )
    ].copy()

    outputs = []
    for scope, count_path, score_path, reference_path in [
        ("all", args.all_counts, args.all_scores, args.reference_all),
        ("direct", args.direct_counts, args.direct_scores, args.reference_direct),
    ]:
        counts = pd.read_csv(count_path, dtype={"island_id": str})
        scores = pd.read_csv(score_path, dtype={"island_id": str})
        reference = pd.read_csv(reference_path)
        full = _fit(counts, scores, cov, cfg)
        _verify_replay(full, reference, scope)
        common = _fit(counts, scores, common_cov, cfg)
        outputs.append(_comparison(full, common, scope))

    comparison = pd.concat(outputs, ignore_index=True)
    summary = _summary(comparison)
    focus = _focus_signals(comparison)

    args.output.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(args.output / "syndrome_common_support_comparison.csv", index=False)
    summary.to_csv(args.output / "syndrome_common_support_summary.csv", index=False)
    focus.to_csv(args.output / "syndrome_common_support_focus_signals.csv", index=False)

    manifest = {
        "contract": "chapter1_syndrome_common_support_diagnostic_v1",
        "role": "post_hoc_diagnostic_not_submission_inference",
        "distance_window_km": [lower, upper],
        "window_definition": "four-region 5th-95th percentile support intersection from final H1 primary union",
        "estimand": "P(raw architecture state | raw colour, architecture trait resolved)",
        "model": "conditional on selfing_core; confirmatory tier highlighted",
        "full_corrected_replay_verified": True,
        "claim_boundary": (
            "Concordance labels do not identify realized pollinators. The diagnostic asks "
            "whether syndrome-associated architecture coupling depends on unequal regional "
            "distance support; it is not a new confirmatory mechanism test."
        ),
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print("=== SYNDROME COMMON-SUPPORT SUMMARY ===")
    print(summary.to_string(index=False))
    print("\n=== FOCUS SIGNALS ===")
    print(focus.to_string(index=False))
    print("\n=== COMMON-SUPPORT FDR-SUPPORTED ===")
    supported = comparison.loc[comparison["common_fdr_supported"].fillna(False)]
    show = [
        "evidence_scope",
        "context",
        "combination",
        "distance_estimate_common",
        "distance_p_common",
        "distance_q_common",
        "n_unique_islands_common",
        "n_clusters_common",
    ]
    print(supported[show].to_string(index=False) if len(supported) else "none")


if __name__ == "__main__":
    main()
