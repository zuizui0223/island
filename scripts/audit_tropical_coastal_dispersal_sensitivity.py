"""Post-hoc tropical far-tail coastal-dispersal sensitivity.

This diagnostic tests whether the >783 km tropical yellow/orange × deep-tube
signal survives exclusion of externally defined maritime/coastal-strand taxa.
It also reports finite-cluster sensitivities and far-tail block leave-one-out
ranges. It is not part of the frozen Chapter 1 submission inference.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import t as student_t

from audit_pollination_syndrome_nonlinearity import (
    _fit_models,
    _load_windows,
    _prepare_focus,
)
from island_v2.chapter1_v13_raw_colour_coupling_audit import (
    _species_raw_states,
    build_colour_conditioned_architecture_counts,
)
from island_v2.flora_status_support import stratum_mask

COMBINATION = "yellow_orange__butterfly_deep_tube_given_colour"
CONTEXT = "tropical"


def _resolve_sets(cfg: dict) -> dict[str, set[str]]:
    raw = cfg["exclusion_sets"]
    resolved: dict[str, set[str]] = {}

    def visit(name: str) -> set[str]:
        if name in resolved:
            return resolved[name]
        spec = raw[name] or {}
        values = {str(x).strip() for x in spec.get("species", []) if str(x).strip()}
        parent = spec.get("inherits")
        if parent:
            values |= visit(str(parent))
        resolved[name] = values
        return values

    for name in raw:
        visit(str(name))
    return resolved


def _eligible_species(species_axis: pd.DataFrame) -> pd.DataFrame:
    states = _species_raw_states(species_axis, evidence_scope="direct")
    out = states.loc[
        states["colour_states"].map(lambda x: "yellow_orange" in x)
        & states["tube_depth_class"].map(bool)
    ].copy()
    out["success"] = out["tube_depth_class"].map(lambda x: "deep" in x).astype(int)
    return out[["accepted_species", "success"]]


def _excluded_island_counts(
    status_flora: pd.DataFrame,
    eligible: pd.DataFrame,
    *,
    stratum: str,
    excluded: set[str],
) -> pd.DataFrame:
    if not excluded:
        return pd.DataFrame(
            columns=["island_id", "excluded_trials", "excluded_successes"]
        )
    flora = status_flora.loc[stratum_mask(status_flora, stratum), [
        "island_id", "accepted_species"
    ]].copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.loc[flora["accepted_species"].isin(excluded)].drop_duplicates(
        ["island_id", "accepted_species"]
    )
    joined = flora.merge(eligible, on="accepted_species", how="inner", validate="many_to_one")
    if joined.empty:
        return pd.DataFrame(
            columns=["island_id", "excluded_trials", "excluded_successes"]
        )
    return joined.groupby("island_id", as_index=False).agg(
        excluded_trials=("success", "size"),
        excluded_successes=("success", "sum"),
    )


def _adjust_counts(
    stratified_counts: pd.DataFrame,
    status_flora: pd.DataFrame,
    eligible: pd.DataFrame,
    *,
    stratum: str,
    excluded: set[str],
) -> pd.DataFrame:
    counts = stratified_counts.loc[
        stratified_counts["stratum"].astype(str).eq(stratum)
        & stratified_counts["combination"].astype(str).eq(COMBINATION)
    ].copy()
    removed = _excluded_island_counts(
        status_flora,
        eligible,
        stratum=stratum,
        excluded=excluded,
    )
    counts = counts.merge(removed, on="island_id", how="left", validate="one_to_one")
    counts[["excluded_trials", "excluded_successes"]] = counts[
        ["excluded_trials", "excluded_successes"]
    ].fillna(0)
    counts["trials"] = counts["trials"] - counts["excluded_trials"]
    counts["successes"] = counts["successes"] - counts["excluded_successes"]
    counts = counts.loc[counts["trials"].gt(0)].copy()
    if (counts["successes"] < 0).any() or (counts["successes"] > counts["trials"]).any():
        raise RuntimeError("invalid adjusted counts after coastal exclusion")
    counts["share"] = counts["successes"] / counts["trials"]
    counts["stratum"] = "all_observed"
    return counts


def _prepare_variant(
    counts: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    *,
    stratum: str,
    lower_km: float,
) -> pd.DataFrame:
    score = scores.loc[scores["stratum"].astype(str).eq(stratum)].copy()
    score["stratum"] = "all_observed"
    return _prepare_focus(
        counts,
        score,
        covariates,
        cfg,
        context=CONTEXT,
        combination=COMBINATION,
        lower_km=lower_km,
    )


def _t_p(beta: float, se: float, df: int) -> float:
    if not (math.isfinite(beta) and math.isfinite(se)) or se <= 0 or df <= 0:
        return float("nan")
    return float(2.0 * student_t.sf(abs(beta / se), df=df))


def _far_tail_removed_summary(
    status_flora: pd.DataFrame,
    eligible: pd.DataFrame,
    cov: pd.DataFrame,
    *,
    stratum: str,
    excluded: set[str],
    knot_km: float,
) -> dict[str, float | int]:
    flora = status_flora.loc[stratum_mask(status_flora, stratum), [
        "island_id", "accepted_species"
    ]].drop_duplicates().copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    joined = flora.merge(eligible, on="accepted_species", how="inner", validate="many_to_one")
    joined = joined.merge(
        cov[["island_id", "distance_to_continent_km"]].drop_duplicates("island_id"),
        on="island_id",
        how="inner",
        validate="many_to_one",
    )
    far = joined.loc[joined["distance_to_continent_km"].gt(knot_km)].copy()
    total_deep = int(far["success"].sum())
    total_trials = int(len(far))
    removed = far.loc[far["accepted_species"].isin(excluded)]
    return {
        "far_tail_trial_occurrences_before": total_trials,
        "far_tail_deep_occurrences_before": total_deep,
        "far_tail_trial_occurrences_removed": int(len(removed)),
        "far_tail_deep_occurrences_removed": int(removed["success"].sum()),
        "share_far_tail_deep_removed": (
            float(removed["success"].sum()) / total_deep if total_deep else np.nan
        ),
        "n_excluded_species_present_far_tail": int(
            removed["accepted_species"].nunique()
        ),
    }


def _fit_one(
    work: pd.DataFrame,
    cfg: dict,
    *,
    knot_km: float,
) -> dict:
    fit = _fit_models(work, cfg, knot_km=knot_km)
    tail_blocks = int(fit["n_blocks_above_knot"])
    return {
        **fit,
        "post_hinge_p_t_tail_df": _t_p(
            float(fit["above_slope_per_log1p_km"]),
            float(fit["above_se"]),
            max(tail_blocks - 1, 1),
        ),
        "tail_df_for_t_sensitivity": max(tail_blocks - 1, 1),
    }


def _block_leaveout(
    work: pd.DataFrame,
    cfg: dict,
    *,
    knot_km: float,
    cluster: str,
    variant: str,
    stratum: str,
) -> pd.DataFrame:
    tail_blocks = sorted(
        work.loc[work["distance_to_continent_km"].gt(knot_km), cluster]
        .astype(str)
        .unique()
    )
    rows = []
    for block in tail_blocks:
        reduced = work.loc[work[cluster].astype(str).ne(block)].copy()
        if reduced["island_id"].nunique() < 50:
            continue
        fit = _fit_one(reduced, cfg, knot_km=knot_km)
        rows.append(
            {
                "variant": variant,
                "stratum": stratum,
                "block_removed": block,
                "post_hinge_slope": fit["above_slope_per_log1p_km"],
                "post_hinge_p_normal": fit["above_p"],
                "post_hinge_p_t_tail_df": fit["post_hinge_p_t_tail_df"],
                "n_far_tail_islands": fit["n_above_knot"],
                "n_far_tail_blocks": fit["n_blocks_above_knot"],
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--status-flora", type=Path, required=True)
    parser.add_argument("--direct-scores", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--pattern-config", type=Path, required=True)
    parser.add_argument("--sensitivity-config", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    species_axis = pd.read_csv(args.species_axis)
    status_flora = pd.read_csv(args.status_flora, dtype={"island_id": str})
    scores = pd.read_csv(args.direct_scores, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    pattern_cfg = yaml.safe_load(args.pattern_config.read_text(encoding="utf-8"))
    sensitivity_cfg = yaml.safe_load(
        args.sensitivity_config.read_text(encoding="utf-8")
    )
    lower_km, _iqr, knot_km = _load_windows(args.windows)
    cluster = str(pattern_cfg["cluster_column"])
    exclusion_sets = _resolve_sets(sensitivity_cfg)
    strata = [str(x) for x in sensitivity_cfg["strata"]]

    eligible = _eligible_species(species_axis)
    stratified_counts = build_colour_conditioned_architecture_counts(
        species_axis,
        status_flora,
        evidence_scope="direct",
        strata=strata,
    )

    headline_rows = []
    jackknife_parts = []
    removed_species_rows = []

    for stratum in strata:
        for variant, excluded in exclusion_sets.items():
            counts = _adjust_counts(
                stratified_counts,
                status_flora,
                eligible,
                stratum=stratum,
                excluded=excluded,
            )
            work = _prepare_variant(
                counts,
                scores,
                cov,
                pattern_cfg,
                stratum=stratum,
                lower_km=lower_km,
            )
            if work["island_id"].nunique() < 50:
                headline_rows.append(
                    {
                        "variant": variant,
                        "stratum": stratum,
                        "status": "fewer_than_50_complete_islands",
                        "n_model_islands": int(work["island_id"].nunique()),
                    }
                )
                continue

            fit = _fit_one(work, pattern_cfg, knot_km=knot_km)
            removed = _far_tail_removed_summary(
                status_flora,
                eligible,
                cov,
                stratum=stratum,
                excluded=excluded,
                knot_km=knot_km,
            )
            headline_rows.append(
                {
                    "variant": variant,
                    "stratum": stratum,
                    "status": "fit",
                    "n_model_islands": int(work["island_id"].nunique()),
                    "n_pre_hinge_islands": fit["n_below_or_at_knot"],
                    "n_far_tail_islands": fit["n_above_knot"],
                    "n_pre_hinge_blocks": fit["n_blocks_below_or_at_knot"],
                    "n_far_tail_blocks": fit["n_blocks_above_knot"],
                    "post_hinge_slope": fit["above_slope_per_log1p_km"],
                    "post_hinge_se": fit["above_se"],
                    "post_hinge_p_normal": fit["above_p"],
                    "post_hinge_p_t_tail_df": fit["post_hinge_p_t_tail_df"],
                    "tail_df_for_t_sensitivity": fit["tail_df_for_t_sensitivity"],
                    "hinge_change": fit["hinge_change"],
                    "hinge_change_p_normal": fit["hinge_change_p"],
                    "delta_aic_hinge_minus_linear": fit[
                        "delta_aic_hinge_minus_linear"
                    ],
                    **removed,
                }
            )
            jackknife = _block_leaveout(
                work,
                pattern_cfg,
                knot_km=knot_km,
                cluster=cluster,
                variant=variant,
                stratum=stratum,
            )
            if not jackknife.empty:
                jackknife_parts.append(jackknife)

            if excluded:
                present = sorted(
                    set(
                        status_flora.loc[
                            status_flora["accepted_species"].isin(excluded),
                            "accepted_species",
                        ].astype(str)
                    )
                )
                for species in present:
                    removed_species_rows.append(
                        {
                            "variant": variant,
                            "stratum": stratum,
                            "accepted_species": species,
                        }
                    )

    headline = pd.DataFrame(headline_rows)
    jackknife = (
        pd.concat(jackknife_parts, ignore_index=True)
        if jackknife_parts
        else pd.DataFrame()
    )
    removed_species = pd.DataFrame(removed_species_rows)

    jack_summary_rows = []
    if not jackknife.empty:
        for (variant, stratum), part in jackknife.groupby(
            ["variant", "stratum"], sort=False
        ):
            slopes = pd.to_numeric(part["post_hinge_slope"], errors="coerce")
            jack_summary_rows.append(
                {
                    "variant": variant,
                    "stratum": stratum,
                    "n_leaveouts": int(len(part)),
                    "n_positive": int(slopes.gt(0).sum()),
                    "min_post_hinge_slope": float(slopes.min()),
                    "max_post_hinge_slope": float(slopes.max()),
                    "median_post_hinge_slope": float(slopes.median()),
                    "max_p_t_tail_df": float(
                        pd.to_numeric(
                            part["post_hinge_p_t_tail_df"], errors="coerce"
                        ).max()
                    ),
                }
            )
    jack_summary = pd.DataFrame(jack_summary_rows)

    args.output.mkdir(parents=True, exist_ok=True)
    headline.to_csv(args.output / "coastal_exclusion_hinge_summary.csv", index=False)
    jackknife.to_csv(args.output / "coastal_exclusion_block_leaveout.csv", index=False)
    jack_summary.to_csv(
        args.output / "coastal_exclusion_block_leaveout_summary.csv", index=False
    )
    removed_species.to_csv(
        args.output / "coastal_exclusion_species_present.csv", index=False
    )

    print("=== headline ===")
    print(headline.to_string(index=False))
    print("\n=== block leaveout summary ===")
    print(jack_summary.to_string(index=False))


if __name__ == "__main__":
    main()
