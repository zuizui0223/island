"""Post-hoc tropical far-tail dispersal and coding sensitivity.

This diagnostic is outside frozen Chapter 1 inference. It tests whether the
>783 km yellow/orange x deep-tube signal survives removal of literature-backed
ocean-dispersed coastal species/lineages and obvious tube-depth coding flags.
It also reports small-cluster sensitivities for the post-hinge slope.
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import t as student_t

from island_v2.chapter1_v13_raw_colour_coupling_audit import _species_raw_states
from island_v2.chapter1_context_analysis import _fit_grouped_binomial_design
from island_v2.flora_status_support import stratum_mask

from audit_pollination_syndrome_nonlinearity import (
    _fit_models,
    _load_windows,
    _prepare_focus,
    _z,
)

COMBINATION = "yellow_orange__butterfly_deep_tube_given_colour"
CONTEXT = "tropical"
STRATA = ("all_native", "native_nonendemic")


def _genus(species: str) -> str:
    text = str(species).strip()
    return text.split()[0] if text else ""


def _scenario_definitions(
    exclusion_csv: Path,
    coding_csv: Path,
) -> dict[str, dict[str, set[str]]]:
    exclusion = pd.read_csv(exclusion_csv)
    strict_species = set(
        exclusion.loc[
            exclusion["evidence_class"].eq("strict_ocean_dispersed"),
            "accepted_species",
        ]
        .dropna()
        .astype(str)
    )
    strict_genera = {_genus(species) for species in strict_species}

    coding = pd.read_csv(coding_csv)
    coding_species = set(
        coding.loc[
            coding["audit_action"].eq("exclude_in_coding_sensitivity"),
            "accepted_species",
        ]
        .dropna()
        .astype(str)
    )

    return {
        "baseline": {"species": set(), "genera": set()},
        "strict_ocean_species": {
            "species": strict_species,
            "genera": set(),
        },
        "strict_ocean_genera": {
            "species": set(),
            "genera": strict_genera,
        },
        "tube_coding_flagged_species": {
            "species": coding_species,
            "genera": set(),
        },
        "strict_ocean_plus_coding": {
            "species": strict_species | coding_species,
            "genera": set(),
        },
    }


def _filter_flora(
    status_flora: pd.DataFrame,
    *,
    species: set[str],
    genera: set[str],
) -> pd.DataFrame:
    frame = status_flora.copy()
    frame["accepted_species"] = frame["accepted_species"].astype(str)
    keep = ~frame["accepted_species"].isin(species)
    if genera:
        keep &= ~frame["accepted_species"].map(_genus).isin(genera)
    return frame.loc[keep].copy()


def _target_occurrences(
    species_axis: pd.DataFrame,
    status_flora: pd.DataFrame,
    *,
    evidence_scope: str,
) -> pd.DataFrame:
    """Species-level rows exactly equivalent to the focal coupling denominator."""
    states = _species_raw_states(species_axis, evidence_scope=evidence_scope)
    eligible = states.loc[
        states["colour_states"].map(lambda values: "yellow_orange" in values)
        & states["tube_depth_class"].map(bool)
    ].copy()
    eligible["success"] = eligible["tube_depth_class"].map(
        lambda values: "deep" in values
    ).astype(int)
    eligible["accepted_species"] = eligible["accepted_species"].astype(str)
    eligible["genus"] = eligible["accepted_species"].map(_genus)

    flora = status_flora.copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    joined = flora.merge(
        eligible[["accepted_species", "genus", "success"]],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    return joined.drop_duplicates(["island_id", "accepted_species"])


def _fast_target_counts(
    occurrences: pd.DataFrame,
    *,
    stratum: str,
    excluded_species: set[str],
    excluded_genera: set[str],
) -> pd.DataFrame:
    subset = occurrences.loc[stratum_mask(occurrences, stratum)].copy()
    keep = ~subset["accepted_species"].isin(excluded_species)
    if excluded_genera:
        keep &= ~subset["genus"].isin(excluded_genera)
    subset = subset.loc[keep].copy()
    if subset.empty:
        return pd.DataFrame(
            columns=[
                "island_id",
                "stratum",
                "combination",
                "successes",
                "trials",
                "share",
            ]
        )
    counts = (
        subset.groupby("island_id", as_index=False)
        .agg(
            successes=("success", "sum"),
            trials=("accepted_species", "size"),
        )
    )
    counts["stratum"] = "all_observed"
    counts["combination"] = COMBINATION
    counts["share"] = counts["successes"] / counts["trials"]
    return counts


def _selfing_occurrences(
    status_flora: pd.DataFrame,
    species_scores: pd.DataFrame,
    *,
    stratum: str,
) -> pd.DataFrame:
    required = {"accepted_species", "syndrome", "syndrome_concordance"}
    missing = required - set(species_scores.columns)
    if missing:
        raise ValueError(f"species syndrome table missing columns: {sorted(missing)}")
    selfing = species_scores.loc[
        species_scores["syndrome"].astype(str).eq("selfing_core"),
        ["accepted_species", "syndrome_concordance"],
    ].drop_duplicates("accepted_species")
    selfing["accepted_species"] = selfing["accepted_species"].astype(str)

    flora = status_flora.loc[stratum_mask(status_flora, stratum)].copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora["genus"] = flora["accepted_species"].map(_genus)
    return (
        flora[["island_id", "accepted_species", "genus"]]
        .drop_duplicates(["island_id", "accepted_species"])
        .merge(selfing, on="accepted_species", how="inner", validate="many_to_one")
    )


def _recomputed_selfing_view(
    selfing_occurrences: pd.DataFrame,
    *,
    excluded_species: set[str],
    excluded_genera: set[str],
) -> pd.DataFrame:
    work = selfing_occurrences.copy()
    keep = ~work["accepted_species"].isin(excluded_species)
    if excluded_genera:
        keep &= ~work["genus"].isin(excluded_genera)
    work = work.loc[keep].copy()
    if work.empty:
        return pd.DataFrame(
            columns=["island_id", "stratum", "syndrome", "syndrome_score"]
        )
    score = (
        work.groupby("island_id", as_index=False)
        .agg(syndrome_score=("syndrome_concordance", "mean"))
    )
    score["stratum"] = "all_observed"
    score["syndrome"] = "selfing_core"
    return score


def _frozen_selfing_view(
    scores: pd.DataFrame,
    *,
    stratum: str,
) -> pd.DataFrame:
    view = scores.loc[
        scores["stratum"].astype(str).eq(stratum)
        & scores["syndrome"].astype(str).eq("selfing_core")
    ].copy()
    view["stratum"] = "all_observed"
    return view


def _prepare_stratum_work(
    occurrences: pd.DataFrame,
    scores: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    *,
    stratum: str,
    excluded_species: set[str],
    excluded_genera: set[str],
    lower_km: float,
    selfing_view: pd.DataFrame | None = None,
) -> pd.DataFrame:
    counts = _fast_target_counts(
        occurrences,
        stratum=stratum,
        excluded_species=excluded_species,
        excluded_genera=excluded_genera,
    )
    score_view = (
        selfing_view.copy()
        if selfing_view is not None
        else _frozen_selfing_view(scores, stratum=stratum)
    )
    if counts.empty or score_view.empty:
        return pd.DataFrame()

    return _prepare_focus(
        counts,
        score_view,
        covariates,
        cfg,
        context=CONTEXT,
        combination=COMBINATION,
        lower_km=lower_km,
    )


def _cluster_t_p(estimate: float, se: float, n_tail_blocks: int) -> float:
    if not math.isfinite(estimate) or not math.isfinite(se) or se <= 0:
        return float("nan")
    df = max(int(n_tail_blocks) - 1, 1)
    stat = estimate / se
    return float(2.0 * student_t.sf(abs(stat), df=df))


def _tail_block_jackknife(
    work: pd.DataFrame,
    cfg: dict,
    *,
    knot_km: float,
    baseline_slope: float,
) -> dict[str, float | int]:
    cluster = str(cfg["cluster_column"])
    tail_blocks = sorted(
        work.loc[
            work["distance_to_continent_km"].gt(knot_km),
            cluster,
        ]
        .dropna()
        .astype(str)
        .unique()
    )
    slopes: list[float] = []
    for block in tail_blocks:
        reduced = work.loc[~work[cluster].astype(str).eq(block)].copy()
        if reduced["island_id"].nunique() < 50:
            continue
        try:
            fit = _fit_models(reduced, cfg, knot_km=knot_km)
        except (ValueError, np.linalg.LinAlgError):
            continue
        slope = float(fit["above_slope_per_log1p_km"])
        if math.isfinite(slope):
            slopes.append(slope)

    out: dict[str, float | int] = {
        "jackknife_tail_blocks_requested": len(tail_blocks),
        "jackknife_tail_blocks_fit": len(slopes),
    }
    if len(slopes) < 3:
        out.update(
            {
                "jackknife_min_slope": float("nan"),
                "jackknife_max_slope": float("nan"),
                "jackknife_median_slope": float("nan"),
                "jackknife_same_sign": 0,
                "jackknife_se": float("nan"),
                "jackknife_t_p": float("nan"),
            }
        )
        return out

    values = np.asarray(slopes, dtype=float)
    mean = float(values.mean())
    g = len(values)
    se = math.sqrt((g - 1.0) / g * float(np.sum((values - mean) ** 2)))
    p = float("nan")
    if se > 0 and math.isfinite(se):
        stat = baseline_slope / se
        p = float(2.0 * student_t.sf(abs(stat), df=max(g - 1, 1)))
    sign = 1 if baseline_slope > 0 else -1 if baseline_slope < 0 else 0
    same = int(np.sum(np.sign(values) == sign))
    out.update(
        {
            "jackknife_min_slope": float(values.min()),
            "jackknife_max_slope": float(values.max()),
            "jackknife_median_slope": float(np.median(values)),
            "jackknife_same_sign": same,
            "jackknife_se": float(se),
            "jackknife_t_p": p,
        }
    )
    return out


def _far_tail_species_table(
    species_axis: pd.DataFrame,
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    *,
    evidence_scope: str,
    stratum: str,
    knot_km: float,
) -> pd.DataFrame:
    states = _species_raw_states(species_axis, evidence_scope=evidence_scope)
    eligible = states.loc[
        states["colour_states"].map(lambda values: "yellow_orange" in values)
        & states["tube_depth_class"].map(bool)
    ].copy()
    eligible["deep"] = eligible["tube_depth_class"].map(
        lambda values: "deep" in values
    )
    eligible["genus"] = eligible["accepted_species"].map(_genus)

    flora = status_flora.loc[stratum_mask(status_flora, stratum)].copy()
    flora = flora[
        [
            "island_id",
            "accepted_species",
            "origin_status",
            "endemic_status",
            "floristic_status",
        ]
    ].drop_duplicates(["island_id", "accepted_species"])
    joined = flora.merge(
        eligible[
            ["accepted_species", "genus", "deep"]
        ],
        on="accepted_species",
        how="inner",
        validate="many_to_one",
    )
    needed = [
        "island_id",
        "analysis_regime",
        "distance_to_continent_km",
        "spatial_block",
    ]
    joined = joined.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="many_to_one",
    )
    joined = joined.loc[
        joined["analysis_regime"].eq(CONTEXT)
        & joined["distance_to_continent_km"].gt(knot_km)
    ].copy()
    if joined.empty:
        return pd.DataFrame()

    out = (
        joined.groupby(
            [
                "accepted_species",
                "genus",
                "origin_status",
                "floristic_status",
                "deep",
            ],
            as_index=False,
        )
        .agg(
            n_islands=("island_id", "nunique"),
            n_occurrences=("island_id", "size"),
            n_blocks=("spatial_block", "nunique"),
        )
        .sort_values(["deep", "n_occurrences"], ascending=[False, False])
    )
    out.insert(0, "stratum", stratum)
    out.insert(0, "evidence_scope", evidence_scope)
    return out


def _ocean_assembly_gradient(
    status_flora: pd.DataFrame,
    covariates: pd.DataFrame,
    cfg: dict,
    *,
    strict_species: set[str],
    support_islands: set[str],
    stratum: str,
) -> pd.DataFrame:
    flora = status_flora.loc[stratum_mask(status_flora, stratum)].copy()
    flora["island_id"] = flora["island_id"].astype(str)
    flora["accepted_species"] = flora["accepted_species"].astype(str)
    flora = flora.loc[flora["island_id"].isin(support_islands)].drop_duplicates(
        ["island_id", "accepted_species"]
    )
    flora["strict_ocean_species"] = flora["accepted_species"].isin(strict_species).astype(int)
    counts = (
        flora.groupby("island_id", as_index=False)
        .agg(
            trials=("accepted_species", "size"),
            successes=("strict_ocean_species", "sum"),
        )
    )

    geography = str(cfg["geography_column"])
    cluster = str(cfg["cluster_column"])
    baseline = [str(x) for x in cfg["baseline_covariates"]]
    context_col = str(cfg["context_column"])
    needed = ["island_id", geography, context_col, cluster, *baseline]
    work = counts.merge(
        covariates[needed].drop_duplicates("island_id"),
        on="island_id",
        how="left",
        validate="one_to_one",
    )
    work = work.loc[work[context_col].eq(CONTEXT)].dropna(
        subset=["successes", "trials", geography, cluster, *baseline]
    )
    names = ["intercept", f"z_{geography}", *[f"z_{x}" for x in baseline]]
    columns = [np.ones(len(work), dtype=float), _z(work[geography])]
    columns.extend(_z(work[x]) for x in baseline)
    coef, fit, _cov = _fit_grouped_binomial_design(
        work["successes"].to_numpy(float),
        work["trials"].to_numpy(float),
        np.column_stack(columns),
        names,
        work[cluster].astype(str).to_numpy(),
    )
    row = coef.loc[coef["predictor"].eq(f"z_{geography}")].iloc[0]
    estimate = float(row["estimate_log_odds"])
    se = float(row["cluster_robust_se"])
    g = int(fit["n_clusters"])
    t_p = (
        float(2.0 * student_t.sf(abs(estimate / se), df=max(g - 1, 1)))
        if se > 0
        else float("nan")
    )
    return pd.DataFrame(
        [
            {
                "stratum": stratum,
                "n_islands": int(fit["n_islands"]),
                "n_blocks": g,
                "strict_ocean_occurrences": int(work["successes"].sum()),
                "total_flora_occurrences": int(work["trials"].sum()),
                "distance_slope_log_odds_per_sd": estimate,
                "cluster_robust_se": se,
                "normal_p": float(row["p_value"]),
                "cluster_t_p": t_p,
            }
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--species-axis", type=Path, required=True)
    parser.add_argument("--status-flora", type=Path, required=True)
    parser.add_argument("--all-scores", type=Path, required=True)
    parser.add_argument("--direct-scores", type=Path, required=True)
    parser.add_argument("--all-species-scores", type=Path, required=True)
    parser.add_argument("--direct-species-scores", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--windows", type=Path, required=True)
    parser.add_argument("--exclusions", type=Path, required=True)
    parser.add_argument("--coding-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    species_axis = pd.read_csv(args.species_axis)
    status_flora = pd.read_csv(args.status_flora, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    lower_km, _iqr_km, knot_km = _load_windows(args.windows)
    scenarios = _scenario_definitions(args.exclusions, args.coding_audit)

    args.output.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    support_rows: list[dict[str, object]] = []
    species_tables: list[pd.DataFrame] = []

    score_paths = {
        "all": (args.all_scores, args.all_species_scores),
        "direct": (args.direct_scores, args.direct_species_scores),
    }
    occurrence_by_scope = {
        scope: _target_occurrences(
            species_axis,
            status_flora,
            evidence_scope=scope,
        )
        for scope in score_paths
    }
    recomputed_rows: list[dict[str, object]] = []
    recompute_manifest: list[dict[str, object]] = []

    for scope, (score_path, species_score_path) in score_paths.items():
        scores = pd.read_csv(score_path, dtype={"island_id": str})
        species_scores = pd.read_csv(species_score_path)
        occurrences = occurrence_by_scope[scope]

        selfing_occ_by_stratum: dict[str, pd.DataFrame] = {}
        if scope == "direct":
            for stratum in STRATA:
                selfing_occ_by_stratum[stratum] = _selfing_occurrences(
                    status_flora,
                    species_scores,
                    stratum=stratum,
                )
                baseline_recomputed = _recomputed_selfing_view(
                    selfing_occ_by_stratum[stratum],
                    excluded_species=set(),
                    excluded_genera=set(),
                )
                frozen = _frozen_selfing_view(scores, stratum=stratum)[
                    ["island_id", "syndrome_score"]
                ].rename(columns={"syndrome_score": "frozen_score"})
                check = frozen.merge(
                    baseline_recomputed[
                        ["island_id", "syndrome_score"]
                    ].rename(columns={"syndrome_score": "recomputed_score"}),
                    on="island_id",
                    how="inner",
                    validate="one_to_one",
                )
                max_delta = float(
                    (check["frozen_score"] - check["recomputed_score"]).abs().max()
                )
                if not math.isfinite(max_delta) or max_delta > 1e-10:
                    raise RuntimeError(
                        f"{stratum}: recomputed selfing baseline mismatch {max_delta}"
                    )
                recompute_manifest.append(
                    {
                        "evidence_scope": scope,
                        "stratum": stratum,
                        "n_selfing_islands_checked": int(len(check)),
                        "max_abs_frozen_vs_recomputed_selfing": max_delta,
                    }
                )

        for stratum in STRATA:
            species_tables.append(
                _far_tail_species_table(
                    species_axis,
                    status_flora,
                    cov,
                    evidence_scope=scope,
                    stratum=stratum,
                    knot_km=knot_km,
                )
            )
            for scenario, exclusion in scenarios.items():
                work = _prepare_stratum_work(
                    occurrences,
                    scores,
                    cov,
                    cfg,
                    stratum=stratum,
                    excluded_species=exclusion["species"],
                    excluded_genera=exclusion["genera"],
                    lower_km=lower_km,
                )
                n_islands = int(work["island_id"].nunique()) if not work.empty else 0
                if n_islands < 50:
                    rows.append(
                        {
                            "evidence_scope": scope,
                            "stratum": stratum,
                            "scenario": scenario,
                            "status": "fewer_than_50_complete_islands",
                            "n_model_islands": n_islands,
                        }
                    )
                    continue

                fit = _fit_models(work, cfg, knot_km=knot_km)
                t_p = _cluster_t_p(
                    float(fit["above_slope_per_log1p_km"]),
                    float(fit["above_se"]),
                    int(fit["n_blocks_above_knot"]),
                )
                jack = _tail_block_jackknife(
                    work,
                    cfg,
                    knot_km=knot_km,
                    baseline_slope=float(fit["above_slope_per_log1p_km"]),
                )
                rows.append(
                    {
                        "evidence_scope": scope,
                        "stratum": stratum,
                        "scenario": scenario,
                        "status": "fit",
                        "n_model_islands": fit["n_islands"],
                        "n_model_blocks": fit["n_clusters"],
                        "n_prehinge_islands": fit["n_below_or_at_knot"],
                        "n_far_tail_islands": fit["n_above_knot"],
                        "n_prehinge_blocks": fit["n_blocks_below_or_at_knot"],
                        "n_far_tail_blocks": fit["n_blocks_above_knot"],
                        "prehinge_slope": fit["below_slope_per_log1p_km"],
                        "prehinge_normal_p": fit["below_p"],
                        "hinge_change": fit["hinge_change"],
                        "hinge_change_normal_p": fit["hinge_change_p"],
                        "post_hinge_slope": fit["above_slope_per_log1p_km"],
                        "post_hinge_se": fit["above_se"],
                        "post_hinge_normal_p": fit["above_p"],
                        "post_hinge_tailblock_t_p": t_p,
                        **jack,
                    }
                )
                support_rows.append(
                    {
                        "evidence_scope": scope,
                        "stratum": stratum,
                        "scenario": scenario,
                        "excluded_species_n": len(exclusion["species"]),
                        "excluded_genera_n": len(exclusion["genera"]),
                        "excluded_species": "|".join(sorted(exclusion["species"])),
                        "excluded_genera": "|".join(sorted(exclusion["genera"])),
                    }
                )

                if (
                    scope == "direct"
                    and scenario
                    in {
                        "baseline",
                        "strict_ocean_species",
                        "strict_ocean_plus_coding",
                    }
                ):
                    recomputed_view = _recomputed_selfing_view(
                        selfing_occ_by_stratum[stratum],
                        excluded_species=exclusion["species"],
                        excluded_genera=exclusion["genera"],
                    )
                    work_recomputed = _prepare_stratum_work(
                        occurrences,
                        scores,
                        cov,
                        cfg,
                        stratum=stratum,
                        excluded_species=exclusion["species"],
                        excluded_genera=exclusion["genera"],
                        lower_km=lower_km,
                        selfing_view=recomputed_view,
                    )
                    if work_recomputed["island_id"].nunique() >= 50:
                        fit_r = _fit_models(
                            work_recomputed,
                            cfg,
                            knot_km=knot_km,
                        )
                        t_r = _cluster_t_p(
                            float(fit_r["above_slope_per_log1p_km"]),
                            float(fit_r["above_se"]),
                            int(fit_r["n_blocks_above_knot"]),
                        )
                        jack_r = _tail_block_jackknife(
                            work_recomputed,
                            cfg,
                            knot_km=knot_km,
                            baseline_slope=float(
                                fit_r["above_slope_per_log1p_km"]
                            ),
                        )
                        recomputed_rows.append(
                            {
                                "evidence_scope": scope,
                                "stratum": stratum,
                                "scenario": scenario,
                                "selfing_mode": "recomputed_after_exclusion",
                                "n_model_islands": fit_r["n_islands"],
                                "n_far_tail_islands": fit_r["n_above_knot"],
                                "n_far_tail_blocks": fit_r["n_blocks_above_knot"],
                                "post_hinge_slope": fit_r[
                                    "above_slope_per_log1p_km"
                                ],
                                "post_hinge_se": fit_r["above_se"],
                                "post_hinge_normal_p": fit_r["above_p"],
                                "post_hinge_tailblock_t_p": t_r,
                                **jack_r,
                            }
                        )

    result = pd.DataFrame(rows)
    support = pd.DataFrame(support_rows)
    recomputed = pd.DataFrame(recomputed_rows)
    recompute_check = pd.DataFrame(recompute_manifest)
    species = (
        pd.concat([x for x in species_tables if not x.empty], ignore_index=True)
        if any(not x.empty for x in species_tables)
        else pd.DataFrame()
    )

    direct_scores = pd.read_csv(args.direct_scores, dtype={"island_id": str})
    baseline_support = _prepare_stratum_work(
        occurrence_by_scope["direct"],
        direct_scores,
        cov,
        cfg,
        stratum="native_nonendemic",
        excluded_species=set(),
        excluded_genera=set(),
        lower_km=lower_km,
    )
    assembly = _ocean_assembly_gradient(
        status_flora,
        cov,
        cfg,
        strict_species=scenarios["strict_ocean_species"]["species"],
        support_islands=set(baseline_support["island_id"].astype(str)),
        stratum="native_nonendemic",
    )

    result.to_csv(args.output / "tropical_far_tail_dispersal_sensitivity.csv", index=False)
    support.to_csv(args.output / "tropical_far_tail_scenario_manifest.csv", index=False)
    species.to_csv(args.output / "tropical_far_tail_native_species.csv", index=False)
    assembly.to_csv(args.output / "tropical_ocean_species_assembly_gradient.csv", index=False)
    recomputed.to_csv(
        args.output / "tropical_far_tail_recomputed_selfing_sensitivity.csv",
        index=False,
    )
    recompute_check.to_csv(
        args.output / "tropical_far_tail_recomputed_selfing_validation.csv",
        index=False,
    )

    primary = result.loc[
        result["evidence_scope"].eq("direct")
        & result["stratum"].eq("native_nonendemic")
    ].copy()
    primary.to_csv(
        args.output / "tropical_far_tail_primary_direct_native_nonendemic.csv",
        index=False,
    )

    print("=== ocean-dispersed species assembly gradient ===")
    print(assembly.to_string(index=False))
    print("\n=== recomputed-selfing validation ===")
    print(recompute_check.to_string(index=False))
    print("\n=== recomputed-selfing sensitivity ===")
    print(recomputed.to_string(index=False))
    print("\n=== direct native-nonendemic sensitivity ===")
    show = [
        "scenario",
        "n_prehinge_islands",
        "n_far_tail_islands",
        "n_far_tail_blocks",
        "post_hinge_slope",
        "post_hinge_normal_p",
        "post_hinge_tailblock_t_p",
        "jackknife_min_slope",
        "jackknife_max_slope",
        "jackknife_same_sign",
        "jackknife_tail_blocks_fit",
        "jackknife_t_p",
    ]
    print(primary[show].to_string(index=False))

    direct_native = species.loc[
        species["evidence_scope"].eq("direct")
        & species["stratum"].eq("native_nonendemic")
        & species["deep"].eq(True)
    ].copy()
    print("\n=== direct native-nonendemic deep species ===")
    print(
        direct_native[
            ["accepted_species", "genus", "n_islands", "n_occurrences", "n_blocks"]
        ]
        .head(50)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
