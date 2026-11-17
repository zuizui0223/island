"""Original-data source-pool regional contrasts; anonymous aggregate outputs only.
Post-hoc diagnostic, never promotion of the frozen H1-H4 submission baseline.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

REGIONS = ("northern_midlatitude", "tropical")
MODES = ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10")
TRAITS = (
    "self_compatibility", "selfing_mating_system", "autonomous_selfing",
    "plain_colour", "generalized_form", "actinomorphic_symmetry",
    "shallow_open_tube",
)
FLORAL = TRAITS[3:]
COMPONENTS = (
    "raw_h1_mean", "source_species_expectation",
    "genus_species_structure_enrichment",
    "total_species_sorting_enrichment",
    "within_genus_species_sorting_enrichment",
)
COVARS = ["log_island_area_km2", *[f"climate_pc{k}" for k in range(1, 5)]]
LIMITS = (0.25149279162763416, 783.0170567911987)


def ols_cluster(part: pd.DataFrame, response: str, *, adjust_sc: bool) -> dict:
    columns = ["log_distance_to_continent_km", *COVARS]
    if adjust_sc:
        columns.append("sc_proxy")
    clean = part.replace([np.inf, -np.inf], np.nan).dropna(
        subset=[response, *columns, "spatial_block"]
    )
    n = len(clean)
    g = clean.spatial_block.astype(str).nunique()
    if n < 50 or g < 10:
        return {"status": "insufficient_support", "n_islands": n, "n_blocks": g}
    design = np.column_stack([np.ones(n), clean[columns].to_numpy(float)])
    y = clean[response].to_numpy(float)
    p = design.shape[1]
    if n <= p or np.linalg.matrix_rank(design) < p:
        return {"status": "singular", "n_islands": n, "n_blocks": g}
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    resid = y - design @ beta
    bread = np.linalg.inv(design.T @ design)
    influence = []
    blocks = clean.spatial_block.astype(str)
    for block in blocks.unique():
        idx = blocks.eq(block).to_numpy()
        influence.append(float((bread @ (design[idx].T @ resid[idx]))[1]))
    raw = np.square(np.asarray(influence))
    correction = g / (g - 1) * (n - 1) / (n - p)
    se = float(np.sqrt(correction * raw.sum()))
    if se <= 0 or not np.isfinite(se):
        return {"status": "variance_unidentified", "n_islands": n, "n_blocks": g}
    return {
        "status": "fit", "n_islands": n, "n_blocks": g,
        "slope_log1p_km": float(beta[1]), "se": se,
        "p_two_sided": float(2 * t.sf(abs(beta[1] / se), df=g - 1)),
        "max_block_variance_share": float(raw.max() / raw.sum()),
    }


def bh(values: np.ndarray) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(len(p), np.nan)
    good = np.where(np.isfinite(p))[0]
    if len(good):
        order = good[np.argsort(p[good])]
        raw = p[order] * len(good) / np.arange(1, len(good) + 1)
        out[order] = np.minimum(1.0, np.minimum.accumulate(raw[::-1])[::-1])
    return out


def prepare(scores: pd.DataFrame, cov: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    required = {
        "island_id", "context", "evidence_scope", "source_mode",
        "outcome", "identity_error", *COMPONENTS,
    }
    if not required.issubset(scores):
        raise ValueError(f"missing score columns: {sorted(required - set(scores))}")
    if scores.identity_error.abs().max() > 1e-12:
        raise ValueError("source decomposition identity does not close")
    keys = ["island_id", "context", "evidence_scope", "source_mode", "outcome"]
    if scores.duplicated(keys).any():
        raise ValueError("duplicate island/outcome/source-mode score rows")
    need_cov = {
        "island_id", "analysis_regime", "spatial_block",
        "log_distance_to_continent_km", *COVARS,
    }
    if not need_cov.issubset(cov):
        raise ValueError(f"missing geography columns: {sorted(need_cov - set(cov))}")
    if cov.island_id.duplicated().any():
        raise ValueError("duplicate corrected geography island IDs")
    signed = assignments.loc[
        assignments.source_mode.isin(MODES),
        ["island_id", "source_mode", "entity_ID"],
    ].drop_duplicates()
    donors = signed.groupby(
        ["island_id", "source_mode"]
    ).entity_ID.apply(
        lambda s: "|".join(sorted(map(str, s)))
    ).rename("donor_signature").reset_index()
    x = scores.loc[
        scores.context.isin(REGIONS) & scores.source_mode.isin(MODES)
    ].copy()
    x = x.merge(
        cov[[
            "island_id", "analysis_regime", "spatial_block",
            "log_distance_to_continent_km", *COVARS,
        ]], on="island_id", validate="many_to_one",
    )
    if not x.context.eq(x.analysis_regime).all():
        raise ValueError("region mismatch between score and corrected geography")
    return x.merge(
        donors, on=["island_id", "source_mode"], validate="many_to_one",
    )


def run(scores: pd.DataFrame, cov: pd.DataFrame, assignments: pd.DataFrame) -> pd.DataFrame:
    data = prepare(scores, cov, assignments)
    rows = []
    for (scope, mode), subset in data.groupby(
        ["evidence_scope", "source_mode"], sort=True
    ):
        sc = subset.loc[
            subset.outcome.eq("self_compatibility"),
            [
                "island_id", "context",
                "total_species_sorting_enrichment",
                "within_genus_species_sorting_enrichment",
            ],
        ].rename(columns={
            "total_species_sorting_enrichment": "sc_total",
            "within_genus_species_sorting_enrichment": "sc_within",
        })
        for window in ("full", "common_05_95"):
            for outcome in TRAITS:
                base = subset.loc[subset.outcome.eq(outcome)].copy()
                if window == "common_05_95":
                    dist = np.expm1(base.log_distance_to_continent_km.astype(float))
                    base = base.loc[dist.between(*LIMITS)].copy()
                for component in COMPONENTS:
                    conditional = (
                        outcome in FLORAL
                        and component in (
                            "total_species_sorting_enrichment",
                            "within_genus_species_sorting_enrichment",
                        )
                    )
                    for adjust_sc in ([False, True] if conditional else [False]):
                        part = base.copy()
                        if adjust_sc:
                            part = part.merge(
                                sc, on=["island_id", "context"],
                                validate="one_to_one",
                            )
                            col = (
                                "sc_within"
                                if component == "within_genus_species_sorting_enrichment"
                                else "sc_total"
                            )
                            part["sc_proxy"] = part[col]
                        fitted = {}
                        donors = {}
                        for region in REGIONS:
                            section = part.loc[part.context.eq(region)]
                            complete_columns = [
                                component, "log_distance_to_continent_km",
                                *COVARS, "spatial_block", "donor_signature",
                            ]
                            if adjust_sc:
                                complete_columns.append("sc_proxy")
                            eligible = section.replace(
                                [np.inf, -np.inf], np.nan
                            ).dropna(subset=complete_columns)
                            donors[region] = int(eligible.donor_signature.nunique())
                            result = ols_cluster(
                                eligible, component, adjust_sc=adjust_sc,
                            )
                            if donors[region] < 5:
                                result["status"] = "insufficient_donors"
                            fitted[region] = result
                        a, b = [fitted[r] for r in REGIONS]
                        record = {
                            "evidence_scope": scope, "source_mode": mode,
                            "window": window, "outcome": outcome,
                            "component": component, "sc_adjusted": adjust_sc,
                            "north_n": a["n_islands"], "north_blocks": a["n_blocks"],
                            "north_donors": donors[REGIONS[0]],
                            "tropical_n": b["n_islands"],
                            "tropical_blocks": b["n_blocks"],
                            "tropical_donors": donors[REGIONS[1]],
                            "north_status": a["status"],
                            "tropical_status": b["status"],
                        }
                        if a["status"] == b["status"] == "fit":
                            diff = a["slope_log1p_km"] - b["slope_log1p_km"]
                            sed = float(np.hypot(a["se"], b["se"]))
                            df = min(a["n_blocks"], b["n_blocks"]) - 1
                            record.update(
                                north_beta=a["slope_log1p_km"],
                                tropical_beta=b["slope_log1p_km"],
                                contrast_north_minus_tropical=diff,
                                contrast_se=sed, contrast_df=df,
                                contrast_p=float(
                                    2 * t.sf(abs(diff / sed), df=df)
                                ),
                                contrast_ci_lo=float(
                                    diff - t.ppf(.975, df) * sed
                                ),
                                contrast_ci_hi=float(
                                    diff + t.ppf(.975, df) * sed
                                ),
                                north_max_block_share=a["max_block_variance_share"],
                                tropical_max_block_share=b["max_block_variance_share"],
                            )
                        rows.append(record)
    out = pd.DataFrame(rows)
    if "contrast_p" not in out:
        out["contrast_p"] = np.nan
    out["contrast_q_bh"] = np.nan
    # BH only within a single sensitivity surface; modes are correlated.
    keys = [
        "evidence_scope", "source_mode", "window", "component", "sc_adjusted",
    ]
    for _, ids in out.groupby(keys).groups.items():
        out.loc[ids, "contrast_q_bh"] = bh(
            out.loc[ids, "contrast_p"].to_numpy(float)
        )
    out["supports_novel_rule"] = False  # channel service is not measured here
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("scores", "covariates", "assignments", "out"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    scores = pd.read_csv(args.scores, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    assignments = pd.read_csv(args.assignments, dtype={"island_id": str})
    result = run(scores, cov, assignments)
    args.out.mkdir(parents=True, exist_ok=True)
    result.to_csv(
        args.out / "anonymous_region_source_rivals.csv", index=False,
    )
    summary = {
        "status": "post_hoc_diagnostic_not_submission_inference",
        "source": "frozen native nonendemic Chapter 1 run 34232450884",
        "scopes": sorted(result.evidence_scope.dropna().unique().tolist()),
        "source_modes": list(MODES), "common_support_km": list(LIMITS),
        "model": "equal-island OLS, raw log1p(km), spatial-block robust",
        "n_panels": len(result),
        "n_fitted_contrasts": int(result.contrast_p.notna().sum()),
        "n_bh_q_below_005": int(result.contrast_q_bh.lt(.05).sum()),
        "rights": "anonymous aggregates only; no island, species or donor IDs",
        "independent_channel_evidence_passed": False,
        "causal_limit": (
            "not within-species evolution, colonization, historical pollinator "
            "loss, effective service, or identified causal mechanism"
        ),
    }
    (args.out / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8",
    )
    print(json.dumps(summary, indent=2), flush=True)
    focus = result.loc[
        result.sc_adjusted
        & result.component.eq("within_genus_species_sorting_enrichment")
        & result.window.eq("common_05_95")
    ]
    print(focus.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
