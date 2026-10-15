"""Exploratory H2 covariate ablation with fixed response-specific island support.

No frozen source labels, result locks, or submission p-values are altered.
The original H2 R coefficient was not an independent prespecified primary test.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy.stats import t

from island_v2.chapter1_v14_h2_decomposition import _bh, run_decomposition

AREA = ["log_island_area_km2"]
CLIMATE = [f"climate_pc{i}" for i in range(1, 5)]
SPECS = {
    "isolation_selfing_only": [],
    "area_only": AREA,
    "climate_only": CLIMATE,
    "area_and_climate": AREA + CLIMATE,
}
CONDITIONAL = {"generalized_accessible", "plain_colour"}


def same_support_covariates(covariates: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Complete cases are selected once with the full original covariate set."""
    req = [
        cfg["geography_column"],
        cfg["context_column"],
        cfg["cluster_column"],
        *AREA,
        *CLIMATE,
    ]
    required = ["island_id", *req]
    missing = set(required) - set(covariates)
    if missing:
        raise ValueError(f"missing canonical covariates: {sorted(missing)}")
    work = covariates.loc[:, required].copy()
    numeric = [cfg["geography_column"], *AREA, *CLIMATE]
    for column in numeric:
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work[cfg["context_column"]] = work[cfg["context_column"]].fillna("")
    work[cfg["cluster_column"]] = work[cfg["cluster_column"]].fillna("")
    work = work.dropna(subset=numeric)
    work = work.loc[
        work[cfg["context_column"]].astype(str).ne("")
        & work[cfg["cluster_column"]].astype(str).ne("")
        & work["island_id"].notna()
    ].copy()
    if work["island_id"].duplicated().any():
        raise ValueError("covariate table has duplicated island_id")
    return work


def finite_cluster_p(estimate: object, se: object, clusters: int) -> float:
    try:
        beta = float(estimate)
        stderr = float(se)
        n_clusters = int(clusters)
    except (TypeError, ValueError):
        return math.nan
    if not math.isfinite(beta) or not math.isfinite(stderr):
        return math.nan
    if stderr <= 0 or n_clusters < 2:
        return math.nan
    return float(2 * t.sf(abs(beta / stderr), n_clusters - 1))


def add_finite_reference(rows: pd.DataFrame) -> pd.DataFrame:
    result = rows.copy()
    result["distance_p_finite_t"] = [
        finite_cluster_p(r.distance_estimate, r.distance_se, r.n_clusters)
        for r in result.itertuples()
    ]
    result["selfing_p_finite_t"] = [
        finite_cluster_p(
            r.selfing_core_estimate, r.selfing_core_se, r.n_clusters
        )
        for r in result.itertuples()
    ]
    result["distance_q_eight_finite_t"] = np.nan
    result["selfing_q_eight_exploratory"] = np.nan
    for (scope, model), indices in result.groupby(
        ["evidence_scope", "controls"]
    ).groups.items():
        family = [
            i for i in indices if result.loc[i, "response"] in CONDITIONAL
        ]
        if len(family) != 8:
            raise ValueError(
                f"Expected fixed eight-test family: {scope} {model} {len(family)}"
            )
        for source, target in [
            ("distance_p_finite_t", "distance_q_eight_finite_t"),
            ("selfing_p_finite_t", "selfing_q_eight_exploratory"),
        ]:
            p = pd.to_numeric(result.loc[family, source], errors="coerce")
            result.loc[family, target] = _bh(p.fillna(1)).to_numpy()
    return result


def provenance_gate(
    fitted: pd.DataFrame, frozen: pd.DataFrame, evidence: str
) -> dict:
    keys = ["context", "response"]
    merged = fitted.merge(
        frozen, on=keys, suffixes=("_new", "_frozen"), validate="one_to_one"
    )
    if len(merged) != 12:
        raise RuntimeError(f"Baseline reproduction has {len(merged)} rows")
    if not (merged.n_islands_new == merged.n_islands_frozen).all():
        raise RuntimeError(f"Baseline island denominators drifted ({evidence})")
    max_deviation = 0.0
    for column in ["distance_estimate", "selfing_core_estimate"]:
        left = pd.to_numeric(merged[f"{column}_new"], errors="coerce")
        right = pd.to_numeric(merged[f"{column}_frozen"], errors="coerce")
        if not np.allclose(left, right, atol=1e-4, rtol=1e-4, equal_nan=True):
            raise RuntimeError(f"Baseline {column} drift ({evidence})")
        max_deviation = max(
            max_deviation, float((left - right).abs().max(skipna=True))
        )
    return {
        "evidence_scope": evidence,
        "baseline_reproduced": True,
        "n_models": len(merged),
        "max_abs_beta_difference": max_deviation,
    }


def partial_r_squared(y: np.ndarray, x: np.ndarray) -> float:
    if x.size == 0 or len(y) < x.shape[1] + 2:
        return math.nan
    z = np.column_stack([np.ones(len(y)), x])
    beta = np.linalg.lstsq(z, y, rcond=None)[0]
    residual = y - z @ beta
    baseline = y - y.mean()
    denominator = float(np.dot(baseline, baseline))
    return (
        1.0 - float(np.dot(residual, residual)) / denominator
        if denominator > 0
        else math.nan
    )


def geographic_overlap(cov: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """Descriptive D~area/climate overlap, not identification of confounding."""
    region_key = cfg["context_column"]
    geographic_key = cfg["geography_column"]
    result = []
    for region, part in cov.groupby(region_key, sort=True):
        y = part[geographic_key].to_numpy(float)
        area = part[AREA].to_numpy(float)
        climate = part[CLIMATE].to_numpy(float)
        joined = np.column_stack([area, climate])
        r = float(np.corrcoef(y, area[:, 0])[0, 1]) if len(part) > 2 else math.nan
        result.append({
            "context": region,
            "n_islands_covariate_complete": int(len(part)),
            "correlation_distance_log_area": r,
            "r2_distance_area": partial_r_squared(y, area),
            "r2_distance_climate": partial_r_squared(y, climate),
            "r2_distance_area_and_climate": partial_r_squared(y, joined),
        })
    return pd.DataFrame(result)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-all", type=Path, required=True)
    parser.add_argument("--counts-direct", type=Path, required=True)
    parser.add_argument("--scores-all", type=Path, required=True)
    parser.add_argument("--scores-direct", type=Path, required=True)
    parser.add_argument("--covariates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if cfg["contract"] != "chapter1_v14_h2_decomposition_v1":
        raise ValueError("unexpected H2 source config")
    if cfg["baseline_covariates"] != AREA + CLIMATE:
        raise ValueError("baseline control family no longer matches canonical H2")
    cov = same_support_covariates(pd.read_csv(args.covariates), cfg)

    assembled = []
    checks = []
    for scope, counts_path, scores_path in [
        ("all_analysis_eligible", args.counts_all, args.scores_all),
        ("direct_only", args.counts_direct, args.scores_direct),
    ]:
        counts = pd.read_csv(counts_path)
        scores = pd.read_csv(scores_path)
        fits = {}
        for name, controls in SPECS.items():
            local_cfg = dict(cfg, baseline_covariates=list(controls))
            fit, _ = run_decomposition(
                counts, scores, cov, local_cfg, evidence_scope=scope
            )
            if len(fit) != 12:
                raise RuntimeError(f"{scope}/{name}: missing H2 fit rows")
            fit["controls"] = name
            fits[name] = fit
            assembled.append(fit)
        baseline = fits["area_and_climate"]
        check = provenance_gate(
            baseline, pd.read_csv(args.frozen / (
                "all" if scope == "all_analysis_eligible" else "direct"
            ) / "h2_decomposition_models.csv"), scope
        )
        checks.append(check)
        reference = baseline.set_index(["context", "response"])["n_islands"]
        for name, fitted in fits.items():
            observed = fitted.set_index(["context", "response"])["n_islands"]
            if not observed.equals(reference):
                raise RuntimeError(f"Sample change between control sets: {scope}/{name}")

    results = add_finite_reference(pd.concat(assembled, ignore_index=True))
    for col in ["distance_estimate", "selfing_core_estimate"]:
        full = results.loc[results.controls.eq("area_and_climate"), [
            "evidence_scope", "context", "response", col
        ]].rename(columns={col: "full_" + col})
        results = results.merge(
            full, on=["evidence_scope", "context", "response"], validate="many_to_one"
        )
        results["delta_" + col] = results[col] - results["full_" + col]
    args.output.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output / "ablation_h2_estimates.csv", index=False)
    geographic_overlap(cov, cfg).to_csv(
        args.output / "area_climate_distance_overlap.csv", index=False
    )
    summary = {
        "status": "posthoc_covariate_sensitivity_not_submission",
        "objective": "Sensitivity to island area and climate; common islands per H2 response",
        "full_case_covariate_islands": len(cov),
        "controls": SPECS,
        "baseline_gates": checks,
        "no_model_data_or_result_lock_changed": True,
        "limitations": [
            "Adjusted conditional associations, not founder effect or causal mediation",
            "Climate PCs have no independent mechanism identification",
            "H2 R coefficient was not the frozen primary estimand",
            "Same-island comparisons do not assess climate/area effect modification",
            "GloPL H2.2/H2.3 have different geographic units and are not refitted",
            "Nominal H1 and frozen H2 testing decisions remain unchanged",
        ],
    }
    (args.output / "manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
