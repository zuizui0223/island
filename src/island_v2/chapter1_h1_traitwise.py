"""Retrospective traitwise H1 replay; no pooled directional or domain scores."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import typer
import yaml
from scipy.stats import t

from island_v2.chapter1_all_data_probability import (
    _fit_single_beta_binomial,
    _prepare,
    _standardize,
    build_broad_counts,
)
from island_v2.chapter1_wcvp_native_compatibility import build_regional_native_compatible_flora


def holm_fixed(values):
    """Missing tests keep their family slot but never receive a reportable q."""
    p = np.asarray(values, dtype=float)
    good = np.isfinite(p)
    safe = np.where(good, p, 1.0)
    order = np.argsort(safe, kind="stable")
    adjusted = np.minimum(1.0, np.maximum.accumulate(safe[order] * np.arange(len(p), 0, -1)))
    out = np.empty(len(p))
    out[order] = adjusted
    out[~good] = np.nan
    return out


def cluster_inference(estimate, influences, n, dimension):
    g = len(influences)
    if g < 2 or n <= dimension or not np.isfinite(influences).all():
        return {"status": "not_testable"}
    correction = g / (g - 1) * (n - 1) / (n - dimension)
    squares = np.square(influences)
    se = float(np.sqrt(correction * squares.sum()))
    if not np.isfinite(estimate) or not np.isfinite(se) or se <= 0:
        return {"status": "not_testable"}
    critical = float(t.ppf(0.975, g - 1))
    return {
        "status": "ok",
        "estimate": float(estimate),
        "se": se,
        "df": g - 1,
        "ci_low": float(estimate - critical * se),
        "ci_high": float(estimate + critical * se),
        "p_two_sided": float(2 * t.sf(abs(estimate / se), g - 1)),
        "n_clusters": g,
        "effective_clusters": float(squares.sum() ** 2 / np.square(squares).sum()),
        "maximum_cluster_variance_share": float(squares.max() / squares.sum()),
    }


def fit_trait(part, outcome, cfg):
    base = {"n_islands": int(part.island_id.nunique()), "n_species_trials": int(part.trials.sum())}
    if base["n_islands"] < cfg["minimum_islands_per_outcome"]:
        return dict(base, status="insufficient_islands")
    predictors = [*cfg["baseline_covariates"], cfg["geography_column"]]
    try:
        design = np.column_stack([np.ones(len(part)), *[_standardize(part[x]) for x in predictors]])
    except (ValueError, typer.BadParameter) as exc:
        return dict(base, status="invalid_design", message=str(exc))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return dict(base, status="rank_deficient_design")
    names = ["intercept", *predictors]
    fit = _fit_single_beta_binomial(
        part.successes.to_numpy(float),
        part.trials.to_numpy(float),
        design,
        names,
        max_iter=cfg["max_iter"],
    )
    retry = False
    if not fit["success"]:
        retry = True
        fit = _fit_single_beta_binomial(
            part.successes.to_numpy(float),
            part.trials.to_numpy(float),
            design,
            names,
            max_iter=cfg["retry_max_iter"],
        )
    base.update(
        optimizer_success=bool(fit["success"]), retry_used=retry, optimizer_message=fit["message"]
    )
    if not fit["success"]:
        return dict(base, status="optimizer_failure")
    slope = len(names) - 1
    labels = part[cfg["cluster_column"]].astype(str).to_numpy()
    influences = np.array(
        [
            fit["bread"][slope] @ fit["score"][labels == label].sum(axis=0)
            for label in np.unique(labels)
        ]
    )
    base.update(
        cluster_inference(float(fit["theta"][slope]), influences, len(part), len(fit["theta"]))
    )
    return base


def run_scope(flora, audit, covariates, cfg, stratum):
    local = dict(cfg, strata=[stratum])
    counts = build_broad_counts(flora, audit, local)
    prepared = _prepare(counts, covariates, local)
    rows = []
    for context in cfg["contexts"]:
        for outcome in cfg["model_outcomes"]:
            part = prepared.loc[
                prepared[cfg["context_column"]].eq(context) & prepared.outcome.eq(outcome)
            ]
            domain = next(
                key for key, members in cfg["response_families"].items() if outcome in members
            )
            result = fit_trait(part, outcome, cfg)
            rows.append(dict(context=context, outcome=outcome, domain=domain, **result))
            print(context, outcome, result["status"], flush=True)
    result = pd.DataFrame(rows)
    if "p_two_sided" not in result:
        result["p_two_sided"] = np.nan
    if cfg.get("multiplicity", "holm") == "none":
        selected_p = result.p_two_sided
    else:
        result["p_holm_28"] = holm_fixed(result.p_two_sided)
        selected_p = result.p_holm_28
    result["interpretation"] = "uncertain_or_not_testable"
    supported = selected_p.lt(cfg["alpha"])
    if "estimate" in result:
        result.loc[supported & result.estimate.gt(0), "interpretation"] = "classic_direction"
        result.loc[supported & result.estimate.lt(0), "interpretation"] = "opposite_direction"
    return result


def main():
    parser = argparse.ArgumentParser(__doc__)
    for arg in [
        "flora",
        "all-audit",
        "direct-audit",
        "covariates",
        "wcvp",
        "mapping",
        "config",
        "output",
    ]:
        parser.add_argument("--" + arg, type=Path, required=True)
    args = parser.parse_args()
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8-sig"))
    args.output.mkdir(parents=True, exist_ok=True)
    paths = {k: v for k, v in vars(args).items() if k != "output"}
    manifest = {
        "design": "retrospective_traitwise_no_scores",
        "inference": cfg["inference"],
        "multiplicity": cfg.get("multiplicity", "holm"),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "inputs": {
            k: {"name": v.name, "sha256": hashlib.sha256(v.read_bytes()).hexdigest()}
            for k, v in paths.items()
        },
    }
    flora = pd.read_csv(args.flora, dtype={"island_id": str})
    cov = pd.read_csv(args.covariates, dtype={"island_id": str})
    native, coverage = build_regional_native_compatible_flora(
        flora, pd.read_csv(args.wcvp), pd.read_csv(args.mapping, dtype={"island_id": str})
    )
    manifest["wcvp_coverage"] = coverage
    frames = []
    for evidence, path in [("all", args.all_audit), ("direct", args.direct_audit)]:
        audit = pd.read_csv(path)
        for scope, source, stratum in [
            ("broad", flora, "all_observed"),
            ("wcvp", native, "all_native"),
        ]:
            result = run_scope(source, audit, cov, cfg, stratum)
            result.insert(0, "flora_scope", scope)
            result.insert(0, "evidence_scope", evidence)
            result.to_csv(args.output / f"{evidence}_{scope}.csv", index=False)
            frames.append(result)
    combined = pd.concat(frames, ignore_index=True)
    combined.to_csv(args.output / "traitwise_results.csv", index=False)
    manifest["status_counts"] = combined.status.value_counts().to_dict()
    manifest["complete"] = bool(combined.status.eq("ok").all() and len(combined) == 112)
    manifest["claim_ceiling"] = (
        "Assemblage associations, not causal within-lineage evolution. WCVP is regional native compatibility, not exact island nativity. Pointwise CIs. See inference field for multiplicity policy; nominal tests do not establish family-wise support."
    )
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
