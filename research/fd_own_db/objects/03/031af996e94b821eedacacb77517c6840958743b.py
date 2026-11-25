"""Privacy-safe scoring of preselected mainland-source negative controls.

Private taxon x island predictions are passed in memory. Each evaluated
taxon name must belong to the same held-out fold in every island, and each
source/placebo probability must have been fitted WITHOUT all names in that
fold. This module validates predictions and produces ONLY group aggregates;
it does NOT select control sources, rebuild licensed GIFT data, or test
evolutionary causation.
"""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd

ALLOWED_TRAITS = frozenset({"self_compatibility", "actinomorphic_symmetry"})
ALLOWED_SCOPES = frozenset({"All", "Direct"})
INPUT_COLUMNS = (
    "scope", "trait", "region", "island_id", "taxon_name",
    "spatial_block", "fold", "state", "p_actual",
    "p_placebo_1", "p_placebo_2", "source_pool_ratio_1",
    "source_pool_ratio_2", "source_coverage_diff_1",
    "source_coverage_diff_2", "source_jaccard_1", "source_jaccard_2",
)


def validate_predictions(frame: pd.DataFrame, *, folds: int = 5) -> pd.DataFrame:
    """Check the declared name-level holdout and outcome-blind control calipers.

    It cannot prove donor fitting excluded evaluated species: this remains an
    input construction/manifest requirement, validated at source replay.
    """
    if folds != 5:
        raise ValueError("The frozen reference uses five species-name folds")
    missing = set(INPUT_COLUMNS) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing source/placebo columns: {sorted(missing)}")
    x = frame[list(INPUT_COLUMNS)].copy()
    if x.empty:
        raise ValueError("Empty source/placebo scoring input")
    for field in ("scope", "trait", "region", "island_id", "taxon_name", "spatial_block"):
        if x[field].isna().any() or x[field].astype(str).str.strip().eq("").any():
            raise ValueError(f"Unknown or blank {field}")
        x[field] = x[field].astype(str).str.strip()
    if not x.scope.isin(ALLOWED_SCOPES).all() or not x.trait.isin(ALLOWED_TRAITS).all():
        raise ValueError("Unknown scope or trait")
    if x.duplicated(["scope", "trait", "island_id", "taxon_name"]).any():
        raise ValueError("Duplicate island x taxon x trait score")
    if x.groupby(["scope", "trait", "island_id"]).region.nunique().gt(1).any():
        raise ValueError("Inconsistent island region")
    if x.groupby(["scope", "trait", "island_id"]).spatial_block.nunique().gt(1).any():
        raise ValueError("Inconsistent spatial block")
    x["fold"] = pd.to_numeric(x["fold"], errors="coerce")
    if x.fold.isna().any() or not x.fold.isin(range(folds)).all():
        raise ValueError("Invalid held-out fold")
    expected = x.taxon_name.map(lambda n: hashlib.sha256(n.encode("utf-8")).digest()[0] % folds)
    if not np.array_equal(expected.to_numpy(), x.fold.to_numpy()):
        raise ValueError("Species-fold leakage/identity mismatch")
    if x.groupby("taxon_name").fold.nunique().gt(1).any():
        raise ValueError("A species has inconsistent folds across islands")
    if x.groupby(["scope", "trait", "taxon_name"]).state.nunique().gt(1).any():
        raise ValueError("Contradictory binary trait states for same taxon")
    numeric = ["state", "p_actual", "p_placebo_1", "p_placebo_2",
               "source_pool_ratio_1", "source_pool_ratio_2",
               "source_coverage_diff_1", "source_coverage_diff_2",
               "source_jaccard_1", "source_jaccard_2"]
    for col in numeric:
        x[col] = pd.to_numeric(x[col], errors="coerce")
        if not np.isfinite(x[col].to_numpy()).all():
            raise ValueError(f"Missing or nonfinite {col}")
    if not x.state.isin((0, 1)).all():
        raise ValueError("State must be an observed binary trait (no imputation)")
    if not ((x[["p_actual", "p_placebo_1", "p_placebo_2"]] >= 0).all().all()
            and (x[["p_actual", "p_placebo_1", "p_placebo_2"]] <= 1).all().all()):
        raise ValueError("Invalid probability")
    for i in (1, 2):
        if not x[f"source_pool_ratio_{i}"].between(.5, 2.0).all():
            raise ValueError("Unmatched source pool size")
        if not x[f"source_coverage_diff_{i}"].between(0, .0200001).all():
            raise ValueError("Unmatched outcome-blind genus donor coverage")
        if not x[f"source_jaccard_{i}"].between(0, .7999999).all():
            raise ValueError("Control source geography too similar")
    return x


def safe_grouped_scores(
    frame: pd.DataFrame, *, bootstrap_draws: int = 1999,
    seed: int = 20261110, min_islands: int = 40, min_blocks: int = 10,
) -> pd.DataFrame:
    """Compare local source vs TWO preselected placebos, with equal block weight.

    Positive means real source predictions are more accurate. Bootstrap is
    conditional on fitted donor models and control assignment; not a causal p.
    """
    if bootstrap_draws < 199 or min_islands < 40 or min_blocks < 10:
        raise ValueError("Privacy/support threshold or bootstrap count too small")
    x = validate_predictions(frame)
    y = x.state.to_numpy(dtype=float)
    p = x.p_actual.to_numpy(dtype=float)
    p1 = x.p_placebo_1.to_numpy(dtype=float)
    p2 = x.p_placebo_2.to_numpy(dtype=float)
    x["brier_advantage"] = .5*((y-p1)**2 + (y-p2)**2)-(y-p)**2
    eps = 1e-6
    actual = np.clip(p, eps, 1-eps)
    c1 = np.clip(p1, eps, 1-eps)
    c2 = np.clip(p2, eps, 1-eps)
    # positive = average alternative loss minus actual loss
    x["log_advantage"] = y*np.log(actual) + (1-y)*np.log1p(-actual) - .5*(
        y*(np.log(c1)+np.log(c2))+(1-y)*(np.log1p(-c1)+np.log1p(-c2))
    )
    out = []
    rng = np.random.default_rng(seed)
    for (scope, trait, region), group in x.groupby(["scope", "trait", "region"], sort=True):
        island = group.groupby(["island_id", "spatial_block"], as_index=False)[
            ["brier_advantage", "log_advantage"]
        ].mean()
        if len(island) < min_islands or island.spatial_block.nunique() < min_blocks:
            continue
        block = island.groupby("spatial_block")[["brier_advantage", "log_advantage"]].mean()
        v = block.brier_advantage.to_numpy()
        w = block.log_advantage.to_numpy()
        n = len(v)
        draws = rng.integers(n, size=(bootstrap_draws, n))
        bootstrap = v[draws].mean(axis=1)
        log_bootstrap = w[draws].mean(axis=1)
        left_out = (v.sum() - v)/(n-1)
        out.append({
            "scope": scope, "trait": trait, "region": region,
            "n_islands": int(len(island)),
            "n_spatial_blocks": int(n),
            "brier_advantage_equal_block": float(v.mean()),
            "brier_ci025": float(np.quantile(bootstrap, .025)),
            "brier_ci975": float(np.quantile(bootstrap, .975)),
            "brier_leave_block_min": float(left_out.min()),
            "brier_leave_block_max": float(left_out.max()),
            "positive_block_fraction": float((v > 0).mean()),
            "log_advantage_equal_block": float(w.mean()),
            "log_ci025": float(np.quantile(log_bootstrap, .025)),
            "log_ci975": float(np.quantile(log_bootstrap, .975)),
        })
    columns = ["scope", "trait", "region", "n_islands", "n_spatial_blocks",
               "brier_advantage_equal_block", "brier_ci025", "brier_ci975",
               "brier_leave_block_min", "brier_leave_block_max",
               "positive_block_fraction", "log_advantage_equal_block",
               "log_ci025", "log_ci975"]
    return pd.DataFrame.from_records(out, columns=columns)
