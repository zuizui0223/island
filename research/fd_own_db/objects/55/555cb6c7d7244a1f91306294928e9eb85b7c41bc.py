"""Paired-ledger native-mainland source-placebo validation (privacy-safe).

PRIVATE input rows describe identical target island x taxon outcomes scored
under two independently archived training ledgers. Source-set signatures and
taxon/island identifiers NEVER appear in the returned aggregate. This audits
predictive sensitivity, not pollination function or evolutionary causation.
"""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd

KEY = ("scope", "trait", "region", "island_id", "control_seed", "taxon_name")
ID_FIELDS = (*KEY, "spatial_block", "ledger", "source_real",
             "source_placebo_1", "source_placebo_2")
NUMERIC = ("fold", "state", "p_real", "p_placebo_1", "p_placebo_2",
           "max_coverage_difference_1", "max_coverage_difference_2",
           "pool_ratio_1", "pool_ratio_2", "source_jaccard_1",
           "source_jaccard_2")
ALLOWED_TRAITS = {"self_incompatibility", "floral_symmetry"}
ALLOWED_SCOPES = {"All", "Direct"}
LEDGERS = {"snapshot", "source"}


def validate_paired_ledgers(source: pd.DataFrame) -> pd.DataFrame:
    """Require identical targets and actual/placebo assignments in both ledgers.

    This cannot check unseen upstream model training. The source construction
    contract must independently prove global exclusion of each heldout fold
    from all native GIFT donor lists.
    """
    needed = set(ID_FIELDS) | set(NUMERIC)
    if not needed.issubset(source.columns):
        raise ValueError(f"Missing private prediction fields: {sorted(needed-set(source.columns))}")
    x = source[list(ID_FIELDS)+list(NUMERIC)].copy()
    if x.empty:
        raise ValueError("No paired-ledger predictions")
    for key in ID_FIELDS:
        if x[key].isna().any() or x[key].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing identity / source field: {key}")
        x[key] = x[key].astype(str).str.strip()
    if not x.scope.isin(ALLOWED_SCOPES).all() or not x.trait.isin(ALLOWED_TRAITS).all():
        raise ValueError("Unknown evidence scope or trait")
    if not x.ledger.isin(LEDGERS).all():
        raise ValueError("Unknown trait ledger")
    if x.duplicated([*KEY, "ledger"]).any():
        raise ValueError("Repeated target name in one ledger")
    if (x.source_real == x.source_placebo_1).any() or (
        x.source_real == x.source_placebo_2).any() or (
        x.source_placebo_1 == x.source_placebo_2).any():
        raise ValueError("Source and placebo signatures must differ")
    for col in NUMERIC:
        x[col] = pd.to_numeric(x[col], errors="coerce")
        if not np.isfinite(x[col].to_numpy(dtype=float)).all():
            raise ValueError(f"Nonfinite private field: {col}")
    if not x.fold.isin(range(5)).all() or not x.state.isin([0, 1]).all():
        raise ValueError("Invalid fold or observed binary state")
    expected = x.taxon_name.map(
        lambda n: hashlib.sha256(n.encode("utf-8")).digest()[0] % 5)
    if not np.array_equal(expected.to_numpy(), x.fold.to_numpy()):
        raise ValueError("Taxon-name heldout-fold mismatch")
    for key in ("p_real", "p_placebo_1", "p_placebo_2"):
        if not x[key].between(0, 1).all():
            raise ValueError(f"Invalid predicted probability: {key}")
    for k in (1, 2):
        if not x[f"max_coverage_difference_{k}"].between(0, .0200001).all():
            raise ValueError("Unmatched donor evidence coverage")
        if not x[f"pool_ratio_{k}"].between(.5, 2).all():
            raise ValueError("Unmatched mainland donor pool size")
        if not x[f"source_jaccard_{k}"].between(0, .7999999).all():
            raise ValueError("Too much source-set overlap")
    if x.groupby(["scope", "trait", "taxon_name"]).state.nunique().gt(1).any():
        raise ValueError("One taxon has conflicting trait states across island records")
    # If any target name occurs in one ledger but not the other, fail.
    if not x.groupby(list(KEY)).ledger.nunique().eq(2).all():
        raise ValueError("Target taxon/island sets differ between ledgers")
    # Verify the same target state, metadata and three source signatures.
    equality = ["state", "fold", "spatial_block", "source_real",
                "source_placebo_1", "source_placebo_2",
                "max_coverage_difference_1", "max_coverage_difference_2",
                "pool_ratio_1", "pool_ratio_2", "source_jaccard_1",
                "source_jaccard_2"]
    if any(x.groupby(list(KEY))[col].nunique().gt(1).any() for col in equality):
        raise ValueError("Different evaluated state or source assignment across ledgers")
    # A source signature is not allowed to change among tested species
    # on the same island, or between the two independently drawn seeds
    # for the REAL source assignment.
    group = ["scope", "trait", "island_id", "control_seed"]
    for col in ("region", "spatial_block", "source_real",
                "source_placebo_1", "source_placebo_2"):
        if x.groupby(group)[col].nunique().gt(1).any():
            raise ValueError("Conflicting island source metadata")
    if x.groupby(["scope", "trait", "island_id"]).source_real.nunique().gt(1).any():
        raise ValueError("Real source differs across seeds")
    return x


def equal_block_paired_summary(
    source: pd.DataFrame, *, min_islands: int = 40, min_blocks: int = 10,
    min_trait_carriers: int = 20, draws: int = 4999, seed: int = 20261110
) -> pd.DataFrame:
    """Aggregate paired Brier score gains: real versus two matched placebos.

    Whole-block bootstrap intervals are conditional on training/reference
    models and observational blocks. They are NOT independent evolution p's.
    """
    if min_islands < 40 or min_blocks < 10 or min_trait_carriers < 20 or draws < 199:
        raise ValueError("Unsafe support gate")
    x = validate_paired_ledgers(source)
    y = x.state.to_numpy(dtype=float)
    real = x.p_real.to_numpy(dtype=float)
    c1 = x.p_placebo_1.to_numpy(dtype=float)
    c2 = x.p_placebo_2.to_numpy(dtype=float)
    x["brier_advantage"] = ((y-c1)**2 + (y-c2)**2)/2 - (y-real)**2
    island = x.groupby(
        ["scope", "trait", "region", "control_seed",
         "spatial_block", "island_id", "ledger"], as_index=False
    ).agg(n_names=("taxon_name", "nunique"),
          advantage=("brier_advantage", "mean"))
    if not island.groupby(
        ["scope", "trait", "region", "control_seed", "island_id"]
    ).n_names.nunique().eq(1).all():
        raise ValueError("Paired ledger carrier counts disagree")
    if (island.n_names < min_trait_carriers).any():
        raise ValueError("Under-supported island trait evidence")
    wide = island.pivot(
        index=["scope", "trait", "region", "control_seed",
               "spatial_block", "island_id"],
        columns="ledger", values="advantage"
    ).reset_index()
    if wide[["snapshot", "source"]].isna().any().any():
        raise ValueError("Missing comparison ledger")
    wide["delta"] = wide["source"] - wide["snapshot"]
    rng = np.random.default_rng(seed)
    result = []
    for (scope, trait, region, control_seed), g in wide.groupby(
        ["scope", "trait", "region", "control_seed"], sort=True
    ):
        blocks = g.groupby("spatial_block")[["snapshot", "source", "delta"]].mean()
        if len(g) < min_islands or len(blocks) < min_blocks:
            continue
        values = blocks.to_numpy(dtype=float)
        r = values[rng.integers(len(blocks), size=(draws, len(blocks)))].mean(axis=1)
        lo, hi = np.quantile(r, [.025, .975], axis=0)
        result.append({
            "scope": scope, "trait": trait, "region": region,
            "control_seed": control_seed, "n_same_target_islands": int(len(g)),
            "n_spatial_blocks": int(len(blocks)),
            "snapshot_gain": float(values[:,0].mean()),
            "snapshot_ci_low": float(lo[0]), "snapshot_ci_high": float(hi[0]),
            "source_gain": float(values[:,1].mean()),
            "source_ci_low": float(lo[1]), "source_ci_high": float(hi[1]),
            "paired_gain_delta": float(values[:,2].mean()),
            "paired_ci_low": float(lo[2]), "paired_ci_high": float(hi[2]),
            "source_positive_block_fraction": float((values[:,1] > 0).mean()),
        })
    return pd.DataFrame.from_records(result, columns=[
        "scope", "trait", "region", "control_seed", "n_same_target_islands",
        "n_spatial_blocks", "snapshot_gain", "snapshot_ci_low",
        "snapshot_ci_high", "source_gain", "source_ci_low", "source_ci_high",
        "paired_gain_delta", "paired_ci_low", "paired_ci_high",
        "source_positive_block_fraction",
    ])
