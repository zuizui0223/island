"""Trait-ledger consistency and paired mainland-source sensitivity.

Runs on rights-authorized PRIVATE per-name prediction tables. No species, genus,
island or spatial-block identities are written to aggregate results. It does
not construct an historical source pool, and bootstrap intervals are conditional
on fitted/selected models rather than independent evolutionary replications.
"""
from __future__ import annotations

import hashlib
import numpy as np
import pandas as pd

STATES = {
    "self_incompatibility": {"SC": 1, "SI": 0},
    "floral_symmetry": {"actinomorphic": 1, "zygomorphic": 0},
}
KEYS = ("scope", "trait", "region", "island_id", "spatial_block", "accepted_species")
REQUIRED_PREDICTION = (*KEYS, "ledger", "fold", "state", "p_global", "p_source")


def strict_states(ledger: pd.DataFrame) -> pd.DataFrame:
    """Keep only explicitly reported single states, not imputed zero states."""
    needed = {"accepted_species", "trait_name", "normalized_value"}
    if not needed.issubset(ledger.columns):
        raise ValueError("Missing original trait evidence fields")
    d = ledger[list(needed)].copy()
    d = d[d.trait_name.isin(STATES)].copy()
    if d.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("Duplicate taxon-trait evidence")
    if d.accepted_species.isna().any():
        raise ValueError("Missing accepted species name")
    d["state"] = [STATES[t].get(v, np.nan) for t, v in
                  zip(d.trait_name, d.normalized_value)]
    d = d.dropna(subset=["state"]).copy()
    return d[["accepted_species", "trait_name", "state"]]


def ledger_overlap(snapshot: pd.DataFrame, source: pd.DataFrame) -> pd.DataFrame:
    """Only grouped support/conflict counts; no taxon IDs in output."""
    a = strict_states(snapshot).rename(columns={"state": "snapshot_state"})
    b = strict_states(source).rename(columns={"state": "source_state"})
    d = a.merge(b, on=["accepted_species", "trait_name"], how="outer",
                validate="one_to_one", indicator=True)
    rows = []
    for trait, g in d.groupby("trait_name"):
        both = g[g._merge.eq("both")]
        rows.append({
            "trait": trait,
            "snapshot_strict": int(g.snapshot_state.notna().sum()),
            "source_strict": int(g.source_state.notna().sum()),
            "both_strict": int(len(both)),
            "only_snapshot": int(g._merge.eq("left_only").sum()),
            "only_source": int(g._merge.eq("right_only").sum()),
            "conflicting_states": int((both.snapshot_state != both.source_state).sum()),
        })
    return pd.DataFrame(rows)


def paired_source_ledger_sensitivity(
    predictions: pd.DataFrame, *, bootstrap_draws: int = 4999,
    min_islands: int = 40, min_blocks: int = 10, seed: int = 20261110
) -> pd.DataFrame:
    """Equal-island then equal-spatial-block paired difference in Brier gain.

    Each same species x island must be present ONCE in each evidence ledger
    (snapshot/source), with the same observed state and SHA-name test fold.
    Positive source gain means local-source Brier < global-genus Brier.
    A positive paired delta means that the source *ledger* gives larger gain.
    This is a sampling/definition sensitivity, NOT geological adaptation.
    """
    if bootstrap_draws < 199 or min_islands < 40 or min_blocks < 10:
        raise ValueError("Invalid minimum reporting support")
    if not set(REQUIRED_PREDICTION).issubset(predictions.columns):
        raise ValueError("Missing private heldout prediction fields")
    d = predictions[list(REQUIRED_PREDICTION)].copy()
    if d.empty:
        raise ValueError("Empty predictions")
    if d[list(KEYS) + ["ledger"]].isna().any().any():
        raise ValueError("Missing identity/ledger field")
    if not d.ledger.isin(["snapshot", "source"]).all():
        raise ValueError("Unknown evidence ledger")
    if not d.trait.isin(STATES).all():
        raise ValueError("Unknown strict binary trait")
    if d.duplicated([*KEYS, "ledger"]).any():
        raise ValueError("Duplicate name/island/ledger outcome")
    for col in ("fold", "state", "p_global", "p_source"):
        d[col] = pd.to_numeric(d[col], errors="coerce")
        if not np.isfinite(d[col].to_numpy(dtype=float)).all():
            raise ValueError("Missing/invalid numerical evidence")
    if not d.state.isin((0, 1)).all():
        raise ValueError("Unobserved nonbinary trait")
    if not (d.fold.isin(range(5))).all():
        raise ValueError("Invalid heldout fold")
    if not d.p_global.between(0, 1).all() or not d.p_source.between(0, 1).all():
        raise ValueError("Invalid probabilities")
    expected = d.accepted_species.map(
        lambda name: hashlib.sha256(str(name).encode("utf-8")).digest()[0] % 5)
    if not np.array_equal(expected.to_numpy(), d.fold.to_numpy()):
        raise ValueError("Species name fold mismatch: possible target leakage")
    # Both ledgers MUST evaluate identical observed species states in each island.
    namekeys = [*KEYS]
    wide_state = d.pivot(index=namekeys, columns="ledger", values="state")
    if wide_state.isna().any().any() or not (wide_state["source"] == wide_state["snapshot"]).all():
        raise ValueError("Paired evaluation species/state differs between ledgers")
    d["gain"] = (d.p_global - d.state)**2 - (d.p_source - d.state)**2
    gain = d.pivot(index=namekeys, columns="ledger", values="gain")
    gain["paired_delta"] = gain["source"] - gain["snapshot"]
    x = gain.reset_index()
    islands = x.groupby(
        ["scope", "trait", "region", "island_id", "spatial_block"],
        as_index=False)[["snapshot", "source", "paired_delta"]].mean()
    block = islands.groupby(
        ["scope", "trait", "region", "spatial_block"],
        as_index=False)[["snapshot", "source", "paired_delta"]].mean()
    rng = np.random.default_rng(seed)
    results = []
    for (scope, trait, region), group in block.groupby(["scope", "trait", "region"], sort=True):
        n_islands = islands.loc[
            (islands.scope == scope) & (islands.trait == trait)
            & (islands.region == region), "island_id"].nunique()
        if n_islands < min_islands or len(group) < min_blocks:
            continue
        v = group.paired_delta.to_numpy()
        resample = v[rng.integers(len(v), size=(bootstrap_draws, len(v)))].mean(axis=1)
        results.append({
            "scope": scope, "trait": trait, "region": region,
            "n_identical_target_islands": int(n_islands),
            "n_spatial_blocks": int(len(group)),
            "snapshot_gain": float(group.snapshot.mean()),
            "source_gain": float(group.source.mean()),
            "source_minus_snapshot_gain": float(v.mean()),
            "block_bootstrap_lower": float(np.quantile(resample, .025)),
            "block_bootstrap_upper": float(np.quantile(resample, .975)),
            "fraction_positive_blocks": float(np.mean(v > 0)),
        })
    return pd.DataFrame(results, columns=[
        "scope", "trait", "region", "n_identical_target_islands",
        "n_spatial_blocks", "snapshot_gain", "source_gain",
        "source_minus_snapshot_gain", "block_bootstrap_lower",
        "block_bootstrap_upper", "fraction_positive_blocks"])
