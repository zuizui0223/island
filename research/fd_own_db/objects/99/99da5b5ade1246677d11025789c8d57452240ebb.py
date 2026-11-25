"""Within-archipelago test of specialist vs partial nectarivorous birds against island floral form.

Estimands, nulls and gates are frozen in docs/RANGE_MAP_NECTARIVORE_CONTRACT_20261125.md.

E1 = beta_specialist - beta_partial on the observed tubular share.
E2 = the same on the share minus its plant genus(-grade) block expectation (primary).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import sparse

from .range_map_nectarivore import OTHER, PARTIAL, SPECIALIST, diet_counts

MIN_INFORMATIVE_GROUPS = 30
MIN_ISLANDS = 200
MIN_SCORED = 5
CONTROLS = ("log_area", "log2_dist", "log2_n", "abs_lat")


def plant_block_expectation(
    plants: pd.DataFrame, block_cols=("genus", "grade"), trait_col: str = "trait"
) -> pd.DataFrame:
    """Observed and block-expected trait share per island.

    ``plants`` has one row per island x scored species with columns ``island``, ``species``, the block
    columns and a 0/1 ``trait``. The block mean is taken over the distinct species of the whole table
    (the region pool), so the expectation is that of exchanging trait labels within blocks.
    """
    block_cols = list(block_cols)
    need = {"island", "species", trait_col, *block_cols}
    if need - set(plants.columns):
        raise KeyError(f"plant table missing {sorted(need - set(plants.columns))}")
    t = plants[trait_col].to_numpy()
    if not np.isin(t, (0, 1)).all():
        raise ValueError("trait must be 0/1")
    if plants.duplicated(["island", "species"]).any():
        raise ValueError("duplicate island x species rows")
    if plants.groupby("species")[trait_col].nunique().max() > 1:
        raise ValueError("a species has conflicting trait values")
    per_species = plants.drop_duplicates("species")
    conflict = plants.groupby("species")[block_cols].nunique().max().max()
    if conflict > 1:
        raise ValueError("a species is assigned to more than one block")
    block_mean = per_species.groupby(block_cols)[trait_col].mean().rename("block_mean")
    d = plants.join(block_mean, on=block_cols)
    out = d.groupby("island").agg(
        n=(trait_col, "size"), y=(trait_col, "sum"), expected_y=("block_mean", "sum")
    )
    out["share"] = out["y"] / out["n"]
    out["expected"] = out["expected_y"] / out["n"]
    out["residual"] = out["share"] - out["expected"]
    return out


def bird_design(counts: np.ndarray, islands: pd.DataFrame, extra=None) -> np.ndarray:
    cols = [
        np.log1p(counts[:, SPECIALIST]),
        np.log1p(counts[:, PARTIAL]),
        np.log1p(counts[:, OTHER]),
    ]
    if extra is not None:
        cols.append(np.log1p(extra))
    cols += [islands[c].to_numpy(float) for c in CONTROLS]
    return np.column_stack(cols)


def add_controls(islands: pd.DataFrame) -> pd.DataFrame:
    """Derive the frozen controls from area_km2, dist_km, n, lat."""
    need = {"area_km2", "dist_km", "n", "lat"}
    if need - set(islands.columns):
        raise KeyError(f"island table missing {sorted(need - set(islands.columns))}")
    d = islands.copy()
    d["log_area"] = np.log(d["area_km2"].astype(float))
    d["log2_dist"] = np.log2(1 + d["dist_km"].astype(float))
    d["log2_n"] = np.log2(d["n"].astype(float))
    d["abs_lat"] = d["lat"].abs().astype(float)
    return d


def group_index(groups) -> list[np.ndarray]:
    g = pd.Series(np.asarray(groups)).reset_index(drop=True)
    return [np.flatnonzero(g.to_numpy() == k) for k in pd.unique(g)]


def within_fit(X: np.ndarray, y: np.ndarray, groups: list[np.ndarray], mult=None):
    """Demean within group; each group gets equal total weight (times ``mult`` for bootstrap copies)."""
    xs, ys, used = [], [], 0
    for k, idx in enumerate(groups):
        m = 1 if mult is None else mult[k]
        if m == 0 or len(idx) < 2:
            continue
        x, v = X[idx], y[idx]
        w = np.sqrt(m / len(idx))
        xs.append((x - x.mean(0)) * w)
        ys.append((v - v.mean()) * w)
        used += 1
    if not xs:
        return None
    x, v = np.vstack(xs), np.concatenate(ys)
    if np.linalg.matrix_rank(x) < x.shape[1]:
        return None
    beta = np.linalg.lstsq(x, v, rcond=None)[0]
    return beta


def contrast(beta) -> float:
    return float(beta[0] - beta[1])


def informative_groups(groups: list[np.ndarray], counts: np.ndarray) -> int:
    return sum(1 for idx in groups if len(idx) >= 2 and np.unique(counts[idx, SPECIALIST]).size > 1)


@dataclass
class NectarTestResult:
    status: str
    n_islands: int
    n_groups: int
    n_informative_groups: int
    estimates: dict
    bootstrap: dict
    permutation: dict

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def permute_within_blocks(labels: np.ndarray, blocks: np.ndarray, rng) -> np.ndarray:
    out = labels.copy()
    for b in np.unique(blocks):
        idx = np.flatnonzero(blocks == b)
        if len(idx) > 1:
            out[idx] = labels[rng.permutation(idx)]
    return out


def run_test(
    islands: pd.DataFrame,
    incidence: sparse.spmatrix,
    diet: np.ndarray,
    bird_blocks: dict[str, np.ndarray],
    n_boot: int = 1999,
    n_perm: int = 1999,
    seed: int = 20261125,
    min_groups: int = MIN_INFORMATIVE_GROUPS,
    min_islands: int = MIN_ISLANDS,
    outcomes=("share", "residual"),
) -> NectarTestResult:
    """Fit E1/E2 with group bootstrap and within-block bird diet permutations.

    ``islands`` rows align with ``incidence`` rows and carry ``group``, ``share``, ``residual`` and the
    raw control columns. ``bird_blocks`` maps a null name (e.g. 'genus', 'family') to per-species ids.
    """
    d = add_controls(islands).reset_index(drop=True)
    if incidence.shape[0] != len(d):
        raise ValueError("incidence rows must align with islands")
    if (d["n"] < MIN_SCORED).any():
        raise ValueError(f"islands with fewer than {MIN_SCORED} scored species must be removed upstream")
    counts = diet_counts(incidence, diet)
    groups = group_index(d["group"])
    n_inf = informative_groups(groups, counts)
    status = "ok" if (n_inf >= min_groups and len(d) >= min_islands) else "insufficient_support"
    X = bird_design(counts, d)
    est, boot, perm = {}, {}, {}
    rng = np.random.default_rng(seed)
    for name in outcomes:
        y = d[name].to_numpy(float)
        within_var = sum(float(np.var(y[idx])) for idx in groups if len(idx) >= 2)
        beta = within_fit(X, y, groups) if within_var > 1e-12 else None
        if beta is None:
            est[name] = None
            boot[name] = None
            perm[name] = {}
            continue
        e = contrast(beta)
        est[name] = {
            "E": e,
            "beta_specialist": float(beta[0]),
            "beta_partial": float(beta[1]),
            "beta_other": float(beta[2]),
        }
        draws, dropped = [], 0
        for _ in range(n_boot):
            pick = rng.integers(0, len(groups), len(groups))
            mult = np.bincount(pick, minlength=len(groups))
            b = within_fit(X, y, groups, mult)
            if b is None:
                dropped += 1
            else:
                draws.append(contrast(b))
        boot[name] = {
            "q025": float(np.quantile(draws, 0.025)) if draws else None,
            "q975": float(np.quantile(draws, 0.975)) if draws else None,
            "valid": len(draws),
            "dropped_rank_deficient": dropped,
        }
        perm[name] = {}
        for block_name, blocks in bird_blocks.items():
            null, failed = [], 0
            for _ in range(n_perm):
                lab = permute_within_blocks(np.asarray(diet), np.asarray(blocks), rng)
                b = within_fit(bird_design(diet_counts(incidence, lab), d), y, groups)
                if b is None:
                    failed += 1
                else:
                    null.append(contrast(b))
            null = np.asarray(null)
            perm[name][block_name] = {
                "p_one_sided": float((1 + np.sum(null >= e)) / (1 + len(null))),
                "null_median": float(np.median(null)) if len(null) else None,
                "valid": int(len(null)),
                "failed": failed,
            }
    return NectarTestResult(status, len(d), len(groups), n_inf, est, boot, perm)


def decision(result: NectarTestResult) -> str:
    """Apply the frozen decision table to the primary outcome."""
    if result.status != "ok":
        return "insufficient_support"
    b1, b2 = result.bootstrap.get("share"), result.bootstrap.get("residual")
    if not b1 or b1["q025"] is None:
        return "not_estimable"
    e1_pos = b1["q025"] > 0
    e2_pos = bool(b2 and b2["q025"] is not None and b2["q025"] > 0)
    p_genus = result.permutation.get("residual", {}).get("genus", {}).get("p_one_sided", 1.0)
    if e2_pos and p_genus < 0.05:
        return "beyond_plant_and_bird_lineage"
    if e2_pos:
        return "beyond_plant_lineage_only"
    if e1_pos:
        return "lineage_cooccurrence"
    if b1["q025"] <= 0 <= b1["q975"]:
        return "no_generalisation"
    return "negative_association"


def run_clade_replication(
    islands: pd.DataFrame,
    incidence: sparse.spmatrix,
    diet: np.ndarray,
    clade: np.ndarray,
    clade_names,
    outcome: str = "residual",
    n_boot: int = 1999,
    seed: int = 20261125,
) -> dict:
    """E for each clade's specialists, other specialists entering as a separate control term."""
    d = add_controls(islands).reset_index(drop=True)
    groups = group_index(d["group"])
    y = d[outcome].to_numpy(float)
    rng = np.random.default_rng(seed)
    out = {}
    for c, name in enumerate(clade_names):
        counts, extra = clade_counts(incidence, diet, clade, c)
        n_inf = informative_groups(groups, counts)
        X = bird_design(counts, d, extra=extra)
        beta = within_fit(X, y, groups)
        if beta is None:
            out[name] = {"informative_groups": n_inf, "E": None}
            continue
        draws = []
        for _ in range(n_boot):
            mult = np.bincount(rng.integers(0, len(groups), len(groups)), minlength=len(groups))
            b = within_fit(X, y, groups, mult)
            if b is not None:
                draws.append(contrast(b))
        out[name] = {
            "informative_groups": n_inf,
            "E": contrast(beta),
            "q025": float(np.quantile(draws, 0.025)) if draws else None,
            "q975": float(np.quantile(draws, 0.975)) if draws else None,
            "valid": len(draws),
        }
    return out


def clade_counts(incidence: sparse.spmatrix, diet: np.ndarray, clade: np.ndarray, focal: int):
    """Counts with specialists of ``focal`` clade kept as specialists and other specialists split off."""
    diet = np.asarray(diet)
    in_focal = (diet == SPECIALIST) & (np.asarray(clade) == focal)
    other_spec = (diet == SPECIALIST) & ~in_focal
    lab = np.where(other_spec, OTHER, diet)
    counts = diet_counts(incidence, lab)
    counts[:, OTHER] -= np.asarray(sparse.csr_matrix(incidence, dtype=float) @ other_spec.astype(float)).ravel()
    extra = np.asarray(sparse.csr_matrix(incidence, dtype=float) @ other_spec.astype(float)).ravel()
    return counts, extra
