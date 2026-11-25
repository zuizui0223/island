"""Synthetic islands for the design check in docs/RANGE_MAP_NECTARIVORE_CONTRACT_20261125.md.

Two data-generating processes give the same positive specialist-tubular association (E1 > 0):

- ``cooccurrence``: islands with more specialist birds hold more species of mostly-tubular plant genera;
  which species of a genus occur is unrelated to the birds.
- ``matching``: plant genera occur independently of the birds; within genera that vary in form, the
  tubular species are more likely to occur where specialists are present.
- ``none``: neither.

Only E2 (share minus plant-genus expectation) should separate the first two.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .nectarivore_flora_test import (
    add_controls,
    bird_design,
    contrast,
    group_index,
    plant_block_expectation,
    within_fit,
)
from .range_map_nectarivore import OTHER, PARTIAL, SPECIALIST

SCENARIOS = ("none", "cooccurrence", "matching")


def plant_pool(rng, n_open=40, n_mixed=15, n_tubular=5, per_genus=8):
    p = np.r_[np.full(n_open, 0.05), np.full(n_mixed, 0.5), np.full(n_tubular, 0.95)]
    rows = []
    for g, pg in enumerate(p):
        for k in range(per_genus):
            rows.append((f"g{g}_s{k}", f"g{g}", int(rng.random() < pg), pg))
    return pd.DataFrame(rows, columns=["species", "genus", "trait", "p_genus"])


def simulate(scenario: str, n_groups: int, islands_per_group: int = 8, strength: float = 1.0, seed: int = 0):
    if scenario not in SCENARIOS:
        raise ValueError(scenario)
    rng = np.random.default_rng(seed)
    pool = plant_pool(rng)
    genera = pool["genus"].unique()
    members = [np.flatnonzero(pool["genus"].to_numpy() == q) for q in genera]
    is_tub_genus = np.array([pool["p_genus"].to_numpy()[m[0]] >= 0.9 for m in members], float)
    trait = pool["trait"].to_numpy()
    isl, plants = [], []
    for g in range(n_groups):
        group_mean = rng.gamma(2.0, 1.5)
        for k in range(islands_per_group):
            spec = rng.poisson(group_mean * rng.uniform(0, 2))
            part = rng.poisson(3)
            oth = rng.poisson(30)
            n = int(rng.integers(15, 40))
            x = np.log1p(spec)
            gw = np.exp(strength * x * is_tub_genus) if scenario == "cooccurrence" else np.ones(len(genera))
            taken = np.zeros(len(pool), bool)
            for q in rng.choice(len(genera), size=n, p=gw / gw.sum()):
                cand = members[q][~taken[members[q]]]
                if len(cand) == 0:
                    continue
                w = np.exp(strength * 1.5 * x * trait[cand]) if scenario == "matching" else np.ones(len(cand))
                taken[cand[rng.choice(len(cand), p=w / w.sum())]] = True
            iid = f"G{g}I{k}"
            plants.append(pool.loc[taken].assign(island=iid, grade="h"))
            isl.append(
                {
                    "island": iid,
                    "group": f"G{g}",
                    "specialist": spec,
                    "partial": part,
                    "other": oth,
                    "area_km2": float(np.exp(rng.normal(3, 1.5))),
                    "dist_km": float(rng.uniform(5, 2000)),
                    "lat": float(rng.uniform(-23, 23)),
                }
            )
    plants = pd.concat(plants, ignore_index=True)
    shares = plant_block_expectation(plants, block_cols=("genus", "grade"))
    d = pd.DataFrame(isl).set_index("island").join(shares, how="inner").reset_index()
    return d[d["n"] >= 5].reset_index(drop=True)


def fit_e1_e2(d: pd.DataFrame, n_boot: int = 199, seed: int = 1):
    d = add_controls(d)
    counts = np.zeros((len(d), 3))
    counts[:, OTHER], counts[:, PARTIAL], counts[:, SPECIALIST] = d["other"], d["partial"], d["specialist"]
    X = bird_design(counts, d)
    groups = group_index(d["group"])
    rng = np.random.default_rng(seed)
    out = {}
    for name in ("share", "residual"):
        y = d[name].to_numpy(float)
        b = within_fit(X, y, groups)
        draws = []
        for _ in range(n_boot):
            mult = np.bincount(rng.integers(0, len(groups), len(groups)), minlength=len(groups))
            bb = within_fit(X, y, groups, mult)
            if bb is not None:
                draws.append(contrast(bb))
        out[name] = (
            contrast(b) if b is not None else np.nan,
            float(np.quantile(draws, 0.025)) if draws else np.nan,
            float(np.quantile(draws, 0.975)) if draws else np.nan,
        )
    return out
