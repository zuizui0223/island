"""Exact finite-population flower functional diversity under nested native-origin strata.

This module computes the EXPECTED n-species rarefied categorical Rao Q
analytically, avoiding stochastic rarefaction noise. It deliberately keeps
marginal trait diversity, joint-profile diversity, and colour/architecture
coupling as distinct endpoints. All inputs must be original species-direct
single-state records. No pollen efficacy or island-local evolution is inferred.

For a population with N distinct accepted species and n sampled without
replacement, E[Q_n] = (n-1)/n * P(two distinct population species differ).
"""

from __future__ import annotations

from collections import Counter
from typing import Iterable

import numpy as np

TRAITS = ("flower_primary_color", "floral_form", "floral_symmetry")
METRICS = (
    "marginal_rao",
    "joint_rao",
    "colour_architecture_pairing_departure",
    "open_radial_share",
    "actinomorphic_share",
)


def _pair_same(keys: Iterable[object], n_total: int) -> float:
    counts = Counter(keys)
    return float(sum(k * (k - 1) for k in counts.values()) /
                 (n_total * (n_total - 1)))


def exact_rarefied_flora_fd(triples: np.ndarray, *, n: int = 30) -> dict[str, float]:
    """Expected n-species categorical Q and colour/architecture coupling.

    Colour labels can be randomly permuted across the ENTIRE species pool,
    while form+symmetry remain an intact observed architecture pair.
    Expected shuffled joint-Q is derived exactly, not from 49 random labels.

    Pairing departure = observed E[Q_joint] minus colour-shuffled E[Q_joint].
    Negative means more same-colour/same-architecture clustering than the
    fixed-marginal colour shuffle; positive means dis-assortative pairing.
    This is NOT mutual information and NOT a lineage-conditioned null.
    """
    a = np.asarray(triples, dtype=object)
    if a.ndim != 2 or a.shape[1] != 3 or not isinstance(n, int) or n < 2:
        raise ValueError("Expected complete N x 3 species-direct data and n>=2")
    N = len(a)
    if N < n or N < 2:
        raise ValueError("Insufficient species for exact equal-n rarefaction")
    rows = []
    for row in a:
        if any(v is None or (isinstance(v, float) and np.isnan(v))
               for v in row):
            raise ValueError("Unknown trait state cannot be interpreted as absent")
        state = tuple(str(v).strip() for v in row)
        if any(not v or v.lower() in ("nan", "none") or "|" in v for v in state):
            raise ValueError("Multi-state/unknown trait cannot form a Direct profile")
        rows.append(state)
    same_color = _pair_same((r[0] for r in rows), N)
    same_form = _pair_same((r[1] for r in rows), N)
    same_symmetry = _pair_same((r[2] for r in rows), N)
    same_arch = _pair_same(((r[1], r[2]) for r in rows), N)
    same_joint = _pair_same(rows, N)
    finite = (n - 1) / n
    marginal = finite * (3 - same_color - same_form - same_symmetry) / 3
    joint = finite * (1 - same_joint)
    shuffled_joint = finite * (1 - same_color * same_arch)
    pairing = joint - shuffled_joint
    return {
        "marginal_rao": float(marginal),
        "joint_rao": float(joint),
        "colour_architecture_pairing_departure": float(pairing),
        "open_radial_share": float(sum(r[1] == "open_radial" for r in rows) / N),
        "actinomorphic_share": float(
            sum(r[2] == "actinomorphic" for r in rows) / N),
    }
