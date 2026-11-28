#!/usr/bin/env python3
"""Same-carrier decomposition of the exploratory tropical symmetry-depth contrast.

The original H1 model, carrier masks and protected source hashes are inherited.
Island-local fixed margins define an exact expectation; the separate species-
coherent randomization is conditional on global genus x joint evidence grade,
NOT on the simultaneous margins of every island.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from chapter1_flower_axis_decoupling import (
    LOCK, TRAITS, DIST, digest, trait_pairs, make_frame, design, estimate,
)

SEED = 20261128
DRAWS = 1999
METRICS = ("observed", "genus_expected", "grade_increment",
           "genus_grade_expected", "within_genus_grade")
PANELS = (("all_recorded", 10), ("all_recorded", 30),
          ("high_only", 10), ("literal_native", 10))


def group_matrix(M, code):
    """One row per nonempty island x stratum, with original species columns."""
    coo = M.tocoo()
    if not np.all(coo.data == 1):
        raise ValueError("nonbinary incidence")
    k = int(code.max()) + 1
    key, ix = np.unique(coo.row.astype(np.int64) * k + code[coo.col],
                        return_inverse=True)
    W = sparse.csr_matrix((coo.data, (ix, coo.col)),
                          shape=(len(key), M.shape[1]))
    n = np.asarray(W.sum(1)).ravel()
    return W, key // k, n


def setup(M, deep, code):
    W, islands, n = group_matrix(M, code)
    d = np.asarray(W @ deep).ravel()
    return W, islands, d / n, (d > 0) & (d < n)


def expected_deep(prep, radial, n_islands):
    W, islands, f, _ = prep
    return np.bincount(islands, weights=(W @ radial) * f,
                       minlength=n_islands)


def components(M, A, D, genus, grade):
    """Exact additive partition of P(A|D)-P(A|not D) on fixed local margins."""
    deep = M @ D
    shallow = M @ (1 - D)
    if np.any(deep <= 0) or np.any(shallow <= 0):
        raise ValueError("empty depth subset")
    total_a = M @ A
    obs_deep = M @ (A * D)
    a = setup(M, D, genus)
    b = setup(M, D, grade)
    eg = expected_deep(a, A, M.shape[0])
    egg = expected_deep(b, A, M.shape[0])
    obs = obs_deep / deep - (total_a - obs_deep) / shallow
    ge = eg / deep - (total_a - eg) / shallow
    gge = egg / deep - (total_a - egg) / shallow
    grade_increment = gge - ge
    within = (obs_deep - egg) * (1 / deep + 1 / shallow)
    y = np.column_stack((obs, ge, grade_increment, gge, within))
    if not np.allclose(obs, ge + grade_increment + within, atol=1e-12):
        raise ValueError("partition mismatch")
    return y, b, 1 / deep + 1 / shallow


def joint_grades(ledger, v, high_only):
    z = ledger.loc[
        ledger.trait_name.isin(TRAITS)
        & ledger.evidence_scope.eq("species_direct")
        & ledger.evidence_quality.isin(["high"] if high_only else ["high", "medium"])
        & ledger.normalized_value.notna(),
        ["accepted_species", "trait_name", "normalized_value", "evidence_quality"],
    ].copy()
    z["normalized_value"] = z.normalized_value.str.strip()
    z = z[(z.trait_name.eq(TRAITS[0])
           & z.normalized_value.isin(["actinomorphic", "zygomorphic"]))
          | (z.trait_name.eq(TRAITS[1])
             & z.normalized_value.isin(["absent_or_open", "shallow", "deep"]))]
    if z.duplicated(["accepted_species", "trait_name"]).any():
        raise ValueError("duplicate grade record")
    grade = z.pivot(index="accepted_species", columns="trait_name",
                    values="evidence_quality").loc[v.index, list(TRAITS)]
    if grade.isna().any().any():
        raise ValueError("missing grade")
    return pd.factorize(v.genus.astype(str) + "|" + grade[TRAITS[0]]
                        + "|" + grade[TRAITS[1]], sort=True)[0]


def perm_test(M, A, D, strata, prep, scale, X, j, observed, draws, seed):
    """Globally shuffle A once per species within strata; repeated islands share it."""
    # The regression slope is linear in its response on the fixed design.
    slope_weights = np.linalg.pinv(X)[j]
    groups = [np.flatnonzero(strata == g) for g in np.unique(strata)]
    movable = [ix for ix in groups if len(ix) > 1 and np.ptp(A[ix]) > 0]
    informative = [ix for ix in movable if np.ptp(D[ix]) > 0]
    # Groups with identical D can never contribute a within-stratum residual.
    rng = np.random.default_rng(seed)
    if not informative or not np.any(prep[3]):
        return {"status": "no_exchangeable_depth_symmetry_support",
                "informative_strata": len(informative),
                "informative_species": sum(map(len, informative)),
                "p": None}
    values = np.empty(draws)
    for k in range(draws):
        aa = A.copy()
        for ix in informative:
            aa[ix] = rng.permutation(A[ix])
        e = expected_deep(prep, aa, M.shape[0])
        y = ((M @ (aa * D)) - e) * scale
        values[k] = slope_weights @ y
    if np.ptp(values) < 1e-14:
        return {"status": "degenerate_null", "p": None,
                "informative_strata": len(informative),
                "informative_species": sum(map(len, informative))}
    p = (1 + int(np.count_nonzero(np.abs(values) >= abs(observed) - 1e-14))) / (draws + 1)
    return {
        "status": "fit", "draws": draws, "p": p,
        "null_mean": float(values.mean()),
        "null_range_025_975": np.quantile(values, [.025, .975]).tolist(),
        "informative_strata": len(informative),
        "informative_species": sum(map(len, informative)),
        "warning": "Species-coherent conditional randomization; assumes exchangeability within genus and joint evidence grade; local margins are recalculated, not simultaneously fixed under randomization.",
    }


def panel(f, ledger, geo, scope, nmin, parent):
    v = trait_pairs(ledger, scope == "high_only")
    grade0 = joint_grades(ledger, v, scope == "high_only")
    d, M, C, g0 = make_frame(f, v, geo, scope)
    d = d.loc[d.analysis_regime.eq("tropical") & d.carrier_n.ge(nmin)
              & d.deep_n.ge(5) & d.shallow_n.ge(5)].sort_values("mi")
    prior = next(x for x in parent["panels"] if
                 x["scope"] == scope and x["region"] == "tropical"
                 and x["endpoint"] == "conditional_difference" and x["nmin"] == nmin)
    b = d.spatial_block.to_numpy(str)
    n, blocks = len(d), len(set(b))
    if n != prior["n"] or blocks != prior["blocks"]:
        raise ValueError("original support not reproduced")
    out = {"scope": scope, "nmin": nmin, "islands": n, "blocks": blocks,
           "parent_status": prior["status"]}
    if n < 50 or blocks < 10:
        return out | {"status": "insufficient_support"}
    X, j = design(d)
    if X is None:
        return out | {"status": "constant_design"}
    m = M[d.mi.to_numpy(int)]
    ix = np.flatnonzero(np.asarray(m.sum(0)).ravel() > 0)
    m = m[:, ix]
    A = v.A.to_numpy(float)[ix]
    D = 1 - v.O.to_numpy(float)[ix]
    genus = pd.factorize(g0[ix], sort=True)[0]
    strata = pd.factorize(grade0[ix], sort=True)[0]
    y, prep, scale = components(m, A, D, genus, strata)
    fit = estimate(X, y, b, j)
    if fit is None:
        return out | {"status": "rank_or_df_failure"}
    betas = {k: fit[str(i)] for i, k in enumerate(METRICS)}
    if prior["status"] != "fit" or not np.isclose(
        betas["observed"]["beta"], prior["beta"], atol=1e-10, rtol=0
    ):
        raise ValueError("original slope not reproduced")
    if not np.isclose(betas["observed"]["beta"],
                      sum(betas[k]["beta"] for k in
                          ("genus_expected", "grade_increment", "within_genus_grade")),
                      atol=1e-10, rtol=0):
        raise ValueError("slope partition mismatch")
    null = perm_test(m, A, D, strata, prep, scale, X, j,
                     betas["within_genus_grade"]["beta"], DRAWS, SEED)
    return out | {
        "status": "fit", "species": len(ix), "genus_strata": int(genus.max()) + 1,
        "joint_grade_strata": int(strata.max()) + 1,
        "islands_with_mixed_depth_stratum": int(np.unique(prep[1][prep[3]]).size),
        "metrics": betas,
        "species_coherent_null": null,
        "interval_note": "Spatial-cluster t intervals are pointwise; repeated species across islands motivate the separate coherent permutation. Neither is an independent validation of the earlier discovery.",
    }


def self_test():
    # Pooled association can arise entirely from group composition.
    m = sparse.csr_matrix(np.ones((1, 8)))
    A = np.array([1, 1, 1, 1, 0, 0, 0, 0], dtype=float)
    D = np.array([1, 1, 1, 0, 1, 0, 0, 0], dtype=float)
    g = np.repeat([0, 1], 4)
    y, _, _ = components(m, A, D, g, g)
    assert np.allclose(y[0], [.5, .5, 0, .5, 0])
    # Same margins, different coupling within the same group.
    A = np.array([1, 1, 0, 0, 1, 1, 0, 0], dtype=float)
    D = A.copy()
    y, prep, scale = components(m, A, D, g, g)
    assert np.allclose(y[0], [1, 0, 0, 0, 1])
    # Joint evidence grade can reproduce an association unexplained by genus.
    y, _, _ = components(m, A, D, np.zeros(8, int), A.astype(int))
    assert np.allclose(y[0], [1, 0, 1, 1, 0])
    # Exact enumeration: the within-stratum residual expectation is zero.
    from itertools import combinations
    yy = []
    for ix in combinations(range(8), 4):
        a = np.zeros(8)
        a[list(ix)] = 1
        y, _, _ = components(m, a, D, np.zeros(8, int), np.zeros(8, int))
        yy.append(y[0, -1])
    assert abs(np.mean(yy)) < 1e-12
    # Fixed seed and shared species assignment in overlapping island rosters.
    M = sparse.csr_matrix([[1,1,1,1,1,0,0,0],
                           [1,0,1,1,1,1,0,0],
                           [1,1,0,1,1,0,1,1]], dtype=float)
    G = np.zeros(8, int)
    y, prep, scale = components(M, A, D, G, G)
    X = np.column_stack([np.ones(3), np.arange(3)])
    obs = float(np.linalg.pinv(X)[1] @ y[:, -1])
    p1 = perm_test(M, A, D, G, prep, scale, X, 1, obs, 99, SEED)
    p2 = perm_test(M, A, D, G, prep, scale, X, 1, obs, 99, SEED)
    assert p1 == p2 and p1["status"] == "fit" and 0 < p1["p"] <= 1
    print("FLORAL_STRATUM_SELF_TEST_PASS")


def main():
    p = argparse.ArgumentParser()
    for k in ("flora", "ledger", "geography", "parent", "out"):
        p.add_argument("--" + k, type=Path)
    p.add_argument("--self-test", action="store_true")
    a = p.parse_args()
    if a.self_test:
        self_test()
        return
    for k, h in LOCK.items():
        if digest(getattr(a, k)) != h:
            raise ValueError("source hash mismatch " + k)
    parent = json.loads(a.parent.read_text())
    if parent["run"] != 45579499069:
        raise ValueError("parent run mismatch")
    f = pd.read_csv(a.flora, dtype=str, usecols=["island_id", "accepted_species", "origin_status"])
    ledger = pd.read_csv(a.ledger, dtype=str)
    geo = pd.read_csv(a.geography, dtype={"island_id": str, "spatial_block": str})
    if len(f) != 1039757 or f.island_id.nunique() != 4505 or len(geo) != 8264:
        raise ValueError("source census mismatch")
    if f.duplicated(["island_id", "accepted_species"]).any() or not geo.island_id.is_unique:
        raise ValueError("duplicate input")
    rows = []
    for scope, nmin in PANELS:
        x = panel(f, ledger, geo, scope, nmin, parent)
        rows.append(x)
        print("FLORAL_STRATUM_PANEL " + json.dumps(x, allow_nan=False), flush=True)
    out = {"owner": "island", "source_hashes": LOCK, "parent_run": 45579499069,
           "primary": "tropical all_recorded n>=10, one two-sided coherent-null test on within-genus-grade conditional difference slope",
           "panels": rows, "metric_order": list(METRICS),
           "warning": "Outcome-motivated diagnosis of the previously observed tropical pattern. Exact local expectation is not causal attribution or phylogenetic correction; grade and native sensitivities are nested, not replication. No submission promotion."}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    print("FLORAL_STRATUM_COMPLETE 4", flush=True)


if __name__ == "__main__":
    main()
