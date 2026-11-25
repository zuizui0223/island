#!/usr/bin/env python3
"""Run the frozen range-map nectarivore x floral-form test inside the protected island CI.

Inputs (all stay inside protected storage):
  --incidence  directory written by scripts/range_map_nectarivore_incidence.py
  --islands    CSV, one row per GIFT island linked to a GSHHG polygon:
               island (GSHHG id), group (archipelago label, or a 10-degree block id), region
               (tropics / north_midlatitude / south_extratropical), dist_km, lat
  --plants     CSV, one row per island x Direct-scored plant species:
               island, species, genus, grade (high/medium), status_scope, and 0/1 endpoint columns

Output: one anonymous JSON summary (no island, plant or bird identifiers). Contract:
docs/RANGE_MAP_NECTARIVORE_CONTRACT_20261125.md. Thresholds are frozen; do not tune them.
"""
import argparse
import json
import pathlib

import numpy as np
import pandas as pd
from scipy import sparse

from fd.nectarivore_flora_test import (
    MIN_SCORED,
    decision,
    plant_block_expectation,
    run_clade_replication,
    run_test,
)
from fd.range_map_nectarivore import CLADES, SPECIALIST

PRIMARY = dict(endpoint="tubular", region="tropics", scope="all_recorded", variant="primary", cover=0.5, grouping="group")
SECONDARY_ENDPOINTS = ("deep", "extended_access", "red_extended")
SECONDARY_REGIONS = ("north_midlatitude", "south_extratropical")
SENSITIVITY = [
    dict(variant="primary", cover=0.1),
    dict(variant="primary", cover=0.0),
    dict(variant="buffer1km", cover=0.5),
    dict(variant="nonbreeding", cover=0.5),
    dict(grouping="block10"),
    dict(scope="native"),
]


def load_incidence(d: pathlib.Path, variant: str, cover: float):
    c = sparse.load_npz(d / f"cover_{variant}.npz").tocsr()
    m = c.copy()
    m.data = (m.data > 0) if cover <= 0 else (m.data >= cover - 1e-12)
    m.eliminate_zeros()
    isl = pd.read_csv(d / "islands.csv")
    sp = pd.read_csv(d / f"species_{variant}.csv")
    if c.shape != (len(isl), len(sp)):
        raise ValueError("incidence shape does not match islands/species tables")
    return m.astype(bool), isl, sp


def analysis_frame(islands, plants, inc_islands, endpoint, region, scope):
    pl = plants[(plants["status_scope"] == scope) if scope != "all_recorded" else slice(None)]
    meta = islands[islands["region"] == region]
    pl = pl[pl["island"].isin(meta["island"]) & pl[endpoint].notna()]
    shares = plant_block_expectation(
        pl.assign(trait=pl[endpoint].astype(int)), block_cols=("genus", "grade")
    ).reset_index()
    d = meta.merge(shares, on="island").merge(inc_islands, on="island", validate="one_to_one")
    d = d[(d["n"] >= MIN_SCORED) & (d["area_km2"] >= 1)].reset_index(drop=True)
    return d


def one_run(args, islands, plants, endpoint, region, scope, variant, cover, grouping, n_boot, n_perm, clades=False):
    m, inc_isl, sp = load_incidence(args.incidence, variant, cover)
    row_of = {k: i for i, k in enumerate(inc_isl["island"])}
    d = analysis_frame(islands, plants, inc_isl, endpoint, region, scope)
    if grouping != "group":
        d = d.assign(group=d[grouping])
    rows = [row_of[k] for k in d["island"]]
    inc = m[rows]
    diet = sp["diet_class"].to_numpy()
    r = run_test(
        d, inc, diet, {"genus": sp["genus_block"].to_numpy(), "family": sp["family_block"].to_numpy()},
        n_boot=n_boot, n_perm=n_perm,
    )
    out = {
        "endpoint": endpoint, "region": region, "scope": scope, "variant": variant, "cover": cover,
        "grouping": grouping, "decision": decision(r), **r.as_dict(),
        "specialist_species_present": int((np.asarray(inc[:, diet == SPECIALIST].sum(axis=0)).ravel() > 0).sum()),
        "bird_genera_with_mixed_diet_present": mixed_genera(inc, diet, sp["genus_block"].to_numpy()),
    }
    if clades and r.status == "ok":
        out["clade_replication"] = run_clade_replication(d, inc, diet, sp["clade"].to_numpy(), CLADES, n_boot=n_boot)
    return out


def mixed_genera(inc, diet, genus):
    present = np.asarray(inc.sum(axis=0)).ravel() > 0
    g = pd.DataFrame({"g": genus[present], "d": diet[present]})
    return int((g.groupby("g")["d"].nunique() > 1).sum())


def holm(ps):
    order = np.argsort(ps)
    adj = np.empty(len(ps))
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (len(ps) - rank) * ps[i])
        adj[i] = min(1.0, running)
    return adj.tolist()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--incidence", type=pathlib.Path, required=True)
    p.add_argument("--islands", type=pathlib.Path, required=True)
    p.add_argument("--plants", type=pathlib.Path, required=True)
    p.add_argument("--n-boot", type=int, default=1999)
    p.add_argument("--n-perm", type=int, default=1999)
    p.add_argument("--out", type=pathlib.Path, required=True)
    a = p.parse_args()
    islands = pd.read_csv(a.islands)
    plants = pd.read_csv(a.plants)
    if "block10" not in islands:
        islands["block10"] = (
            (np.floor(islands["lat"] / 10)).astype(int).astype(str) + "_" + (np.floor(islands["lon"] / 10)).astype(int).astype(str)
            if "lon" in islands else islands["group"]
        )
    kw = dict(n_boot=a.n_boot, n_perm=a.n_perm)
    primary = one_run(a, islands, plants, **PRIMARY, **kw, clades=True)
    secondary = []
    for e in SECONDARY_ENDPOINTS:
        secondary.append(one_run(a, islands, plants, **{**PRIMARY, "endpoint": e}, **kw))
    for rg in SECONDARY_REGIONS:
        secondary.append(one_run(a, islands, plants, **{**PRIMARY, "region": rg}, **kw))
    ps = [
        s["permutation"].get("residual", {}).get("genus", {}).get("p_one_sided", 1.0) if s["status"] == "ok" else 1.0
        for s in secondary
    ]
    for s, q in zip(secondary, holm(ps)):
        s["holm_p_residual_genus"] = q
    sensitivity = [one_run(a, islands, plants, **{**PRIMARY, **s}, **kw) for s in SENSITIVITY]
    summary = {
        "contract": "docs/RANGE_MAP_NECTARIVORE_CONTRACT_20261125.md",
        "primary": primary,
        "secondary": secondary,
        "sensitivity": sensitivity,
        "note": "Range-map co-occurrence, species-wide floral records. Not visitation or pollen transfer.",
    }
    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(summary, indent=2, allow_nan=False))
    print("PRIMARY_DECISION", primary["decision"], flush=True)


if __name__ == "__main__":
    main()
