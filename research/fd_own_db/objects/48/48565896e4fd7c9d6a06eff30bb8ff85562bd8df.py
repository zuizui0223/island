"""Outcome-blind maximally geography-disjoint island origin pair matching.

Motivation: original 37 oceanic/continental pairs belong to only five broad
spatial connected components. Re-select pairs from *all* geological-typed
All+Direct complete-floral islands, maximizing number of matches under the
SAME area, distance, climate, name-support and coverage calipers, while
preventing any shared geographic spatial_block and/or source archipelago
across separate matched pairs.

This evaluates an alternate and smaller population, NOT independent replication
of the original pair cohort. Lack of an archived source-archipelago label is
treated as unknown per island, never as shared ancestry. Geological origin is
not randomized and real lineage/observation dependence may remain.

Only aggregates, audit parameters, and input hashes are released. No matched
island or source-region identifiers leave this local pipeline.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
import json
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import lil_matrix

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .pipeline import REF_INPUT, _write_csv, file_sha256
from .repertoire import (
    REGIONS, COVARIATES, MATCH_LOG_RATIO, MAX_COVERAGE_DIFFERENCE,
    MAX_CLIMATE_Z_DISTANCE, MAX_CLIMATE_PC3_DIFFERENCE,
    _features, _paired_swap_statistic, _upper_tail,
)
from .g1_fixed_margins import _bh

GROUPINGS = ("spatial_block", "weigelt_archipelago", "both", "both_verified")
ENDPOINTS = ("exact_profile_L2", "categorical_similarity_kernel")
SCOPES = ("all", "direct")
TEST_FIELDS = (
    "grouping", "scope", "endpoint", "n_pairs",
    "n_northern_midlatitude", "n_tropical",
    "observed", "null_median", "null_q025", "null_q975",
    "p_exact_upper", "n_exact_sign_configurations",
    "minimum_attainable_p", "q_bh_posthoc16", "q_by_posthoc16",
)

AUDIT_FIELDS = (
    "grouping", "n_candidate_edges", "n_selected_pairs",
    "n_continental_islands", "n_oceanic_islands",
    "n_northern_midlatitude", "n_tropical",
    "n_distinct_spatial_blocks_used",
    "n_distinct_labelled_archipelagos_used", "n_unknown_archipelago_endpoints",
    "max_pair_cost", "median_pair_cost", "solver_status",
)


def _missing_label(value) -> bool:
    return not isinstance(value, str) or value.strip().lower() in {
        "", "na", "n/a", "nan", "null", "none", "unknown", "unresolved", "-"
    }


def _geographic_resources(edge: dict, grouping: str) -> set[str]:
    if grouping not in GROUPINGS:
        raise ValueError("Unsupported geographic independence condition")
    keys = {
        f"island:{edge['continental_id']}",
        f"island:{edge['oceanic_id']}",
    }
    if grouping in ("spatial_block", "both", "both_verified"):
        for side in ("c", "o"):
            b = edge[f"{side}_spatial_block"]
            if _missing_label(b):
                raise ValueError("Missing mandatory spatial block")
            keys.add(f"block:{b}")
    if grouping in ("weigelt_archipelago", "both", "both_verified"):
        for side, name in (("c", "continental"), ("o", "oceanic")):
            label = edge[f"{side}_archipelago"]
            if _missing_label(label):
                if grouping == "both_verified":
                    raise ValueError("Cannot verify source archipelago for selected island")
                # A missing source label is NOT a known independent archipelago.
                # Its own island cannot be reused, but it creates no spurious
                # shared grouping across other islands.
                keys.add(f"unknown_archipelago_for:{edge[name+'_id']}")
            else:
                keys.add(f"archipelago:{label}")
    return keys


def _build_edges(g0_all: list[dict], g0_direct: list[dict],
                 geography: dict[str, dict]) -> tuple[list[dict], dict]:
    """Recreate EXACT original PR4 geographic/caliper eligibility.

    Do not select or inspect any flower-colour/form outcome in this function.
    """
    direct_ids = {r["island_id"] for r in g0_direct}
    candidates = []
    for r in g0_all:
        island = r["island_id"]
        geo = geography.get(island, {})
        origin = geo.get("geological_origin")
        if origin not in ("continental", "oceanic") or island not in direct_ids:
            continue
        if r["analysis_regime"] not in REGIONS:
            continue
        features = _features(r, geo)
        if not np.isfinite(features).all():
            continue
        candidates.append({
            "island_id": island,
            "regime": r["analysis_regime"],
            "origin": origin, "features": features,
            "spatial_block": geo.get("spatial_block"),
            "archipelago": geo.get("weigelt_archipelago"),
        })
    edges = []
    supports = []
    for regime in REGIONS:
        population = [x for x in candidates if x["regime"] == regime]
        continental = sorted(
            (x for x in population if x["origin"] == "continental"),
            key=lambda x: x["island_id"],
        )
        oceanic = sorted(
            (x for x in population if x["origin"] == "oceanic"),
            key=lambda x: x["island_id"],
        )
        if not continental or not oceanic:
            supports.append({"region": regime, "n_edges": 0})
            continue
        a = np.vstack([x["features"] for x in continental])
        b = np.vstack([x["features"] for x in oceanic])
        diff = a[:, None, :] - b[None, :, :]
        sd = np.maximum(np.std(np.vstack((a, b)), axis=0),
                        [0.35, 0.35, 0.25, 0.04, 0.1, 0.1, 0.1, 0.1])
        standardized = diff / sd
        costs = np.linalg.norm(standardized * np.asarray(
            [1, 1, 1, 1, 1, 1, 1.6, 1]), axis=2)
        allowed = (
            (np.abs(diff[:, :, 0]) <= MATCH_LOG_RATIO)
            & (np.abs(diff[:, :, 1]) <= MATCH_LOG_RATIO)
            & (np.abs(diff[:, :, 2]) <= MATCH_LOG_RATIO)
            & (np.abs(diff[:, :, 3]) <= MAX_COVERAGE_DIFFERENCE)
            & (np.linalg.norm(standardized[:, :, 4:], axis=2) <=
               MAX_CLIMATE_Z_DISTANCE)
            & (np.abs(diff[:, :, 6]) <= MAX_CLIMATE_PC3_DIFFERENCE)
        )
        for i, j in zip(*np.where(allowed)):
            ca, oc = continental[int(i)], oceanic[int(j)]
            edges.append({
                "continental_id": ca["island_id"],
                "oceanic_id": oc["island_id"],
                "region": regime,
                "c_spatial_block": ca["spatial_block"],
                "o_spatial_block": oc["spatial_block"],
                "c_archipelago": ca["archipelago"],
                "o_archipelago": oc["archipelago"],
                "cost": float(costs[i, j]),
            })
        supports.append({
            "region": regime,
            "n_eligible_continental": len(continental),
            "n_eligible_oceanic": len(oceanic),
            "n_edges": int(allowed.sum()),
        })
    if not edges:
        raise ValueError("No valid island origin matching edges")
    edges.sort(key=lambda x: (x["region"], x["continental_id"], x["oceanic_id"]))
    return edges, {
        "source_candidates_total": len(candidates),
        "per_region": supports,
        "calipers": {
            "max_abs_log_area_ratio": float(MATCH_LOG_RATIO),
            "max_abs_log_isolation_ratio": float(MATCH_LOG_RATIO),
            "max_abs_log_complete_name_ratio": float(MATCH_LOG_RATIO),
            "max_abs_trait_coverage_difference": MAX_COVERAGE_DIFFERENCE,
            "max_climate_PC_z_norm": MAX_CLIMATE_Z_DISTANCE,
            "max_abs_raw_climate_PC3": MAX_CLIMATE_PC3_DIFFERENCE,
            "covariate_weighting": "original fd PR4 climate PC3 weight 1.6",
        },
    }


def independent_pair_matching(
    edges: list[dict], grouping: str, *,
    assert_optimal: bool = True,
) -> tuple[list[dict], dict]:
    """Binary set-packing MILP; cardinality dominates covariate distance.

    Each island and every source geographic resource may occur in <=1 match.
    This delivers maximal observed geographic disjointness under the given
    masks, NOT a guarantee of independent island lineage history.
    """
    if not edges:
        raise ValueError("No candidate origin-contrasting pairs")
    if grouping == "both_verified":
        edges = [e for e in edges
                 if not _missing_label(e["c_archipelago"])
                 and not _missing_label(e["o_archipelago"])]
    if not edges:
        raise ValueError("No complete-archipelago source pairs")
    resources = [_geographic_resources(e, grouping) for e in edges]
    resource_names = sorted(set().union(*resources))
    index = {k: i for i, k in enumerate(resource_names)}
    matrix = lil_matrix((len(index), len(edges)), dtype=float)
    for j, used in enumerate(resources):
        for label in used:
            matrix[index[label], j] = 1.
    costs = np.asarray([e["cost"] for e in edges], dtype=float)
    if np.any(~np.isfinite(costs)):
        raise ValueError("Invalid candidate covariate dissimilarities")
    # Strong >100 point cardinality bonus. The source cost is < few units.
    objective = -1000.0 + costs / 10
    solution = milp(
        c=objective, integrality=np.ones(len(edges), dtype=np.int32),
        bounds=Bounds(0, 1),
        constraints=LinearConstraint(matrix.tocsr(), lb=0, ub=1),
        options={"mip_rel_gap": 0.001, "time_limit": 70},
    )
    if solution.x is None:
        raise ValueError("No geographic set-packing solution found")
    if assert_optimal and solution.status != 0:
        raise ValueError("Matching optimality not verified")
    selected = [edge for edge, flag in zip(edges, solution.x) if flag > .5]
    occupied = set()
    for e in selected:
        labels = _geographic_resources(e, grouping)
        if labels & occupied:
            raise AssertionError("Two independent matches reuse a geographic source")
        occupied.update(labels)
    spatial = set()
    arch = set()
    unknown = 0
    for e in selected:
        spatial.update([e["c_spatial_block"], e["o_spatial_block"]])
        for name in ("c_archipelago", "o_archipelago"):
            a = e[name]
            if _missing_label(a):
                unknown += 1
            else:
                arch.add(a)
    report = {
        "grouping": grouping, "n_candidate_edges": len(edges),
        "n_selected_pairs": len(selected),
        "n_continental_islands": len({e["continental_id"] for e in selected}),
        "n_oceanic_islands": len({e["oceanic_id"] for e in selected}),
        "n_northern_midlatitude": sum(
            e["region"] == "northern_midlatitude" for e in selected),
        "n_tropical": sum(e["region"] == "tropical" for e in selected),
        "n_distinct_spatial_blocks_used": len(spatial),
        "n_distinct_labelled_archipelagos_used": len(arch),
        "n_unknown_archipelago_endpoints": unknown,
        "max_pair_cost": float(max(e["cost"] for e in selected)) if selected else None,
        "median_pair_cost": float(np.median([e["cost"] for e in selected]))
                            if selected else None,
        "solver_status": int(solution.status),
    }
    return selected, report


def _exact_pair_null(delta: np.ndarray, kernel: np.ndarray) -> dict:
    n = len(delta)
    if not 8 <= n <= 16:
        raise ValueError("Exact independently matched pair count outside [8,16]")
    gram = delta @ kernel @ delta.T
    obs = float(max(0, gram.sum()/n**2))
    null = []
    for signs in product((-1.0, 1.0), repeat=n):
        v = np.asarray(signs)
        null.append(float(max(0, v @ gram @ v / n**2)))
    values = np.asarray(null)
    tol = 128*np.finfo(float).eps * max(1, abs(obs))
    return {
        "observed": obs,
        "null_median": float(np.median(values)),
        "null_q025": float(np.quantile(values, .025)),
        "null_q975": float(np.quantile(values, .975)),
        "p_exact_upper": float(np.mean(values >= obs-tol)),
        "n_exact_sign_configurations": len(values),
        "minimum_attainable_p": 2/len(values),
    }


def _profile_delta(edges: list[dict], carrier: dict[str,dict[int,int]],
                   space: ProfileSpace) -> np.ndarray:
    cache = {}
    for e in edges:
        for origin in ("continental", "oceanic"):
            key = e[origin+"_id"]
            if key not in carrier:
                raise ValueError("Missing island profile carrier input")
            if key not in cache:
                selected = carrier[key]
                if not selected:
                    raise ValueError("No recorded floral profile for matched island")
                vec = np.zeros(len(space.ids), dtype=float)
                vec[[space.lookup[p] for p in selected]] = 1/len(selected)
                cache[key] = vec
    return np.vstack([
        cache[e["oceanic_id"]] - cache[e["continental_id"]]
        for e in edges
    ])


def run(input_root: Path, g0_root: Path, crosswalk: Path, output: Path) -> dict:
    input_root, g0_root, crosswalk, output = map(Path,
        (input_root, g0_root, crosswalk, output))
    output.mkdir(parents=True, exist_ok=True)
    all_g0 = list(csv_rows(g0_root/"all_islands.csv"))
    direct_g0 = list(csv_rows(g0_root/"direct_islands.csv"))
    geo = {r["island_id"]:r for r in csv_rows(crosswalk)}
    edges, candidate_audit = _build_edges(all_g0,direct_g0,geo)
    space = ProfileSpace(load_profiles(input_root/REF_INPUT))
    kernel = 1-space.matrix
    carr = {scope: load_carriers(input_root/"floral_space_20261106"/scope/
                                "island_profile_carriers.csv.gz")
            for scope in SCOPES}
    report, tests = [], []
    for grouping in GROUPINGS:
        selected, diag = independent_pair_matching(edges, grouping)
        if not 8 <= len(selected) <= 16:
            raise ValueError("Expected independent match count 8..16; revise exact enumeration")
        report.append(diag)
        for scope in SCOPES:
            delta = _profile_delta(selected,carr[scope],space)
            for endpoint, matrix in (
                ("exact_profile_L2",np.eye(len(space.ids))),
                ("categorical_similarity_kernel",kernel),
            ):
                stat = _exact_pair_null(delta,matrix)
                tests.append({
                    "grouping":grouping, "scope":scope,"endpoint":endpoint,
                    "n_pairs":len(selected),
                    "n_northern_midlatitude":diag["n_northern_midlatitude"],
                    "n_tropical":diag["n_tropical"],
                    **stat,
                })
    # One explicitly post-hoc 16-comparison family: 4 correlated pairing
    # definitions × 2 evidence scopes × 2 category endpoints.
    q = _bh([r["p_exact_upper"] for r in tests])
    harmonic = sum(1/i for i in range(1,len(tests)+1))
    for row, adj in zip(tests,q):
        row["q_bh_posthoc16"]=float(adj)
        row["q_by_posthoc16"]=float(min(1,adj*harmonic))
    if len(report)!=4 or len(tests)!=16:
        raise AssertionError("Unexpected geographical matching comparison family")
    _write_csv(output/"independent_matching_support.csv",report,AUDIT_FIELDS)
    _write_csv(output/"independent_origin_repertoire_tests.csv",tests,TEST_FIELDS)
    manifest={
        "schema":"fd_geography_disjoint_matching_v1",
        "status":"posthoc_geographic_disjointness_stress_test",
        "source_population":"All+Direct common G0 complete colour/form coverage, independent geological origin labels",
        "geographic_calipers_identical_to":"fd PR4, no floral profile outcome in pair assignment",
        "edge_audit":candidate_audit,
        "matching":"global source-resource-disjoint binary MILP; cardinality prioritized over standardized covariate distance",
        "grouping_types":list(GROUPINGS),
        "n_test_family":len(tests),
        "verified_archipelago_sensitivity": "both_verified requires both island endpoints have nonmissing published archipelago labels; missing never means verified independence",
        "p_exact":"all 2^n pair sign flips, no Monte Carlo correction; squared statistic has ± symmetry",
        "interpretation":"alternate, geographically less dependent matched subsample, NOT validation of lineage/phylogenetic independence or true geological treatment randomization",
        "data_publication":"NO island IDs, source block IDs, pair membership, original taxon names or restricted genus evidence released",
        "source_sha256":{
            "all_g0":file_sha256(g0_root/"all_islands.csv"),
            "direct_g0":file_sha256(g0_root/"direct_islands.csv"),
            "geological_crosswalk":file_sha256(crosswalk),
            "profile_reference":file_sha256(input_root/REF_INPUT),
        },
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"status":manifest["status"],
                      "matching_support":report,
                      "tests":[{"grouping":r["grouping"],"scope":r["scope"],
                                "endpoint":r["endpoint"],"n":r["n_pairs"],
                                "p":r["p_exact_upper"],"q":r["q_bh_posthoc16"]}
                               for r in tests]},ensure_ascii=False))
    return manifest


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    for label in ("input", "g0", "crosswalk", "output"):
        p.add_argument("--"+label,type=Path,required=True)
    args=p.parse_args()
    run(args.input,args.g0,args.crosswalk,args.output)


if __name__=="__main__":
    main()
