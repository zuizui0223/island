"""Geological-origin floral repertoire with outcome-blind geographic matching.

Public inputs: source-pinned island x profile carrier counts and independently
audited geological origin crosswalk. No species x trait provenance is copied.

Primary matched contrast is exact within geographic regime, with calipers
on area, distance, complete-name support, coverage and climate. It answers
whether documented floral-profile repertoires differ in matched recorded
floras; it does NOT condition on true colonist source pools or genus identity.

Paired sign swaps are a conditional/exploratory reference, not randomization of
the actual geological history. Profile-specific exact tests use McNemar's
discordant-pair binomial reference; FDR includes all 496 profiles.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .g1_fixed_margins import _bh
from .pipeline import REF_INPUT, _write_csv, file_sha256

ORIGINS = ("continental", "oceanic")
MIN_REGION_PAIRS = 8
MATCH_LOG_RATIO = math.log(3.0)
MAX_COVERAGE_DIFFERENCE = 0.25
MAX_CLIMATE_Z_DISTANCE = 3.0
MAX_CLIMATE_PC3_DIFFERENCE = 2.0
PERMUTATIONS = 9999
PROFILE_FAMILY = 496
REGIONS = ("northern_midlatitude", "tropical")
COVARIATES = ("log_area", "log_isolation", "log_names", "joint_coverage",
              "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4")


def _numeric(value) -> float:
    try:
        x = float(value)
        return x if math.isfinite(x) else float("nan")
    except (ValueError, TypeError):
        return float("nan")


def _features(g0: dict, geology: dict) -> np.ndarray:
    area = _numeric(g0.get("area_km2"))
    distance = _numeric(g0.get("distance_to_continent_km"))
    names = _numeric(g0.get("n_names"))
    cover = _numeric(g0.get("joint_coverage"))
    if area <= 0 or distance < 0 or names <= 0 or not np.isfinite(cover):
        return np.full(len(COVARIATES), np.nan)
    return np.asarray(
        [math.log(area), math.log1p(distance), math.log(names), cover]
        + [_numeric(geology.get(f"climate_pc{i}")) for i in range(1, 5)],
        dtype=float,
    )


def match_origin_islands(g0_all: list[dict], g0_direct: list[dict],
                         crosswalk: dict[str, dict], *,
                         coverage_max: float = MAX_COVERAGE_DIFFERENCE,
                         climate_pc3_max: float = MAX_CLIMATE_PC3_DIFFERENCE,
                         max_pair_geographic_km: float | None = None
                         ) -> tuple[list[dict], dict]:
    """Match once using All covariates, retain only islands with Direct support.

    The same fixed pairs are then analyzed for All and Direct; nested scopes
    are NOT counted as independent replications.
    """
    direct_ids = {row["island_id"] for row in g0_direct}
    eligible = []
    for row in g0_all:
        island = row["island_id"]
        geo = crosswalk.get(island, {})
        if island not in direct_ids or geo.get("geological_origin") not in ORIGINS:
            continue
        features = _features(row, geo)
        if np.isfinite(features).all():
            eligible.append({
                "island_id": island,
                "analysis_regime": row["analysis_regime"],
                "origin": geo["geological_origin"],
                "features": features,
            })
    pairs, regional = [], []
    for region in REGIONS:
        cohort = [row for row in eligible if row["analysis_regime"] == region]
        continental = sorted((r for r in cohort if r["origin"] == "continental"),
                             key=lambda x: x["island_id"])
        oceanic = sorted((r for r in cohort if r["origin"] == "oceanic"),
                         key=lambda x: x["island_id"])
        if not continental or not oceanic:
            regional.append({"region": region, "n_matched": 0, "status": "no_overlap"})
            continue
        a = np.vstack([row["features"] for row in continental])
        b = np.vstack([row["features"] for row in oceanic])
        pooled = np.vstack([a, b])
        sd = np.maximum(np.std(pooled, axis=0), [0.35, 0.35, 0.25, 0.04, 0.1, 0.1, 0.1, 0.1])
        delta = a[:, None, :] - b[None, :, :]
        standardized = delta / sd
        cost = np.linalg.norm(standardized * np.array([1, 1, 1, 1, 1, 1, 1.6, 1]), axis=2)
        valid = (
            (np.abs(delta[:, :, 0]) <= MATCH_LOG_RATIO)
            & (np.abs(delta[:, :, 1]) <= MATCH_LOG_RATIO)
            & (np.abs(delta[:, :, 2]) <= MATCH_LOG_RATIO)
            & (np.abs(delta[:, :, 3]) <= coverage_max)
            & (np.linalg.norm(standardized[:, :, 4:], axis=2) <= MAX_CLIMATE_Z_DISTANCE)
            & (np.abs(delta[:, :, 6]) <= climate_pc3_max)
        )
        if max_pair_geographic_km is not None:
            # Geographic source-neighborhood sensitivity: reduce comparisons
            # spanning completely different botanical continents, without
            # claiming this identifies actual historical colonist pools.
            if max_pair_geographic_km <= 0:
                raise ValueError("Geographic pairing radius must be positive")
            lat_c = np.radians([_numeric(crosswalk[r["island_id"]]["island_latitude"])
                                for r in continental])
            lat_o = np.radians([_numeric(crosswalk[r["island_id"]]["island_latitude"])
                                for r in oceanic])
            lon_c = np.radians([_numeric(crosswalk[r["island_id"]]["island_longitude"])
                                for r in continental])
            lon_o = np.radians([_numeric(crosswalk[r["island_id"]]["island_longitude"])
                                for r in oceanic])
            cos_angle = (np.sin(lat_c[:, None]) * np.sin(lat_o[None, :])
                         + np.cos(lat_c[:, None]) * np.cos(lat_o[None, :])
                         * np.cos(lon_c[:, None] - lon_o[None, :]))
            km = 6371.0088 * np.arccos(np.clip(cos_angle, -1.0, 1.0))
            valid &= np.isfinite(km) & (km <= max_pair_geographic_km)
        # One distinct dummy per continental island: maximize valid matches
        # before minimizing their standardized covariate dissimilarity.
        matrix = np.hstack([np.where(valid, cost, 1e6),
                            np.full((len(a), len(a)), 1e3)])
        ii, jj = linear_sum_assignment(matrix)
        chosen = [(i, j) for i, j in zip(ii, jj) if j < len(b) and valid[i, j]]
        for i, j in chosen:
            pairs.append({
                "region": region,
                "continental_id": continental[i]["island_id"],
                "oceanic_id": oceanic[j]["island_id"],
                "distance": float(cost[i, j]),
            })
        after = np.asarray([b[j] - a[i] for i, j in chosen]) if chosen else np.zeros((0, 8))
        regional.append({
            "region": region, "n_available_continental": len(a),
            "n_available_oceanic": len(b), "n_matched": len(chosen),
            "n_eligible_edges": int(valid.sum()),
            "mean_oceanic_minus_continental": after.mean(axis=0).tolist() if chosen else [],
            "max_abs_smd_pooled": float(np.max(np.abs(after.mean(axis=0) / sd))) if chosen else None,
            "smd_covariates": {c: float(v) for c, v in zip(
                COVARIATES, after.mean(axis=0) / sd)} if chosen else {},
            "status": "eligible_for_main_inference" if len(chosen) >= MIN_REGION_PAIRS
                      else "insufficient_matched_pairs",
        })
    supported = [r["region"] for r in regional
                 if r.get("n_matched", 0) >= MIN_REGION_PAIRS]
    main_pairs = [p for p in pairs if p["region"] in supported]
    if not main_pairs:
        raise ValueError("No geographic regime has sufficient origin-matched island support")
    audit = {
        "n_geological_origin_eligible_with_complete_all_and_direct": len(eligible),
        "n_main_pairs": len(main_pairs),
        "n_distinct_islands": 2 * len(main_pairs),
        "source_regimes_in_main": supported,
        "regional_diagnostics": regional,
        "distance_calipers": {
            "log_area_max": MATCH_LOG_RATIO, "log_isolation_max": MATCH_LOG_RATIO,
            "log_names_max": MATCH_LOG_RATIO,
            "joint_coverage_abs_max": coverage_max,
            "standardized_climate_euclidean_max": MAX_CLIMATE_Z_DISTANCE,
            "climate_pc3_raw_abs_max": climate_pc3_max,
            "max_pair_geographic_km": max_pair_geographic_km,
        },
        "matching": "maximum valid cardinality, minimum weighted standardized distance; no replacement",
        "selection_blind_to_flower_profile_identity": True,
        "not_matched_on": ["genus", "historical source pool", "island species identity"],
    }
    return main_pairs, audit


def _paired_swap_statistic(delta: np.ndarray, kernel: np.ndarray,
                           signs: np.ndarray) -> tuple[float, np.ndarray]:
    """Kernel distance between the two mean profile-frequency repertoires."""
    gram = delta @ kernel @ delta.T
    obs = float(gram.sum() / (len(delta) ** 2))
    null = np.einsum("bi,ij,bj->b", signs, gram, signs, optimize=True) / (len(delta) ** 2)
    return max(obs, 0.0), np.maximum(null, 0.0)


def _upper_tail(observed: float, reference: np.ndarray) -> float:
    tol = 128 * np.finfo(float).eps * max(1.0, abs(observed))
    return float((1 + np.count_nonzero(reference >= observed - tol)) / (len(reference) + 1))


def _two_sided(observed: float, reference: np.ndarray) -> float:
    tol = 128 * np.finfo(float).eps * max(1.0, abs(observed))
    p_lower = (1 + np.count_nonzero(reference <= observed + tol)) / (len(reference) + 1)
    p_upper = (1 + np.count_nonzero(reference >= observed - tol)) / (len(reference) + 1)
    return float(min(1, 2 * min(p_lower, p_upper)))


def _binary_profile_tests(oceanic: np.ndarray, continental: np.ndarray,
                          labels: list[tuple[int, str, str]], scope: str
                          ) -> list[dict]:
    records = []
    for column, (profile, color, form) in enumerate(labels):
        o, c = oceanic[:, column], continental[:, column]
        o_only = int(np.count_nonzero(o & ~c))
        c_only = int(np.count_nonzero(c & ~o))
        n_discordant = o_only + c_only
        p = float(binomtest(o_only, n_discordant, 0.5).pvalue) if n_discordant else 1.0
        records.append({
            "scope": scope, "profile_id": profile,
            "flower_primary_color": color, "floral_form": form,
            "matched_pairs": len(o),
            "oceanic_incidence_fraction": float(o.mean()),
            "continental_incidence_fraction": float(c.mean()),
            "risk_difference_o_minus_c": float(o.mean() - c.mean()),
            "n_discordant": n_discordant,
            "n_oceanic_only": o_only,
            "n_continental_only": c_only,
            "p_exact_pair_swap": p,
        })
    if len(records) != PROFILE_FAMILY:
        raise ValueError("Unexpected profile catalogue: cannot change fixed test family size")
    q = _bh([row["p_exact_pair_swap"] for row in records])
    harmonic = float(sum(1 / i for i in range(1, PROFILE_FAMILY + 1)))
    for i, row in enumerate(records):
        row["q_bh_496"] = float(q[i])
        row["q_by_496"] = float(min(1, q[i] * harmonic))
    return records


def _state_contrasts(oceanic: np.ndarray, continental: np.ndarray,
                     profiles: dict, profile_ids: list[int], scope: str) -> list[dict]:
    records = []
    for dimension, index, key in (
        ("color", 0, "flower_primary_color"),
        ("form", 1, "floral_form"),
    ):
        states = sorted({s for p in profiles.values() for s in p[index]})
        for state in states:
            cols = [j for j, profile in enumerate(profile_ids)
                    if state in profiles[profile][index]]
            o = oceanic[:, cols].any(axis=1)
            c = continental[:, cols].any(axis=1)
            pos = int(np.count_nonzero(o & ~c))
            neg = int(np.count_nonzero(c & ~o))
            n = pos + neg
            records.append({
                "scope": scope, "dimension": dimension, "state": state,
                "n_pairs": len(o), "oceanic_fraction": float(o.mean()),
                "continental_fraction": float(c.mean()),
                "risk_difference": float(o.mean() - c.mean()),
                "discordant_pairs": n, "oceanic_only": pos, "continental_only": neg,
                "p_exact_pair_swap": float(binomtest(pos, n).pvalue) if n else 1.0,
            })
    q = _bh([r["p_exact_pair_swap"] for r in records])
    harmonic = sum(1 / x for x in range(1, len(records) + 1))
    for i, row in enumerate(records):
        row["q_bh_states"] = float(q[i])
        row["q_by_states"] = float(min(1, q[i] * harmonic))
    return records


def run(input_root: Path, g0_root: Path, crosswalk_path: Path, output: Path,
        *, permutations: int = PERMUTATIONS, seed: int = 20261108) -> dict:
    if permutations < 999:
        raise ValueError("At least 999 random swaps required")
    input_root, g0_root, crosswalk_path, output = map(
        Path, (input_root, g0_root, crosswalk_path, output))
    output.mkdir(parents=True, exist_ok=True)
    crosswalk = {r["island_id"]: r for r in csv_rows(crosswalk_path)}
    g0_all = list(csv_rows(g0_root / "all_islands.csv"))
    g0_direct = list(csv_rows(g0_root / "direct_islands.csv"))
    matched, audit = match_origin_islands(g0_all, g0_direct, crosswalk)
    profiles = load_profiles(input_root / REF_INPUT)
    space = ProfileSpace(profiles)
    labels = [(p, "|".join(sorted(profiles[p][0])), "|".join(sorted(profiles[p][1])))
              for p in space.ids]
    kernel = 1.0 - space.matrix
    eigen = np.linalg.eigvalsh(kernel)
    # Average state-set Jaccard similarities yield a positive Tanimoto kernel.
    # Reject an unsupported indefinite kernel rather than silently fake distances.
    if eigen.min() < -1e-7:
        raise ValueError("Floral Jaccard similarity kernel not positive semidefinite")
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(permutations, len(matched)))
    group_records, profile_records, state_records, q_records = [], [], [], []
    source_files = {"profile_reference": input_root / REF_INPUT,
                    "island_geology": crosswalk_path,
                    "all_islands": g0_root / "all_islands.csv",
                    "direct_islands": g0_root / "direct_islands.csv"}
    paired_ids = {
        "continental": [m["continental_id"] for m in matched],
        "oceanic": [m["oceanic_id"] for m in matched],
    }
    for scope in ("all", "direct"):
        carrier_path = input_root / "floral_space_20261106" / scope / "island_profile_carriers.csv.gz"
        source_files[f"{scope}_carriers"] = carrier_path
        carrier = load_carriers(carrier_path)
        oceanic = np.zeros((len(matched), len(space.ids)), dtype=bool)
        continental = np.zeros_like(oceanic)
        qdiff = []
        for i, pair in enumerate(matched):
            c, o = pair["continental_id"], pair["oceanic_id"]
            for key, matrix in ((c, continental), (o, oceanic)):
                if key not in carrier:
                    raise ValueError(f"{scope}: missing matched island {key}")
                indices = [space.lookup[p] for p in carrier[key]]
                matrix[i, indices] = True
            qdiff.append(space.island(carrier[o])["rao_q_profile_equal"]
                         - space.island(carrier[c])["rao_q_profile_equal"])
        qdiff = np.asarray(qdiff, dtype=float)
        sampled = (signs * qdiff[None, :]).mean(axis=1)
        q_records.append({
            "scope": scope, "n_matched_pairs": len(matched),
            "mean_oceanic_minus_continental_rao_q": float(qdiff.mean()),
            "null_q025": float(np.quantile(sampled, 0.025)),
            "null_q975": float(np.quantile(sampled, 0.975)),
            "p_swap_two_sided": _two_sided(float(qdiff.mean()), sampled),
        })
        # Normalize by the number of profiles found on each island, so the
        # multivariate statistic describes repertoire POSITION, not richness.
        po = oceanic.astype(float) / oceanic.sum(axis=1, keepdims=True)
        pc = continental.astype(float) / continental.sum(axis=1, keepdims=True)
        delta = po - pc
        for name, k in (("categorical_similarity_kernel", kernel),
                        ("exact_profile_L2", np.eye(len(space.ids)))):
            observed, ref = _paired_swap_statistic(delta, k, signs)
            group_records.append({
                "scope": scope, "endpoint": name, "n_matched_pairs": len(matched),
                "observed_statistic": observed,
                "null_median": float(np.median(ref)),
                "null_q025": float(np.quantile(ref, 0.025)),
                "null_q975": float(np.quantile(ref, 0.975)),
                "p_pair_swap_upper": _upper_tail(observed, ref),
            })
        profile_records.extend(_binary_profile_tests(oceanic, continental, labels, scope))
        state_records.extend(_state_contrasts(oceanic, continental, profiles, space.ids, scope))
    # Preserve the entire geology-labelled eligible population as an
    # explicitly descriptive comparison. These are NOT the matched-pair
    # inferential cohorts, and no p-values are assigned to this table.
    full_records = []
    for scope, eligible_rows in (("all", g0_all), ("direct", g0_direct)):
        carriers = load_carriers(input_root / "floral_space_20261106" / scope /
                                 "island_profile_carriers.csv.gz")
        for region in ("all",) + REGIONS:
            members = {}
            for origin in ORIGINS:
                members[origin] = sorted(r["island_id"] for r in eligible_rows
                    if r["island_id"] in crosswalk
                    and crosswalk[r["island_id"]].get("geological_origin") == origin
                    and (region == "all" or r["analysis_regime"] == region))
            for profile, color, form in labels:
                a = members["continental"]
                b = members["oceanic"]
                fa = (sum(profile in carriers[i] for i in a) / len(a)) if a else None
                fb = (sum(profile in carriers[i] for i in b) / len(b)) if b else None
                full_records.append({
                    "scope": scope, "region": region, "profile_id": profile,
                    "flower_primary_color": color, "floral_form": form,
                    "n_continental_islands": len(a),
                    "n_oceanic_islands": len(b),
                    "continental_incidence_fraction": fa,
                    "oceanic_incidence_fraction": fb,
                    "difference_oceanic_minus_continental":
                        fb-fa if fa is not None and fb is not None else None,
                    "interpretation": "descriptive_unmatched_not_significance",
                })
    _write_csv(output / "all_eligible_profile_prevalence.csv", full_records,
               cols(full_records) if False else (
                   "scope", "region", "profile_id", "flower_primary_color",
                   "floral_form", "n_continental_islands", "n_oceanic_islands",
                   "continental_incidence_fraction", "oceanic_incidence_fraction",
                   "difference_oceanic_minus_continental", "interpretation",
               ))

    q_global = _bh([r["p_pair_swap_upper"] for r in group_records])
    for i, row in enumerate(group_records):
        row["q_bh_4"] = float(q_global[i])
        row["q_by_4"] = float(min(1, q_global[i] * (1 + 1/2 + 1/3 + 1/4)))
    cols = lambda rows: tuple(dict.fromkeys(k for row in rows for k in row))
    _write_csv(output / "global_repertoire_tests.csv", group_records, cols(group_records))
    _write_csv(output / "paired_rao_q.csv", q_records, cols(q_records))
    _write_csv(output / "profile_specific_results.csv", profile_records, cols(profile_records))
    _write_csv(output / "flower_state_results.csv", state_records, cols(state_records))
    manifest = {
        "status": "exploratory_matched_origin_repertoire",
        "upstream_island_commit": "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c",
        "matching": audit,
        "n_permutations": permutations,
        "seed": seed,
        "exact_profiles_test_family": PROFILE_FAMILY,
        "global_multivariate_family": 4,
        "source_columns": "region, area, distance, climate PCs 1-4, complete name count, trait coverage",
        "not_conditioned_on": ["genus", "lineage", "true source-region colonist pool",
                               "native-only flora", "spatially verified geological identity"],
        "analysis_unit": "matched island pair, not independent species",
        "method": "paired sign swap on richness-normalized recorded profile frequency differences",
        "claim_limits": "matching is observational with imperfect climate balance; category sets are recorded descriptions, not demonstrated realized pollinator service",
        "source_sha256": {name: file_sha256(path) for name, path in source_files.items()},
        "n_profile_bh_supported": {scope: sum(1 for x in profile_records
                                              if x["scope"] == scope and x["q_bh_496"] < 0.05)
                                   for scope in ("all", "direct")},
        "n_global_bh_supported": sum(x["q_bh_4"] < 0.05 for x in group_records),
        "n_global_by_supported": sum(x["q_by_4"] < 0.05 for x in group_records),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": manifest["status"], "n_pairs": len(matched),
        "region_matches": [(x["region"], x["n_matched"]) for x in audit["regional_diagnostics"]],
        "group_results": group_records,
        "q": q_records,
        "profile_bh_supported": manifest["n_profile_bh_supported"],
        "states_bh": {s: sum(r["q_bh_states"] < 0.05 for r in state_records
                             if r["scope"] == s) for s in ("all", "direct")},
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("input", "g0", "crosswalk", "output"):
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=PERMUTATIONS)
    args = parser.parse_args()
    run(args.input, args.g0, args.crosswalk, args.output,
        permutations=args.permutations)


if __name__ == "__main__":
    main()
