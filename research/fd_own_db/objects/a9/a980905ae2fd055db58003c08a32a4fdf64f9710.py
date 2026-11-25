"""Distance-native floral packing diagnostics against fixed-richness support nulls.

This deliberately measures RECORDED state-set phenotypes, not latent functional
hypervolume, biological absence, selection, fitness or pollination services.
N0 draws from all 496 publicly recorded reference profiles; N1 draws from
other eligible islands in the same broad analysis regime (leave-one-island-out).
The latter is an observational comparison, NOT a native-mainland source or
genus-conditioned null. All emitted outputs are aggregate-only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .pipeline import REF_INPUT, _write_csv, file_sha256
from .repertoire import match_origin_islands
from .geography_independent import _build_edges, independent_pair_matching
from .spatial_clusters import connected_block_components
from .g1_fixed_margins import _bh

SCOPES = ("all", "direct")
ENDPOINTS = ("mean_pairwise_distance", "mean_nearest_profile_distance")
NULLS = ("N0_global496", "N1_regional_observed_LOO",
         "N1b_regional_islandfreq_LOO")
PAIRED_NULLS = NULLS[:2]  # frequency-weighted N1b is descriptive falsification only
MODES = ("full", "equal30")
COHORTS = ("original_block_components", "geography_both_verified")
SUMMARY_FIELDS = (
    "scope", "null", "mode", "endpoint", "region", "geological_origin",
    "n_islands", "status", "median_n_profiles", "median_observed",
    "median_null_mean", "median_observed_minus_null", "median_SES",
    "mean_SES", "fraction_SES_below_zero", "median_null_sd",
    "median_thinning_q025", "median_thinning_q975",
)
PAIR_FIELDS = (
    "cohort", "scope", "null", "mode", "endpoint", "status",
    "n_available_pairs", "n_independent_sign_groups",
    "observed_mean_SES_oceanic_minus_continental",
    "median_pair_SES_oceanic_minus_continental",
    "exact_two_sided_sign_p", "sign_null_q025", "sign_null_q975",
    "n_sign_configurations", "q_bh_posthoc32", "q_by_posthoc32",
)


def _seed(base: int, *labels: object) -> int:
    payload = ":".join(map(str, (base, *labels))).encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def distance_statistics(distance: np.ndarray, selected: np.ndarray) -> np.ndarray:
    """Mean distinct-pair distance and mean within-support nearest-neighbor D.

    Distinct profiles are equally weighted, diagonals are excluded. Both
    endpoints are in original categorical set-Jaccard D, not PCoA.
    """
    ix = np.asarray(selected, dtype=int)
    if ix.ndim != 1 or len(ix) < 2 or len(np.unique(ix)) != len(ix):
        raise ValueError("At least two DISTINCT reference profile indices required")
    if np.any(ix < 0) or np.any(ix >= len(distance)):
        raise ValueError("Profile index outside reference")
    d = distance[np.ix_(ix, ix)].copy()
    pair = float(d.sum() / (len(ix) * (len(ix) - 1)))
    np.fill_diagonal(d, np.inf)
    near = float(np.min(d, axis=1).mean())
    return np.asarray([pair, near])


def sample_null(distance: np.ndarray, pool: np.ndarray, n: int, *,
                draws: int, seed: int,
                weights: np.ndarray | None = None) -> np.ndarray | None:
    """Draw n distinct profiles without replacement, optional fixed incidence weights."""
    available = np.asarray(pool, dtype=int)
    if n < 2 or draws < 20 or len(np.unique(available)) != len(available):
        raise ValueError("Invalid reference pool, sample size or null draws")
    if len(available) < n:
        return None  # do not replace unavailable states or fabricate absences
    probability = None
    if weights is not None:
        w = np.asarray(weights, dtype=float)
        if w.shape != available.shape or not np.isfinite(w).all() or np.any(w <= 0):
            raise ValueError("Incidence weights must be positive and pool-aligned")
        probability = w / w.sum()
    rng = np.random.default_rng(seed)
    out = np.empty((draws, len(ENDPOINTS)), dtype=float)
    for j in range(draws):
        take = rng.choice(available, size=n, replace=False, p=probability)
        out[j] = distance_statistics(distance, take)
    return out


def regional_profile_counts(rows: list[dict],
                            carriers: dict[str, dict[int, int]]) -> dict[str, Counter]:
    """Counts islands, not taxa, bearing each profile in a broad regime."""
    counts: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        island = r["island_id"]
        if island not in carriers:
            raise ValueError("G0 eligible island absent from public carrier input")
        counts[r["analysis_regime"]].update(carriers[island].keys())
    return dict(counts)


def regional_leave_one_out_pool(regional: dict[str, Counter], regime: str,
                                focal_profiles: set[int],
                                lookup: dict[int, int]) -> np.ndarray:
    """Exclude the focal island from construction of the reference support."""
    freq = regional.get(regime, Counter())
    ids = sorted(p for p, n in freq.items()
                 if n > int(p in focal_profiles))
    return np.asarray([lookup[p] for p in ids], dtype=int)


def _score(observed: np.ndarray, null: np.ndarray,
           *, thinning: np.ndarray | None = None) -> list[dict | None]:
    result = []
    for j, endpoint in enumerate(ENDPOINTS):
        mu = float(null[:, j].mean())
        sd = float(null[:, j].std(ddof=1))
        if not np.isfinite(sd) or sd <= 1e-12:
            result.append(None)
            continue
        value = float(observed[j])
        if thinning is None:
            lo, hi = value, value
        else:
            lo, hi = map(float, np.quantile(thinning[:, j], [.025, .975]))
        result.append({
            "endpoint": endpoint, "observed": value,
            "null_mean": mu, "null_sd": sd,
            "delta": value - mu, "ses": (value - mu) / sd,
            "thinning_q025": (lo - mu) / sd,
            "thinning_q975": (hi - mu) / sd,
        })
    return result


def exact_sign_flip(differences: np.ndarray,
                    components: list[list[int]]) -> dict:
    """Two-sided mean difference using all geographically grouped sign flips."""
    delta = np.asarray(differences, dtype=float)
    if delta.ndim != 1 or len(delta) < 8:
        raise ValueError("Insufficient matched differences")
    if len(components) > 18 or sorted(i for c in components for i in c) != list(range(len(delta))):
        raise ValueError("Invalid geographic component partition")
    observed = float(delta.mean())
    null = []
    for signs in itertools.product((-1.0, 1.0), repeat=len(components)):
        flipped = delta.copy()
        for component, sign in zip(components, signs):
            flipped[component] *= sign
        null.append(float(flipped.mean()))
    null = np.asarray(null)
    tol = 128 * np.finfo(float).eps * max(1.0, abs(observed))
    return {
        "observed_mean_SES_oceanic_minus_continental": observed,
        "median_pair_SES_oceanic_minus_continental": float(np.median(delta)),
        "exact_two_sided_sign_p": float(np.mean(np.abs(null) >= abs(observed) - tol)),
        "sign_null_q025": float(np.quantile(null, .025)),
        "sign_null_q975": float(np.quantile(null, .975)),
        "n_sign_configurations": len(null),
    }


def _summary(scores: list[dict], *, min_release_n: int = 5) -> list[dict]:
    groups: dict[tuple, list[dict]] = defaultdict(list)
    for r in scores:
        for regime in ("all", r["region"]):
            groups[tuple(r[k] for k in ("scope", "null", "mode", "endpoint")) +
                   (regime, r["geological_origin"])].append(r)
    summaries = []
    for key, values in sorted(groups.items()):
        row = dict(zip(("scope", "null", "mode", "endpoint", "region",
                        "geological_origin"), key))
        row["n_islands"] = len(values)
        row["status"] = "observational_reference_only" if len(values) >= min_release_n else "small_cell_suppressed"
        if len(values) >= min_release_n:
            for source, dest, fn in (
                ("n_profiles", "median_n_profiles", np.median),
                ("observed", "median_observed", np.median),
                ("null_mean", "median_null_mean", np.median),
                ("delta", "median_observed_minus_null", np.median),
                ("ses", "median_SES", np.median),
                ("ses", "mean_SES", np.mean),
                ("null_sd", "median_null_sd", np.median),
                ("thinning_q025", "median_thinning_q025", np.median),
                ("thinning_q975", "median_thinning_q975", np.median),
            ):
                row[dest] = float(fn([v[source] for v in values]))
            row["fraction_SES_below_zero"] = float(
                np.mean([v["ses"] < 0 for v in values]))
        summaries.append(row)
    return summaries


def _paired_results(scores: list[dict], original: list[dict],
                    independent: list[dict], geography: dict[str, dict]
                    ) -> list[dict]:
    lookup = {
        (r["scope"], r["null"], r["mode"], r["endpoint"], r["island_id"]): r["ses"]
        for r in scores
    }
    rows = []
    for cohort, pairs in ((COHORTS[0], original), (COHORTS[1], independent)):
        for scope in SCOPES:
            for null in PAIRED_NULLS:
                for mode in MODES:
                    for endpoint in ENDPOINTS:
                        key = (scope, null, mode, endpoint)
                        kept = [
                            pair for pair in pairs
                            if key + (pair["continental_id"],) in lookup
                            and key + (pair["oceanic_id"],) in lookup
                        ]
                        row = dict(zip(("cohort", "scope", "null", "mode", "endpoint"),
                                       (cohort, *key)))
                        row["n_available_pairs"] = len(kept)
                        if len(kept) < 8:
                            row["status"] = "insufficient_complete_pairs"
                            rows.append(row)
                            continue
                        diff = np.asarray([
                            lookup[key + (p["oceanic_id"],)] -
                            lookup[key + (p["continental_id"],)]
                            for p in kept
                        ])
                        groups = (connected_block_components(kept, geography)
                                  if cohort == COHORTS[0]
                                  else [[i] for i in range(len(kept))])
                        row.update(exact_sign_flip(diff, groups))
                        row["n_independent_sign_groups"] = len(groups)
                        row["status"] = ("conservative_spatial_component_stress_test"
                                         if cohort == COHORTS[0] else
                                         "posthoc_geographic_disjoint_pair_sensitivity")
                        rows.append(row)
    if len(rows) != 32:
        raise AssertionError("Predeclared family is 2 cohorts x 2 scopes x 2 nulls x 2 modes x 2 endpoints")
    # Missing slots remain p=1 in the 32-test exploratory family.
    q = _bh([r.get("exact_two_sided_sign_p", 1.0) for r in rows])
    harmonic = sum(1 / j for j in range(1, len(rows) + 1))
    for row, adj in zip(rows, q):
        if "exact_two_sided_sign_p" in row:
            row["q_bh_posthoc32"] = float(adj)
            row["q_by_posthoc32"] = float(min(1., adj * harmonic))
    return rows


def run(input_root: Path, g0_root: Path, crosswalk: Path, output: Path,
        *, draws: int = 999, thinnings: int = 99, seed: int = 20261109
        ) -> dict:
    if draws < 20 or thinnings < 20:
        raise ValueError("Require >=20 independent reference and thinning draws")
    input_root, g0_root, crosswalk, output = map(
        Path, (input_root, g0_root, crosswalk, output))
    output.mkdir(parents=True, exist_ok=True)
    space = ProfileSpace(load_profiles(input_root / REF_INPUT))
    distance = space.matrix
    if len(space.ids) != 496:
        raise ValueError("Frozen N0 global 496-profile universe changed")
    geo = {r["island_id"]: r for r in csv_rows(crosswalk)}
    g0 = {scope: list(csv_rows(g0_root / f"{scope}_islands.csv"))
          for scope in SCOPES}
    original, _ = match_origin_islands(g0["all"], g0["direct"], geo)
    edges, _ = _build_edges(g0["all"], g0["direct"], geo)
    strict, match_audit = independent_pair_matching(edges, "both_verified")
    global_reference = np.arange(len(space.ids), dtype=int)
    global_cache: dict[int, np.ndarray] = {}
    rows = []
    carrier_sources = {}
    for scope in SCOPES:
        path = input_root / "floral_space_20261106" / scope / "island_profile_carriers.csv.gz"
        carrier_sources[scope] = path
        carriers = load_carriers(path)
        regional = regional_profile_counts(g0[scope], carriers)
        for row in g0[scope]:
            island = row["island_id"]
            origin = geo.get(island, {}).get("geological_origin")
            if origin not in ("continental", "oceanic"):
                continue
            profiles = set(carriers[island])
            selected = np.asarray([space.lookup[p] for p in sorted(profiles)], dtype=int)
            n = len(selected)
            if n < 2:
                continue
            thin = None
            if n >= 30:
                rng = np.random.default_rng(_seed(seed, scope, island, "observed_thinning"))
                thin = np.asarray([
                    distance_statistics(distance, rng.choice(selected, 30, replace=False))
                    for _ in range(thinnings)
                ])
            regional_pool = regional_leave_one_out_pool(
                regional, row["analysis_regime"], profiles, space.lookup)
            # Each occurrence is ONE other eligible island's recorded support,
            # not a carrier-name count or inferred biological abundance.
            regional_weights = np.asarray([
                regional[row["analysis_regime"]][space.ids[int(idx)]]
                - int(space.ids[int(idx)] in profiles)
                for idx in regional_pool
            ], dtype=float)
            if np.any(regional_weights <= 0):
                raise AssertionError("Leave-one-out regional incidence must be positive")
            for mode in MODES:
                if mode == "equal30" and thin is None:
                    continue
                n_draw = n if mode == "full" else 30
                observed_sample = (
                    distance_statistics(distance, selected)
                    if mode == "full" else thin.mean(axis=0)
                )
                for null in NULLS:
                    if null == "N0_global496":
                        if n_draw not in global_cache:
                            global_cache[n_draw] = sample_null(
                                distance, global_reference, n_draw,
                                draws=draws, seed=_seed(seed, "N0", n_draw))
                        simulated = global_cache[n_draw]
                    elif null == "N1_regional_observed_LOO":
                        simulated = sample_null(
                            distance, regional_pool, n_draw, draws=draws,
                            seed=_seed(seed, "N1", scope, island, mode))
                    else:
                        simulated = sample_null(
                            distance, regional_pool, n_draw, draws=draws,
                            seed=_seed(seed, "N1b", scope, island, mode),
                            weights=regional_weights)
                    if simulated is None:
                        continue
                    for value in _score(observed_sample, simulated,
                                        thinning=thin if mode == "equal30" else None):
                        if value is not None:
                            rows.append({
                                "scope": scope, "null": null, "mode": mode,
                                "region": row["analysis_regime"],
                                "geological_origin": origin,
                                "island_id": island, "n_profiles": n,
                                **value,
                            })
    summaries = _summary(rows)
    paired = _paired_results(rows, original, strict, geo)
    _write_csv(output / "distance_ses_origin_summary.csv", summaries, SUMMARY_FIELDS)
    _write_csv(output / "distance_ses_paired_spatial_sensitivity.csv", paired, PAIR_FIELDS)
    manifest = {
        "schema": "fd_original_distance_fixed_support_SES_v1",
        "status": "observational_posthoc_distance_packing_sensitivity",
        "reference": "original colour/form state-set Jaccard D, no projection",
        "global_null": "uniform n distinct profiles from all 496 frozen reference states",
        "regional_null": "same-regime observed profile union, excluding entire focal island; neither mainland donor nor genus conditional",
        "regional_frequency_sensitivity": "N1b selects distinct profiles without replacement in proportion to number of OTHER eligible regional islands bearing each profile; descriptive only, excluded from 32 geology tests",
        "full": "sample size = number of distinct observed island profiles",
        "equal30": f"mean over {thinnings} repeated observed 30-distinct-profile subsamples, compared with fixed n=30 reference",
        "reference_draws": draws, "seed": seed,
        "endpoints": list(ENDPOINTS), "nested_scopes": list(SCOPES),
        "origin_pair_cohorts": {
            "original_matched_pairs": len(original),
            "original_spatial_block_components": len(connected_block_components(original, geo)),
            "geography_verified_disjoint_pairs": len(strict),
            "strict_unknown_source_archipelago_endpoints": match_audit["n_unknown_archipelago_endpoints"],
        },
        "paired_tests": "32 fixed exploratory slots for N0 and N1 only; N1b is descriptive falsification; two-sided exact component sign flips or geographic-disjoint pair flips; BH/BY across all slots with missing set to p=1",
        "release": "aggregate summaries only, suppress summaries n<5; no island, taxon, paired identity, genus or source names",
        "interpretation_limit": "SES conditions on recorded support only; no ecological absence, evolutionary cause, source/genus null or phylogenetic independence inferred",
        "source_sha256": {
            "profile_reference": file_sha256(input_root / REF_INPUT),
            "all_g0": file_sha256(g0_root / "all_islands.csv"),
            "direct_g0": file_sha256(g0_root / "direct_islands.csv"),
            "geology_crosswalk": file_sha256(crosswalk),
            **{f"{scope}_carrier_aggregate": file_sha256(p)
               for scope, p in carrier_sources.items()},
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({
        "schema": manifest["schema"], "n_island_endpoint_ses_internal": len(rows),
        "n_public_summary_rows": len(summaries),
        "n_public_pair_test_slots": len(paired),
        "n_verified_disjoint_pairs": len(strict),
        "min_pair_p": min((r["exact_two_sided_sign_p"] for r in paired
                          if "exact_two_sided_sign_p" in r), default=None),
        "note": "Interpret only after spatial multiplicity and source/lineage controls",
    }))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("input", "g0", "crosswalk", "output"):
        parser.add_argument("--" + arg, type=Path, required=True)
    parser.add_argument("--draws", type=int, default=999)
    parser.add_argument("--thinnings", type=int, default=99)
    args = parser.parse_args()
    run(args.input, args.g0, args.crosswalk, args.output,
        draws=args.draws, thinnings=args.thinnings)


if __name__ == "__main__":
    main()
