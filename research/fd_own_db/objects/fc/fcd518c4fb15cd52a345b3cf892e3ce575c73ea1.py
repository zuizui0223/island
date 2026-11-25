"""G1 fixed-margins exploratory null for distance-based floral profiles.

This deliberately tests equal-PROFILE-weight Rao Q and nearest-profile
distance.  Binary island/profile Curveball trades preserve observed island
profile richness and across-island profile prevalence.  Species carrier
abundances, shared-taxon removal, source/genus assignments and native status
are NOT preserved and cannot be inferred from these aggregate inputs.

The selected island and pair cohorts are chosen using only IDs from the
already frozen PR321 matched pairs, before looking at phenotype distances.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .pipeline import REF_INPUT, SCOPES, _write_csv, file_sha256

REGIMES = (
    "northern_midlatitude", "northern_high_latitude",
    "tropical", "southern_extratropical",
)
METRICS = ("rao_q_equal_profile", "pair_nn_equal_profile")
FIELDS = (
    "metric", "scope", "analysis_regime", "encoding", "n_islands_pool",
    "n_islands_scored", "n_pairs_scored", "observed", "null_median",
    "difference", "null_q025", "null_q975", "p_two_sided_approx",
    "q_bh_32", "q_by_32", "rhat", "max_abs_lag1", "n_distinct",
    "mixing_pass", "direction",
)
DRAW_FIELDS = ("metric", "scope", "analysis_regime", "chain", "draw", "statistic")


def _stable_seed(seed: int, *labels: str) -> int:
    data = ":".join([str(seed), *labels]).encode("utf-8")
    return int.from_bytes(hashlib.sha256(data).digest()[:8], "little")


def _bh(values: list[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    n = len(p)
    order = np.argsort(p, kind="stable")
    ranks = np.arange(1, n + 1)
    adjusted = np.minimum.accumulate((p[order] * n / ranks)[::-1])[::-1]
    result = np.empty(n, dtype=float)
    result[order] = np.minimum(1.0, adjusted)
    return result


def _curveball_once(rows: list[set[int]], rng: np.random.Generator) -> None:
    n = len(rows)
    if n < 2:
        return
    i, j = rng.choice(n, size=2, replace=False)
    ai, bj = rows[int(i)], rows[int(j)]
    only_a, only_b = ai - bj, bj - ai
    if not only_a or not only_b:
        return
    exclusive = sorted(only_a | only_b)
    rng.shuffle(exclusive)
    shared = ai & bj
    split = len(only_a)
    rows[int(i)] = shared | set(exclusive[:split])
    rows[int(j)] = shared | set(exclusive[split:])


def _margins(rows: list[set[int]]) -> tuple[list[int], dict[int, int]]:
    sizes = [len(r) for r in rows]
    sums: dict[int, int] = {}
    for r in rows:
        for p in r:
            sums[p] = sums.get(p, 0) + 1
    return sizes, sums


def _local_q(rows: list[set[int]], selected: list[int], matrix: np.ndarray) -> float:
    values = []
    for i in selected:
        ids = np.fromiter(sorted(rows[i]), dtype=np.intp)
        d = matrix[np.ix_(ids, ids)]
        values.append(float(d.mean()))
    return float(np.median(values))


def _pair_nn(rows: list[set[int]], pairs: list[tuple[int, int]], matrix: np.ndarray) -> float:
    values = []
    for a, b in pairs:
        ia = np.fromiter(sorted(rows[a]), dtype=np.intp)
        ib = np.fromiter(sorted(rows[b]), dtype=np.intp)
        d = matrix[np.ix_(ia, ib)]
        values.append(float((d.min(axis=1).mean() + d.min(axis=0).mean()) / 2.0))
    return float(np.median(values))


def _mixing(chains: list[np.ndarray]) -> tuple[float, float, int]:
    if len(chains) != 2 or len(chains[0]) != len(chains[1]):
        raise ValueError("Expected two equally long chains")
    samples = np.vstack(chains)
    n = samples.shape[1]
    within = samples.var(axis=1, ddof=1)
    w = float(within.mean())
    if w <= 1e-24:
        rhat = 1.0 if np.ptp(samples) <= 1e-12 else float("inf")
    else:
        between = n * float(samples.mean(axis=1).var(ddof=1))
        rhat = float(np.sqrt(((n - 1) / n * w + between / n) / w))
    acfs = []
    for chain in samples:
        if chain.var() <= 1e-24:
            acfs.append(0.0)
        else:
            acfs.append(float(np.corrcoef(chain[:-1], chain[1:])[0, 1]))
    distinct = int(len(np.unique(np.round(samples.reshape(-1), decimals=12))))
    return rhat, max(abs(v) if np.isfinite(v) else 1.0 for v in acfs), distinct


def _p_value(obs: float, values: np.ndarray) -> float:
    tolerance = 100.0 * np.finfo(float).eps * max(1.0, abs(obs))
    lower = (1 + int((values <= obs + tolerance).sum())) / (len(values) + 1)
    upper = (1 + int((values >= obs - tolerance).sum())) / (len(values) + 1)
    return min(1.0, 2.0 * min(lower, upper))


def _read_pairs(g0: Path) -> dict[tuple[str, str], list[tuple[str, str]]]:
    cohort: dict[tuple[str, str], list[tuple[str, str]]] = {}
    for scope in SCOPES:
        for row in csv_rows(g0 / f"{scope}_pairs.csv"):
            regime = row["analysis_regime"]
            if regime not in REGIMES:
                raise ValueError(f"Unknown analysis regime: {regime}")
            pair = tuple(sorted((row["island_a"], row["island_b"])))
            if pair[0] == pair[1]:
                raise ValueError("Self pairs are disallowed")
            cohort.setdefault((scope, regime), []).append(pair)
    # Upstream PR321 matched pairs are identical for All/Direct in each encoding.
    for one, two in (("all", "direct"), ("all_single_state", "direct_single_state")):
        for regime in REGIMES:
            if sorted(cohort.get((one, regime), [])) != sorted(cohort.get((two, regime), [])):
                raise ValueError(f"Matched All/Direct cohorts do not agree: {regime}/{one}/{two}")
    for key, pairs in cohort.items():
        if len(pairs) != len(set(pairs)):
            raise ValueError(f"Duplicate pair in {key}")
    return cohort


def _draw_cohort(pairs: list[tuple[str, str]], *, seed: int, local_n: int,
                 pair_n: int) -> tuple[list[str], list[tuple[str, str]], list[str]]:
    pool = sorted({x for pair in pairs for x in pair})
    if not pool or not pairs:
        return pool, [], []
    pair_gen = np.random.default_rng(_stable_seed(seed, "pairs"))
    island_gen = np.random.default_rng(_stable_seed(seed, "islands"))
    n_pairs = min(pair_n, len(pairs))
    ordered = sorted(pairs)
    sample_pairs = [ordered[k] for k in sorted(pair_gen.choice(len(ordered), size=n_pairs, replace=False))]
    n_local = min(local_n, len(pool))
    sample_islands = [pool[k] for k in sorted(island_gen.choice(len(pool), size=n_local, replace=False))]
    return pool, sample_pairs, sample_islands


def run_g1(input_root: Path, g0_root: Path, output_root: Path, *,
           local_n: int = 128, pair_n: int = 128, draws_per_chain: int = 200,
           burn_trades_per_island: int = 12, thin_trades_per_island: int = 3,
           seed: int = 20261108) -> dict:
    for name, val in (("local_n", local_n), ("pair_n", pair_n),
                      ("draws_per_chain", draws_per_chain),
                      ("burn_trades_per_island", burn_trades_per_island),
                      ("thin_trades_per_island", thin_trades_per_island)):
        if isinstance(val, bool) or not isinstance(val, int) or val < 1:
            raise ValueError(f"{name} must be positive integer")
    if draws_per_chain < 8:
        raise ValueError("At least eight MCMC draws per chain required")
    input_root, g0_root, output_root = map(Path, (input_root, g0_root, output_root))
    output_root.mkdir(parents=True, exist_ok=True)
    ref = input_root / REF_INPUT
    space = ProfileSpace(load_profiles(ref))
    cohorts = _read_pairs(g0_root)
    summaries, draws, cohorts_audit = [], [], []
    carrier_paths = {s: input_root / "floral_space_20261106" / s /
                     "island_profile_carriers.csv.gz" for s in SCOPES}
    carriers = {scope: load_carriers(path) for scope, path in carrier_paths.items()}
    for scope in SCOPES:
        encoding = "single_state" if "single_state" in scope else "full_state_sets"
        for regime in REGIMES:
            original_pairs = cohorts.get((scope, regime), [])
            matched_seed = _stable_seed(seed, encoding, regime)
            pool, sampled_pairs, local_islands = _draw_cohort(
                original_pairs, seed=matched_seed, local_n=local_n, pair_n=pair_n,
            )
            if not pool:
                for metric in METRICS:
                    summaries.append({"metric": metric, "scope": scope,
                                      "analysis_regime": regime, "encoding": encoding,
                                      "n_islands_pool": 0, "n_islands_scored": 0,
                                      "n_pairs_scored": 0, "observed": float("nan"),
                                      "null_median": float("nan"), "difference": float("nan"),
                                      "null_q025": float("nan"), "null_q975": float("nan"),
                                      "p_two_sided_approx": float("nan"), "rhat": float("nan"),
                                      "max_abs_lag1": float("nan"), "n_distinct": 0,
                                      "mixing_pass": False, "direction": "undefined"})
                continue
            rows = []
            for island in pool:
                if island not in carriers[scope]:
                    raise ValueError(f"{scope}: Missing profile incidence for paired island {island}")
                rows.append({space.lookup[p] for p in carriers[scope][island]})
            row_lookup = {island: i for i, island in enumerate(pool)}
            local_ix = [row_lookup[i] for i in local_islands]
            pair_ix = [(row_lookup[a], row_lookup[b]) for a, b in sampled_pairs]
            margins = _margins(rows)
            observed = {
                METRICS[0]: _local_q(rows, local_ix, space.matrix),
                METRICS[1]: _pair_nn(rows, pair_ix, space.matrix),
            }
            chains: dict[str, list[np.ndarray]] = {metric: [] for metric in METRICS}
            for chain in range(2):
                rng = np.random.default_rng(_stable_seed(seed, scope, regime, "chain", str(chain)))
                current = [s.copy() for s in rows]
                for _ in range(burn_trades_per_island * len(current)):
                    _curveball_once(current, rng)
                trace = {metric: np.zeros(draws_per_chain, dtype=float) for metric in METRICS}
                for draw in range(draws_per_chain):
                    for _ in range(thin_trades_per_island * len(current)):
                        _curveball_once(current, rng)
                    trace[METRICS[0]][draw] = _local_q(current, local_ix, space.matrix)
                    trace[METRICS[1]][draw] = _pair_nn(current, pair_ix, space.matrix)
                if _margins(current) != margins:
                    raise AssertionError(f"Curveball did not preserve exact row/column margins: {scope}/{regime}")
                for metric in METRICS:
                    chains[metric].append(trace[metric])
                    draws.extend({"metric": metric, "scope": scope,
                                  "analysis_regime": regime, "chain": chain,
                                  "draw": t, "statistic": float(v)}
                                 for t, v in enumerate(trace[metric]))
            cohorts_audit.append({"scope": scope, "analysis_regime": regime,
                                  "encoding": encoding, "n_islands_pool": len(pool),
                                  "n_islands_scored": len(local_ix),
                                  "n_pairs_original": len(original_pairs),
                                  "n_pairs_scored": len(pair_ix),
                                  "local_selection_hash": hashlib.sha256(
                                      "|".join(local_islands).encode()).hexdigest(),
                                  "pair_selection_hash": hashlib.sha256(
                                      "|".join(a + ":" + b for a,b in sampled_pairs).encode()).hexdigest()})
            for metric in METRICS:
                sampled = np.concatenate(chains[metric])
                rhat, acf, distinct = _mixing(chains[metric])
                mixing_pass = bool(np.isfinite(rhat) and rhat <= 1.05
                                   and acf <= 0.30 and distinct >= 8)
                expected = float(np.median(sampled))
                summaries.append({
                    "metric": metric, "scope": scope, "analysis_regime": regime,
                    "encoding": encoding, "n_islands_pool": len(pool),
                    "n_islands_scored": len(local_ix), "n_pairs_scored": len(pair_ix),
                    "observed": observed[metric], "null_median": expected,
                    "difference": observed[metric] - expected,
                    "null_q025": float(np.quantile(sampled, 0.025)),
                    "null_q975": float(np.quantile(sampled, 0.975)),
                    "p_two_sided_approx": _p_value(observed[metric], sampled)
                    if mixing_pass else float("nan"),
                    "rhat": rhat, "max_abs_lag1": acf, "n_distinct": distinct,
                    "mixing_pass": mixing_pass,
                    "direction": "greater" if observed[metric] > expected
                    else "lower" if observed[metric] < expected else "equal",
                })
    # Two endpoints x four scopes x four regions: a single declared 32-family.
    if len(summaries) != 32:
        raise AssertionError(f"Expected 32 test-family slots: {len(summaries)}")
    probs = [float(x["p_two_sided_approx"]) if np.isfinite(x["p_two_sided_approx"])
             else 1.0 for x in summaries]
    bh = _bh(probs)
    harmonic = np.sum(1.0 / np.arange(1, len(probs) + 1))
    for i, row in enumerate(summaries):
        row["q_bh_32"] = float(bh[i]) if row["mixing_pass"] else float("nan")
        row["q_by_32"] = float(min(1.0, bh[i] * harmonic)) if row["mixing_pass"] else float("nan")
    _write_csv(output_root / "g1_summary.csv", summaries, FIELDS)
    _write_csv(output_root / "g1_draws.csv", draws, DRAW_FIELDS)
    _write_csv(output_root / "cohort_audit.csv", cohorts_audit,
               ("scope","analysis_regime","encoding","n_islands_pool","n_islands_scored",
                "n_pairs_original","n_pairs_scored","local_selection_hash","pair_selection_hash"))
    manifest = {
        "schema": "fd_g1_fixed_binary_profile_margin_v1",
        "status": "exploratory_mcmc_approximate",
        "upstream_commit": "d4a5dafab5f0737ea6c1c1c0263fc20decf6035c",
        "parameters": {"local_n": local_n, "pair_n": pair_n,
                       "draws_per_chain": draws_per_chain, "n_chains": 2,
                       "burn_trades_per_island": burn_trades_per_island,
                       "thin_trades_per_island": thin_trades_per_island, "seed": seed},
        "conditioned_on": ["island profile count (row sum)",
                          "number of sampled islands per profile (column sum)",
                          "fixed PR321 matched pairs and preregistered ID-only subsamples"],
        "not_conditioned_on": ["carrier taxon counts", "species identity",
                               "genus/source pool", "native status", "traits availability process"],
        "endpoints": list(METRICS),
        "weights": "equal for unique recorded profiles; not carrier-weighted",
        "inference": "MCMC approximate 2-sided inclusive-tail p; BH/BY over fixed 32 slots",
        "caution": "MCMC draws and overlapping pairs/scopes are not independent replications",
        "observed_test_slots": len(summaries),
        "mixing_pass_slots": sum(bool(x["mixing_pass"]) for x in summaries),
        "source_sha256": {"reference": file_sha256(ref),
                          "g0_pair_csv": {scope: file_sha256(g0_root / f"{scope}_pairs.csv")
                                          for scope in SCOPES},
                          "carrier": {scope: file_sha256(path) for scope,path in carrier_paths.items()}},
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"status": manifest["status"],
                      "mixing_pass": manifest["mixing_pass_slots"],
                      "min_bh": float(min((x["q_bh_32"] for x in summaries
                                             if np.isfinite(x["q_bh_32"])), default=float("nan"))),
                      "n_positive_bh": sum(1 for x in summaries if x["mixing_pass"]
                                           and x["q_bh_32"] < 0.05)},
                     ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--g0", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=200, help="draws PER chain")
    parser.add_argument("--local-n", type=int, default=128)
    parser.add_argument("--pair-n", type=int, default=128)
    args = parser.parse_args()
    run_g1(args.input, args.g0, args.output, local_n=args.local_n,
           pair_n=args.pair_n, draws_per_chain=args.draws)


if __name__ == "__main__":
    main()
