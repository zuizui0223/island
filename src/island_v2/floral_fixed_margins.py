"""Fixed row/column margin sensitivity for island floral-profile sharing.

A *complementary estimand* to the genus-conditioned, shared-taxon-excluded
comparison in PR #321. This uses the public island x profile incidence, fixes
both per-island profile richness and regional profile occupancy across islands,
and retains the same sampled pair IDs. It scores whole island floral profiles:
it DOES NOT remove shared taxa or fix genus composition, source pools, or range
size. Do not interpret a departure as floral convergence or selection.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter, defaultdict
from pathlib import Path
import random
import statistics

SCOPES = ("all", "direct", "all_single_state", "direct_single_state")
REGIMES = ("northern_midlatitude", "northern_high_latitude", "tropical", "southern_extratropical")
ENCODINGS = (("all", "direct"), ("all_single_state", "direct_single_state"))


def read_rows(path: Path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def write_csv(path: Path, fields: list[str], rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        with path.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as out:
                    writer = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
                    writer.writeheader()
                    writer.writerows(rows)
    else:
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


def jaccard(a: set[int], b: set[int]) -> float:
    union = len(a | b)
    if not union:
        raise ValueError("empty union: absent profiles cannot be treated as zeros")
    return len(a & b) / union


def curveball_trade(rows: list[set[int]], rng: random.Random) -> bool:
    """One symmetric Curveball trade, exactly preserving row and column sums."""
    if len(rows) < 2:
        return False
    i = rng.randrange(len(rows))
    j = rng.randrange(len(rows) - 1)
    if j >= i:
        j += 1
    a, b = rows[i], rows[j]
    own_a, own_b = a - b, b - a
    if not own_a or not own_b:
        return False
    pooled = sorted(own_a | own_b)
    chosen = set(rng.sample(pooled, len(own_a)))
    common = a & b
    rows[i] = common | chosen
    rows[j] = common | (set(pooled) - chosen)
    return True


def margins(rows: list[set[int]]):
    return tuple(map(len, rows)), Counter(code for row in rows for code in row)


def summary_stat(rows: list[set[int]], pairs: list[tuple[int, int]]) -> float:
    return statistics.median(jaccard(rows[a], rows[b]) for a, b in pairs)


def _percentile(v: list[float], pct: float) -> float:
    v = sorted(v)
    if not v:
        raise ValueError("No null draws")
    t = (len(v) - 1) * pct
    i = int(t)
    return v[i] if i + 1 == len(v) else v[i] + (v[i + 1] - v[i]) * (t - i)


def _rhat(chains: list[list[float]]) -> float | None:
    if len(chains) < 2 or min(map(len, chains)) < 3 or len(set(map(len, chains))) != 1:
        return None
    n = len(chains[0])
    variances = [statistics.variance(chain) for chain in chains]
    w = statistics.mean(variances)
    if w <= 1e-20:
        return None
    b = n * statistics.variance([statistics.mean(chain) for chain in chains])
    return (((n - 1) / n * w + b / n) / w) ** 0.5


def _lag1(chains: list[list[float]]) -> float | None:
    values = []
    for chain in chains:
        if len(chain) < 3:
            continue
        mean = statistics.mean(chain)
        d = sum((x - mean) ** 2 for x in chain)
        if d > 1e-20:
            values.append(sum((chain[i] - mean) * (chain[i - 1] - mean)
                              for i in range(1, len(chain))) / d)
    return max(values) if values else None


def run_comparison(rows: list[set[int]], pair_indices: list[tuple[int, int]], *,
                   seed: int, per_chain: int = 200, burn_multiplier: int = 12,
                   gap_multiplier: int = 3) -> tuple[dict, list[dict]]:
    if len(rows) < 2 or not pair_indices or per_chain < 3:
        raise ValueError("At least 2 islands, 1 pair, 3 draws per chain required")
    original_margins = margins(rows)
    if any(not row for row in rows):
        raise ValueError("Every sampled island needs at least one profile")
    for a, b in pair_indices:
        if a == b or a < 0 or b < 0 or a >= len(rows) or b >= len(rows):
            raise ValueError("Invalid pair indices")
    observed = summary_stat(rows, pair_indices)
    chains, output = [], []
    for chain_id in range(2):
        state = [set(s) for s in rows]
        rng = random.Random(seed + 10000019 * chain_id)
        for _ in range(burn_multiplier * len(state)):
            curveball_trade(state, rng)
        chain = []
        for step in range(per_chain):
            for _ in range(gap_multiplier * len(state)):
                curveball_trade(state, rng)
            score = summary_stat(state, pair_indices)
            chain.append(score)
            output.append({"chain": chain_id, "draw": step + 1, "median_jaccard": score})
        if margins(state) != original_margins:
            raise AssertionError("Curveball margin invariant violated")
        chains.append(chain)
    flat = [v for c in chains for v in c]
    rhat, lag1 = _rhat(chains), _lag1(chains)
    diagnostics_ok = (rhat is not None and rhat <= 1.05
                      and lag1 is not None and lag1 <= 0.30
                      and len(set(flat)) >= 8)
    # Autocorrelated MCMC samples are not iid random permutations. Tail values
    # below are EXPLORATORY approximate frequencies; undefined if mixing fails.
    eps = 1e-12
    p_lower = (1 + sum(x <= observed + eps for x in flat)) / (1 + len(flat))
    p_upper = (1 + sum(x >= observed - eps for x in flat)) / (1 + len(flat))
    p = min(1., 2 * min(p_lower, p_upper)) if diagnostics_ok else None
    return {
        "observed_median_jaccard": observed,
        "null_median": statistics.median(flat),
        "delta_jaccard": observed - statistics.median(flat),
        "null_q025": _percentile(flat, 0.025),
        "null_q975": _percentile(flat, 0.975),
        "mcmc_rhat": rhat,
        "mcmc_max_lag1": lag1,
        "null_distinct_values": len(set(flat)),
        "null_draws": len(flat),
        "p_two_sided_exploratory": p,
        "mixing_status": "pass" if diagnostics_ok else "diagnostic_fail_no_p_value",
    }, output


def _adjust(pvals: list[float | None], *, by: bool = False) -> list[float | None]:
    m = len(pvals)
    ordered = sorted((1 if p is None else p, i) for i, p in enumerate(pvals))
    out = [None] * m
    prev = 1.
    harmonic = sum(1 / j for j in range(1, m + 1)) if by else 1.
    for rank in range(m, 0, -1):
        p, i = ordered[rank - 1]
        prev = min(prev, p * m * harmonic / rank, 1.)
        if pvals[i] is not None:
            out[i] = prev
    return out


def _filehash(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def _incidence(root: Path, scope: str, needed: set[str]) -> dict[str, set[int]]:
    out: dict[str, set[int]] = defaultdict(set)
    for row in read_rows(root / scope / "island_profile_carriers.csv.gz"):
        island = row["island_id"]
        if island not in needed:
            continue
        if int(row["n_species"]) < 1:
            raise ValueError("Nonpositive profile carrier count")
        value = int(row["profile_id"])
        if value in out[island]:
            raise ValueError("Duplicate island/profile incidence")
        out[island].add(value)
    if needed - out.keys():
        raise ValueError(f"Missing profile incidence for {len(needed - out.keys())} islands")
    metrics = {}
    for row in read_rows(root / scope / "island_metrics.csv"):
        if row["island_id"] in needed:
            metrics[row["island_id"]] = int(row["n_profiles"])
    for island in needed:
        if metrics.get(island) != len(out[island]):
            raise ValueError("Profile incidence and saved island richness disagree")
    return out


def run_suite(previous: Path, root: Path, output: Path, *, per_chain: int = 200,
              seed: int = 20261107, burn_multiplier: int = 12,
              gap_multiplier: int = 3) -> dict:
    previous, root, output = map(Path, (previous, root, output))
    if per_chain < 3 or seed < 0 or burn_multiplier < 1 or gap_multiplier < 1:
        raise ValueError("Invalid simulation arguments")
    by_scope = defaultdict(list)
    for row in read_rows(previous / "pairwise_observed_and_null.csv.gz"):
        if row["evidence_scope"] not in SCOPES or row["analysis_regime"] not in REGIMES:
            raise ValueError("Unexpected evidence scope or regime")
        by_scope[row["evidence_scope"], row["analysis_regime"]].append(row)
    pair_sets = {}
    for a, b in ENCODINGS:
        for regime in REGIMES:
            pair_a = [(r["island_a"], r["island_b"]) for r in by_scope[a, regime]]
            pair_b = [(r["island_a"], r["island_b"]) for r in by_scope[b, regime]]
            if pair_a != pair_b or not pair_a:
                raise ValueError("The original All and Direct pair cohorts are not matched")
            if len(set(pair_a)) != len(pair_a):
                raise ValueError("Duplicate sampled pair")
            pair_sets[a, regime] = pair_sets[b, regime] = pair_a
    results, draws = [], []
    for scope_i, scope in enumerate(SCOPES):
        needed = set(r["island_a"] for regime in REGIMES for r in by_scope[scope, regime])
        needed.update(r["island_b"] for regime in REGIMES for r in by_scope[scope, regime])
        lookup = _incidence(root, scope, needed)
        for regime_i, regime in enumerate(REGIMES):
            pairs = pair_sets[scope, regime]
            ids = sorted(set(x for pair in pairs for x in pair))
            index = {island: k for k, island in enumerate(ids)}
            incidence = [lookup[name] for name in ids]
            code_pairs = [(index[a], index[b]) for a, b in pairs]
            result, drawrows = run_comparison(
                incidence, code_pairs, per_chain=per_chain,
                seed=seed + scope_i * 100003 + regime_i * 101,
                burn_multiplier=burn_multiplier, gap_multiplier=gap_multiplier)
            no_shared = [i for i, row in enumerate(by_scope[scope, regime])
                         if int(row["n_shared_taxa"]) == 0]
            # A separate no-shared-names diagnostic can be interpreted against
            # the same *whole-island* incidence null. No additional p-value.
            result.update({
                "evidence_scope": scope, "analysis_regime": regime,
                "n_pairs": len(pairs), "n_islands": len(ids),
                "n_profile_codes": len(set().union(*incidence)),
                "n_pairs_zero_shared_taxa": len(no_shared),
                "observed_zero_shared_median_jaccard": (
                    statistics.median(jaccard(incidence[code_pairs[k][0]], incidence[code_pairs[k][1]])
                                      for k in no_shared) if no_shared else None),
            })
            results.append(result)
            draws.extend({"evidence_scope": scope, "analysis_regime": regime, **r} for r in drawrows)
    if len(results) != 16:
        raise AssertionError("Expected 16 fixed comparisons")
    q_bh = _adjust([r["p_two_sided_exploratory"] for r in results])
    q_by = _adjust([r["p_two_sided_exploratory"] for r in results], by=True)
    for i, r in enumerate(results):
        r["q_bh_16_exploratory"] = q_bh[i]
        r["q_by_16_exploratory"] = q_by[i]
    output.mkdir(parents=True, exist_ok=True)
    summary_file = output / "permutation_summary.csv"
    draws_file = output / "permutation_draws.csv.gz"
    write_csv(summary_file, list(results[0]), results)
    write_csv(draws_file, list(draws[0]), draws)
    manifest = {
        "status": "complete", "design": "exploratory fixed row/column margins; separate from genus null",
        "endpoint": "whole-island recorded profile sets, NOT shared-name-exclusive",
        "conditional_margins": "profile richness per sampled island; profile prevalence across sampled islands in each regime/scope",
        "pair_ids": "identical to PR321 genus-conditioned pair cohort, All/Direct matched",
        "not_preserved": ["genus composition", "taxon co-occurrence", "shared-name-excluded profile sets",
                          "native status", "historic source pools", "sampling detection"],
        "mcmc": "two Curveball Markov chains; approximate p only if diagnostics pass; not exact independent randomization",
        "params": dict(per_chain=per_chain, seed=seed, burn_multiplier=burn_multiplier,
                       gap_multiplier=gap_multiplier),
        "inputs_sha256": {
            "prior_pair_table": _filehash(previous / "pairwise_observed_and_null.csv.gz"),
            **{scope + "_incidence": _filehash(root / scope / "island_profile_carriers.csv.gz") for scope in SCOPES},
            **{scope + "_metrics": _filehash(root / scope / "island_metrics.csv") for scope in SCOPES},
        },
        "outputs_sha256": {name: _filehash(output / name) for name in
                           (summary_file.name, draws_file.name)},
        "comparisons": len(results),
        "significant_bh_0_05_exploratory": sum(r["q_bh_16_exploratory"] is not None
                                               and r["q_bh_16_exploratory"] < 0.05 for r in results),
        "mixing_pass_count": sum(r["mixing_status"] == "pass" for r in results),
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = ["# Fixed-margin floral-profile sharing sensitivity", "",
             "This **exploratory complementary null** uses the same island pairs as PR #321 but retains all taxon",
             "names in the whole-island incidence. It holds each island's recorded number of floral profiles",
             "and each profile's number of sampled islands fixed, then makes symmetric Curveball trades.",
             "It **does not** condition on genus/source pool or remove shared taxa, so it is not a stronger",
             "version of the original conditional test. A significant result would not establish evolution or",
             "pollinator replacement. Markov-chain tail frequencies are approximate, not exact p-values.", "",
             "| Scope | Region | Pairs | Islands | Observed J | Null median | ΔJ | Approx. BH q | Mixing |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for r in results:
        f = lambda x: "NA" if x is None else f"{x:.4f}"
        lines.append(f"| {r['evidence_scope']} | {r['analysis_regime']} | {r['n_pairs']} | {r['n_islands']} | "
                     f"{f(r['observed_median_jaccard'])} | {f(r['null_median'])} | "
                     f"{f(r['delta_jaccard'])} | {f(r['q_bh_16_exploratory'])} | {r['mixing_status']} |")
    lines.extend(["", "The observed profile set of any island is reused across its pairs; pairs are not independent.",
                  "Rows with zero shared taxon names are enumerated only as a descriptive sensitivity,",
                  "because their numbers vary and the permuted incidence cannot preserve taxon identity.",
                  "Null 2.5–97.5% quantiles are randomization reference ranges, not effect confidence intervals.",
                  "See [summary](permutation_summary.csv), [draws](permutation_draws.csv.gz) and [manifest](manifest.json)."])
    (output / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--previous", type=Path, default=Path("results/floral_space_pairs_20261106"))
    parser.add_argument("--root", type=Path, default=Path("results/floral_space_20261106"))
    parser.add_argument("--output", type=Path, default=Path("results/floral_space_fixed_margins_20261107"))
    parser.add_argument("--per-chain", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261107)
    parser.add_argument("--burn-multiplier", type=int, default=12)
    parser.add_argument("--gap-multiplier", type=int, default=3)
    args = parser.parse_args()
    manifest = run_suite(args.previous, args.root, args.output, per_chain=args.per_chain,
                         seed=args.seed, burn_multiplier=args.burn_multiplier,
                         gap_multiplier=args.gap_multiplier)
    print(json.dumps({k: manifest[k] for k in ("status", "comparisons", "mixing_pass_count",
                                                        "significant_bh_0_05_exploratory")}, sort_keys=True))


if __name__ == "__main__":
    main()
