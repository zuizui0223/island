"""Observed turnover and nestedness in shared-name-excluded island floral profiles.

Consumes only aggregate pair observations from PR321. No inference, profile
identity, species/source rows, or new null distributions are generated.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

SCOPES = ("all", "direct", "all_single_state", "direct_single_state")
REGIMES = ("northern_midlatitude", "northern_high_latitude",
           "tropical", "southern_extratropical")
ENCODINGS = (("all", "direct"), ("all_single_state", "direct_single_state"))
OUT_COLS = ("evidence_scope", "analysis_regime", "island_a", "island_b",
            "n_shared_taxa", "n_profiles_a", "n_profiles_b",
            "profile_intersection", "profile_union", "beta_jaccard",
            "beta_turnover", "beta_nestedness", "turnover_share",
            "max_possible_jaccard", "richness_normalized_jaccard")


def partition(k_a: int, k_b: int, overlap: int) -> dict:
    """Additive Jaccard beta partition with richness-bound similarity."""
    if any(type(n) is not int for n in (k_a, k_b, overlap)):
        raise ValueError("Profile counts must be integers")
    if min(k_a, k_b) < 1 or overlap < 0 or overlap > min(k_a, k_b):
        raise ValueError("Invalid nonempty profile set/intersection counts")
    b, c = k_a - overlap, k_b - overlap
    union = overlap + b + c
    beta = (b + c) / union
    turnover = (2 * min(b, c) / (overlap + 2 * min(b, c))) if min(b, c) else 0.0
    nesting = max(0.0, beta - turnover)
    max_j = min(k_a, k_b) / max(k_a, k_b)
    return dict(beta_jaccard=beta, beta_turnover=turnover,
                beta_nestedness=nesting, turnover_share=turnover / beta if beta else None,
                max_possible_jaccard=max_j,
                richness_normalized_jaccard=(overlap / union) / max_j)


def read_csv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def write_csv(path: Path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
    else:
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


def analyze(source: Path, output: Path) -> list[dict]:
    """Require the original 16 matched cohorts and rebuild observed beta exactly."""
    source, output = Path(source), Path(output)
    required = {"evidence_scope", "analysis_regime", "island_a", "island_b",
                "n_profiles_a", "n_profiles_b", "profile_intersection",
                "profile_union", "profile_jaccard", "n_shared_taxa",
                "n_unique_a", "n_unique_b"}
    groups, ids = defaultdict(list), defaultdict(list)
    for row in read_csv(source):
        if not required <= row.keys():
            raise ValueError(f"Missing required columns: {sorted(required - row.keys())}")
        scope, regime = row["evidence_scope"], row["analysis_regime"]
        if scope not in SCOPES or regime not in REGIMES:
            raise ValueError("Unknown evidence scope or region")
        a, b = row["island_a"], row["island_b"]
        if not a or not b or a >= b:
            raise ValueError("Invalid unordered island pair")
        k_a, k_b, intersect = (int(row[x]) for x in
                               ("n_profiles_a", "n_profiles_b", "profile_intersection"))
        union, n_shared = int(row["profile_union"]), int(row["n_shared_taxa"])
        if (union != k_a + k_b - intersect or n_shared < 0
                or not 1 <= k_a <= int(row["n_unique_a"])
                or not 1 <= k_b <= int(row["n_unique_b"])):
            raise ValueError("Inconsistent stored profile/retained-taxon counts")
        old_j = float(row["profile_jaccard"])
        if not math.isfinite(old_j) or abs(old_j - intersect / union) > 1e-12:
            raise ValueError("Inconsistent stored profile Jaccard")
        derived = partition(k_a, k_b, intersect)
        key = (scope, regime)
        groups[key].append({col: row.get(col) for col in OUT_COLS if col not in derived} | derived)
        ids[key].append((a, b))
    if set(groups) != {(s, r) for s in SCOPES for r in REGIMES}:
        raise ValueError("Exactly 16 nonempty comparison groups required")
    for key, pair_ids in ids.items():
        if len(pair_ids) != len(set(pair_ids)):
            raise ValueError(f"Duplicate unordered pair: {key}")
    for all_scope, direct_scope in ENCODINGS:
        for regime in REGIMES:
            if ids[all_scope, regime] != ids[direct_scope, regime]:
                raise ValueError("All/Direct pair IDs are not matched")

    summaries, derived_pairs = [], []
    for scope in SCOPES:
        for regime in REGIMES:
            group = groups[scope, regime]
            endpoints = defaultdict(list)
            for r in group:
                endpoints[r["island_a"]].append(r)
                endpoints[r["island_b"]].append(r)
            def island_mean(metric):
                return mean(mean(r[metric] for r in adjacent) for adjacent in endpoints.values())
            positive = [r["turnover_share"] for r in group if r["turnover_share"] is not None]
            summary = dict(
                evidence_scope=scope, analysis_regime=regime, n_pairs=len(group),
                n_islands=len(endpoints),
                n_pairs_no_shared_taxa=sum(int(r["n_shared_taxa"]) == 0 for r in group),
                median_beta_jaccard=median(r["beta_jaccard"] for r in group),
                median_beta_turnover=median(r["beta_turnover"] for r in group),
                median_beta_nestedness=median(r["beta_nestedness"] for r in group),
                median_turnover_share=median(positive) if positive else None,
                fraction_turnover_dominant=mean(r["beta_turnover"] > r["beta_nestedness"] for r in group),
                median_richness_normalized_jaccard=median(r["richness_normalized_jaccard"] for r in group),
                island_balanced_mean_turnover=island_mean("beta_turnover"),
                island_balanced_mean_nestedness=island_mean("beta_nestedness"),
                island_balanced_mean_richness_normalized_jaccard=island_mean(
                    "richness_normalized_jaccard"),
            )
            summaries.append(summary)
            derived_pairs.extend(group)

    write_csv(output / "pairwise_beta_partition.csv.gz", OUT_COLS, derived_pairs)
    write_csv(output / "beta_partition_summary.csv", list(summaries[0]), summaries)
    checksum = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    manifest = {
        "status": "complete", "n_scope_region": 16, "n_pair_rows": len(derived_pairs),
        "source": str(source), "source_sha256": checksum(source),
        "output_sha256": {n: checksum(output / n) for n in
                          ("pairwise_beta_partition.csv.gz", "beta_partition_summary.csv")},
        "endpoint": "shared-taxon-name-exclusive observed profile sets",
        "identity": "beta_jaccard = beta_turnover + beta_nestedness PER PAIR",
        "richness_normalized_jaccard": "observed J / min(KA,KB)*max(KA,KB); not a null",
        "weights": "pair medians plus degree-balanced equal-island means",
        "inference": "descriptive only; no new significance tests or causal identification",
        "limitations": ["All/Direct and island pairs overlap", "Coarse categorical profiles",
                        "No genus/source/observation-process-conditioned null for these components",
                        "No hypervolume or evolutionary-convergence inference"]
    }
    # The algebraic definition above is J / (min(KA,KB) / max(KA,KB)).
    manifest["richness_normalized_jaccard"] = "observed J divided by [min(KA,KB)/max(KA,KB)]; not a null"
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = [
        "# Shared-taxon-exclusive floral dissimilarity: turnover or nestedness?",
        "",
        "This analysis partitions **observed** Jaccard dissimilarity of floral-profile",
        "sets after removing shared taxon names in each sampled island pair.",
        "Turnover is replacement of profile identities; nestedness-resultant dissimilarity",
        "reflects unequal profile richness. No new p-values are calculated.",
        "",
        "| Scope | Region | Pairs | Median turnover | Median nestedness | Median turnover share | Median relative Jaccard |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for x in summaries:
        share = x["median_turnover_share"]
        lines.append(f"| {x['evidence_scope']} | {x['analysis_regime']} | {x['n_pairs']:,} | "
                     f"{x['median_beta_turnover']:.3f} | {x['median_beta_nestedness']:.3f} | "
                     f"{share:.3f} | {x['median_richness_normalized_jaccard']:.3f} |"
                     if share is not None else
                     f"| {x['evidence_scope']} | {x['analysis_regime']} | {x['n_pairs']:,} | "
                     f"{x['median_beta_turnover']:.3f} | {x['median_beta_nestedness']:.3f} | NA | "
                     f"{x['median_richness_normalized_jaccard']:.3f} |")
    lines += [
        "", "Per pair: a=shared, b=only A, c=only B. Jaccard beta=(b+c)/(a+b+c),",
        "turnover=2min(b,c)/(a+2min(b,c)); nestedness-resultant=Jaccard beta−turnover.",
        "Turnover share is omitted for identical profile sets (beta=0).",
        "Relative Jaccard=observed J / [min(KA,KB)/max(KA,KB)] standardizes the",
        "**attainable maximum only**; it does not fix profile prevalence or genus composition.",
        "Component medians do **not** add to the median total; pair-level identities do.",
        "Pairs sharing an island are dependent. Equal-island means are separately in CSV.",
        "The [whole-island fixed-margin reference](../floral_space_fixed_margins_20261107/README.md)",
        "includes shared taxa and is NOT the null distribution of this pair-exclusive result.",
    ]
    (output / "README.md").write_text("\n".join(lines) + "\n")
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(
        "results/floral_space_pairs_20261106/pairwise_observed_and_null.csv.gz"))
    parser.add_argument("--output", type=Path, default=Path(
        "results/floral_space_beta_partition_20261108"))
    args = parser.parse_args()
    results = analyze(args.source, args.output)
    print(json.dumps({"status": "complete", "comparisons": len(results),
                      "rows": sum(r["n_pairs"] for r in results)}))


if __name__ == "__main__":
    main()
