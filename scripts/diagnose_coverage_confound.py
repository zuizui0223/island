#!/usr/bin/env python3
"""Reproduce the coverage-confound diagnosis from frozen repository files.

Three questions, in order:

1. Is the readiness gate itself the binding constraint, or the evidence?
2. If a partial flora sample is used, how wrong is the composition estimate?
3. Is the coverage shortfall correlated with island isolation -- that is, with
   the predictor the study is testing?

Inputs are the committed strict coverage ledger from the PR #132 branch plus the
frozen GBIF flora artifacts on main. No network access, no CI artifacts.

    python3 scripts/diagnose_coverage_confound.py --ledger <strict_coverage.csv.gz>

The ledger lives on `codex/open-web-evidence-pilot` at
data/v2/staging/traits/direct_llm_pilot/
  20260807_europe_pmc_reproduction_integrated/strict_species_axis_coverage.csv.gz
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import re
import statistics
from pathlib import Path

MASTER = Path("data/v2/staging/gbif/collected/island_taxa.csv")
OCCURRENCES = Path("data/v2/staging/gbif/collected/island_species_occurrences.csv.gz")

COLOUR = "flower_colour"
STRUCTURE = "floral_structural_complexity"
REPRODUCTION = "reproductive_assurance"
AXES = (COLOUR, STRUCTURE, REPRODUCTION)

DIRECT_QUALITIES = {"high", "medium"}


def load_master() -> dict[str, tuple[int, int, str]]:
    """species -> (n_islands, n_gbif_records, family)."""
    out: dict[str, tuple[int, int, str]] = {}
    with MASTER.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                out[row["accepted_species"]] = (
                    int(row["n_islands"]),
                    int(row["n_records"]),
                    row["family"],
                )
            except (KeyError, ValueError):
                continue
    return out


def load_islands() -> dict[str, set[str]]:
    islands: dict[str, set[str]] = collections.defaultdict(set)
    with gzip.open(OCCURRENCES, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            islands[row["island_id"]].add(row["species"])
    return islands


def load_ledger(path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]], dict[str, str]]:
    """Return direct-evidence species per axis, unresolved species per axis, colour values."""
    direct: dict[str, set[str]] = collections.defaultdict(set)
    unresolved: dict[str, set[str]] = collections.defaultdict(set)
    colour_value: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            species, axis = row["accepted_species"], row["axis"]
            quality = (row.get("quality") or "").strip()
            if quality in DIRECT_QUALITIES:
                direct[axis].add(species)
                if axis == COLOUR:
                    match = re.search(r'\["([^"]+)"', row.get("trait_composition") or "")
                    if match:
                        colour_value[species] = match.group(1)
            if not quality:
                unresolved[axis].add(species)
    return direct, unresolved, colour_value


def ready_islands(
    islands: dict[str, set[str]], covered: set[str], min_species: int, min_fraction: float
) -> set[str]:
    out = set()
    for island_id, flora in islands.items():
        hits = len(flora & covered)
        if hits >= min_species and hits / len(flora) >= min_fraction:
            out.add(island_id)
    return out


def report_gate_sensitivity(islands, direct) -> None:
    print("1) GATE SENSITIVITY (direct evidence only)\n")
    print(f'{"min species / min fraction":<28}{"colour":>9}{"structure":>11}{"reprod":>9}{"ALL 3":>9}')
    for min_species, min_fraction in (
        (30, 0.5), (30, 0.4), (30, 0.3), (30, 0.2), (30, 0.1), (30, 0.0),
        (50, 0.3), (50, 0.2), (100, 0.2),
    ):
        sets = [ready_islands(islands, direct[a], min_species, min_fraction) for a in AXES]
        both = sets[0] & sets[1] & sets[2]
        label = f"{min_species} / {min_fraction}"
        print(f"{label:<28}{len(sets[0]):>9,}{len(sets[1]):>11,}{len(sets[2]):>9,}{len(both):>9,}")


def report_bias(islands, colour_value, master) -> None:
    """Subsample high-coverage islands the way real coverage accrues, and measure drift."""
    truth_islands = []
    for island_id, flora in islands.items():
        covered = [s for s in flora if s in colour_value]
        if len(covered) >= 50 and len(covered) / len(flora) >= 0.8:
            truth_islands.append(covered)

    def composition(species: list[str]) -> dict[str, float]:
        counts = collections.Counter(colour_value[s] for s in species)
        total = sum(counts.values())
        return {k: v / total for k, v in counts.items()}

    print(f"\n2) BIAS OF A PARTIAL SAMPLE ({len(truth_islands):,} ground-truth islands)\n")
    print("   Selection mimics how coverage actually accrues: most-widespread species first.\n")
    print(f'{"subsample":<12}{"median n":>10}{"median max shift":>20}{"median TVD":>13}{"p90 TVD":>10}')
    for fraction in (0.5, 0.3, 0.2, 0.1):
        shifts, tvds, counts = [], [], []
        for covered in truth_islands:
            truth = composition(covered)
            keep = max(1, int(len(covered) * fraction))
            sample = sorted(covered, key=lambda s: -master.get(s, (0, 0, ""))[0])[:keep]
            estimate = composition(sample)
            counts.append(len(sample))
            deltas = {
                k: abs(estimate.get(k, 0.0) - truth.get(k, 0.0))
                for k in set(truth) | set(estimate)
            }
            shifts.append(max(deltas.values()))
            tvds.append(sum(deltas.values()) / 2)
        p90 = sorted(tvds)[int(0.9 * len(tvds))]
        print(
            f'{f"{int(fraction * 100)}%":<12}{statistics.median(counts):>10.0f}'
            f"{statistics.median(shifts):>20.3f}{statistics.median(tvds):>13.3f}{p90:>10.3f}"
        )


def _spearman(xs: list[float], ys: list[float]) -> float:
    def ranks(values: list[float]) -> list[int]:
        order = sorted(range(len(values)), key=lambda i: values[i])
        out = [0] * len(values)
        for position, index in enumerate(order):
            out[index] = position
        return out

    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    numerator = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    denominator = (
        sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)
    ) ** 0.5
    return numerator / denominator


def report_confound(islands, direct, master) -> None:
    """Endemism share is the isolation proxy: it needs no external covariate file."""
    rows = []
    for flora in islands.values():
        if len(flora) < 30:
            continue
        endemic = sum(1 for s in flora if master.get(s, (0, 0, ""))[0] == 1) / len(flora)
        coverage = len(flora & direct[COLOUR]) / len(flora)
        rows.append((endemic, coverage))

    print(f"\n3) COVERAGE vs ISOLATION ({len(rows):,} islands with >=30 species)\n")
    print(f'{"endemic share of flora":<26}{"islands":>9}{"median colour coverage":>25}')
    for low, high, label in (
        (0.0, 0.02, "0-2%"), (0.02, 0.05, "2-5%"), (0.05, 0.10, "5-10%"),
        (0.10, 0.20, "10-20%"), (0.20, 0.40, "20-40%"), (0.40, 1.01, "40%+"),
    ):
        selected = [c for e, c in rows if low <= e < high]
        if selected:
            print(f"{label:<26}{len(selected):>9,}{statistics.median(selected) * 100:>24.1f}%")

    rho = _spearman([r[0] for r in rows], [r[1] for r in rows])
    print(f"\n   Spearman rho (endemism vs coverage) = {rho:+.3f}")
    print("   Negative means the coverage shortfall tracks the study's own predictor.")


def report_missing_kind(unresolved, master) -> None:
    """Split the unresolved set into 'no data exists' versus 'method failed'."""
    everything = set().union(*unresolved.values())
    print("\n4) IS THE UNRESOLVED SET MISSING DATA, OR MISSING A METHOD?\n")
    print(f'{"group":<38}{"species":>9}{"median GBIF records":>22}{"<20 recs":>10}')
    groups = (
        ("unresolved, 1 island", lambda n: n == 1),
        ("unresolved, 2-3 islands", lambda n: 2 <= n <= 3),
        ("unresolved, 4-10 islands", lambda n: 4 <= n <= 10),
        ("unresolved, 11+ islands", lambda n: n >= 11),
    )
    for label, predicate in groups:
        selected = [s for s in everything if s in master and predicate(master[s][0])]
        records = [master[s][1] for s in selected]
        if records:
            thin = 100 * sum(1 for r in records if r < 20) / len(records)
            print(
                f"{label:<38}{len(selected):>9,}"
                f"{statistics.median(records):>22.0f}{thin:>9.1f}%"
            )

    tail = {s for s in everything if s in master and master[s][0] <= 3}
    families = collections.Counter(master[s][2] for s in tail)
    print(f"\n   Tail (unresolved, <=3 islands) = {len(tail):,} species. Top families:")
    for family, count in families.most_common(6):
        print(f"      {family:<20}{count:>7,}  {100 * count / len(tail):5.1f}%")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, required=True)
    args = parser.parse_args()

    master = load_master()
    islands = load_islands()
    direct, unresolved, colour_value = load_ledger(args.ledger)

    print(f"islands: {len(islands):,}   master species: {len(master):,}")
    print(f"species with a direct colour value: {len(colour_value):,}\n")

    report_gate_sensitivity(islands, direct)
    report_bias(islands, colour_value, master)
    report_confound(islands, direct, master)
    report_missing_kind(unresolved, master)


if __name__ == "__main__":
    main()
