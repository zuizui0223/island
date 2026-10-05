"""Run the first population-level within-lineage island pilot.

This route is deliberately separate from the Chapter 1 assemblage analysis. It uses
published population-level multilocus outcrossing rates (tm) for Arabidopsis lyrata and
an independent paired island-mainland result for Phormium tenax.

The Arabidopsis test asks only whether the three explicitly named Great Lakes island
sites have different tm from the remaining 15 sampled populations. It is an exact
label-permutation test conditional on the observed 18 populations. It is not equivalent
to the Chapter 1 oceanic isolation gradient and cannot by itself identify heritable
evolution rather than environmental effects.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def exact_island_permutation(table: pd.DataFrame) -> dict[str, float | int]:
    tm = pd.to_numeric(table["tm"], errors="raise").to_numpy(float)
    island = table["island_status"].astype(str).eq("island").to_numpy()
    n_island = int(island.sum())
    if len(tm) != 18 or n_island != 3:
        raise ValueError(f"expected 18 populations and 3 island sites, got {len(tm)} / {n_island}")
    observed = float(tm[island].mean() - tm[~island].mean())

    diffs: list[float] = []
    for combo in itertools.combinations(range(len(tm)), n_island):
        mask = np.zeros(len(tm), dtype=bool)
        mask[list(combo)] = True
        diffs.append(float(tm[mask].mean() - tm[~mask].mean()))
    null = np.asarray(diffs, dtype=float)
    tol = 1e-12
    return {
        "n_populations": int(len(tm)),
        "n_island_populations": n_island,
        "n_nonisland_populations": int((~island).sum()),
        "island_mean_tm": float(tm[island].mean()),
        "nonisland_mean_tm": float(tm[~island].mean()),
        "island_minus_nonisland_tm": observed,
        "island_minus_nonisland_selfing_fraction": float(-observed),
        "exact_two_sided_p": float(np.mean(np.abs(null) >= abs(observed) - tol)),
        "exact_one_sided_lower_tm_p": float(np.mean(null <= observed + tol)),
        "n_exact_labelings": int(len(null)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--arabidopsis-csv", type=Path, required=True)
    parser.add_argument("--phormium-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    arabidopsis = pd.read_csv(args.arabidopsis_csv)
    phormium = pd.read_csv(args.phormium_csv)
    result = exact_island_permutation(arabidopsis)

    if len(phormium) != 2 or not phormium["source_result"].eq("no_significant_difference").all():
        raise ValueError("expected two source-reported Phormium island-mainland null comparisons")

    summary = {
        "contract": "within_lineage_population_pilot_20261005_v1",
        "arabidopsis": result,
        "phormium": {
            "n_matched_island_mainland_pairs": 2,
            "source_reported_seedling_outcrossing_difference": "no significant difference in either pair",
        },
        "interpretation": (
            "The currently recoverable population-level evidence does not support a simple "
            "within-species rule that island populations are uniformly more selfing. "
            "Arabidopsis contains one strongly selfing island population alongside two "
            "strongly outcrossing island populations, and the paired Phormium study reports "
            "no significant seedling-outcrossing difference in either island-mainland pair."
        ),
        "claim_boundary": (
            "This is a two-species feasibility pilot. Arabidopsis sites are freshwater Great "
            "Lakes populations rather than the Chapter 1 global oceanic isolation estimand; "
            "Phormium values are retained as the source-reported paired result. Neither route "
            "separates genetic evolution from plasticity without additional evidence."
        ),
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "RESULT.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    pd.DataFrame([result]).to_csv(
        args.output_dir / "arabidopsis_exact_island_test.csv", index=False
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
