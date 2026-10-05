"""Population-level within-lineage island pilot.

The pilot uses published population-level multilocus outcrossing rates (tm) for
Arabidopsis lyrata and source-reported paired island-mainland inference for
Phormium tenax.

It is intentionally separate from the Chapter 1 assemblage analysis. The
Arabidopsis test is an exact label-permutation test conditional on the 18 sampled
populations; it is not equivalent to the global oceanic-isolation estimand.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd


def exact_island_permutation(table: pd.DataFrame) -> dict[str, float | int]:
    """Compare mean tm of explicitly named island sites with the other sites."""
    tm = pd.to_numeric(table["tm"], errors="raise").to_numpy(float)
    island = table["island_status"].astype(str).eq("island").to_numpy()
    n_island = int(island.sum())
    if len(tm) != 18 or n_island != 3:
        raise ValueError(
            f"expected 18 populations and 3 island sites, got {len(tm)} / {n_island}"
        )

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
        "exact_two_sided_p": float(
            np.mean(np.abs(null) >= abs(observed) - tol)
        ),
        "exact_one_sided_lower_tm_p": float(np.mean(null <= observed + tol)),
        "n_exact_labelings": int(len(null)),
    }


def run_pilot(
    arabidopsis: pd.DataFrame,
    phormium: pd.DataFrame,
) -> dict[str, object]:
    """Run the two-species feasibility pilot and return a serializable summary."""
    result = exact_island_permutation(arabidopsis)

    if len(phormium) != 2 or not phormium["source_result"].eq(
        "no_significant_difference"
    ).all():
        raise ValueError(
            "expected two source-reported Phormium island-mainland null comparisons"
        )

    return {
        "contract": "within_lineage_population_pilot_20261005_v1",
        "arabidopsis": result,
        "phormium": {
            "n_matched_island_mainland_pairs": 2,
            "source_reported_seedling_outcrossing_difference": (
                "no significant difference in either pair"
            ),
        },
        "interpretation": (
            "The currently recoverable population-level evidence does not support a "
            "simple within-species rule that island populations are uniformly more "
            "selfing. Arabidopsis contains one strongly selfing island population "
            "alongside two strongly outcrossing island populations, and the paired "
            "Phormium study reports no significant seedling-outcrossing difference "
            "in either island-mainland pair."
        ),
        "claim_boundary": (
            "This is a two-species feasibility pilot. Arabidopsis sites are "
            "freshwater Great Lakes populations rather than the Chapter 1 global "
            "oceanic-isolation estimand; Phormium values are retained as the "
            "source-reported paired result. Neither route separates genetic "
            "evolution from plasticity without additional evidence."
        ),
    }


def write_outputs(
    arabidopsis_csv: Path,
    phormium_csv: Path,
    output_dir: Path,
) -> dict[str, object]:
    """Read curated inputs, run the pilot, and write deterministic outputs."""
    arabidopsis = pd.read_csv(arabidopsis_csv)
    phormium = pd.read_csv(phormium_csv)
    summary = run_pilot(arabidopsis, phormium)

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "RESULT.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    pd.DataFrame([summary["arabidopsis"]]).to_csv(
        output_dir / "arabidopsis_exact_island_test.csv", index=False
    )
    return summary
