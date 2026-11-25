"""Five-fold conservative sensitivity for marginal-fixed floral source scores.

Input is NON-IDENTIFYING 10-row aggregate with one fold-level mean for each
All/Direct scope. This deliberately never accesses source species/genus,
island IDs, geography coordinates or private donor records.

Exact 2^5 fold sign configurations are a STRESS TEST only: crossfit folds
reuse much of the training pool and are not established exchangeable
independent ecological replicates. This is NOT a confirmatory p-value.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

FIELDS = (
    "scope", "endpoint", "n_folds", "n_islands", "n_positive_fold_means",
    "island_weighted_mean", "leave_one_fold_out_min",
    "leave_one_fold_out_max", "sign_configurations",
    "two_sided_sign_tail_if_exchangeable", "status",
)
ENDPOINTS = (
    "mean_joint_before", "mean_joint_marginal_fixed",
    "mean_brier_marginal_fixed",
)


def fold_sign_sensitivity(table: pd.DataFrame) -> list[dict]:
    """Exact 32 sign configurations and leave-one-geography-fold stability."""
    required = {"scope", "fold", "n_islands", *ENDPOINTS}
    if not required.issubset(table):
        raise ValueError("Missing nonidentifying geographic fold summary columns")
    out = []
    for scope, data in table.groupby("scope", sort=True):
        part = data.sort_values("fold")
        if (len(part) != 5 or part.fold.nunique() != 5 or
                set(part.fold.astype(int)) != set(range(5))):
            raise ValueError("Exactly five geographic fold summary rows expected")
        ns = part.n_islands.to_numpy(dtype=float)
        if not np.isfinite(ns).all() or np.any(ns <= 0):
            raise ValueError("Invalid fold island support")
        for endpoint in ENDPOINTS:
            diffs = part[endpoint].to_numpy(dtype=float)
            if not np.isfinite(diffs).all():
                raise ValueError("Nonfinite heldout fold score")
            observed = float(np.average(diffs, weights=ns))
            signs = np.asarray(list(itertools.product((-1., 1.), repeat=5)))
            null = (signs @ (ns * diffs)) / ns.sum()
            eps = np.finfo(float).eps * 128 * max(1.0, abs(observed))
            tail = float(np.mean(np.abs(null) >= abs(observed)-eps))
            leave_out = [(float((ns@diffs-ns[i]*diffs[i])/(ns.sum()-ns[i])))
                         for i in range(5)]
            out.append({
                "scope": scope, "endpoint": endpoint,
                "n_folds": 5, "n_islands": int(ns.sum()),
                "n_positive_fold_means": int(np.sum(diffs > 0)),
                "island_weighted_mean": observed,
                "leave_one_fold_out_min": min(leave_out),
                "leave_one_fold_out_max": max(leave_out),
                "sign_configurations": len(null),
                "two_sided_sign_tail_if_exchangeable": tail,
                "status": "coarse_diagnostic_not_confirmatory_folds_share_donors",
            })
    return out


def run(summary: Path, output: Path) -> dict:
    summary, output = Path(summary), Path(output)
    observed = pd.read_csv(summary)
    scores = fold_sign_sensitivity(observed)
    output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(scores, columns=FIELDS).to_csv(
        output/"fold_coarse_sign_sensitivity.csv", index=False)
    with summary.open("rb") as handle:
        digest = hashlib.file_digest(handle, "sha256").hexdigest()
    manifest = {
        "schema": "fd_marginal_fixed_five_fold_sign_sensitivity_v1",
        "input_summary_sha256": digest,
        "status": "posthoc_coarse_geographic_fold_robustness",
        "test_design": "5 geographic heldout folds, 32 sign configurations, n-island weighted means, leave-one-fold-out influence",
        "critical_warning": "Fold outcomes share overlapping training donor sets and source taxonomy; exact signs are only a stress test, not randomization p-values or independent evolutionary replicates",
        "no_individual_taxon_island_source_data": True,
    }
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({"n_scopes": len(set(x["scope"] for x in scores)),
                      "n_diagnostic_rows": len(scores),
                      "interpretation": "not independent hypothesis significance"}))
    return manifest


def main() -> None:
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args=p.parse_args()
    run(args.summary,args.output)


if __name__=="__main__":
    main()
