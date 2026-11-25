"""Post-hoc robustness of matched floral repertoire to climate/coverage balance.

Primary frozen matching is unchanged. Report stricter calipers and each main
biogeographic regime separately; results are sensitivity checks, not new
independent confirmation or additional discoveries.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .pipeline import REF_INPUT, _write_csv
from .repertoire import (
    PERMUTATIONS, REGIONS, _paired_swap_statistic, _two_sided,
    _upper_tail, match_origin_islands,
)

COLUMNS = ("matching", "region", "scope", "n_pairs", "n_distinct_islands",
           "rao_q_diff", "p_rao_q_two_sided",
           "kernel_statistic", "kernel_reference_median", "kernel_p_upper",
           "exact_statistic", "exact_reference_median", "exact_p_upper",
           "max_abs_smd_pooled", "status")


def analyze(input_root: Path, g0_root: Path, crosswalk_path: Path, output: Path,
            *, permutations: int = PERMUTATIONS, seed: int = 20261108) -> dict:
    input_root, g0_root, crosswalk_path, output = map(
        Path, (input_root, g0_root, crosswalk_path, output))
    output.mkdir(parents=True, exist_ok=True)
    crosswalk = {r["island_id"]: r for r in csv_rows(crosswalk_path)}
    all_rows = list(csv_rows(g0_root / "all_islands.csv"))
    direct_rows = list(csv_rows(g0_root / "direct_islands.csv"))
    main, main_audit = match_origin_islands(all_rows, direct_rows, crosswalk)
    # The stricter thresholds are post-hoc diagnostic sensitivity choices,
    # made to reduce residual PC3 and trait-coverage imbalance.
    strict, strict_audit = match_origin_islands(
        all_rows, direct_rows, crosswalk,
        coverage_max=0.12, climate_pc3_max=1.0)
    geographic, geographic_audit = match_origin_islands(
        all_rows, direct_rows, crosswalk,
        max_pair_geographic_km=5000.0)
    near_region, near_region_audit = match_origin_islands(
        all_rows, direct_rows, crosswalk,
        max_pair_geographic_km=3000.0)
    profiles = load_profiles(input_root / REF_INPUT)
    space = ProfileSpace(profiles)
    kernel = 1.0 - space.matrix
    carriers = {
        scope: load_carriers(input_root / "floral_space_20261106" /
                             scope / "island_profile_carriers.csv.gz")
        for scope in ("all", "direct")
    }
    records = []
    for cohort, pairs, audit in (
        ("original", main, main_audit),
        ("stricter_covariate_balance", strict, strict_audit),
        ("within_5000km_source_proxy", geographic, geographic_audit),
        ("within_3000km_source_proxy", near_region, near_region_audit),
    ):
        for region in ("all",) + REGIONS:
            chosen = [p for p in pairs if region == "all" or p["region"] == region]
            for scope in ("all", "direct"):
                base = {
                    "matching": cohort, "region": region,
                    "scope": scope, "n_pairs": len(chosen),
                    "n_distinct_islands": 2 * len(chosen),
                    "max_abs_smd_pooled": max(
                        (row["max_abs_smd_pooled"] for row in audit["regional_diagnostics"]
                         if row["region"] in (REGIONS if region == "all" else (region,))
                         and row["max_abs_smd_pooled"] is not None),
                        default=None),
                }
                if len(chosen) < 8:
                    records.append({**base, "status": "insufficient_pairs"})
                    continue
                source = carriers[scope]
                o = np.zeros((len(chosen), len(space.ids)), dtype=float)
                c = np.zeros_like(o)
                qdiff = np.zeros(len(chosen), dtype=float)
                for i, pair in enumerate(chosen):
                    for name, matrix in ((pair["continental_id"], c),
                                         (pair["oceanic_id"], o)):
                        ids = [space.lookup[p] for p in source[name]]
                        matrix[i, ids] = 1.0 / len(ids)
                    qdiff[i] = (
                        space.island(source[pair["oceanic_id"]])["rao_q_profile_equal"]
                        - space.island(source[pair["continental_id"]])["rao_q_profile_equal"]
                    )
                rng = np.random.default_rng(
                    seed + (12 if cohort.startswith("stricter") else 0)
                    + (22 if cohort.startswith("within_5000") else
                       32 if cohort.startswith("within_3000") else 0)
                    + (3 if region == "tropical" else 1 if region == "northern_midlatitude" else 0)
                )
                signs = rng.choice([-1.0, 1.0], size=(permutations, len(chosen)))
                delta = o - c
                q_obs = float(qdiff.mean())
                q_null = (signs * qdiff).mean(axis=1)
                kern_obs, kern_null = _paired_swap_statistic(delta, kernel, signs)
                exact_obs, exact_null = _paired_swap_statistic(
                    delta, np.eye(len(space.ids)), signs)
                records.append({
                    **base, "status": "exploratory_posthoc_sensitivity",
                    "rao_q_diff": q_obs,
                    "p_rao_q_two_sided": _two_sided(q_obs, q_null),
                    "kernel_statistic": kern_obs,
                    "kernel_reference_median": float(np.median(kern_null)),
                    "kernel_p_upper": _upper_tail(kern_obs, kern_null),
                    "exact_statistic": exact_obs,
                    "exact_reference_median": float(np.median(exact_null)),
                    "exact_p_upper": _upper_tail(exact_obs, exact_null),
                })
    _write_csv(output / "sensitivity_repertoire.csv", records, COLUMNS)
    manifest = {
        "schema": "fd_origin_repertoire_sensitivity_20261108",
        "status": "exploratory_posthoc_not_independent_confirmation",
        "original_matching": main_audit,
        "strict_matching": strict_audit,
        "geographic_5000km_matching": geographic_audit,
        "geographic_3000km_matching": near_region_audit,
        "n_records": len(records),
        "n_strict_pairs": len(strict),
        "n_geographic_5000km_pairs": len(geographic),
        "n_geographic_3000km_pairs": len(near_region),
        "different_repertoire_tests_are_nested_sensitivities": True,
        "not_conditioned_on": ["genus", "taxon identities", "source pool",
                               "site-specific pollinator identity"],
        "not_multiplicity_corrected": "secondary sensitivity raw p values only; primary hypothesis family belongs to repertoire.py",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_strict_pairs": len(strict), "records": records},
                     ensure_ascii=False))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("input", "g0", "crosswalk", "output"):
        parser.add_argument("--" + arg, type=Path, required=True)
    args = parser.parse_args()
    analyze(args.input, args.g0, args.crosswalk, args.output)


if __name__ == "__main__":
    main()
