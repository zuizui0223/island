"""Origin-aware exploratory floral functional-diversity contrasts.

'geological_origin' = source-verified Zizka 2022 supplementary labels only.
'LGM_connectivity' = Weigelt 2013 GMMC historical mainland connection
and MUST NOT be relabeled oceanic geological origin if unconnected.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path

import numpy as np

from .core import csv_rows
from .origin import run_crosswalk
from .pipeline import SCOPES, _write_csv, file_sha256

GROUPS = (
    ("geological_origin", "continental", "oceanic"),
    ("LGM_connectivity", "lgm_connected", "lgm_unconnected"),
)
LOCAL_STATS = ("rao_q_profile_equal", "rao_q", "kde90_area_pc12_equal_h035")
PAIR_STATS = ("functional_nn_dissimilarity", "mean_cross_profile_distance")


def _num(value):
    try:
        x = float(value)
        return x if np.isfinite(x) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def _median(rows: list[dict], key: str) -> float | None:
    data = np.asarray([_num(x.get(key)) for x in rows], dtype=float)
    data = data[np.isfinite(data)]
    return float(np.median(data)) if len(data) else None


def descriptive(islands: list[dict], pairs: list[dict], *, key: str,
                first: str, second: str, scope: str) -> tuple[list[dict], list[dict]]:
    rows, pair_rows = [], []
    loc = {i["island_id"]: i for i in islands}
    groups = ("all",) + tuple(sorted({r.get("analysis_regime", "") for r in islands}))
    for regime in groups:
        eligible = [r for r in islands if regime == "all" or r["analysis_regime"] == regime]
        for label in (first, second, "mixed", "unresolved"):
            sample = [r for r in eligible if r.get(key, "unresolved") == label]
            rows.append({
                "scope": scope, "axis": key, "regime": regime, "group": label,
                "n_islands": len(sample),
                "n_unique_spatial_blocks": len({r.get("spatial_block") for r in sample}),
                "median_n_names": _median(sample, "n_names"),
                "median_joint_coverage": _median(sample, "joint_coverage"),
                "median_log_isolation": _median(
                    [{"value": math.log1p(_num(r["distance_to_continent_km"]))}
                     for r in sample if _num(r["distance_to_continent_km"]) >= 0], "value"),
                **{f"median_{stat}": _median(sample, stat) for stat in LOCAL_STATS},
            })
        for label in (first, second):
            use = []
            endpoints = set()
            for r in pairs:
                if regime != "all" and r["analysis_regime"] != regime:
                    continue
                a, b = loc.get(r["island_a"]), loc.get(r["island_b"])
                if a is None or b is None:
                    continue
                if a.get(key) == label and b.get(key) == label:
                    use.append(r)
                    endpoints.update((r["island_a"], r["island_b"]))
            pair_rows.append({"scope": scope, "axis": key, "regime": regime,
                              "group": label, "n_pairs": len(use),
                              "n_distinct_endpoints": len(endpoints),
                              **{f"median_{stat}": _median(use, stat)
                                 for stat in PAIR_STATS}})
    return rows, pair_rows


def _ols_origin_effect(data: list[dict], class_key: str, one: str, two: str,
                       stat: str, *, seed: int) -> dict:
    """Fixed-design wild cluster bootstrap; associations, not causal effects."""
    items = []
    for r in data:
        if r.get(class_key) not in (one, two):
            continue
        vals = [_num(r.get(col)) for col in (
            stat, "area_km2", "distance_to_continent_km", "n_names", "joint_coverage",
            "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4"
        )]
        if not all(np.isfinite(vals)) or min(vals[1:4]) < 0 or vals[1] == 0 or vals[3] == 0:
            continue
        items.append((r, vals))
    counts = Counter(r[class_key] for r, vals in items)
    blocks = {r["spatial_block"] for r, vals in items}
    result = {"origin_axis": class_key, "response": stat, "n_records": len(items),
              "n_first": counts[one], "n_second": counts[two],
              "n_spatial_blocks": len(blocks), "contrast": f"{two} minus {one}",
              "status": "insufficient_matched_support"}
    if min(counts[one], counts[two]) < 35 or len(blocks) < 8:
        return result
    regions = sorted({r["analysis_regime"] for r, vals in items})
    if len({r["analysis_regime"] for r,vals in items if r[class_key] == one} &
           {r["analysis_regime"] for r,vals in items if r[class_key] == two}) == 0:
        result["status"] = "no_regional_overlap"
        return result
    y = np.array([vals[0] for r, vals in items])
    x = np.array([
        [1, int(r[class_key] == two),
         math.log(vals[1]), math.log1p(vals[2]), math.log(vals[3]), vals[4],
         vals[5], vals[6], vals[7], vals[8],
         *[int(r["analysis_regime"] == reg) for reg in regions[1:]]]
        for r, vals in items], dtype=float)
    # Center and standardize numeric covariates; binary group indicator untouched.
    x[:, 2:10] = (x[:, 2:10] - x[:, 2:10].mean(axis=0)) / np.maximum(
        x[:, 2:10].std(axis=0), 1.0e-9)
    # Preserve intercept and origin contrast. Retain a nuisance covariate
    # only if it adds numerical rank; otherwise constant indicators/collinear
    # climate features make the model non-identifiable for no scientific gain.
    keep = [0, 1]
    if np.linalg.matrix_rank(x[:, keep]) != 2:
        result["status"] = "nonidentifiable_origin_contrast"
        return result
    for column in range(2, x.shape[1]):
        if np.linalg.matrix_rank(x[:, keep + [column]]) > len(keep):
            keep.append(column)
    result["used_design_columns"] = len(keep)
    result["dropped_collinear_columns"] = x.shape[1] - len(keep)
    x = x[:, keep]
    estimate = float(np.linalg.lstsq(x, y, rcond=None)[0][1])
    rng = np.random.default_rng(seed)
    byblock = {}
    for i, (r, vals) in enumerate(items):
        byblock.setdefault(r["spatial_block"], []).append(i)
    group_ids = sorted(byblock)
    # Wild-cluster residual multiplier bootstrap preserves ALL original rows,
    # origin representation and fixed design rank at each iteration. This is
    # preferable to resampling blocks when rare island types cause a high
    # proportion of non-identifiable bootstrap replicates.
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    fitted = x @ beta
    residual = y - fitted
    block_ix = np.array([group_ids.index(r["spatial_block"]) for r, vals in items])
    coefficients = []
    for repeat in range(999):
        signs = rng.choice([-1.0, 1.0], size=len(group_ids), replace=True)
        generated = fitted + residual * signs[block_ix]
        coefficients.append(float(np.linalg.lstsq(x, generated, rcond=None)[0][1]))
    result["coefficient"] = estimate
    result["bootstrap_method"] = "wild_cluster_Rademacher_999_fixed_design"
    result["n_bootstraps_usable"] = len(coefficients)
    result["status"] = "exploratory_adjusted_wild_cluster_bootstrap"
    result["bootstrap_ci_025"], result["bootstrap_ci_975"] = [
        float(v) for v in np.quantile(coefficients, [0.025, 0.975])]
    return result


def run_analysis(g0_root: Path, hv_root: Path, output: Path, crosswalk_root: Path,
                 *, fetch: bool = True) -> dict:
    g0_root, hv_root, output, crosswalk_root = map(
        Path, (g0_root, hv_root, output, crosswalk_root))
    output.mkdir(parents=True, exist_ok=True)
    if fetch:
        source = run_crosswalk(crosswalk_root)
    else:
        source = json.loads((crosswalk_root / "crosswalk_audit.json").read_text())
    origin = {r["island_id"]: r for r in csv_rows(crosswalk_root / "island_origin_crosswalk.csv")}
    if len(origin) != source["target_gshhg"]:
        raise ValueError("Crosswalk count differs from geography manifest")
    summary_rows, pair_rows, adjusted_rows, cohort_counts = [], [], [], []
    for scope in SCOPES:
        g0 = list(csv_rows(g0_root / f"{scope}_islands.csv"))
        hv = {r["island_id"]: r for r in csv_rows(hv_root / f"{scope}_hypervolumes.csv")}
        pairs = list(csv_rows(g0_root / f"{scope}_pairs.csv"))
        if set(hv) != {r["island_id"] for r in g0}:
            raise ValueError(f"{scope}: KDE and G0 island supports do not match")
        local = []
        for g in g0:
            island = g["island_id"]
            if island not in origin:
                raise ValueError("Missing geographic island from origin crosswalk")
            local.append({**g, **{k: v for k,v in origin[island].items() if k not in g},
                          **{k:v for k,v in hv[island].items()
                             if k.startswith("kde90_area_")}})
        cohort_counts.append({
            "scope": scope, "eligible_islands": len(local),
            "known_geological_class": sum(r["geological_origin"] in ("oceanic", "continental") for r in local),
            "known_LGM_connection": sum(r["LGM_connectivity"] in ("lgm_connected", "lgm_unconnected") for r in local),
            "n_pair_rows": len(pairs),
        })
        for key, one, two in GROUPS:
            a, b = descriptive(local, pairs, key=key, first=one, second=two, scope=scope)
            summary_rows.extend(a)
            pair_rows.extend(b)
            for response in LOCAL_STATS:
                fitted = _ols_origin_effect(local, key, one, two, response,
                                            seed=20261108 + len(adjusted_rows))
                adjusted_rows.append({"scope": scope, **fitted})
    if not summary_rows:
        raise ValueError("No classified islands in G0")
    _write_csv(output / "stratified_island_medians.csv", summary_rows, tuple(summary_rows[0]))
    _write_csv(output / "stratified_pair_medians.csv", pair_rows, tuple(pair_rows[0]))
    # Some exploratory adjusted contrasts have missing CI; use stable union fields.
    reg_fields = list(dict.fromkeys(k for row in adjusted_rows for k in row))
    _write_csv(output / "adjusted_origin_contrasts.csv", adjusted_rows, tuple(reg_fields))
    report = {
        "status": "exploratory_origin_split_with_unresolved_retained",
        "geological_origin_evidence": "Zizka 2022 manually sourced 83 missing GIFT labels ONLY",
        "historical_connectivity_evidence": "Weigelt et al. 2013 GMMC; 0 unconnected is NOT synonymous with oceanic origin",
        "classified_source_audit": source,
        "cohorts": cohort_counts,
        "regression": "four overlapping scope sensitivities; log(area)+log1p(distance)+log(name richness)+coverage+climate PCs+regime; 999-draw Rademacher wild-spatial-block bootstrap CI, not multiplicity corrected",
        "pair_dependence": "Pairs reuse islands, no independent-pair p values",
        "test_scope": "Descriptive comparison of recorded floral phenotype space only; source/lineage and native status not matched",
        "files": ["stratified_island_medians.csv", "stratified_pair_medians.csv",
                  "adjusted_origin_contrasts.csv"],
    }
    (output / "manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": report["status"], "cohorts": cohort_counts,
                      "contrast_status": [{"axis":r["origin_axis"],"response":r["response"],
                                          "status":r["status"],
                                          "estimate":r.get("coefficient")}
                                         for r in adjusted_rows]}, ensure_ascii=False))
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--g0", type=Path, required=True)
    ap.add_argument("--hypervolume", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--crosswalk", type=Path, required=True)
    ap.add_argument("--no-fetch", action="store_true")
    args = ap.parse_args()
    run_analysis(args.g0, args.hypervolume, args.output,
                 args.crosswalk, fetch=not args.no_fetch)


if __name__ == "__main__":
    main()
