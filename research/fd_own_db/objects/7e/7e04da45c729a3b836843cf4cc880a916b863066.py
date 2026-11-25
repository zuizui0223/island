"""External GBIF-effort robustness audit for geological-origin floral repertoires.

A pinned, earlier island-v2 GBIF collection generated PUBLIC *Tracheophyta*
island-level effort metadata. These are not a calibrated detection history,
not the same sampling frame as the later reviewed angiosperm floral checklist,
and never justify treating an unrecorded genus as truly absent.

This explicitly POST-HOC screening only asks: do the previously reported
geology-matched profile composition differences persist after excluding
island pairs whose recording opportunity metadata are sparse/unequal?

Published artifacts contain ONLY aggregate metrics, not island-level
geological mappings, taxon/genus names, or exact matched pair IDs.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import urllib.request

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_carriers, load_profiles
from .pipeline import REF_INPUT
from .repertoire import _paired_swap_statistic, _upper_tail, match_origin_islands
from .g1_fixed_margins import _bh

ISLAND_SOURCE_SHA = "3bfa818ffbc7a0dea4319933c4b0a0eb41cfabb5"
EFFORT_RELATIVE_PATH = "data/v2/staging/gbif/collected/island_observation_effort.csv"
EFFORT_URL = (
    "https://raw.githubusercontent.com/zuizui0223/island/"
    + ISLAND_SOURCE_SHA + "/" + EFFORT_RELATIVE_PATH
)
METRICS = (
    "n_unique_gbif_ids", "n_datasets", "n_species", "n_preserved_specimen",
    "n_human_observation", "n_uncertainty_reported",
    "median_coordinate_uncertainty_m",
    "n_coordinate_uncertainty_gt_1000m",
)
CHECKS = (
    "all_original_pairs",
    "both_effort_available",
    "min30records_2datasets",
    "min30_balanced4x_uncertainty10km",
)
ENDPOINTS = ("exact_profile_L2", "categorical_similarity_kernel")
SCOPES = ("all", "direct")


def download_pinned_effort(destination: Path) -> str:
    """Fetch only pinned public aggregated effort metadata, with byte limit."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(
        EFFORT_URL, headers={"User-Agent": "fd-effort-sensitivity/1.0"})
    sha = hashlib.sha256()
    size = 0
    with urllib.request.urlopen(req, timeout=90) as source, destination.open("wb") as target:
        while data := source.read(1 << 20):
            size += len(data)
            if size > 12_000_000:
                raise ValueError("Unexpectedly large GBIF aggregate effort source")
            sha.update(data)
            target.write(data)
    if size < 50:
        raise ValueError("Empty source effort input")
    return sha.hexdigest()


def read_effort(path: Path) -> dict[str, dict]:
    df = pd.read_csv(path, dtype={"island_id": str})
    required = {"island_id", *METRICS}
    if not required.issubset(df.columns):
        raise ValueError(f"Missing effort fields: {sorted(required - set(df.columns))}")
    if df.island_id.isna().any() or df.island_id.duplicated().any():
        raise ValueError("Invalid or repeated effort island identifiers")
    for name in METRICS:
        df[name] = pd.to_numeric(df[name], errors="coerce")
    for name in ("n_unique_gbif_ids", "n_datasets", "n_species"):
        if (df[name].dropna() < 0).any():
            raise ValueError("Negative effort count")
    return df.set_index("island_id").to_dict("index")


def pair_effort_label(pair: dict, effort: dict[str, dict]) -> tuple[str, dict]:
    """Classify a pair under predetermined nested screening thresholds."""
    c, o = effort.get(pair["continental_id"]), effort.get(pair["oceanic_id"])
    if c is None or o is None:
        return "no_effort_metadata", {"has_both": False}
    def numeric(row: dict, key: str) -> float:
        value = row.get(key)
        return float(value) if value is not None else float("nan")

    record_c, record_o = numeric(c, "n_unique_gbif_ids"), numeric(o, "n_unique_gbif_ids")
    datasets_c, datasets_o = numeric(c, "n_datasets"), numeric(o, "n_datasets")
    uncertainty_c = numeric(c, "median_coordinate_uncertainty_m")
    uncertainty_o = numeric(o, "median_coordinate_uncertainty_m")
    available = bool(np.isfinite([record_c, record_o, datasets_c, datasets_o]).all())
    moderate = bool(
        available and min(record_c, record_o) >= 30
        and min(datasets_c, datasets_o) >= 2)
    balanced = bool(
        moderate
        and abs(math.log1p(record_c) - math.log1p(record_o)) <= math.log(4)
        and np.isfinite([uncertainty_c, uncertainty_o]).all()
        and max(uncertainty_c, uncertainty_o) <= 10_000
    )
    return "available", {
        "has_both": available,
        "moderate": moderate,
        "balanced": balanced,
        "log_record_ratio_o_minus_c": (
            float(math.log1p(record_o) - math.log1p(record_c))
            if available else None),
        "log_dataset_ratio_o_minus_c": (
            float(math.log1p(datasets_o) - math.log1p(datasets_c))
            if available else None),
    }


def eligible_pairs(pairs: list[dict], effort: dict[str, dict], selection: str
                   ) -> tuple[list[dict], dict]:
    if selection not in CHECKS:
        raise ValueError("Unrecognized effort screening check")
    included = []
    records = []
    for pair in pairs:
        label, diag = pair_effort_label(pair, effort)
        records.append(diag)
        yes = (
            selection == "all_original_pairs"
            or (selection == "both_effort_available" and diag.get("has_both", False))
            or (selection == "min30records_2datasets" and diag.get("moderate", False))
            or (selection == "min30_balanced4x_uncertainty10km" and diag.get("balanced", False))
        )
        if yes:
            included.append(pair)
    ratios = [r["log_record_ratio_o_minus_c"]
              for r in records if r.get("log_record_ratio_o_minus_c") is not None]
    audit = {
        "selection": selection,
        "n_original_pairs": len(pairs),
        "n_pairs_retained": len(included),
        "n_effort_metadata_both": sum(bool(r.get("has_both")) for r in records),
        "n_minimum_effort_both": sum(bool(r.get("moderate")) for r in records),
        "n_balanced_strict_both": sum(bool(r.get("balanced")) for r in records),
        "median_log_effort_o_minus_c": float(np.median(ratios)) if ratios else None,
        "n_northern_midlatitude": sum(p["region"] == "northern_midlatitude" for p in included),
        "n_tropical": sum(p["region"] == "tropical" for p in included),
    }
    return included, audit


def _repertoire_for_island(carriers: dict, identifier: str,
                           space: ProfileSpace) -> np.ndarray:
    if identifier not in carriers:
        raise ValueError("Matched island absent from public floral carrier aggregate")
    profile_ids = [space.lookup[p] for p in carriers[identifier]]
    if not profile_ids:
        raise ValueError("Empty floral profile support")
    v = np.zeros(len(space.ids), dtype=float)
    v[profile_ids] = 1.0 / len(profile_ids)
    return v


def score_subset(pairs: list[dict], carriers: dict, space: ProfileSpace,
                 *, draws: int, seed: int) -> list[dict]:
    if len(pairs) < 8:
        return []
    if draws < 999:
        raise ValueError("At least 999 paired swaps required")
    v = {}
    for p in pairs:
        for key in ("continental_id", "oceanic_id"):
            i = p[key]
            if i not in v:
                v[i] = _repertoire_for_island(carriers, i, space)
    delta = np.stack([
        v[p["oceanic_id"]] - v[p["continental_id"]] for p in pairs
    ])
    rng = np.random.default_rng(seed)
    signs = rng.choice((-1.0, 1.0), size=(draws, len(pairs)))
    result = []
    for endpoint, kernel in (
        ("exact_profile_L2", np.eye(len(space.ids))),
        ("categorical_similarity_kernel", 1 - space.matrix),
    ):
        observed, reference = _paired_swap_statistic(delta, kernel, signs)
        result.append({
            "endpoint": endpoint, "n_pairs": len(pairs),
            "n_islands": 2 * len(pairs),
            "observed_statistic": observed,
            "null_median": float(np.median(reference)),
            "null_q025": float(np.quantile(reference, .025)),
            "null_q975": float(np.quantile(reference, .975)),
            "p_swap_upper": _upper_tail(observed, reference),
        })
    return result


def run(input_root: Path, g0_root: Path, crosswalk_path: Path,
        effort_path: Path, output: Path, *, draws: int = 9999,
        seed: int = 20261108) -> dict:
    input_root, g0_root, crosswalk_path, effort_path, output = map(
        Path, (input_root, g0_root, crosswalk_path, effort_path, output))
    output.mkdir(parents=True, exist_ok=True)
    from .core import csv_rows
    g0_all = list(csv_rows(g0_root / "all_islands.csv"))
    g0_direct = list(csv_rows(g0_root / "direct_islands.csv"))
    crosswalk = {r["island_id"]: r for r in csv_rows(crosswalk_path)}
    matched, original_audit = match_origin_islands(g0_all, g0_direct, crosswalk)
    effort = read_effort(effort_path)
    space = ProfileSpace(load_profiles(input_root / REF_INPUT))
    results, groups = [], []
    for selection in CHECKS:
        pairs, audit = eligible_pairs(matched, effort, selection)
        groups.append(audit)
        for scope in SCOPES:
            carriers = load_carriers(
                input_root / "floral_space_20261106" / scope /
                "island_profile_carriers.csv.gz")
            scored = score_subset(pairs, carriers, space, draws=draws,
                                  seed=seed + CHECKS.index(selection) * 17)
            for row in scored:
                results.append({
                    "selection": selection, "scope": scope, **row,
                    "status": "posthoc_effort_screen_not_independent_confirmation",
                })
            if not scored:
                for endpoint in ENDPOINTS:
                    results.append({
                        "selection": selection, "scope": scope,
                        "endpoint": endpoint,
                        "status": "insufficient_pairs_n_less_than_8",
                        "n_pairs": len(pairs),
                    })
    expected_tests = len(CHECKS) * len(SCOPES) * len(ENDPOINTS)
    if len(results) != expected_tests:
        raise ValueError("Unexpected number of declared effort screen test slots")
    probs = [float(r["p_swap_upper"]) if "p_swap_upper" in r else 1.
             for r in results]
    bh = _bh(probs)
    H = sum(1/i for i in range(1, expected_tests+1))
    for j, row in enumerate(results):
        if "p_swap_upper" in row:
            row["q_bh_posthoc16"] = float(bh[j])
            row["q_by_posthoc16"] = float(min(1., bh[j] * H))
    pd.DataFrame(groups).to_csv(
        output / "effort_pair_support_summary.csv", index=False)
    pd.DataFrame(results).to_csv(
        output / "effort_screen_repertoire_tests.csv", index=False)
    # In the public G0 all cohort, proportion of geological islands with
    # this EARLIER effort table (not a completeness estimate).
    checks = []
    for row in g0_all:
        geo = crosswalk.get(row["island_id"], {})
        if geo.get("geological_origin") not in ("continental", "oceanic"):
            continue
        item = effort.get(row["island_id"])
        checks.append({
            "origin": geo["geological_origin"],
            "has_earlier_v2_effort": bool(item),
            "n_gbif_records": (
                item.get("n_unique_gbif_ids") if item else np.nan),
            "n_datasets": item.get("n_datasets") if item else np.nan,
        })
    table = pd.DataFrame(checks)
    coverage = []
    for origin, group in table.groupby("origin"):
        available = group.loc[group.has_earlier_v2_effort]
        coverage.append({
            "origin": origin,
            "n_flora_support_islands": len(group),
            "n_with_prior_effort_metadata": len(available),
            "effort_metadata_coverage_fraction": len(available) / len(group),
            "median_gbif_unique_ids": (
                float(available.n_gbif_records.median()) if len(available) else None),
            "median_gbif_datasets": (
                float(available.n_datasets.median()) if len(available) else None),
        })
    pd.DataFrame(coverage).to_csv(
        output / "effort_coverage_by_origin.csv", index=False)
    source_hash = hashlib.sha256(effort_path.read_bytes()).hexdigest()
    manifest = {
        "schema": "fd_posthoc_earlier_gbif_effort_screen_v1",
        "status": "descriptive_effort_screen_not_occupancy_or_detection_model",
        "source_island_commit": ISLAND_SOURCE_SHA,
        "source_public_effort_url": EFFORT_URL,
        "source_effort_sha256": source_hash,
        "sampled_previous_GBIF_scope": "Tracheophyta GBIF collection",
        "floral_outcome_scope": "reviewed complete-colour/form angiosperm taxon-name profiles",
        "n_original_geology_matched_pairs": len(matched),
        "n_prior_effort_islands": len(effort),
        "screens": list(CHECKS),
        "n_label_swap_draws": draws,
        "posthoc_declared_family_size": expected_tests,
        "posthoc_multiplicity": "BH/BY over 16 correlated tests, untestable slots assigned p=1",
        "n_matched_p_value_tests_computed": sum(
            "p_swap_upper" in r for r in results),
        "nonclaims": [
            "No source-resolved true absence can be inferred from an unrecorded taxon",
            "Effort metadata from an earlier campaign is NOT aligned per-taxon detection history",
            "Two nested All/Direct scopes are not independent confirmations",
            "Different post-hoc quality subsamples are correlated and selection-sensitive",
            "No calibrated arrival or establishment probability was estimated",
        ],
        "privacy": "published aggregate only; no island-by-origin identities or restricted names",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "status": manifest["status"],
        "support": groups,
        "coverage": coverage,
        "comparisons": [
            {"subset": r["selection"], "scope": r["scope"],
             "endpoint": r.get("endpoint"), "n_pairs": r.get("n_pairs"),
             "p": r.get("p_swap_upper")}
            for r in results
        ],
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for x in ("input", "g0", "crosswalk", "output"):
        parser.add_argument("--" + x, type=Path, required=True)
    parser.add_argument("--effort", type=Path)
    parser.add_argument("--draws", type=int, default=9999)
    args = parser.parse_args()
    effort = args.effort or args.output.parent / "pinned_prior_gbif_effort.csv"
    if args.effort is None:
        download_pinned_effort(effort)
    run(args.input, args.g0, args.crosswalk, effort, args.output,
        draws=args.draws)


if __name__ == "__main__":
    main()
