"""Independent privacy/score audit of original GIFT same-species floral×reproduction models.

Reads ONLY original-anonymous GitHub Actions JSON/ZIP artifacts. Never accepts
botanical accepted-species rosters or island/source identification tables.
"""
from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path
import zipfile

REGIONS = {"northern_midlatitude", "tropical", "southern_extratropical"}
SCOPES = {"all_recorded", "original_native"}
ARCH = ("floral_form", "floral_symmetry", "tube_depth_class")
ASSURANCE = ("self_incompatibility", "autonomous_selfing_capacity")
ENDPOINTS = {a + "__" + b for a in ARCH for b in ASSURANCE}
GROUPS = ("physical_block", "first_source_group")
PRIVATE_KEYS = {"accepted_species", "island_id", "entity_ID", "genus",
                "genus_name", "latitude", "longitude", "spatial_block",
                "first_source", "source_genus", "island_latitude",
                "island_longitude"}


def check_private_keys(obj):
    """Private source ids are prohibited even in deeply nested dict payloads."""
    if isinstance(obj, dict):
        if PRIVATE_KEYS & set(obj):
            raise ValueError("Original protected source identity key exported")
        for v in obj.values():
            check_private_keys(v)
    elif isinstance(obj, list):
        for v in obj:
            check_private_keys(v)


def read_anonymous(path: Path) -> dict:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            names = [name for name in archive.namelist()
                     if name.endswith(".json") and not name.startswith("__MACOSX")]
            if len(names) != 1:
                raise ValueError("Exactly one anonymous JSON expected in source ZIP")
            content = archive.read(names[0])
    else:
        content = path.read_bytes()
    data = json.loads(content)
    if not isinstance(data, dict):
        raise ValueError("Expected an original source JSON object")
    check_private_keys(data)
    return data


def sign(row: dict, metric: str) -> str:
    """Require BOTH proper scores under BOTH independently chosen groupings."""
    if not row.get("support_both_independent_source_groups"):
        return "unsupported"
    bounds = []
    for grp in GROUPS:
        frame = row.get(grp, {})
        if frame.get("supported_for_inference") is not True:
            raise ValueError("Original grouped support flag contradicts source panel")
        for loss in ("Brier", "logloss"):
            item = frame.get(metric + "_" + loss, {})
            interval = item.get("bootstrap_95")
            if (not isinstance(interval, list) or len(interval) != 2
                    or not all(isinstance(x, (int, float)) and math.isfinite(x)
                               for x in interval)
                    or interval[0] > interval[1]):
                raise ValueError("Original paired cluster interval missing or malformed")
            bounds.append(interval)
    if all(lo > 0 for lo, hi in bounds):
        return "positive"
    if all(hi < 0 for lo, hi in bounds):
        return "negative"
    return "ambiguous"


def audit_original_joint(data: dict) -> dict:
    if data.get("original_parent_GIFT_botanical_floras") != 4505:
        raise ValueError("Original 4505-island botanical parent changed")
    endpoints = data.get("frozen_pilot_endpoints")
    if not isinstance(endpoints, list) or set(endpoints) != ENDPOINTS or len(endpoints) != 6:
        raise ValueError("Original six same-species morphology×reproduction endpoints changed")
    rows = data.get("regional_source_heldout_results")
    if not isinstance(rows, list) or len(rows) != 36:
        raise ValueError("Missing original 6×3×2 source outcome panels")
    seen = set()
    outcome = []
    counts = Counter()
    conditional = bool(data.get("genus_conditional_independence_retains_between_genus_compositional_covariance"))
    for row in rows:
        if (row.get("original_parent_botanical_floras") != 4505
                or row.get("original_full_accepted_name_group_ban_on_both_models") is not True
                or row.get("no_pollinator_guild_or_within_species_island_adaptation_measured") is not True):
            raise ValueError("Original source protection and noncausality controls absent")
        region, scope, trait = row["region"], row["scope"], row["trait_endpoint"]
        key = region, scope, trait
        if key in seen or region not in REGIONS or scope not in SCOPES or trait not in ENDPOINTS:
            raise ValueError("Duplicate or unregistered source region×scope×true pair")
        seen.add(key)
        n = row.get("n_original_target_islands_scored", -1)
        if not isinstance(n, int) or n < 0:
            raise ValueError("Invalid original GIFT scored island count")
        metrics = {
            "joint_vs_global_independence": (
                "predicted_joint_gain_vs_own_independence",
                "external_genus_joint_dependency_over_own_marginals_supported"),
            "local_prior_joint_vs_independence": (
                "regional_blind_joint_gain_vs_own_independence",
                "regional_blind_joint_dependency_over_own_marginals_supported"),
        }
        if conditional:
            metrics["joint_vs_within_genus_independence"] = (
                "predicted_joint_gain_vs_genus_conditional_independence",
                "external_genus_joint_dependency_over_genus_conditional_marginals_supported")
        doc = {"region": region, "scope": scope, "trait": trait,
               "n_original_islands": n}
        for label, (stem, expected_field) in metrics.items():
            v = sign(row, stem)
            if row.get(expected_field) is not (v == "positive"):
                raise ValueError("Source stored positive flag disagrees with proper-loss gate: "+stem)
            doc[label] = v
            counts[(scope, region, label, v)] += 1
        outcome.append(doc)
    if len(seen) != 36:
        raise ValueError("Original panels incomplete")
    return {
        "status": data.get("status"),
        "original_parent_islands": 4505,
        "n_real_same_species_joint_endpoints": 6,
        "n_source_panels": 36,
        "stronger_genus_conditional_null_included": conditional,
        "source_verdict_counts": {"|".join(k): v for k, v in sorted(counts.items())},
        "panels": sorted(outcome, key=lambda v: (v["scope"], v["region"], v["trait"])),
        "nonclaims": ("Not causal animal pollination, within-species island evolution, "
                      "phylogenetic selection or independent mainland ancestry. "
                      "Intervals exploratory and not multiplicity-adjusted."),
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("anonymous_original_zip_or_json", type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = audit_original_joint(read_anonymous(args.anonymous_original_zip_or_json))
    s = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(s, encoding="utf-8")
    print(s)


if __name__ == "__main__":
    main()
