"""Independently audit source-anonymous GIFT stage2 proper-loss summaries.

Accept GitHub Actions ZIP with a single JSON or raw JSON. Never read protected
plant/island name ledgers. Output complete signed regional decisions rather than
only successes; do not relabel empirical prediction as ecological selection.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import zipfile

REGIONS = {"northern_midlatitude", "tropical", "southern_extratropical"}
SCOPES = {"all_recorded", "original_native"}
CLUSTERS = ("physical_block", "first_source_group")
PRIVATE_KEYS = (
    "island_id", "accepted_species", "entity_ID", "genus_name", "source_genus",
    "spatial_block", "first_source", "island_latitude", "island_longitude",
    "latitude", "longitude",
)


def load_anon(path: Path) -> dict:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            items = [name for name in archive.namelist()
                     if name.lower().endswith(".json")]
            if len(items) != 1:
                raise ValueError("Expected exactly one anonymous JSON in ZIP")
            content = archive.read(items[0]).decode("utf-8")
    else:
        content = path.read_text(encoding="utf-8")
    for name in PRIVATE_KEYS:
        if '"' + name + '":' in content:
            raise ValueError("Private original source identifier key in JSON")
    data = json.loads(content)
    if not isinstance(data, dict):
        raise ValueError("Original anonymous source root must be an object")
    return data


def direction(row: dict, *, model: str, baseline: str) -> str:
    """95% exploratory interval must agree on BOTH scores and BOTH clusters."""
    if not all(row.get(k, {}).get("supported_for_inference") is True
               for k in CLUSTERS):
        return "insufficient"
    lows, highs = [], []
    for cluster in CLUSTERS:
        for loss in ("Brier", "logloss"):
            col = f"{model}_gain_vs_{baseline}_{loss}"
            item = row[cluster].get(col, {})
            ci = item.get("bootstrap_95")
            if (not isinstance(ci, list) or len(ci) != 2
                    or any(not isinstance(x, (int, float)) for x in ci)
                    or ci[0] > ci[1]):
                raise ValueError("Missing or invalid cluster confidence interval " + col)
            lows.append(ci[0])
            highs.append(ci[1])
    if min(lows) > 0:
        return "positive"
    if max(highs) < 0:
        return "negative"
    return "ambiguous"


def verify(data: dict) -> dict:
    if data.get("original_parent_GIFT_botanical_floras") != 4505:
        raise ValueError("Original GIFT population not 4505")
    eps = data.get("frozen_pilot_endpoints", [])
    if not isinstance(eps, list) or len(eps) not in (2, 13) or len(set(eps)) != len(eps):
        raise ValueError("Expected two endpoints or full thirteen")
    rows = data.get("regional_source_heldout_results", [])
    if not isinstance(rows, list) or len(rows) != len(eps) * 3 * 2:
        raise ValueError("Missing original region×endpoint×status panels")
    identity = [(r["region"], r["scope"], r["trait_endpoint"]) for r in rows]
    if (len(set(identity)) != len(identity)
            or {x[0] for x in identity} != REGIONS
            or {x[1] for x in identity} != SCOPES
            or {x[2] for x in identity} != set(eps)):
        raise ValueError("Repeated or unexpected region/endpoints/status identity")
    result, counts = [], Counter()
    for row in rows:
        if (row.get("original_parent_botanical_floras") != 4505
                or row.get("original_full_accepted_name_group_ban_on_both_models") is not True):
            raise ValueError("Original full-botanical name exclusion/parent missing")
        if row["n_original_target_islands_scored"] < 0:
            raise ValueError("Negative original source cohort")
        d = {
            "region": row["region"],
            "scope": row["scope"],
            "endpoint": row["trait_endpoint"],
            "n_islands": row["n_original_target_islands_scored"],
        }
        checks = [
            ("local_blind", "environment_source_nearest", "regional_blind",
             "positive_vs_local_blind_both_scores_and_groups"),
            ("external_blind", "environment_source_nearest", "external_blind",
             "positive_vs_external_blind_both_scores_and_groups"),
            ("pooled_genus", "environment_source_nearest", "regional_genus_prior",
             None),
        ]
        if "expanded_positive_over_constrained_both_scores_and_groups" in row:
            checks += [
                ("expanded_local_blind", "expanded_environment_source_nearest",
                 "regional_blind", "expanded_positive_vs_local_blind_both_scores_and_groups"),
                ("expanded_vs_constrained", "expanded_environment_source_nearest",
                 "constrained", "expanded_positive_over_constrained_both_scores_and_groups"),
            ]
        for key, model, baseline, expected_flag in checks:
            observed = direction(row, model=model, baseline=baseline)
            if expected_flag is not None and row.get(expected_flag) is not (observed == "positive"):
                raise ValueError("Stored source positive gate disagrees with independent audit: "
                                 + expected_flag)
            d[key] = observed
            counts[(row["scope"], key, observed)] += 1
        result.append(d)
    columns = sorted({k for row in result for k in row}
                     - {"region", "scope", "endpoint", "n_islands"})
    common = {}
    for scope in sorted(SCOPES):
        for col in columns:
            common[f"{scope}:{col}"] = sorted(
                ep for ep in eps
                if all(any(r["scope"] == scope and r["region"] == reg
                           and r["endpoint"] == ep and r.get(col) == "positive"
                           for r in result) for reg in REGIONS))
    return {
        "source_status": data.get("status"),
        "original_parent": 4505,
        "n_endpoints": len(eps),
        "n_panels": len(rows),
        "status_counts": {"|".join(key): val for key, val in sorted(counts.items())},
        "three_region_supported_endpoints": common,
        "panels": sorted(result, key=lambda x: (x["scope"], x["endpoint"], x["region"])),
        "limits": "Predictive botanical Direct category score, exploratory group intervals; not guild-mediated adaptation or ancestral mainland state.",
    }


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("source_anonymous_zip_or_json", type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = verify(load_anon(args.source_anonymous_zip_or_json))
    content = json.dumps(result, sort_keys=True, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(content, encoding="utf-8")
    print(content)


if __name__ == "__main__":
    main()
