"""Auditable full-78 original GIFT Direct trait-category vs nominal Rao-Q crossgate.

Input: ONLY the 2 anonymous original GIFT Actions JSON/ZIP reports. Preserve
all 4,505 original botanical island floras and full 13 x 3 x 2 cells. The
valid source stage2 and source Rao Q models must share exactly the same focal
plant-carrier cohort and proper scoring metrics on each cell.

Output is an inspectable anonymous 78-panel table with double geographical
grouping and double proper-loss sign gates, not a pollination mechanism claim.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
import zipfile

REGIONS = ("northern_midlatitude", "tropical", "southern_extratropical")
SCOPES = ("all_recorded", "original_native")
GROUPS = ("physical_block", "first_source_group")
SOURCE_MODEL = "environment_source_nearest"
REFERENCES = ("regional_blind", "external_blind", "regional_genus_prior")
Q_ESTIMATORS = ("plugin", "unbiased")
PROTECTED_KEYS = frozenset({
    "island_id", "accepted_species", "entity_ID", "genus_name", "source_genus",
    "spatial_block", "first_source", "island_latitude", "island_longitude",
    "latitude", "longitude",
})
EXPECTED_PARENT = 4505
EXPECTED_ENDPOINTS = 13
EXPECTED_CELLS = 78


def read_source(path: Path, expected_sha256: str | None = None) -> dict:
    if expected_sha256 and hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError("Source-original artifact ZIP digest drift")
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as archive:
            entries = [s for s in archive.namelist() if s.endswith(".json")]
            if len(entries) != 1:
                raise ValueError("Expected exactly one anonymous original GIFT JSON")
            obj = json.loads(archive.read(entries[0]))
    else:
        obj = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError("Expected JSON root object")
    check_privacy(obj)
    return obj


def check_privacy(obj):
    if isinstance(obj, dict):
        if PROTECTED_KEYS.intersection(obj):
            raise ValueError("Original plant/island/source identity key exported")
        for child in obj.values():
            check_privacy(child)
    elif isinstance(obj, list):
        for child in obj:
            check_privacy(child)


def indexed_panels(data: dict) -> tuple[set[str], dict[tuple, dict]]:
    if data.get("original_parent_GIFT_botanical_floras") != EXPECTED_PARENT:
        raise ValueError("Full 4505 botanical-island parent changed")
    eps = data.get("frozen_pilot_endpoints")
    if (not isinstance(eps, list) or len(eps) != EXPECTED_ENDPOINTS
            or len(set(eps)) != EXPECTED_ENDPOINTS):
        raise ValueError("Full 13 Direct original endpoints missing")
    cells = data.get("regional_source_heldout_results")
    if not isinstance(cells, list) or len(cells) != EXPECTED_CELLS:
        raise ValueError("Missing original 13×3×2 source panels")
    result = {}
    for r in cells:
        region, scope, trait = r.get("region"), r.get("scope"), r.get("trait_endpoint")
        if (region not in REGIONS or scope not in SCOPES or trait not in eps
                or r.get("original_parent_botanical_floras") != EXPECTED_PARENT
                or r.get("original_full_accepted_name_group_ban_on_both_models") is not True):
            raise ValueError("Unexpected original botanical source cohort/embargo")
        key = region, scope, trait
        if key in result:
            raise ValueError("Duplicate original endpoint×status×region panel")
        result[key] = r
    if len(result) != EXPECTED_CELLS:
        raise ValueError("Original 78 independent target cohort keys missing")
    return set(eps), result


def four_way_sign(row: dict, field: str) -> str:
    if not row.get("support_both_independent_source_groups"):
        return "unsupported"
    intervals = []
    for grp in GROUPS:
        summary = row.get(grp, {})
        if not summary.get("supported_for_inference"):
            raise ValueError("Contradictory source geographic group support")
        for loss in ("Brier", "logloss") if not field.endswith("_abs_error") else ("",):
            key = field+"_"+loss if loss else field
            val = summary.get(key, {}).get("bootstrap_95")
            if (not isinstance(val, list) or len(val) != 2
                    or not all(isinstance(v, (int, float)) and math.isfinite(v)
                               for v in val)
                    or val[0] > val[1]):
                raise ValueError("Incomplete original source paired confidence intervals "+key)
            intervals.append(val)
    if all(lo > 0 for lo, hi in intervals):
        return "positive"
    if all(hi < 0 for lo, hi in intervals):
        return "negative"
    return "ambiguous"


def rao_sign(row: dict, estimator: str) -> str:
    if estimator not in Q_ESTIMATORS:
        raise ValueError("Unregistered Gini–Simpson source estimator")
    return four_way_sign(
        row, f"{SOURCE_MODEL}_gain_vs_regional_blind_categorical_simpson_{estimator}_abs_error")


def evaluate(category_data: dict, rao_data: dict, *, tol: float = 1e-10) -> dict:
    endpoints, first = indexed_panels(category_data)
    rao_endpoints, second = indexed_panels(rao_data)
    if endpoints != rao_endpoints or set(first) != set(second):
        raise ValueError("Original category and diversity endpoints are not identical")
    if rao_data.get("categorical_simpson_secondary_not_multivariate_hypervolume") is not True:
        raise ValueError("Expected secondary nominal Rao-Q, not continuous trait volume")
    all_rows = []
    counts = Counter()
    for key in sorted(first):
        c, q = first[key], second[key]
        if (c["n_original_target_islands_scored"] != q["n_original_target_islands_scored"]
                or c["n_stage1_original_region_geography_eligible"] !=
                   q["n_stage1_original_region_geography_eligible"]
                or c["n_original_region_flora"] != q["n_original_region_flora"]
                or c["support_both_independent_source_groups"] !=
                   q["support_both_independent_source_groups"]):
            raise ValueError("Original scored Direct botanical carrier cohort or source support drift")
        for grp in GROUPS:
            for reference in REFERENCES:
                for loss in ("Brier", "logloss"):
                    col = f"{SOURCE_MODEL}_gain_vs_{reference}_{loss}"
                    left = c.get(grp, {}).get(col)
                    right = q.get(grp, {}).get(col)
                    if (left is None) != (right is None):
                        raise ValueError("Inconsistent source-model proper score presence")
                    if left is not None:
                        for val_key in ("group_equal_mean",):
                            if not math.isclose(left[val_key],right[val_key],
                                                rel_tol=0.0,abs_tol=tol):
                                raise ValueError("Original proper-score model changed between source runs")
        score = {}
        for baseline in REFERENCES:
            label = four_way_sign(c, f"{SOURCE_MODEL}_gain_vs_{baseline}")
            score[baseline] = label
        dplugin = rao_sign(q, "plugin")
        dunbiased = rao_sign(q, "unbiased")
        for estimator, outcome in (("plugin",dplugin),("unbiased",dunbiased)):
            flag = q.get(f"categorical_simpson_{estimator}_abs_error_improves_both_groups")
            if flag is not (outcome == "positive"):
                raise ValueError("Source Q stored positive flag contradicts signed grouped interval")
        if c.get("positive_vs_local_blind_both_scores_and_groups") is not (
                score["regional_blind"] == "positive"):
            raise ValueError("Source category local-blind positive flag inconsistent")
        if c.get("positive_vs_external_blind_both_scores_and_groups") is not (
                score["external_blind"] == "positive"):
            raise ValueError("Source category external-blind positive flag inconsistent")
        category_gain = score["regional_blind"] == "positive"
        diversity_gain = dunbiased == "positive"
        triple = all(score[base] == "positive" for base in REFERENCES)
        out = dict(
            region=key[0], scope=key[1], trait_endpoint=key[2],
            n_original_islands=c["n_original_target_islands_scored"],
            both_source_groups_qualified=c["support_both_independent_source_groups"],
            category_vs_local=score["regional_blind"],
            category_vs_external=score["external_blind"],
            category_vs_pooled_genus=score["regional_genus_prior"],
            nominal_RaoQ_plugin_vs_local=dplugin,
            nominal_RaoQ_unbiased_vs_local=dunbiased,
            category_and_unbiased_Q_both_positive=category_gain and diversity_gain,
            strict_3_category_baseline_and_unbiased_Q=triple and diversity_gain,
        )
        all_rows.append(out)
        counts[(key[1], "category_gain", category_gain)] += 1
        counts[(key[1], "diversity_gain", diversity_gain)] += 1
        counts[(key[1], "joint_support", category_gain and diversity_gain)] += 1
        counts[(key[1], "triple_and_diversity", triple and diversity_gain)] += 1
    return dict(
        original_GIFT_botanical_parent=EXPECTED_PARENT,
        n_genuine_Direct_endpoints=EXPECTED_ENDPOINTS,
        n_original_source_panels=len(all_rows),
        summary_counts={"|".join(map(str,k)):v for k,v in sorted(counts.items())},
        panels=all_rows,
        limitations=(
            "Exploratory unadjusted 95% interval gates across dependent botanical "
            "source traits, source-trained species-wide phenotype categories; "
            "nominal categorical Q not continuous multivariate morphology volume; "
            "no observed within-species island trait evolution or independent "
            "mainland ancestry/pollinator-guild/fitness transition."),
    )


def main():
    p=argparse.ArgumentParser()
    p.add_argument("category_zip_or_json", type=Path)
    p.add_argument("raoq_zip_or_json", type=Path)
    p.add_argument("--category-sha256")
    p.add_argument("--raoq-sha256")
    p.add_argument("--output-json", type=Path)
    p.add_argument("--output-csv", type=Path)
    opts=p.parse_args()
    cat=read_source(opts.category_zip_or_json,opts.category_sha256)
    q=read_source(opts.raoq_zip_or_json,opts.raoq_sha256)
    result=evaluate(cat,q)
    if opts.output_json:
        opts.output_json.parent.mkdir(parents=True,exist_ok=True)
        opts.output_json.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    if opts.output_csv:
        opts.output_csv.parent.mkdir(parents=True,exist_ok=True)
        with opts.output_csv.open("w",newline="",encoding="utf-8") as fh:
            writer=csv.DictWriter(fh,fieldnames=list(result["panels"][0]))
            writer.writeheader()
            writer.writerows(result["panels"])
    print(json.dumps(dict(parent=result["original_GIFT_botanical_parent"],
                          cells=result["n_original_source_panels"],
                          summary=result["summary_counts"]),indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
