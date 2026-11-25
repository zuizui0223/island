"""Cross-validated island floral repertoire imprint of named archipelagos.

Public-data analysis, deliberately NOT geological origin treatment assignment.
For each island with >=3 other same-archipelago, same-regime islands, predict
the held-out island's 496-profile recorded set distribution using:
  - island-blind SAME-REGIME, OTHER-ARCHIPELAGO mean (background);
  - 3 nearest SAME-ARCHIPELAGO island donors, smoothed 50% to background;
  - distance-caliper-matched 3 OTHER-ARCHIPELAGO donors, smoothed identically.

The target island is excluded from every donor pool. The archipelago group,
not the island, is the inference cluster. Missing/NA source archipelago labels
never form shared groups. This is predictive recorded phenotype composition,
not arrival, geological filtering, functional niche equivalence or a test of
independent taxonomic lineages. No protected taxon-level records are used.

Published outputs contain ONLY aggregate statistics and input hashes; source
archipelago names, island IDs, donor identities are never emitted.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .g1_fixed_margins import _bh
from .pipeline import REF_INPUT, _write_csv, file_sha256

SCOPES = ("all", "direct")
ENDPOINTS = ("exact_profile_L2", "categorical_similarity_kernel")
CONTRASTS = (
    "within_vs_region_outside",
    "within_vs_geography_only_nearest",
    "within_vs_distance_matched_outside",
)
REFERENCE_NEIGHBORS = 3
SOURCE_BLEND = .50
DISTANCE_RATIO_LIMIT = 3.0
WITHIN_RADIUS_KM = 1500.0
PERMUTATIONS = 9999
BOOTSTRAP = 4999
REGIME_ORDER = (
    "northern_high_latitude", "northern_midlatitude",
    "southern_extratropical", "tropical",
)
GROUP_FIELDS = (
    "scope", "analysis_regime", "n_named_archipelago_islands",
    "n_archipelago_regime_groups_at_least_4",
    "n_targets_with_source_predictions",
    "n_groups_with_source_predictions",
    "n_targets_with_distance_matched_controls",
    "n_groups_with_distance_matched_controls",
    "median_same_archipelago_third_donor_km",
    "median_matched_distance_gap_log",
    "n_three_nearest_geo_donors_all_same_archipelago",
    "fraction_three_nearest_geo_donors_all_same_archipelago",
    "mean_source_geo_donor_overlap_of_three",
)
RESULT_FIELDS = (
    "scope", "endpoint", "comparison", "n_islands",
    "n_archipelagos", "n_regional_subgroups", "mean_group_improvement",
    "group_median_improvement", "cluster_bootstrap_q025",
    "cluster_bootstrap_q975", "p_cluster_signflip_one_sided",
    "q_bh_12", "q_by_12",
    "n_randomizations", "n_group_bootstrap",
    "status",
)
PRED_FIELDS = (
    "scope", "endpoint", "n_islands", "n_archipelagos",
    "mean_region_outside_error", "mean_within_group_error",
    "mean_distance_matched_other_error", "n_distance_matched",
    "mean_geography_only_error",
)


def missing_archipelago(value) -> bool:
    return not isinstance(value, str) or value.strip().lower() in {
        "", "na", "n/a", "nan", "null", "none", "unknown", "unresolved", "-"
    }


def haversine_km(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    """Pairwise great-circle distances for actual WGS84 island centroids."""
    lat = np.radians(np.asarray(latitude, dtype=np.float64))
    lon = np.radians(np.asarray(longitude, dtype=np.float64))
    if (len(lat) != len(lon) or not len(lat) or
            np.any(~np.isfinite(lat)) or np.any(~np.isfinite(lon)) or
            np.any(np.abs(lat) > math.pi / 2) or
            np.any(np.abs(lon) > math.pi)):
        raise ValueError("Invalid geolocated island centroid population")
    angle = (
        np.sin(lat[:, None]) * np.sin(lat[None, :])
        + np.cos(lat[:, None]) * np.cos(lat[None, :])
        * np.cos(lon[:, None] - lon[None, :])
    )
    return 6371.0088 * np.arccos(np.clip(angle, -1., 1.))


def nearest_geographic_donors(
    target_index: int, distances: np.ndarray,
) -> np.ndarray:
    """Select exactly three nearest same-regime islands using location only.

    No archipelago identity or floral phenotype enters donor selection.
    Excludes the target even when another centroid has zero distance.
    """
    n = len(distances)
    if distances.shape != (n,n) or n < REFERENCE_NEIGHBORS + 1:
        raise ValueError("Too few geocoded reference islands")
    if not 0 <= target_index < n:
        raise ValueError("Missing target island index")
    other = np.flatnonzero(np.arange(n) != target_index)
    return other[np.argsort(distances[target_index, other])[:REFERENCE_NEIGHBORS]]


def donor_selection(
    target_index: int,
    group: np.ndarray,
    distances: np.ndarray,
    *, max_distance_ratio: float = DISTANCE_RATIO_LIMIT,
    within_radius_km: float = WITHIN_RADIUS_KM,
) -> tuple[np.ndarray, np.ndarray | None, float, float | None]:
    """Outcome-blind same-archipelago nearest three and outside controls.

    Other-archipelago donors are selected by a one-to-one Hungarian assignment
    to match their three target geographic distances to the within donors'
    three distances. All three log1p-distance gaps must be <= ln(3).
    The 1500-km upper bound prevents very broad source labels from counting
    as locally shared archipelago evidence in the strict sensitivity.
    """
    n = len(group)
    if not 0 <= target_index < n or distances.shape != (n, n):
        raise ValueError("Distances and target index disagree")
    same = np.flatnonzero(group == group[target_index])
    same = same[same != target_index]
    outside = np.flatnonzero(group != group[target_index])
    if len(same) < REFERENCE_NEIGHBORS or len(outside) < REFERENCE_NEIGHBORS:
        raise ValueError("Insufficient other-island donors")
    nearest = same[np.argsort(distances[target_index, same])[:REFERENCE_NEIGHBORS]]
    near_max = float(distances[target_index, nearest].max())
    if near_max > within_radius_km:
        return nearest, None, near_max, None
    target_logdist = np.log1p(distances[target_index, nearest])
    outside_logdist = np.log1p(distances[target_index, outside])
    cost = np.abs(target_logdist[:, None] - outside_logdist[None, :])
    row, col = linear_sum_assignment(cost)
    if len(row) != REFERENCE_NEIGHBORS:
        raise ValueError("Incomplete geographical donor matching")
    gap = float(cost[row, col].max())
    if gap > math.log(max_distance_ratio):
        return nearest, None, near_max, gap
    return nearest, outside[col], near_max, gap


def _kernel_score(observed: np.ndarray, predicted: np.ndarray,
                  kernel: np.ndarray) -> float:
    difference = observed - predicted
    value = float(difference @ kernel @ difference)
    if not np.isfinite(value):
        raise ValueError("Nonfinite profile-prediction squared error")
    return max(value, 0.)


def _cluster_test(archipelago_values: np.ndarray, *, permutations: int,
                  bootstrap: int, seed: int) -> dict:
    """One archipelago = one sign/one bootstrap unit.

    This dependence adjustment does not make archipelagos phylogenetically
    independent or correct shared donors across groups. It is a more honest
    unit than treating hundreds of islands as independent causal replicates.
    """
    values = np.asarray(archipelago_values, dtype=float)
    if values.ndim != 1 or len(values) < 8 or not np.isfinite(values).all():
        raise ValueError("At least eight finite independent source groups required")
    if permutations < 999 or bootstrap < 999:
        raise ValueError("Insufficient randomization or bootstrap draws")
    rng = np.random.default_rng(seed)
    observed = float(values.mean())
    signs = rng.choice((-1., 1.), size=(permutations, len(values)))
    shuffled = (signs * values[None, :]).mean(axis=1)
    p_one_sided = (1 + int(np.count_nonzero(
        shuffled >= observed - 1e-12 * max(1., abs(observed))
    ))) / (permutations + 1)
    indices = rng.integers(len(values), size=(bootstrap, len(values)))
    ci = np.quantile(values[indices].mean(axis=1), [.025, .975])
    return {
        "mean_group_improvement": observed,
        "group_median_improvement": float(np.median(values)),
        "cluster_bootstrap_q025": float(ci[0]),
        "cluster_bootstrap_q975": float(ci[1]),
        "p_cluster_signflip_one_sided": float(p_one_sided),
    }


def run(input_root: Path, g0_root: Path, crosswalk_path: Path,
        output: Path, *, seed: int = 20261109,
        permutations: int = PERMUTATIONS,
        bootstrap: int = BOOTSTRAP) -> dict:
    input_root, g0_root, crosswalk_path, output = map(
        Path, (input_root, g0_root, crosswalk_path, output))
    output.mkdir(parents=True, exist_ok=True)
    common = {r["island_id"] for r in csv_rows(g0_root / "all_islands.csv")}
    common &= {r["island_id"] for r in csv_rows(g0_root / "direct_islands.csv")}
    geog = {r["island_id"]: r for r in csv_rows(crosswalk_path)}
    chosen = []
    for island in sorted(common):
        g = geog.get(island)
        if g is None or missing_archipelago(g.get("weigelt_archipelago")):
            continue
        try:
            lat = float(g["island_latitude"])
            lon = float(g["island_longitude"])
        except (TypeError, ValueError, KeyError):
            continue
        if not np.isfinite([lat, lon]).all():
            continue
        chosen.append({
            "island_id": island,
            "archipelago": g["weigelt_archipelago"].strip(),
            "regime": g["analysis_regime"],
            "lat": lat,
            "lon": lon,
        })
    if len(chosen) < 20:
        raise ValueError("No globally useful source-archipelago inventory")
    space = ProfileSpace(load_profiles(input_root / REF_INPUT))
    kernels = {
        "exact_profile_L2": np.eye(len(space.ids)),
        "categorical_similarity_kernel": 1. - space.matrix,
    }
    # Every pair of All and Direct evaluations uses identical source/target
    # islands, source labels and (pre-score) geographic donor choices.
    by_region: dict[str, list[dict]] = {}
    for r in chosen:
        by_region.setdefault(r["regime"], []).append(r)
    carriers = {
        s: load_carriers(input_root / "floral_space_20261106" / s /
                         "island_profile_carriers.csv.gz")
        for s in SCOPES
    }
    target = {s: [] for s in SCOPES}
    group_audits = []
    global_audits = defaultdict(lambda: {
        "n_named_archipelago_islands": 0,
        "n_archipelago_regime_groups_at_least_4": 0,
        "n_targets_with_source_predictions": 0,
        "n_targets_with_distance_matched_controls": 0,
        "third_distances": [],
        "match_gaps": [],
        "n_geo_nearest_identical": 0,
        "shared_donor_counts": [],
        "main_groups": set(),
        "caliper_groups": set(),
    })
    for region in REGIME_ORDER:
        members = by_region.get(region, [])
        if not members:
            continue
        n = len(members)
        names = [x["island_id"] for x in members]
        groups = np.asarray([x["archipelago"] for x in members])
        counts = Counter(groups)
        dist = haversine_km(
            np.asarray([x["lat"] for x in members]),
            np.asarray([x["lon"] for x in members]),
        )
        support = np.asarray([
            counts[groups[i]] >= REFERENCE_NEIGHBORS + 1
            and n - counts[groups[i]] >= REFERENCE_NEIGHBORS
            for i in range(n)
        ], dtype=bool)
        audit = global_audits[region]
        audit["n_named_archipelago_islands"] = n
        audit["n_archipelago_regime_groups_at_least_4"] = sum(
            k >= REFERENCE_NEIGHBORS + 1 and n-k >= REFERENCE_NEIGHBORS
            for k in counts.values()
        )
        candidates = np.flatnonzero(support)
        if not len(candidates):
            continue
        # Common cohort and source-identity coverage from the SAME frozen
        # public profiles. A name carrier counts as one *profile type* here.
        data = {}
        for scope in SCOPES:
            matrix = np.zeros((n, len(space.ids)), dtype=float)
            for i, name in enumerate(names):
                raw = carriers[scope].get(name)
                if not raw:
                    raise ValueError("Eligible island absent from public profile carriers")
                ids = [space.lookup[p] for p in raw]
                matrix[i, ids] = 1. / len(ids)
            data[scope] = matrix
        for i in candidates:
            outside = np.flatnonzero(groups != groups[i])
            within, matched, maxwithin, gap = donor_selection(i, groups, dist)
            # Geographic-only negative control: nearest three in the SAME
            # originally declared regime, WITHOUT consulting archipelago
            # labels to select them. Compare identical three-donor and
            # 50%-reference mixture structures.
            nearest_geo = nearest_geographic_donors(i, dist)
            n_shared_geo = len(set(within) & set(nearest_geo))
            audit["n_geo_nearest_identical"] += int(n_shared_geo==REFERENCE_NEIGHBORS)
            audit["shared_donor_counts"].append(n_shared_geo)
            audit["n_targets_with_source_predictions"] += 1
            audit["main_groups"].add(groups[i])
            audit["third_distances"].append(maxwithin)
            if matched is not None:
                audit["n_targets_with_distance_matched_controls"] += 1
                audit["caliper_groups"].add(groups[i])
                if gap is not None:
                    audit["match_gaps"].append(gap)
            for scope in SCOPES:
                matrix = data[scope]
                y = matrix[i]
                bg = matrix[outside].mean(axis=0)  # Entire archipelago held out
                same = SOURCE_BLEND * matrix[within].mean(axis=0) + (1-SOURCE_BLEND)*bg
                other = (SOURCE_BLEND * matrix[matched].mean(axis=0)
                         + (1-SOURCE_BLEND)*bg) if matched is not None else None
                geo_only = (SOURCE_BLEND * matrix[nearest_geo].mean(axis=0)
                            + (1-SOURCE_BLEND)*bg)
                if not np.allclose([bg.sum(),same.sum()],1.,atol=1e-9):
                    raise ValueError("Non-normalized holdout predicted distributions")
                for endpoint, kernel in kernels.items():
                    baseline_error = _kernel_score(y, bg, kernel)
                    within_error = _kernel_score(y, same, kernel)
                    outside_error = (
                        _kernel_score(y, other, kernel) if other is not None
                        else None)
                    geo_error = _kernel_score(y, geo_only, kernel)
                    target[scope].append({
                        "archipelago": groups[i], "regime": region,
                        "endpoint": endpoint,
                        "within_vs_region_outside": baseline_error-within_error,
                        "within_vs_geography_only_nearest": geo_error-within_error,
                        "within_vs_distance_matched_outside": (
                            outside_error-within_error if outside_error is not None
                            else None),
                        "baseline_error": baseline_error,
                        "within_error": within_error,
                        "outside_error": outside_error,
                        "geo_error": geo_error,
                    })
    summary = []
    region_counts = []
    for region in REGIME_ORDER:
        a = global_audits.get(region)
        if not a:
            continue
        region_counts.append({
            "analysis_regime": region,
            "n_named_archipelago_islands": a["n_named_archipelago_islands"],
            "n_archipelago_regime_groups_at_least_4":
                a["n_archipelago_regime_groups_at_least_4"],
            "n_targets_with_source_predictions":
                a["n_targets_with_source_predictions"],
            "n_groups_with_source_predictions": len(a["main_groups"]),
            "n_targets_with_distance_matched_controls":
                a["n_targets_with_distance_matched_controls"],
            "n_groups_with_distance_matched_controls": len(a["caliper_groups"]),
            "median_same_archipelago_third_donor_km":
                float(np.median(a["third_distances"])) if a["third_distances"] else None,
            "median_matched_distance_gap_log":
                float(np.median(a["match_gaps"])) if a["match_gaps"] else None,
            "n_three_nearest_geo_donors_all_same_archipelago":
                a["n_geo_nearest_identical"],
            "fraction_three_nearest_geo_donors_all_same_archipelago":
                a["n_geo_nearest_identical"] /
                a["n_targets_with_source_predictions"] if
                a["n_targets_with_source_predictions"] else None,
            "mean_source_geo_donor_overlap_of_three":
                float(np.mean(a["shared_donor_counts"])) if
                a["shared_donor_counts"] else None,
        })
    results = []
    predictive = []
    for scope in SCOPES:
        rows = target[scope]
        for endpoint in ENDPOINTS:
            same_endpoint = [r for r in rows if r["endpoint"] == endpoint]
            dist_matched = [r for r in same_endpoint if r["outside_error"] is not None]
            predictive.append({
                "scope": scope, "endpoint": endpoint,
                "n_islands": len(same_endpoint),
                "n_archipelagos": len({r["archipelago"] for r in same_endpoint}),
                "mean_region_outside_error": float(np.mean([
                    r["baseline_error"] for r in same_endpoint
                ])) if same_endpoint else None,
                "mean_within_group_error": float(np.mean([
                    r["within_error"] for r in same_endpoint
                ])) if same_endpoint else None,
                "mean_distance_matched_other_error": float(np.mean([
                    r["outside_error"] for r in dist_matched
                ])) if dist_matched else None,
                "n_distance_matched": len(dist_matched),
                "mean_geography_only_error": float(np.mean([
                    r["geo_error"] for r in same_endpoint
                ])) if same_endpoint else None,
            })
            for comparison in CONTRASTS:
                eligible = [r for r in same_endpoint
                            if r[comparison] is not None]
                # Predefined: one average per named source archipelago, NOT
                # per island. A label spanning regions is still one unit.
                values = defaultdict(list)
                subgroups = set()
                for r in eligible:
                    values[r["archipelago"]].append(r[comparison])
                    subgroups.add((r["archipelago"], r["regime"]))
                result = {
                    "scope": scope, "endpoint": endpoint,
                    "comparison": comparison,
                    "n_islands": len(eligible),
                    "n_archipelagos": len(values),
                    "n_regional_subgroups": len(subgroups),
                    "n_randomizations": permutations,
                    "n_group_bootstrap": bootstrap,
                }
                if len(values) < 8:
                    result["status"] = "insufficient_archipelago_replication"
                else:
                    one_per_archipelago = np.asarray([
                        np.mean(v) for v in values.values()
                    ])
                    stat = _cluster_test(
                        one_per_archipelago,
                        permutations=permutations,
                        bootstrap=bootstrap,
                        seed=seed + (0 if scope=="all" else 100) +
                             (0 if endpoint==ENDPOINTS[0] else 10) +
                             CONTRASTS.index(comparison),
                    )
                    result.update(stat)
                    result["status"] = "posthoc_archipelago_heldout_predictive_comparison"
                results.append(result)
    # Keep every predeclared comparison in the FDR family, even unsupported.
    if len(results) != 12:
        raise ValueError("Expected exactly twelve scope×metric×contrast slots")
    q = _bh([r.get("p_cluster_signflip_one_sided",1.) for r in results])
    H = sum(1/i for i in range(1,len(results)+1))
    for i,r in enumerate(results):
        if "p_cluster_signflip_one_sided" in r:
            r["q_bh_12"] = float(q[i])
            r["q_by_12"] = float(min(1., q[i]*H))
    _write_csv(output/"archipelago_holdout_tests.csv", results, RESULT_FIELDS)
    _write_csv(output/"archipelago_support_by_region.csv",region_counts,
               ("analysis_regime",) + GROUP_FIELDS[2:])
    _write_csv(output/"predictive_errors.csv",predictive,PRED_FIELDS)
    manifest={
        "schema":"fd_leave_one_island_out_archipelago_profile_v1",
        "status":"independent_archipelago_units_observational_prediction_not_geological_effect",
        "unit_of_training":"other archipelago islands same geographic regime",
        "unit_of_evaluation":"held-out distinct island floral category-profile distribution",
        "unit_of_inference":"named Weigelt source archipelago; source group average, not independent island",
        "common_evidence_islands":len(chosen),
        "n_named_source_archipelagos":len({r["archipelago"] for r in chosen}),
        "n_regime_archipelago_subgroups_at_least4":sum(
            x["n_archipelago_regime_groups_at_least_4"] for x in region_counts),
        "n_targets_tested":len(target["all"])//2,
        "n_distance_matched_targets":sum(
            x["n_targets_with_distance_matched_controls"] for x in region_counts),
        "n_posthoc_comparison_slots":len(results),
        "settings":{
            "n_nearest_source_donors":REFERENCE_NEIGHBORS,
            "geography_only_negative_control":
                "three nearest eligible islands same regime, ignoring archipelago identity for source selection; 50% blend with SAME outside-group background",
            "negative_control_design":"added as post-hoc follow-up after first 8-test source-vs-distance audit; full revised 12-slot FDR family",
            "blend_source_to_region_reference":SOURCE_BLEND,
            "other_archipelago_distance_ratio_limit":DISTANCE_RATIO_LIMIT,
            "within_donor_radius_limit_km":WITHIN_RADIUS_KM,
            "reference_population":"all and direct both >=50 complete colour/form name units",
            "source_label_missing_policy":"NA/blank/unresolved omitted, not a shared archipelago",
            "regional_background":"all known other-archipelago islands in same regime, excludes every island of test archipelago",
            "original_profile_space_dim":len(space.ids),
            "n_sign_flips_per_comparison":permutations,
            "n_archipelago_bootstrap_draws":bootstrap,
        },
        "statistical_limitations":[
            "Source archipelago labels may be broad and not phylogenetic colonist source definitions",
            "Training islands may share taxon identities with test island; this is leave-island-out, NOT heldout species or evolutionary validation",
            "Island richness and trait reporting mechanisms remain observational",
            "Source predictions reflect ordinary spatial proximity and environmental similarities not independently identified evolutionary persistence",
            "Geographical matching on centroid distance does not guarantee identical environments, archipelago age, source plant composition or field sampling quality",
            "All and Direct evidence scopes overlap; within and distance-controlled tests share islands/archipelagos",
            "Cluster sign flips assume exchangeability of geographic-group predictive gains; between-group shared donor baselines can still correlate",
            "No geological class contrast, true absence or dispersal/establishment probability is estimated",
        ],
        "source_sha256":{
            "profile_reference":file_sha256(input_root / REF_INPUT),
            "crosswalk":file_sha256(crosswalk_path),
            **{f"{s}_g0":file_sha256(g0_root/f"{s}_islands.csv")
               for s in SCOPES},
            **{f"{s}_carriers":file_sha256(
                input_root/"floral_space_20261106"/s/
                "island_profile_carriers.csv.gz") for s in SCOPES},
        },
        "privacy":"Only aggregate results and source checksums; no island IDs, named archipelagos, geographic source IDs, species, genus or trait evidence rows emitted",
    }
    (output/"manifest.json").write_text(
        json.dumps(manifest,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps({
        "status":manifest["status"],
        "n_common_islands":manifest["common_evidence_islands"],
        "n_archipelagos":manifest["n_named_source_archipelagos"],
        "region_support":region_counts,
        "tests":[{k:r.get(k) for k in (
            "scope","endpoint","comparison","n_islands","n_archipelagos",
            "mean_group_improvement","p_cluster_signflip_one_sided",
            "q_bh_12","status"
        )} for r in results],
    },ensure_ascii=False))
    return manifest


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("input","g0","crosswalk","output"):
        p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--permutations",type=int,default=PERMUTATIONS)
    p.add_argument("--bootstrap",type=int,default=BOOTSTRAP)
    args=p.parse_args()
    run(args.input,args.g0,args.crosswalk,args.output,
        permutations=args.permutations,bootstrap=args.bootstrap)


if __name__=="__main__":
    main()
