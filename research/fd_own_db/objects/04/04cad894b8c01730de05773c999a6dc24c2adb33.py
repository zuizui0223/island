"""Fixed-absolute-distance donor radius falsification of recorded island floral packing.

N1c uses ten nearest islands, which gives very different geographic reach for
oceanic and continental focal islands. N1d replaces a fixed *number* of
neighbors with identical absolute 1000/2000 km reach for both geological
classes. Recorded source-island profile incidence, not lineages, is modeled.
Post-hoc descriptive only; no newly optimized radius or geological p-values.
No island IDs, donor maps, or raw coordinates leave this module.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from .core import ProfileSpace, csv_rows, load_carriers, load_profiles
from .distance_ses import (
    MODES, SUMMARY_FIELDS, _score, _seed, _summary, distance_statistics,
    regional_profile_counts, regional_leave_one_out_pool, sample_null,
)
from .local_geography_ses import great_circle_distance_km, local_pool
from .pipeline import REF_INPUT, _write_csv, file_sha256

SCOPES = ("all", "direct")
RADII_KM = (1000, 2000)  # post-hoc fixed design; do NOT select after inspecting SES
MIN_DONOR_ISLANDS = 3
REGIONAL_NULL = "N1b_regional_islandfreq_LOO"
RADIUS_NULL = {radius: f"N1d_radius{radius}km_geo_LOO" for radius in RADII_KM}

SUPPORT_FIELDS = (
    "scope", "geological_origin", "radius_km", "n_geology_islands",
    "n_with_at_least_3_regional_donor_islands", "n_equal30_estimable",
    "n_full_estimable", "median_donor_island_count",
    "median_donor_distinct_profiles", "median_farthest_donor_distance_km",
    "status",
)
MATCHED_FIELDS = (
    "scope", "radius_km", "mode", "endpoint", "region", "geological_origin",
    "n_islands_both_nulls", "status",
    "median_N1b_SES", "median_N1d_SES", "median_within_island_SES_change",
    "fraction_negative_N1b", "fraction_negative_N1d",
    "median_N1b_observed_minus_null_D", "median_N1d_observed_minus_null_D",
    "median_N1d_minus_N1b_expected_D",
)


def within_radius_donors(focal: str, regime: str,
                         eligible: dict[str, dict], geo: dict[str, dict],
                         radius_km: float) -> tuple[list[str], float]:
    """All OTHER eligible same-regime islands in a preselected absolute radius."""
    if focal not in eligible or focal not in geo or not (
        np.isfinite(radius_km) and 0 < radius_km <= 20020
    ):
        raise ValueError("Invalid focal island or fixed geographic radius")
    try:
        lat = float(geo[focal]["island_latitude"])
        lon = float(geo[focal]["island_longitude"])
    except (KeyError, ValueError, TypeError) as exc:
        raise ValueError("Focal coordinate unavailable") from exc
    keys, coords = [], []
    for island, row in eligible.items():
        if island == focal or row["analysis_regime"] != regime:
            continue
        match = geo.get(island)
        if not match:
            continue
        try:
            y, x = float(match["island_latitude"]), float(match["island_longitude"])
        except (KeyError, ValueError, TypeError):
            continue
        if not np.isfinite([y, x]).all() or abs(y) > 90 or abs(x) > 180:
            continue
        keys.append(island)
        coords.append((y, x))
    if not keys:
        return [], float("nan")
    xy = np.asarray(coords, dtype=float)
    dd = great_circle_distance_km(lat, lon, xy[:, 0], xy[:, 1])
    included = sorted(
        ((float(d), island) for d, island in zip(dd, keys)
         if float(d) <= radius_km + 1e-9),
        key=lambda t: (t[0], t[1])
    )
    return [island for d, island in included], (
        included[-1][0] if included else float("nan"))


def fixed_radius_support_summary(rows: list[dict]) -> list[dict]:
    groups = defaultdict(list)
    for row in rows:
        groups[row["scope"], row["geological_origin"], row["radius_km"]].append(row)
    result = []
    for (scope, origin, radius), data in sorted(groups.items()):
        supported = [r for r in data if r["n_donor_islands"] >= MIN_DONOR_ISLANDS]
        record = {
            "scope": scope, "geological_origin": origin, "radius_km": radius,
            "n_geology_islands": len(data),
            "n_with_at_least_3_regional_donor_islands": len(supported),
            "n_equal30_estimable": sum(r["has_equal30"] for r in data),
            "n_full_estimable": sum(r["has_full"] for r in data),
            "status": "fixed_distance_posthoc_observation_reference",
        }
        if len(supported) < 5:
            record["status"] = "small_cell_suppressed"
        else:
            record.update({
                "median_donor_island_count": float(np.median(
                    [r["n_donor_islands"] for r in supported])),
                "median_donor_distinct_profiles": float(np.median(
                    [r["n_donor_profiles"] for r in supported])),
                "median_farthest_donor_distance_km": float(np.median(
                    [r["farthest_km"] for r in supported])),
            })
        result.append(record)
    return result


def matched_radius_summary(scores: list[dict], *, min_release_n: int = 5
                           ) -> list[dict]:
    """Use ONLY islands scored by N1b and the specified fixed-radius null."""
    by_island = defaultdict(dict)
    for row in scores:
        key = tuple(row[col] for col in (
            "scope", "mode", "endpoint", "region",
            "geological_origin", "island_id"))
        if row["null"] in (REGIONAL_NULL, *RADIUS_NULL.values()):
            if row["null"] in by_island[key]:
                raise ValueError("Duplicate island-null-reference metric")
            by_island[key][row["null"]] = row
    groups = defaultdict(list)
    for key, variants in by_island.items():
        b = variants.get(REGIONAL_NULL)
        if b is None:
            continue
        for radius, label in RADIUS_NULL.items():
            local = variants.get(label)
            if local is None:
                continue
            if not np.isclose(b["observed"], local["observed"],
                              atol=1e-10, rtol=1e-10):
                raise AssertionError("Cannot compare different observations as null change")
            g = (key[0], radius, key[1], key[2], key[3], key[4])
            groups[g].append((b, local))
    folded = defaultdict(list)
    for key, records in groups.items():
        folded[key].extend(records)
        if key[4] != "all":
            folded[(key[0], key[1], key[2], key[3], "all", key[5])].extend(records)
    output = []
    for key, pair_values in sorted(folded.items()):
        row = dict(zip(
            ("scope", "radius_km", "mode", "endpoint",
             "region", "geological_origin"), key))
        row["n_islands_both_nulls"] = len(pair_values)
        row["status"] = ("paired_same_island_descriptive"
                         if len(pair_values) >= min_release_n
                         else "small_cell_suppressed")
        if len(pair_values) >= min_release_n:
            a = np.asarray([r[0]["ses"] for r in pair_values])
            b = np.asarray([r[1]["ses"] for r in pair_values])
            row.update({
                "median_N1b_SES": float(np.median(a)),
                "median_N1d_SES": float(np.median(b)),
                "median_within_island_SES_change": float(np.median(b-a)),
                "fraction_negative_N1b": float(np.mean(a < 0)),
                "fraction_negative_N1d": float(np.mean(b < 0)),
                "median_N1b_observed_minus_null_D": float(np.median(
                    [p[0]["delta"] for p in pair_values])),
                "median_N1d_observed_minus_null_D": float(np.median(
                    [p[1]["delta"] for p in pair_values])),
                "median_N1d_minus_N1b_expected_D": float(np.median(
                    [p[1]["null_mean"]-p[0]["null_mean"] for p in pair_values])),
            })
        output.append(row)
    return output


def run(input_root: Path, g0_root: Path, crosswalk: Path, output: Path, *,
        draws: int = 999, thinnings: int = 99, seed: int = 20261109) -> dict:
    """Generate only aggregate N1b vs fixed-radius N1d reference summaries."""
    if draws < 20 or thinnings < 20:
        raise ValueError("Insufficient reference or observed thinning")
    input_root, g0_root, crosswalk, output = map(
        Path, (input_root, g0_root, crosswalk, output))
    space = ProfileSpace(load_profiles(input_root / REF_INPUT))
    if len(space.ids) != 496:
        raise ValueError("Frozen public reference changed")
    geo = {r["island_id"]: r for r in csv_rows(crosswalk)}
    scores, supports, carrier_files = [], [], {}
    for scope in SCOPES:
        g0 = list(csv_rows(g0_root / f"{scope}_islands.csv"))
        elig = {r["island_id"]: r for r in g0}
        path = input_root / "floral_space_20261106" / scope / "island_profile_carriers.csv.gz"
        carrier_files[scope] = path
        carriers = load_carriers(path)
        if set(elig) - set(carriers):
            raise ValueError("Missing eligible public floral carriers")
        regional = regional_profile_counts(g0, carriers)
        for row in g0:
            island = row["island_id"]
            origin = geo.get(island, {}).get("geological_origin")
            if origin not in ("continental", "oceanic"):
                continue
            profiles = set(carriers[island])
            ids = np.asarray([space.lookup[p] for p in sorted(profiles)], dtype=int)
            n = len(ids)
            if n < 2:
                continue
            regional_pool = regional_leave_one_out_pool(
                regional, row["analysis_regime"], profiles, space.lookup)
            regional_weights = np.asarray([
                regional[row["analysis_regime"]][space.ids[int(idx)]]
                - int(space.ids[int(idx)] in profiles)
                for idx in regional_pool], dtype=float)
            if np.any(regional_weights <= 0):
                raise AssertionError("Focal-island-excluded regional weights invalid")
            locations = {}
            for radius in RADII_KM:
                neighbors, farthest = within_radius_donors(
                    island, row["analysis_regime"], elig, geo, radius)
                enough = len(neighbors) >= MIN_DONOR_ISLANDS
                pool, weights = (local_pool(neighbors, carriers, space.lookup)
                                 if enough else
                                 (np.asarray([], dtype=int), np.asarray([], dtype=float)))
                locations[radius] = (pool, weights)
                supports.append({
                    "scope": scope, "geological_origin": origin,
                    "radius_km": radius, "n_donor_islands": len(neighbors),
                    "n_donor_profiles": len(pool), "farthest_km": farthest,
                    "has_full": enough and len(pool) >= n,
                    "has_equal30": enough and n >= 30 and len(pool) >= 30,
                })
            thinning = None
            if n >= 30:
                rng = np.random.default_rng(_seed(seed, scope, island, "N1d_30"))
                thinning = np.asarray([
                    distance_statistics(space.matrix, rng.choice(ids, 30, replace=False))
                    for _ in range(thinnings)
                ])
            for mode in MODES:
                if mode == "equal30" and thinning is None:
                    continue
                size = n if mode == "full" else 30
                observed = (distance_statistics(space.matrix, ids)
                            if mode == "full" else thinning.mean(axis=0))
                references = [(REGIONAL_NULL, regional_pool, regional_weights)]
                references += [(RADIUS_NULL[r], *locations[r]) for r in RADII_KM]
                for null, pool, weights in references:
                    if len(pool) < size:
                        continue
                    simulated = sample_null(
                        space.matrix, pool, size, draws=draws,
                        seed=_seed(seed, scope, island, mode, null),
                        weights=weights,
                    )
                    if simulated is None:
                        continue
                    for score in _score(
                            observed, simulated,
                            thinning=thinning if mode == "equal30" else None):
                        if score is not None:
                            scores.append({
                                "scope": scope, "island_id": island,
                                "region": row["analysis_regime"],
                                "geological_origin": origin,
                                "mode": mode, "null": null,
                                "n_profiles": n, **score,
                            })
    summary = _summary(scores)
    support = fixed_radius_support_summary(supports)
    paired = matched_radius_summary(scores)
    output.mkdir(parents=True, exist_ok=True)
    _write_csv(output / "fixed_radius_ses_summary.csv", summary, SUMMARY_FIELDS)
    _write_csv(output / "fixed_radius_ses_support.csv", support, SUPPORT_FIELDS)
    _write_csv(output / "fixed_radius_matched_reference.csv", paired, MATCHED_FIELDS)
    manifest = {
        "schema": "fd_fixed_1000_2000km_regional_flower_reference_v1",
        "status": "posthoc_descriptive_fixed_geographic_reach_sensitivity",
        "radii_km": list(RADII_KM),
        "minimum_other_eligible_donor_islands": MIN_DONOR_ISLANDS,
        "reference": "original equal-colour/form state-set Jaccard D",
        "nulls": {"N1b": "other regional island-incidence weighted profile pool",
                  "N1d": "other eligible same-regime islands within fixed km radius; pool incidence by donor island"},
        "support": "at least three eligible other islands AND enough distinct donor profiles; never pad with replacement",
        "comparison": "N1b and each N1d compared on identical physical focal islands; SES and original-D differences",
        "endpoints": ["mean off-diagonal D", "mean within-support nearest-other-profile D"],
        "draws": draws, "thinnings": thinnings, "seed": seed,
        "no_new_hypothesis_tests": True,
        "limitations": "regional observational profile records not genus/source/detection-controlled; radius cutoff may select geographically and geologically unbalanced reference availability; no adaptation inference",
        "privacy": "grouped rows with n<5 suppressed; no island IDs, actual source coordinates, or donor lists written",
        "source_sha256": {
            "profile_reference": file_sha256(input_root / REF_INPUT),
            "geological_crosswalk": file_sha256(crosswalk),
            **{scope+"_carrier": file_sha256(path) for scope, path in carrier_files.items()},
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({
        "schema": manifest["schema"],
        "n_aggregate_summary_rows": len(summary),
        "n_radius_support_rows": len(support),
        "n_same_focal_baseline_comparisons": len(paired),
        "n_internal_scores": len(scores),
        "interpretation": "No geological origin p-value or selection inference",
    }))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for field in ("input", "g0", "crosswalk", "output"):
        p.add_argument("--"+field, type=Path, required=True)
    p.add_argument("--draws", type=int, default=999)
    p.add_argument("--thinnings", type=int, default=99)
    a = p.parse_args()
    run(a.input, a.g0, a.crosswalk, a.output,
        draws=a.draws, thinnings=a.thinnings)


if __name__ == "__main__":
    main()
