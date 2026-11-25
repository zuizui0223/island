"""Projection-limited 17-dimensional floral KDE hypervolume.

From the frozen categorical set-Jaccard distance D across 496 profiles,
classical PCoA on sqrt(D) yields a nonnegative gram matrix. The first 17
axes retain 80.49% of the reference positive inertia. Monte Carlo estimates
90% Gaussian KDE highest-density-region volume (log units) in these 17 axes,
for a common Euclidean scaling and bandwidth.  It is NOT full categorical
trait-space volume or measured ecological pollination function.

All output is aggregate-only; no island identifications or private joins
are emitted to the results/ directories.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp

from .core import ProfileSpace, load_profiles, load_carriers, csv_rows
from .g1_fixed_margins import _bh
from .pipeline import REF_INPUT, file_sha256
from .repertoire import REGIONS, match_origin_islands, _two_sided

DIMENSIONS = 17
KEEP_MASS = 0.90
BANDWIDTHS = (0.55, 0.85)
SCOPES = ("all", "direct")


def pcoa_embedding_17(space: ProfileSpace) -> tuple[np.ndarray, dict]:
    """Eigen-embed sqrt(original distance), with ONE common scalar scaling.

    Whitening every low-inertia axis separately would artificially elevate
    minor floral dimensions, so no axis-wise rescaling is performed.
    """
    d = space.matrix
    center = d - d.mean(axis=0)[None, :] - d.mean(axis=1)[:, None] + d.mean()
    values, vectors = np.linalg.eigh(-0.5 * center)
    order = np.argsort(values)[::-1]
    values, vectors = values[order], vectors[:, order]
    positive = values > 1e-10
    if len(values[positive]) < DIMENSIONS or np.any(values < -1e-7):
        raise ValueError("Insufficient positive inertia or non-Euclidean sqrt-D")
    axes = vectors[:, :DIMENSIONS] * np.sqrt(values[:DIMENSIONS])
    reference_sd = float(np.sqrt(np.sum(axes.var(axis=0)) / DIMENSIONS))
    if not np.isfinite(reference_sd) or reference_sd <= 0:
        raise ValueError("Reference axes have no scale")
    pc = axes / reference_sd
    audit = {
        "positive_inertia_2axes": float(values[:2].sum() / values[positive].sum()),
        "positive_inertia_17axes": float(values[:DIMENSIONS].sum() / values[positive].sum()),
        "scale_factor": reference_sd,
        "n_ref_profiles": len(space.ids),
        "n_positive_axes": int(positive.sum()),
    }
    return pc, audit


def kde_hdr_17(points: np.ndarray, *, bandwidth: float,
               n_draws: int, seed: int) -> tuple[float, float]:
    """Estimate volume via importance identity E_f[1_HDR/f] with MC sampling.

    HDR threshold is 10th percentile of log f at KDE-sampled points. Uses a
    midpoint empirical quantile; MC SE is based on importance weights and
    includes the estimated threshold uncertainty only approximately.
    """
    if points.ndim != 2 or points.shape[1] != DIMENSIONS:
        raise ValueError("Expected exactly 17 PCoA axes")
    if len(points) < DIMENSIONS or n_draws < 1000:
        raise ValueError("Too few observed profiles or MC draws")
    if not np.isfinite(points).all() or bandwidth <= 0:
        raise ValueError("Bad KDE inputs/bandwidth")
    rng = np.random.default_rng(seed)
    centers = points[rng.integers(len(points), size=n_draws)]
    samples = centers + rng.normal(0, bandwidth, size=(n_draws, DIMENSIONS))
    dsq = np.maximum(
        0.0, (samples ** 2).sum(axis=1)[:, None]
        + (points ** 2).sum(axis=1)[None, :]
        - 2 * samples @ points.T,
    )
    logdens = (logsumexp(-dsq / (2 * bandwidth * bandwidth), axis=1)
               - math.log(len(points))
               - DIMENSIONS / 2 * math.log(2 * math.pi * bandwidth ** 2))
    in_hdr = logdens >= float(np.quantile(logdens, 1 - KEEP_MASS))
    logvolume = float(logsumexp(-logdens[in_hdr]) - math.log(n_draws))
    logw = np.full(n_draws, -np.inf)
    logw[in_hdr] = -logdens[in_hdr]
    weights = np.exp(logw - np.max(logw))
    relative_error = float(
        weights.std(ddof=1) / np.sqrt(n_draws) / weights.mean())
    if not np.isfinite(logvolume) or not np.isfinite(relative_error):
        raise ValueError("Numerically unstable high-dimensional integration")
    return logvolume, relative_error


def _paired_origin_stats(scores: pd.DataFrame, pairs: list[dict],
                         *, scope: str, bandwidth: float, seed: int) -> list[dict]:
    lookup = scores.set_index("island_id")
    results = []
    for region in ("global", *REGIONS):
        selected = pairs if region == "global" else [p for p in pairs if p["region"] == region]
        selected = [p for p in selected
                    if p["continental_id"] in lookup.index and
                    p["oceanic_id"] in lookup.index]
        if len(selected) < 8:
            continue
        difference = np.asarray([
            lookup.at[p["oceanic_id"], "log_volume"] -
            lookup.at[p["continental_id"], "log_volume"]
            for p in selected
        ], dtype=float)
        rng = np.random.default_rng(seed)
        permutation = (rng.choice((-1., 1.), size=(9999, len(difference))) *
                       difference[None, :]).mean(axis=1)
        observed = float(difference.mean())
        p_two = _two_sided(observed, permutation)
        ii = rng.integers(len(difference), size=(9999, len(difference)))
        qlo, qhi = np.quantile(difference[ii].mean(axis=1), [.025, .975])
        # Test whether the projected volume difference is more than
        # the difference in the number of recorded distinct profiles.
        profile_log_ratio = np.asarray([
            np.log(lookup.at[p["oceanic_id"], "n_profiles"] /
                   lookup.at[p["continental_id"], "n_profiles"])
            for p in selected
        ])
        design = np.column_stack([np.ones(len(difference)), profile_log_ratio])
        beta = np.linalg.lstsq(design, difference, rcond=None)[0]
        intercepts = []
        for bootstrap in ii:
            if np.linalg.matrix_rank(design[bootstrap]) == design.shape[1]:
                intercepts.append(float(np.linalg.lstsq(
                    design[bootstrap], difference[bootstrap], rcond=None)[0][0]))
        adjlo, adjhi = (
            np.quantile(intercepts, [.025, .975]) if len(intercepts) >= 8000
            else (np.nan, np.nan)
        )
        results.append({
            "profile_richness_adjusted_intercept": float(beta[0]),
            "profile_richness_coefficient": float(beta[1]),
            "profile_richness_adjusted_bootstrap_q025": float(adjlo),
            "profile_richness_adjusted_bootstrap_q975": float(adjhi),
            "scope": scope, "bandwidth": bandwidth, "region": region,
            "n_distinct_pairs": len(selected),
            "mean_log_volume_diff_o_minus_c": observed,
            "paired_bootstrap_q025": float(qlo),
            "paired_bootstrap_q975": float(qhi),
            "p_pair_swap_two_sided": p_two,
            "status": "exploratory_projection_only",
        })
    return results


def run(input_root: Path, g0_root: Path, crosswalk: Path, output: Path, *,
        scopes: tuple[str, ...] = SCOPES, n_draws: int = 12000,
        seed: int = 20261108) -> dict:
    input_root, g0_root, crosswalk, output = map(
        Path, (input_root, g0_root, crosswalk, output))
    if any(s not in SCOPES for s in scopes) or not scopes:
        raise ValueError("Select All, Direct or both")
    output.mkdir(parents=True, exist_ok=True)
    space = ProfileSpace(load_profiles(input_root / REF_INPUT))
    coords, pcoa = pcoa_embedding_17(space)
    if not .79 < pcoa["positive_inertia_17axes"] < .82:
        raise ValueError("Frozen profile-reference PCoA variance changed")
    origins = {r["island_id"]: r for r in csv_rows(crosswalk)
               if r["geological_origin"] in ("continental", "oceanic")}
    all_rows = list(csv_rows(g0_root / "all_islands.csv"))
    direct_rows = list(csv_rows(g0_root / "direct_islands.csv"))
    matches, match_audit = match_origin_islands(
        all_rows, direct_rows,
        {r["island_id"]: r for r in csv_rows(crosswalk)},
    )
    records, summaries, paired, quality = [], [], [], []
    for scope_index, scope in enumerate(scopes):
        g0 = list(csv_rows(g0_root / f"{scope}_islands.csv"))
        ids = sorted(r["island_id"] for r in g0 if r["island_id"] in origins)
        carrier = load_carriers(input_root / "floral_space_20261106" /
                                scope / "island_profile_carriers.csv.gz")
        group = {r["island_id"]: r["analysis_regime"] for r in g0}
        for index, island in enumerate(ids):
            if island not in carrier:
                raise ValueError("Missing public island-profile record")
            indices = [space.lookup[p] for p in sorted(carrier[island])]
            if len(indices) < DIMENSIONS:
                raise ValueError("Geology-supported eligible island has too few profiles")
            for bw_index, bandwidth in enumerate(BANDWIDTHS):
                values, errors = [], []
                for replicate in range(2):
                    value, err = kde_hdr_17(
                        coords[indices], bandwidth=bandwidth, n_draws=n_draws,
                        seed=seed + scope_index * 100000000 + index * 101
                             + bw_index * 13 + replicate)
                    values.append(value)
                    errors.append(err)
                records.append({
                    "scope": scope, "island_id": island,
                    "origin": origins[island]["geological_origin"],
                    "regime": group[island], "bandwidth": bandwidth,
                    "n_profiles": len(indices),
                    "log_volume": float(np.mean(values)),
                    "replicate_abs_log_delta": float(abs(values[0] - values[1])),
                    "mean_MC_relative_error": float(np.mean(errors)),
                })
        records_scope = pd.DataFrame([r for r in records if r["scope"] == scope])
        for bw_index, bw in enumerate(BANDWIDTHS):
            bw_rows = records_scope.loc[records_scope.bandwidth == bw]
            quality.append({
                "scope": scope, "bandwidth": bw, "n_islands": len(bw_rows),
                "median_MC_relative_error": float(
                    bw_rows.mean_MC_relative_error.median()),
                "median_replicate_abs_log_delta": float(
                    bw_rows.replicate_abs_log_delta.median()),
                "p95_replicate_abs_log_delta": float(
                    bw_rows.replicate_abs_log_delta.quantile(.95)),
            })
            for origin, subgroup in bw_rows.groupby("origin"):
                summaries.append({
                    "scope": scope, "bandwidth": bw,
                    "geological_origin": origin,
                    "n_islands": len(subgroup),
                    "median_17D_projected_log_hypervolume": float(
                        subgroup.log_volume.median()),
                    "median_n_profiles": float(subgroup.n_profiles.median()),
                    "median_mc_relative_error": float(
                        subgroup.mean_MC_relative_error.median()),
                })
            paired += _paired_origin_stats(
                bw_rows, matches, scope=scope, bandwidth=bw,
                seed=seed + bw_index + scope_index * 23,
            )
    # Post-hoc paired 17D projection: 2 bandwidths × two nested scope analyses.
    primary = [r for r in paired if r["region"] == "global"]
    if len(primary) != 2 * len(scopes):
        raise ValueError("Expected two bandwidth tests for each evidence scope")
    qs = _bh([r["p_pair_swap_two_sided"] for r in primary])
    harmonic = sum(1/i for i in range(1, len(primary)+1))
    for row, q in zip(primary, qs):
        row["q_bh_projection_family"] = float(q)
        row["q_by_projection_family"] = float(min(1, q*harmonic))
    pd.DataFrame(summaries).to_csv(output / "kde17_origin_aggregate.csv", index=False)
    pd.DataFrame(paired).to_csv(output / "kde17_paired_sensitivity.csv", index=False)
    pd.DataFrame(quality).to_csv(output / "kde17_monte_carlo_quality.csv", index=False)
    manifest = {
        "schema": "fd_pcoa17_kde90_monte_carlo_v1",
        "status": "projected_17D_KDE_not_full_categorical_volume",
        "method": "90%-mass 17D Gaussian-KDE HDR volume estimated by sampling directly from fitted mixture and E_f(1_HDR/f)",
        "projection": pcoa,
        "bandwidths": list(BANDWIDTHS),
        "bandwidth_units": "common scalar-ref-standardized original PCoA coordinates",
        "n_draws_per_replicate": n_draws, "replicates": 2,
        "scopes": list(scopes), "n_source_classified_matched_islands": len(matches)*2,
        "n_all_origin_island_records": len(records),
        "n_global_paired_tests": len(primary),
        "interpretation": "17D retains only 80.49% of positive PCoA inertia; smoothing occupies unobserved state-set configurations; remaining ~19.51% of positive inertia is omitted",
        "not_shown": "full-space hypervolume, causal ecological selection, pollination fitness, independent origin replicate",
        "source_hashes_sha256": {
            "profile_reference": file_sha256(input_root/REF_INPUT),
            "geology_crosswalk": file_sha256(crosswalk),
        },
        "data_rights": "aggregate-only outputs; no island IDs, taxon identities, source IDs",
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({
        "n_island_records": len(records), "inertia17": pcoa["positive_inertia_17axes"],
        "summaries": summaries, "global_paired_tests": primary,
        "monte_carlo_quality": quality,
    }, ensure_ascii=False))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("input", "g0", "crosswalk", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--scopes", default="all,direct")
    p.add_argument("--samples", type=int, default=12000)
    args = p.parse_args()
    run(args.input, args.g0, args.crosswalk, args.output,
        scopes=tuple(args.scopes.split(",")), n_draws=args.samples)


if __name__ == "__main__":
    main()
