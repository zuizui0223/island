"""Measured but projection-limited floral hypervolumes.

The complete trait-space is a finite categorical set, not a continuous Euclidean
space.  Report (1) observed occupied categorical reference cells, and
(2) a sensitivity-only 90%-mass Gaussian KDE area in a common **two-axis**
PCoA projection (15.35% positive inertia), never full floral hypervolume.

KDE uses equal weights across PRESENT unique profiles and a common, prechosen
isotropic bandwidth after scaling each reference PCoA axis by the global
496-reference standard deviation. Species carrier weights are a separate
sensitivity. Numerical integration uses a fixed common grid.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from .core import csv_rows, load_carriers
from .pipeline import REF_INPUT, SCOPES, _write_csv, file_sha256

BANDWIDTHS = (0.35, 0.55)
GRID_N = 88
MASS = 0.90
FIELDS = (
    "scope", "island_id", "analysis_regime", "n_names", "n_profiles",
    "observed_profile_cells", "observed_reference_occupancy_fraction",
    "kde90_area_pc12_equal_h035", "kde90_area_pc12_carrier_h035",
    "kde90_area_pc12_equal_h055", "kde90_area_pc12_carrier_h055",
    "min_grid_captured_mass",
)
SUMMARY_FIELDS = (
    "scope", "analysis_regime", "n_islands", "median_occupied_reference_cells",
    "median_kde90_pc12_equal_h035", "median_kde90_pc12_equal_h055",
    "median_kde90_pc12_carrier_h035",
)


def _reference_pc12(path: Path) -> tuple[list[int], np.ndarray, np.ndarray]:
    ids, pcs = [], []
    for row in csv_rows(path):
        ids.append(int(row["profile_id"]))
        pcs.append([float(row["PCoA1"]), float(row["PCoA2"])])
    order = np.argsort(ids)
    ordered = [ids[i] for i in order]
    if len(ordered) < 2 or len(set(ordered)) != len(ordered):
        raise ValueError("Reference requires distinct PCoA positions keyed by unique IDs")
    coords = np.asarray(pcs, dtype=float)[order]
    if not np.isfinite(coords).all():
        raise ValueError("Invalid PCoA reference coordinates")
    sd = coords.std(axis=0)
    if np.any(sd <= 0):
        raise ValueError("Both reference axes must have nonzero standard deviation")
    return ordered, coords / sd, sd


def _grid_kernels(coords: np.ndarray, h: float, n_grid: int) -> tuple[np.ndarray, float]:
    if h <= 0 or n_grid < 12:
        raise ValueError("Positive bandwidth and at least 12 grid points per axis required")
    # Support of ALL 496 profiles + four kernel SD at the boundary.
    low = coords.min(axis=0) - 4.0 * h
    high = coords.max(axis=0) + 4.0 * h
    delta = (high - low) / n_grid
    xx = low[0] + (np.arange(n_grid) + 0.5) * delta[0]
    yy = low[1] + (np.arange(n_grid) + 0.5) * delta[1]
    mesh_x, mesh_y = np.meshgrid(xx, yy, indexing="ij")
    points = np.column_stack((mesh_x.ravel(), mesh_y.ravel()))
    kernels = np.empty((len(coords), len(points)), dtype=np.float64)
    coeff = 1.0 / (2.0 * np.pi * h * h)
    # Keep temporary allocation bounded for the 496 x n_grid^2 density kernels.
    for start in range(0, len(coords), 32):
        stop = min(start + 32, len(coords))
        diff = coords[start:stop, None, :] - points[None, :, :]
        kernels[start:stop] = coeff * np.exp(-0.5 * np.sum(diff * diff, axis=2) / (h*h))
    return kernels, float(delta[0] * delta[1])


def _area90(pdf: np.ndarray, grid_cell_area: float) -> tuple[np.ndarray, np.ndarray]:
    """90%-mass HDR via grid cells sorted by mixture density, with no rescaling.

    A cell containing the mass threshold is included fully. Grid mass must
    cover essentially the whole mixture, otherwise the result is rejected.
    """
    if pdf.ndim != 2 or not np.isfinite(pdf).all() or np.any(pdf < 0):
        raise ValueError("Expected finite nonnegative densities per island")
    total_mass = pdf.sum(axis=1) * grid_cell_area
    if np.any(total_mass < 0.994) or np.any(total_mass > 1.006):
        raise ValueError("Fixed grid fails KDE mass capture; expand grid/resolution")
    sorted_density = np.sort(pdf, axis=1)[:, ::-1]
    cumulative = np.cumsum(sorted_density, axis=1) * grid_cell_area
    first_reached = (cumulative >= MASS).argmax(axis=1)
    if np.any(cumulative[:, -1] < MASS):
        raise ValueError("90% mass lies outside grid")
    return (first_reached + 1) * grid_cell_area, total_mass


def _normalize_reference_counts(rows: list[dict[int, int]],
                                lookup: dict[int, int], n_profiles: int
                                ) -> tuple[np.ndarray, np.ndarray]:
    equal = np.zeros((len(rows), n_profiles), dtype=np.float64)
    carriers = np.zeros_like(equal)
    for i, row in enumerate(rows):
        unknown = set(row) - lookup.keys()
        if unknown:
            raise ValueError(f"Unknown profiles: {sorted(unknown)[:5]}")
        inds = [lookup[key] for key in row]
        values = np.asarray(list(row.values()), dtype=float)
        if np.any(values <= 0):
            raise ValueError("Nonpositive carrier count")
        equal[i, inds] = 1.0 / len(row)
        carriers[i, inds] = values / values.sum()
    return equal, carriers


def run_hypervolume(input_root: Path, g0_root: Path, output_root: Path, *,
                    bandwidths: tuple[float, float] = BANDWIDTHS,
                    grid_n: int = GRID_N, batch_size: int = 40) -> dict:
    input_root, g0_root, output_root = map(Path, (input_root, g0_root, output_root))
    output_root.mkdir(parents=True, exist_ok=True)
    ref = input_root / REF_INPUT
    ids, coords, sd = _reference_pc12(ref)
    lookup = {v: i for i, v in enumerate(ids)}
    kernels = {h: _grid_kernels(coords, h, grid_n) for h in bandwidths}
    carriers = {
        scope: load_carriers(
            input_root / "floral_space_20261106" / scope / "island_profile_carriers.csv.gz"
        ) for scope in SCOPES
    }
    summaries, all_rows, source_paths = [], [], {}
    for scope in SCOPES:
        g0 = g0_root / f"{scope}_islands.csv"
        source_paths[scope] = g0
        eligible = list(csv_rows(g0))
        output = []
        for begin in range(0, len(eligible), batch_size):
            block = eligible[begin:begin + batch_size]
            profile_rows = [carriers[scope][record["island_id"]] for record in block]
            equal, weighted = _normalize_reference_counts(profile_rows, lookup, len(ids))
            volumes, captured = {}, []
            for h in bandwidths:
                kernel, cell_area = kernels[h]
                volume_eq, mass_eq = _area90(equal @ kernel, cell_area)
                volume_weight, mass_weight = _area90(weighted @ kernel, cell_area)
                volumes[h] = (volume_eq, volume_weight)
                captured.extend([mass_eq, mass_weight])
            min_mass = np.min(np.vstack(captured), axis=0)
            for i, record in enumerate(block):
                n_names = int(record["n_names"])
                n_profiles = int(record["n_profiles"])
                if n_profiles != len(profile_rows[i]) or n_names != sum(profile_rows[i].values()):
                    raise ValueError("G0 island richness differs from public carrier counts")
                output.append({
                    "scope": scope, "island_id": record["island_id"],
                    "analysis_regime": record["analysis_regime"],
                    "n_names": n_names, "n_profiles": n_profiles,
                    "observed_profile_cells": n_profiles,
                    "observed_reference_occupancy_fraction": n_profiles / len(ids),
                    "kde90_area_pc12_equal_h035": float(volumes[bandwidths[0]][0][i]),
                    "kde90_area_pc12_carrier_h035": float(volumes[bandwidths[0]][1][i]),
                    "kde90_area_pc12_equal_h055": float(volumes[bandwidths[1]][0][i]),
                    "kde90_area_pc12_carrier_h055": float(volumes[bandwidths[1]][1][i]),
                    "min_grid_captured_mass": float(min_mass[i]),
                })
        _write_csv(output_root / f"{scope}_hypervolumes.csv", output, FIELDS)
        all_rows.extend(output)
        for region in sorted({r["analysis_regime"] for r in output}):
            selected = [r for r in output if r["analysis_regime"] == region]
            summaries.append({
                "scope": scope, "analysis_regime": region, "n_islands": len(selected),
                "median_occupied_reference_cells": float(
                    np.median([r["observed_profile_cells"] for r in selected])),
                "median_kde90_pc12_equal_h035": float(
                    np.median([r["kde90_area_pc12_equal_h035"] for r in selected])),
                "median_kde90_pc12_equal_h055": float(
                    np.median([r["kde90_area_pc12_equal_h055"] for r in selected])),
                "median_kde90_pc12_carrier_h035": float(
                    np.median([r["kde90_area_pc12_carrier_h035"] for r in selected])),
            })
    _write_csv(output_root / "regional_hypervolume_summary.csv", summaries, SUMMARY_FIELDS)
    manifest = {
        "schema": "fd_projection_limited_hypervolume_v1",
        "status": "computed_projected_density_breadth_descriptive",
        "categorical_hypervolume": "distinct recorded reference profiles per island; a cardinality, not Lebesgue volume",
        "kde_hypervolume": "90% HDR area of common PCoA1-PCoA2 Gaussian mixture",
        "do_not_claim": "full multidimensional floral hypervolume or actual pollinator functional niche volume",
        "reference_profiles": len(ids), "projection_positive_inertia": 0.15353075806441305,
        "two_axis_projection_loss": 1.0 - 0.15353075806441305,
        "bandwidth_std_reference_units": list(bandwidths), "axis_sd_unstandardized": sd.tolist(),
        "grid_n_each_axis": grid_n, "mass": MASS, "smoothing": "isotropic Gaussian kernel with FIXED bandwidth globally",
        "primary_weighting": "each distinct observed profile equal; species-carrier weighting is a sensitivity",
        "integrator": "deterministic rectangular midpoint grid; include complete threshold cell",
        "note": "equal-profile KDE area can be positive for a single recorded category due solely to smoothing",
        "source_sha256": {"profile_reference": file_sha256(ref),
                          **{scope: file_sha256(path) for scope, path in source_paths.items()}},
        "n_islands": len(all_rows), "n_summaries": len(summaries),
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print(json.dumps({"status": manifest["status"], "n_islands": len(all_rows),
                      "regions": summaries[:4]}, ensure_ascii=False))
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--g0", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--grid", type=int, default=GRID_N)
    args = parser.parse_args()
    run_hypervolume(args.input, args.g0, args.output, grid_n=args.grid)


if __name__ == "__main__":
    main()
