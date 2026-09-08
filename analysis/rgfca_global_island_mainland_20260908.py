#!/usr/bin/env python3
"""Global RGFCA island-versus-mainland/reference repeated colour contrast.

Exploratory bridge between the FCP RGFCA image measurements and the island-v2
geographic universe.  This script does *not* search for a shared geographic
boundary.  It classifies each RGFCA photograph by exact point-in-polygon
intersection with the frozen island-v2 polygon universe and then asks whether
flower colour differs between island and non-island observations of the same
species under repeated equal-depth resampling.

Important scope note: island-v2's canonical polygon universe contains landmasses
>5 km2 and <=7,000,000 km2.  Therefore the complement is called
``nonisland_reference`` in outputs rather than unqualified ``mainland``: it is
dominated by continents but can contain islands <=5 km2.  This is retained as a
visible sensitivity issue rather than silently relabelled.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

PALETTE = (
    "white", "yellow", "orange", "red", "pink", "magenta", "purple",
    "blue", "bronze", "green", "brown", "black",
)
CHROMATIC = ("yellow", "orange", "red", "pink", "magenta", "purple", "blue")
MUTED = ("white", "bronze", "green", "brown", "black")
OUTCOMES = (
    "flower_white_share",
    "flower_chromatic_share",
    "flower_minus_background_white_share",
    "flower_minus_background_chromatic_share",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--rgfca-csv", type=Path, required=True)
    p.add_argument("--islands-gpkg", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--draws", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260908)
    p.add_argument("--depth-grid", nargs="+", type=int, default=[3, 5, 10])
    p.add_argument("--lat-bin-grid", nargs="+", type=int, default=[10, 5])
    return p.parse_args()


def classify_islands(frame: pd.DataFrame, gpkg: Path) -> pd.DataFrame:
    required = {"measurement_id", "latitude", "longitude"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"RGFCA table missing coordinate columns: {sorted(missing)}")
    islands = gpd.read_file(gpkg)[["island_id", "area_km2", "geometry"]].copy()
    if islands.crs is None:
        raise ValueError("island polygon CRS is missing")
    points = gpd.GeoDataFrame(
        frame[["measurement_id", "latitude", "longitude"]].copy(),
        geometry=gpd.points_from_xy(frame["longitude"], frame["latitude"]),
        crs="EPSG:4326",
    ).to_crs(islands.crs)
    joined = gpd.sjoin(points, islands, how="left", predicate="intersects")
    # Overlapping/nested polygons can yield >1 hit.  Deterministically assign the
    # smallest eligible intersected landmass; island/not-island status itself is
    # unchanged by this tie rule.
    joined = joined.sort_values(["measurement_id", "area_km2"], na_position="last")
    one = joined.drop_duplicates("measurement_id", keep="first")
    out = frame.merge(
        one[["measurement_id", "island_id", "area_km2"]],
        on="measurement_id",
        how="left",
        validate="one_to_one",
    )
    out["is_island_5kmplus"] = out["island_id"].notna()
    return out


def add_palette_outcomes(frame: pd.DataFrame) -> pd.DataFrame:
    d = frame.copy()
    for prefix in ("palette_count_", "background_palette_count_"):
        cols = [prefix + c for c in PALETTE]
        missing = set(cols) - set(d.columns)
        if missing:
            raise ValueError(f"RGFCA table missing palette columns: {sorted(missing)}")
        denominator = d[cols].sum(axis=1)
        if (denominator <= 0).any():
            raise ValueError("classifiable rows contain zero palette denominator")
        d[prefix + "white_share"] = d[prefix + "white"] / denominator
        d[prefix + "chromatic_share"] = d[[prefix + c for c in CHROMATIC]].sum(axis=1) / denominator
        d[prefix + "muted_share"] = d[[prefix + c for c in MUTED]].sum(axis=1) / denominator
    for suffix in ("white_share", "chromatic_share", "muted_share"):
        d["flower_" + suffix] = d["palette_count_" + suffix]
        d["background_" + suffix] = d["background_palette_count_" + suffix]
        d["flower_minus_background_" + suffix] = d["flower_" + suffix] - d["background_" + suffix]
    return d


def signflip_p(delta: pd.Series, rng: np.random.Generator, draws: int = 999) -> float:
    values = delta.dropna().to_numpy(float)
    observed = float(values.mean())
    signs = rng.choice((-1.0, 1.0), size=(draws, len(values)))
    null = (signs * values).mean(axis=1)
    return float((1 + np.sum(np.abs(null) >= abs(observed))) / (draws + 1))


def equal_depth(frame: pd.DataFrame, k: int, draws: int, rng: np.random.Generator) -> tuple[list[dict], list[dict]]:
    support = frame.groupby(["species", "is_island_5kmplus"]).size().unstack(fill_value=0)
    support = support.rename(columns={False: "nonisland", True: "island"})
    eligible = support.loc[(support["island"] >= k) & (support["nonisland"] >= k)].index.tolist()
    work = frame.loc[frame["species"].isin(eligible)].copy()
    full = work.groupby(["species", "is_island_5kmplus"])[list(OUTCOMES)].mean().unstack()
    full_rows: list[dict] = []
    for outcome in OUTCOMES:
        delta = full[(outcome, True)] - full[(outcome, False)]
        full_rows.append({
            "analysis": "full_species_equal",
            "k_min_each_side": k,
            "outcome": outcome,
            "n_species": int(len(delta)),
            "estimate_island_minus_nonisland": float(delta.mean()),
            "median_species_delta": float(delta.median()),
            "species_delta_positive_fraction": float((delta > 0).mean()),
            "signflip_p_two_sided": signflip_p(delta, rng),
        })

    aggregate = np.zeros((draws, len(OUTCOMES)), dtype=float)
    for species in eligible:
        group = work.loc[work["species"].eq(species)]
        island = group.loc[group["is_island_5kmplus"], list(OUTCOMES)].to_numpy(float)
        nonisland = group.loc[~group["is_island_5kmplus"], list(OUTCOMES)].to_numpy(float)
        i_idx = np.argpartition(rng.random((draws, len(island))), k - 1, axis=1)[:, :k]
        n_idx = np.argpartition(rng.random((draws, len(nonisland))), k - 1, axis=1)[:, :k]
        aggregate += island[i_idx].mean(axis=1) - nonisland[n_idx].mean(axis=1)
    aggregate /= len(eligible)
    repeated_rows: list[dict] = []
    for j, outcome in enumerate(OUTCOMES):
        x = aggregate[:, j]
        repeated_rows.append({
            "analysis": "repeated_equal_depth",
            "k_each_side": k,
            "outcome": outcome,
            "n_species": int(len(eligible)),
            "draws": int(draws),
            "mean": float(x.mean()),
            "median": float(np.median(x)),
            "q025": float(np.quantile(x, 0.025)),
            "q975": float(np.quantile(x, 0.975)),
            "positive_fraction": float(np.mean(x > 0)),
        })
    return full_rows, repeated_rows


def lat_quarter_match(frame: pd.DataFrame, width: int, draws: int, rng: np.random.Generator) -> tuple[list[dict], list[dict]]:
    d = frame.copy()
    observed = pd.to_datetime(d["observed_on"], errors="coerce")
    d["quarter"] = observed.dt.quarter.astype("Int64")
    d["lat_bin"] = (np.floor((d["latitude"] + 90) / width) * width - 90).astype(int)
    keys = ["species", "lat_bin", "quarter"]
    counts = d.groupby(keys + ["is_island_5kmplus"]).size().unstack(fill_value=0)
    both = counts.loc[(counts.get(True, 0) > 0) & (counts.get(False, 0) > 0)].reset_index()[keys]
    work = d.merge(both, on=keys, how="inner", validate="many_to_many")
    means = work.groupby(keys + ["is_island_5kmplus"])[list(OUTCOMES)].mean().unstack()

    full_rows: list[dict] = []
    for outcome in OUTCOMES:
        cell_delta = means[(outcome, True)] - means[(outcome, False)]
        species_delta = cell_delta.groupby(level="species").mean()
        full_rows.append({
            "analysis": "lat_quarter_matched_full",
            "lat_bin_degrees": width,
            "outcome": outcome,
            "n_species": int(len(species_delta)),
            "n_cells": int(len(cell_delta)),
            "estimate_island_minus_nonisland": float(species_delta.mean()),
            "median_species_delta": float(species_delta.median()),
            "species_delta_positive_fraction": float((species_delta > 0).mean()),
            "signflip_p_two_sided": signflip_p(species_delta, rng),
        })

    cell_groups = []
    for key, group in work.groupby(keys, sort=False):
        island = group.loc[group["is_island_5kmplus"], list(OUTCOMES)].to_numpy(float)
        nonisland = group.loc[~group["is_island_5kmplus"], list(OUTCOMES)].to_numpy(float)
        if len(island) and len(nonisland):
            cell_groups.append((key[0], island, nonisland))
    species = sorted({s for s, _, _ in cell_groups})
    index = {s: i for i, s in enumerate(species)}
    sums = np.zeros((draws, len(species), len(OUTCOMES)), dtype=float)
    n_cells = np.zeros(len(species), dtype=float)
    for species_name, island, nonisland in cell_groups:
        i = index[species_name]
        sums[:, i, :] += island[rng.integers(len(island), size=draws)] - nonisland[rng.integers(len(nonisland), size=draws)]
        n_cells[i] += 1
    aggregate = (sums / n_cells[None, :, None]).mean(axis=1)
    repeated_rows: list[dict] = []
    for j, outcome in enumerate(OUTCOMES):
        x = aggregate[:, j]
        repeated_rows.append({
            "analysis": "lat_quarter_matched_repeated",
            "lat_bin_degrees": width,
            "outcome": outcome,
            "n_species": int(len(species)),
            "n_cells": int(len(cell_groups)),
            "draws": int(draws),
            "mean": float(x.mean()),
            "median": float(np.median(x)),
            "q025": float(np.quantile(x, 0.025)),
            "q975": float(np.quantile(x, 0.975)),
            "positive_fraction": float(np.mean(x > 0)),
        })
    return full_rows, repeated_rows


def main() -> None:
    args = parse_args()
    rng = np.random.default_rng(args.seed)
    source = pd.read_csv(args.rgfca_csv)
    if "global_classifiable" not in source:
        raise ValueError("RGFCA input lacks global_classifiable")
    classified = classify_islands(source, args.islands_gpkg)
    classified = classified.loc[classified["global_classifiable"].eq(True)].copy()
    data = add_palette_outcomes(classified)

    support = data.groupby(["species", "is_island_5kmplus"]).size().unstack(fill_value=0)
    support = support.rename(columns={False: "nonisland", True: "island"})
    support_rows = [
        {"k_each_side": k, "n_species": int(((support["island"] >= k) & (support["nonisland"] >= k)).sum())}
        for k in (1, 2, 3, 5, 10, 20)
    ]

    full_rows: list[dict] = []
    repeated_rows: list[dict] = []
    for k in args.depth_grid:
        a, b = equal_depth(data, k, args.draws, rng)
        full_rows.extend(a)
        repeated_rows.extend(b)

    matched_full: list[dict] = []
    matched_repeated: list[dict] = []
    for width in args.lat_bin_grid:
        a, b = lat_quarter_match(data, width, args.draws, rng)
        matched_full.extend(a)
        matched_repeated.extend(b)

    assemblage_rows = []
    for side in (False, True):
        means = data.loc[data["is_island_5kmplus"].eq(side)].groupby("species")[list(OUTCOMES)].mean()
        for outcome in OUTCOMES:
            assemblage_rows.append({
                "side": "island_5kmplus" if side else "nonisland_reference",
                "outcome": outcome,
                "n_species": int(len(means)),
                "species_equal_mean": float(means[outcome].mean()),
                "species_equal_median": float(means[outcome].median()),
            })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(support_rows).to_csv(args.output_dir / "support.csv", index=False)
    pd.DataFrame(full_rows).to_csv(args.output_dir / "full_species_equal.csv", index=False)
    pd.DataFrame(repeated_rows).to_csv(args.output_dir / "repeated_equal_depth.csv", index=False)
    pd.DataFrame(matched_full).to_csv(args.output_dir / "matched_full.csv", index=False)
    pd.DataFrame(matched_repeated).to_csv(args.output_dir / "matched_repeated.csv", index=False)
    pd.DataFrame(assemblage_rows).to_csv(args.output_dir / "assemblage_species_equal.csv", index=False)

    manifest = {
        "contract": "rgfca_global_island_mainland_exploratory_v1",
        "status": "exploratory_not_canonical_chapter1_inference",
        "classifiable_rows": int(len(data)),
        "island_rows": int(data["is_island_5kmplus"].sum()),
        "nonisland_rows": int((~data["is_island_5kmplus"]).sum()),
        "species": int(data["species"].nunique()),
        "draws": int(args.draws),
        "seed": int(args.seed),
        "island_definition": "exact point-in-polygon against island-v2 universe; >5 km2 and <=7,000,000 km2",
        "nonisland_definition": "complement of the >=5-km2 eligible island polygons; continent-dominated reference that can include <=5-km2 islands",
        "primary_estimand": "equal-species island-minus-nonisland within-species colour difference",
        "matching_sensitivity": "same species x latitude bin x calendar quarter",
        "background_control": "flower palette share minus matched local-background palette share",
        "interpretation_guards": [
            "no causal island-evolution claim from this opened exploratory cohort",
            "repeated-resampling quantiles are stability intervals, not population confidence intervals",
            "within-species contrast is separate from island-flora compositional filtering",
            "nonisland_reference is not silently relabelled strict mainland until <=5-km2 islands are explicitly separated",
            "palette shares are image-derived colour summaries, not calibrated reflectance spectra",
        ],
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
