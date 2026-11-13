"""Exploratory island floral profiles in one shared categorical trait space.

The regional pools below are unions of recorded floras, not reconstructed
ancestral or mainland source pools. Occupancy expectations are descriptive
conditional baselines, not p values or evidence of evolutionary convergence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln

from island_v2.floral_space_geometry import (
    pcoa,
    profile_distances,
    profile_metrics,
    rarefied_profile_richness,
)

CORE_TRAITS = ("flower_primary_color", "floral_form")
PAIR_COLUMNS = [
    "island_a", "island_b", "region_a", "region_b", "n_species_a", "n_species_b",
    "n_shared_species", "species_jaccard", "profile_jaccard",
    "n_nonshared_species_a", "n_nonshared_species_b",
    "profile_jaccard_without_shared_species",
]


def _log_choose(n: np.ndarray | float, k: int) -> np.ndarray:
    n = np.asarray(n, dtype=float)
    out = np.full(n.shape, -np.inf)
    valid = n >= k
    out[valid] = gammaln(n[valid] + 1) - gammaln(k + 1) - gammaln(n[valid] - k + 1)
    return out


def _occupancy_probabilities(counts: np.ndarray, size: int) -> tuple[np.ndarray, np.ndarray]:
    """Hypergeometric probabilities of zero and exactly one carrier of each type."""
    counts = np.asarray(counts, dtype=int)
    total = int(counts.sum())
    if not 0 <= size <= total:
        raise ValueError("Sample size must lie within the finite species pool")
    if size == 0:
        return np.ones(counts.size), np.zeros(counts.size)
    denominator = _log_choose(total, size)
    p0 = np.exp(_log_choose(total - counts, size) - denominator)
    p1 = counts * np.exp(_log_choose(total - counts, size - 1) - denominator)
    return np.clip(p0, 0, 1), np.clip(p1, 0, 1)


def _genus_expectation(
    genus_counts: Counter,
    pool: dict[str, np.ndarray],
    n_profiles: int,
    cache: dict,
    region: str,
) -> tuple[float, float]:
    p0 = np.ones(n_profiles)
    p1 = np.zeros(n_profiles)
    for genus, size in sorted(genus_counts.items()):
        key = (region, genus, size)
        if key not in cache:
            cache[key] = _occupancy_probabilities(pool[genus], size)
        g0, g1 = cache[key]
        p1 = p1 * g0 + p0 * g1
        p0 *= g0
    return float(np.sum(1 - p0)), float(np.sum(p1))


def _complete_profiles(frame: pd.DataFrame, traits: tuple[str, ...]) -> pd.DataFrame:
    if frame["accepted_species"].isna().any() or frame["accepted_species"].duplicated().any():
        raise ValueError("Profiles require one nonmissing row per accepted species")
    missing = set(traits) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing trait columns: {sorted(missing)}")
    out = frame.copy()
    for trait in traits:
        out[trait] = out[trait].fillna("").astype(str).str.strip()
    out = out.loc[out[list(traits)].ne("").all(axis=1)].copy()
    out["profile_key"] = [
        tuple(tuple(sorted(set(value.split("|")))) for value in values)
        for values in out[list(traits)].itertuples(index=False, name=None)
    ]
    out["genus"] = out["accepted_species"].str.split().str[0]
    return out


def _jaccard(a: set, b: set) -> float:
    union = a | b
    return len(a & b) / len(union) if union else float("nan")


def _sample_pairs(ids: list[str], count: int, seed: int) -> list[tuple[str, str]]:
    n = len(ids)
    total = n * (n - 1) // 2
    count = min(count, total)
    if count <= 0:
        return []
    if total <= count:
        return [(ids[i], ids[j]) for i in range(n) for j in range(i + 1, n)]
    rng = np.random.default_rng(seed)
    selected: set[tuple[int, int]] = set()
    while len(selected) < count:
        i, j = sorted(rng.choice(n, size=2, replace=False).tolist())
        selected.add((i, j))
    return [(ids[i], ids[j]) for i, j in sorted(selected)]


def analyze(
    profiles: pd.DataFrame,
    flora: pd.DataFrame,
    geography: pd.DataFrame,
    output: Path,
    *,
    space_profiles: pd.DataFrame | None = None,
    traits: tuple[str, ...] = CORE_TRAITS,
    rarefaction_size: int = 50,
    min_species: int = 50,
    pair_count: int = 2000,
    seed: int = 20261106,
) -> dict:
    """Write one evidence scope; optionally fit coordinates on a shared parent scope."""
    if flora[["island_id", "accepted_species"]].isna().any().any():
        raise ValueError("Flora identifiers must not be missing")
    if flora.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Flora must be unique at island x species grain")
    if geography["island_id"].isna().any() or geography["island_id"].duplicated().any():
        raise ValueError("Geography must be unique at island grain")
    if "analysis_regime" not in geography or geography["analysis_regime"].isna().any():
        raise ValueError("Every geographic unit needs a regional pool label")
    output.mkdir(parents=True, exist_ok=True)
    # The upstream field is called accepted_species but includes a few genus-only
    # and unidentified records. Keep hybrids as named taxon units, without
    # pretending that this lexical gate resolves their accepted species rank.
    original_population = len(profiles)
    uninomial = profiles.accepted_species.str.split().str.len().lt(2)
    unresolved = profiles.accepted_species.str.contains(
        r"(?:^|\s)(?:sp|spp|cf|aff)\.?(?:\s|$)", regex=True, na=False
    )
    exclusions = profiles.loc[uninomial | unresolved, ["accepted_species"]].copy()
    exclusions["reason"] = np.where(uninomial.loc[exclusions.index], "uninomial", "unresolved_name")
    exclusions.to_csv(output / "excluded_taxon_names.csv", index=False)
    profiles = profiles.loc[~(uninomial | unresolved)].copy()
    if space_profiles is not None:
        keep = space_profiles.accepted_species.str.split().str.len().ge(2)
        keep &= ~space_profiles.accepted_species.str.contains(
            r"(?:^|\s)(?:sp|spp|cf|aff)\.?(?:\s|$)", regex=True, na=False
        )
        space_profiles = space_profiles.loc[keep].copy()
    complete = _complete_profiles(profiles, traits)
    reference = _complete_profiles(profiles if space_profiles is None else space_profiles, traits)
    unique = sorted(reference["profile_key"].unique())
    if len(unique) < 2:
        raise ValueError("At least two distinct complete profiles are required for ordination")
    lookup = {key: i for i, key in enumerate(unique)}
    if not set(complete["profile_key"]).issubset(lookup):
        raise ValueError("Analysis profiles must occur in the common reference space")
    complete["profile_id"] = complete["profile_key"].map(lookup)
    distances = profile_distances(unique)
    # The mean set-Jaccard dissimilarity has substantial negative inertia here.
    # Its square root is an explicit distance transformation, not an imputed trait.
    space = pcoa(distances, correction="sqrt")
    diagnostic = pd.DataFrame(space.diagnostics)
    diagnostic.to_csv(output / "ordination_diagnostics.csv", index=False)
    lingoes_space = pcoa(distances, correction="lingoes")
    pd.DataFrame(lingoes_space.diagnostics).to_csv(
        output / "lingoes_ordination_diagnostics.csv", index=False
    )
    profile_table = pd.DataFrame(
        [{"profile_id": i, **{t: "|".join(v) for t, v in zip(traits, key)}}
         for i, key in enumerate(unique)]
    )
    coordinate_table = pd.DataFrame(
        space.coordinates,
        columns=[f"PCoA{i + 1}" for i in range(space.coordinates.shape[1])],
    )
    profile_table = pd.concat([profile_table, coordinate_table], axis=1)
    global_counts = complete["profile_id"].value_counts()
    profile_table["n_species_scope"] = profile_table["profile_id"].map(global_counts).fillna(0).astype(int)
    profile_table.to_csv(output / "profile_coordinates.csv", index=False)
    complete.drop(columns="profile_key").to_csv(output / "species_profiles.csv.gz", index=False)

    geo_columns = [c for c in ["island_id", "analysis_regime", "island_latitude", "island_longitude",
                              "area_km2", "distance_to_continent_km", "spatial_block"]
                   if c in geography]
    valid_flora = flora.merge(geography[geo_columns], on="island_id", how="inner", validate="many_to_one")
    excluded_rows = len(flora) - len(valid_flora)
    applicable = valid_flora.loc[valid_flora["accepted_species"].isin(profiles["accepted_species"])].copy()
    membership = applicable.merge(
        complete[["accepted_species", "profile_id", "genus"]],
        on="accepted_species", how="inner", validate="many_to_one",
    )
    if membership.empty:
        raise ValueError("No complete species match the corrected island flora")
    species_lookup = complete.set_index("accepted_species")["profile_id"].to_dict()
    raw_sizes = valid_flora.groupby("island_id").size()
    applicable_sizes = applicable.groupby("island_id").size()
    full_species_sets = applicable.groupby("island_id")["accepted_species"].agg(set).to_dict()
    island_sets = membership.groupby("island_id")["accepted_species"].agg(set).to_dict()

    regional_pool = membership[["analysis_regime", "accepted_species", "profile_id", "genus"]].drop_duplicates()
    genus_pools: dict[str, dict[str, np.ndarray]] = {}
    regional_counts: dict[str, np.ndarray] = {}
    for region, part in regional_pool.groupby("analysis_regime", sort=True):
        regional_counts[region] = np.bincount(part.profile_id, minlength=len(unique))
        genus_pools[region] = {
            genus: np.bincount(g.profile_id, minlength=len(unique))
            for genus, g in part.groupby("genus", sort=True)
        }
    expectation_cache: dict = {}
    region_expectation_cache: dict = {}
    metric_rows = []
    carrier_rows = []
    for island, group in membership.groupby("island_id", sort=True):
        profile_ids = group.profile_id.to_numpy()
        metrics = profile_metrics(profile_ids)
        region = group.analysis_regime.iloc[0]
        n = len(group)
        rkey = (region, n)
        if rkey not in region_expectation_cache:
            r0, r1 = _occupancy_probabilities(regional_counts[region], n)
            region_expectation_cache[rkey] = (float(np.sum(1 - r0)), float(np.sum(r1)))
        expected_region, _ = region_expectation_cache[rkey]
        expected_genus, singleton_genus = _genus_expectation(
            Counter(group.genus), genus_pools[region], len(unique), expectation_cache, region,
        )
        row = {
            "island_id": island,
            **{c: group[c].iloc[0] for c in geo_columns if c != "island_id"},
            "n_flora_species": int(raw_sizes[island]),
            "n_applicable_species": int(applicable_sizes[island]),
            **metrics,
            "joint_coverage": n / applicable_sizes[island],
            "rarefaction_size": rarefaction_size,
            "rarefied_profiles": rarefied_profile_richness(profile_ids, rarefaction_size)
            if n >= rarefaction_size else np.nan,
            "expected_profiles_regional_pool": expected_region,
            "expected_profiles_genus_pool": expected_genus,
            "expected_singleton_profiles_genus_pool": singleton_genus,
            "profile_ratio_to_regional_pool": metrics["n_profiles"] / expected_region,
            "profile_ratio_to_genus_pool": metrics["n_profiles"] / expected_genus,
            "profile_difference_from_genus_pool": metrics["n_profiles"] - expected_genus,
        }
        metric_rows.append(row)
        carrier_rows.extend({"island_id": island, "profile_id": p, "n_species": count}
                            for p, count in sorted(Counter(profile_ids).items()))
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(output / "island_metrics.csv", index=False)
    pd.DataFrame(carrier_rows).to_csv(output / "island_profile_carriers.csv.gz", index=False)
    # Explicitly retain excluded/no-data islands as coverage rows, not phenotype zeros.
    coverage = geography[geo_columns].merge(
        pd.concat([raw_sizes.rename("n_flora_species"), applicable_sizes.rename("n_applicable_species")], axis=1),
        on="island_id", how="left", validate="one_to_one",
    ).merge(metrics[["island_id", "n_species", "joint_coverage"]],
            on="island_id", how="left", validate="one_to_one")
    coverage["data_status"] = np.where(coverage.n_flora_species.isna(), "no_flora_records",
                                      np.where(coverage.n_species.isna(), "no_complete_profiles", "observed"))
    coverage.to_csv(output / "island_coverage.csv", index=False)

    eligible = metrics.loc[metrics.n_species >= min_species].sort_values("island_id")
    region_lookup = metrics.set_index("island_id")["analysis_regime"].to_dict()
    pair_rows = []
    for a, b in _sample_pairs(eligible.island_id.tolist(), pair_count, seed):
        sa, sb = island_sets[a], island_sets[b]
        common = sa & sb
        pa, pb = {species_lookup[s] for s in sa}, {species_lookup[s] for s in sb}
        na, nb = sa - common, sb - common
        ua, ub = {species_lookup[s] for s in na}, {species_lookup[s] for s in nb}
        pair_rows.append({
            "island_a": a, "island_b": b, "region_a": region_lookup[a], "region_b": region_lookup[b],
            "n_species_a": len(sa), "n_species_b": len(sb), "n_shared_species": len(common),
            "species_jaccard": _jaccard(sa, sb), "profile_jaccard": _jaccard(pa, pb),
            "all_applicable_species_jaccard": _jaccard(full_species_sets[a], full_species_sets[b]),
            "n_nonshared_species_a": len(na), "n_nonshared_species_b": len(nb),
            "profile_jaccard_without_shared_species": _jaccard(ua, ub)
            if na and nb else np.nan,
        })
    pairs = pd.DataFrame(pair_rows, columns=PAIR_COLUMNS + ["all_applicable_species_jaccard"])
    pairs.to_csv(output / "pairwise_comparisons.csv", index=False)
    regional_rows = []
    for region, part in eligible.groupby("analysis_regime", sort=True):
        regional_rows.append({
            "analysis_regime": region, "n_islands": len(part),
            **{f"median_{c}": float(part[c].median()) for c in [
                "n_species", "n_profiles", "singleton_fraction", "mean_redundancy", "joint_coverage",
                "rarefied_profiles", "profile_ratio_to_genus_pool", "profile_ratio_to_regional_pool",
            ]},
        })
    pd.DataFrame(regional_rows).to_csv(output / "regional_summary.csv", index=False)
    summary = {
        "study_status": "exploratory_descriptive", "traits": list(traits),
        "upstream_population_name_units": original_population,
        "n_excluded_uninomial_names": int(uninomial.sum()),
        "n_excluded_unresolved_names": int((unresolved & ~uninomial).sum()),
        "n_profile_population_species": len(profiles),
        "n_complete_species": len(complete), "n_complete_species_in_flora": membership.accepted_species.nunique(),
        "n_profiles": int(complete.profile_id.nunique()),
        "n_reference_profiles": len(unique), "n_islands_with_complete_profiles": len(metrics),
        "n_geographic_units": len(geography), "geography_excluded_flora_rows": excluded_rows,
        "n_islands_min_species": len(eligible), "min_species": min_species,
        "rarefaction_size": rarefaction_size, "pair_count": len(pairs), "seed": seed,
        "pcoa": {
            "correction": "sqrt", "distance": "equal_weight_mean_set_jaccard",
            "negative_inertia_fraction": space.raw_negative_inertia_fraction,
            "lingoes_constant": space.lingoes_constant, "positive_dimensions": len(space.eigenvalues),
            "two_dimensional_raw_rmse": float(diagnostic.loc[diagnostic.dimension.eq(2), "raw_rmse"].iloc[0])
            if len(diagnostic) >= 2 else None,
            "minimum_raw_rmse_dimension": int(diagnostic.loc[diagnostic.raw_rmse.idxmin(), "dimension"]),
            "two_dimensional_inertia_fraction": float(
                diagnostic.loc[diagnostic.dimension.eq(2), "positive_inertia_fraction"].iloc[0]
            ) if len(diagnostic) >= 2 else 1.0,
            "dimensions_for_80_percent_inertia": int(
                diagnostic.loc[diagnostic.positive_inertia_fraction.ge(0.8), "dimension"].iloc[0]
            ),
            "dimensions_for_95_percent_inertia": int(
                diagnostic.loc[diagnostic.positive_inertia_fraction.ge(0.95), "dimension"].iloc[0]
            ),
        },
        "baseline": "finite sampling without replacement from the recorded species union within each broad latitudinal regime; genus counts fixed in genus baseline",
        "limitations": [
            "Opportunistic realized floras, not censuses or absence data",
            "Recorded species state sets are not joint flower-level morph observations",
            "Regional recorded pools are not historical colonist source pools",
            "Profile overlap is descriptive; categories mechanically increase similarity",
            "Pair comparisons share islands and are not independent replicates",
            "All includes validated genus-level fill and uneven direct-evidence coverage",
            "Accepted-species keys are taxon-name units including hybrid formulas; exact species rank was not re-resolved",
        ],
    }
    if len(eligible):
        summary["eligible_medians"] = {
            c: float(eligible[c].median()) for c in [
                "n_species", "n_profiles", "singleton_fraction", "mean_redundancy", "joint_coverage",
                "rarefied_profiles", "profile_ratio_to_genus_pool", "profile_ratio_to_regional_pool",
            ]
        }
    if len(pairs):
        summary["pair_medians"] = {c: float(pairs[c].median()) for c in [
            "species_jaccard", "profile_jaccard", "profile_jaccard_without_shared_species",
        ]}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--space-profiles", type=Path)
    parser.add_argument("--flora", type=Path, required=True)
    parser.add_argument("--geography", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--traits", nargs="+", default=list(CORE_TRAITS))
    parser.add_argument("--rarefaction-size", type=int, default=50)
    parser.add_argument("--min-species", type=int, default=50)
    parser.add_argument("--pair-count", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    summary = analyze(
        pd.read_csv(args.profiles), pd.read_csv(args.flora), pd.read_csv(args.geography), args.output,
        space_profiles=pd.read_csv(args.space_profiles) if args.space_profiles else None,
        traits=tuple(args.traits), rarefaction_size=args.rarefaction_size,
        min_species=args.min_species, pair_count=args.pair_count, seed=args.seed,
    )
    input_files = {key: value for key, value in vars(args).items()
                   if isinstance(value, Path) and key != "output" and value.is_file()}
    manifest = {key: {"path": str(value), "sha256": hashlib.sha256(value.read_bytes()).hexdigest()}
                for key, value in input_files.items()}
    (args.output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(summary, ensure_ascii=False))


if __name__ == "__main__":
    main()
