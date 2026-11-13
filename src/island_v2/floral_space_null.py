"""Joint, genus-conditioned randomization of recorded floral profiles.

The observed taxon-by-island incidence is fixed. Within each latitude regime,
profile labels are permuted among the complete recorded taxa within a genus.
One taxon receives one assignment across every island in that regime. Regional
assignments are separate conditional references, not reconstructed global traits.

This calibrates the exploratory median profile-count/expectation ratio against
the same operational genus pools used by ``floral_space``. It does not test a
historical colonization model, evolutionary convergence or causal selection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from island_v2.floral_space import _genus_expectation

SCOPES = ("all", "direct", "all_single_state", "direct_single_state")
REGIMES = ("northern_midlatitude", "northern_high_latitude", "tropical", "southern_extratropical")


@dataclass(frozen=True)
class RegionPermutationDesign:
    analysis_regime: str
    taxon_names: tuple[str, ...]
    profile_ids: np.ndarray
    genus_indices: tuple[np.ndarray, ...]
    permutable_genus_indices: tuple[np.ndarray, ...]
    island_ids: tuple[str, ...]
    member_island_indices: np.ndarray
    member_taxon_indices: np.ndarray
    expected_profiles: np.ndarray
    observed_profile_counts: np.ndarray
    n_profile_codes: int
    maximum_expectation_difference: float


def build_membership(
    species_profiles: pd.DataFrame, flora: pd.DataFrame, geography: pd.DataFrame,
) -> pd.DataFrame:
    """Join the complete scope to corrected geography without restricting islands."""
    columns = ["accepted_species", "genus", "profile_id"]
    if species_profiles[columns].isna().any().any():
        raise ValueError("Complete profiles must have taxon, genus and profile identifiers")
    if species_profiles.accepted_species.duplicated().any():
        raise ValueError("Species profiles must have a unique taxon key")
    if flora[["island_id", "accepted_species"]].isna().any().any():
        raise ValueError("Flora identifiers must not be missing")
    if flora.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Flora membership must be unique at island by taxon grain")
    if geography[["island_id", "analysis_regime"]].isna().any().any():
        raise ValueError("Geography must have nonmissing island and regional identifiers")
    if geography.island_id.duplicated().any():
        raise ValueError("Geography must have a unique island key")
    unknown_regimes = set(geography.analysis_regime) - set(REGIMES)
    if unknown_regimes:
        raise ValueError(f"Unexpected analysis regime labels: {sorted(unknown_regimes)}")
    profiles = species_profiles[columns].copy()
    numeric = pd.to_numeric(profiles.profile_id, errors="raise").to_numpy(float)
    if not np.isfinite(numeric).all() or (numeric < 0).any() or (numeric != np.floor(numeric)).any():
        raise ValueError("Profile identifiers must be nonnegative integers")
    profiles["profile_id"] = numeric.astype(np.int64)
    return flora[["island_id", "accepted_species"]].merge(
        geography[["island_id", "analysis_regime"]], on="island_id", how="inner",
        validate="many_to_one",
    ).merge(profiles, on="accepted_species", how="inner", validate="many_to_one")


def prepare_region(
    membership: pd.DataFrame, metrics: pd.DataFrame, regime: str, *, min_species: int = 50,
) -> RegionPermutationDesign:
    """Prepare a regional joint permutation and verify its frozen denominators."""
    if metrics.island_id.duplicated().any():
        raise ValueError("Island metrics must have unique island identifiers")
    part = membership.loc[membership.analysis_regime.eq(regime)].copy()
    if part.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Regional membership must remain unique at island by taxon grain")
    # All recorded complete taxa enter the pool, including those found only on
    # islands that do not qualify for the statistic below.
    pool = part[["accepted_species", "genus", "profile_id"]].drop_duplicates()
    if pool.accepted_species.duplicated().any():
        raise ValueError("A regional taxon must have one genus and one profile")
    pool = pool.sort_values(["genus", "accepted_species"]).reset_index(drop=True)
    selected = metrics.loc[
        metrics.analysis_regime.eq(regime) & metrics.n_species.ge(min_species)
    ].sort_values("island_id")
    taxon_names = tuple(pool.accepted_species)
    island_ids = tuple(selected.island_id)
    taxon_lookup = {name: i for i, name in enumerate(taxon_names)}
    island_lookup = {name: i for i, name in enumerate(island_ids)}
    eligible_membership = part.loc[part.island_id.isin(island_ids)].sort_values(
        ["island_id", "accepted_species"]
    )
    member_islands = eligible_membership.island_id.map(island_lookup).to_numpy(np.int64)
    member_taxa = eligible_membership.accepted_species.map(taxon_lookup).to_numpy(np.int64)
    profile_ids = pool.profile_id.to_numpy(np.int64)
    n_codes = int(profile_ids.max()) + 1 if len(profile_ids) else 0
    genus_indices = tuple(np.asarray(index, dtype=np.int64)
                          for index in pool.groupby("genus", sort=True).indices.values())
    variable_groups = tuple(index for index in genus_indices
                            if np.unique(profile_ids[index]).size > 1)
    expected = selected.expected_profiles_genus_pool.to_numpy(float)
    if not np.isfinite(expected).all() or (expected <= 0).any():
        raise ValueError("Analytic expectations must be positive and finite")
    observed_species = np.bincount(member_islands, minlength=len(island_ids))
    if not np.array_equal(observed_species, selected.n_species.to_numpy()):
        raise ValueError("Observed membership species counts disagree with saved metrics")

    # Rebuild the analytic baseline once, not on every permutation. Within-genus
    # profile multisets are invariant under the randomization, so E stays fixed.
    genus_pools = {
        genus: np.bincount(group.profile_id, minlength=n_codes)
        for genus, group in pool.groupby("genus", sort=True)
    }
    recomputed = []
    cache: dict = {}
    for _, group in eligible_membership.groupby("island_id", sort=True):
        value, _ = _genus_expectation(Counter(group.genus), genus_pools, n_codes, cache, regime)
        recomputed.append(value)
    recomputed_array = np.asarray(recomputed)
    if not np.allclose(recomputed_array, expected, atol=1e-8, rtol=1e-8):
        raise ValueError("Saved analytic genus expectation disagrees with this regional pool")
    maximum_error = float(np.max(np.abs(recomputed_array - expected))) if len(expected) else 0.0
    design = RegionPermutationDesign(
        analysis_regime=regime, taxon_names=taxon_names, profile_ids=profile_ids,
        genus_indices=genus_indices, permutable_genus_indices=variable_groups,
        island_ids=island_ids, member_island_indices=member_islands,
        member_taxon_indices=member_taxa, expected_profiles=expected,
        observed_profile_counts=selected.n_profiles.to_numpy(np.int64),
        n_profile_codes=n_codes, maximum_expectation_difference=maximum_error,
    )
    if not np.array_equal(island_profile_counts(design, profile_ids), design.observed_profile_counts):
        raise ValueError("Observed profile counts disagree with saved metrics")
    return design


def draw_profile_assignment(design: RegionPermutationDesign, rng: np.random.Generator) -> np.ndarray:
    """Draw one regional assignment, retaining identities and duplicate draws."""
    assigned = design.profile_ids.copy()
    for index in design.permutable_genus_indices:
        assigned[index] = rng.permutation(design.profile_ids[index])
    return assigned


def island_profile_counts(design: RegionPermutationDesign, assignment: np.ndarray) -> np.ndarray:
    """Count unique profiles with one taxon assignment shared across islands."""
    values = np.asarray(assignment)
    if values.shape != design.profile_ids.shape or not np.issubdtype(values.dtype, np.integer):
        raise ValueError("Assignment must contain one integer profile ID per regional taxon")
    if (values < 0).any() or (values >= design.n_profile_codes).any():
        raise ValueError("Assignment contains a profile outside the prepared reference codes")
    occupied = np.zeros((len(design.island_ids), design.n_profile_codes), dtype=bool)
    occupied[design.member_island_indices, values[design.member_taxon_indices]] = True
    return np.count_nonzero(occupied, axis=1)


def median_profile_ratio(profile_counts, expected_profiles) -> float:
    counts, expected = np.asarray(profile_counts, dtype=float), np.asarray(expected_profiles, dtype=float)
    if counts.ndim != 1 or counts.shape != expected.shape or not len(counts):
        raise ValueError("Observed and expected profile counts must be aligned nonempty vectors")
    if not np.isfinite(counts).all() or not np.isfinite(expected).all() or (expected <= 0).any():
        raise ValueError("Profile counts must be finite and expectations positive")
    return float(np.median(counts / expected))


def finite_two_sided_p(observed: float, null_draws) -> dict[str, float]:
    """Twice the smaller inclusive finite-corrected tail, capped at one.

    Numerical ties within 100 machine epsilons times max(1, abs(observed)) count
    in both tails. A range of null draws is not a confidence interval.
    """
    draws = np.asarray(null_draws, dtype=float)
    if draws.ndim != 1 or not len(draws) or not np.isfinite(draws).all() or not np.isfinite(observed):
        raise ValueError("A finite observed statistic and nonempty finite null draws are required")
    tolerance = 100 * np.finfo(float).eps * max(1.0, abs(observed))
    lower = float((1 + np.count_nonzero(draws <= observed + tolerance)) / (len(draws) + 1))
    upper = float((1 + np.count_nonzero(draws >= observed - tolerance)) / (len(draws) + 1))
    return {
        "p_lower": lower, "p_upper": upper, "p_two_sided": min(1.0, 2 * min(lower, upper)),
        "null_median": float(np.median(draws)), "null_q025": float(np.quantile(draws, 0.025)),
        "null_q975": float(np.quantile(draws, 0.975)), "null_sd": float(np.std(draws)),
    }


def bh_adjust_fixed(pvalues) -> np.ndarray:
    """BH adjustment with missing comparisons retained in the family as p=1."""
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1 or not len(p) or np.any((p < 0) | (p > 1)):
        raise ValueError("P values must be a nonempty vector between zero and one, or NaN")
    finite = np.isfinite(p)
    safe = np.where(finite, p, 1.0)
    order = np.argsort(safe, kind="stable")
    adjusted = safe[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum(1.0, np.minimum.accumulate(adjusted[::-1])[::-1])
    result = np.empty(len(p))
    result[order] = adjusted
    result[~finite] = np.nan
    return result


def run_comparison(
    design: RegionPermutationDesign, *, permutations: int = 999, seed=20261106,
) -> tuple[dict, np.ndarray]:
    if not isinstance(permutations, int) or permutations < 1:
        raise ValueError("A positive integer permutation count is required")
    base = {
        "analysis_regime": design.analysis_regime, "n_islands": len(design.island_ids),
        "n_taxa": len(design.taxon_names), "n_genera": len(design.genus_indices),
        "n_permutable_genera": len(design.permutable_genus_indices),
        "n_pool_only_taxa": len(design.taxon_names) - len(np.unique(design.member_taxon_indices)),
        "maximum_expectation_difference": design.maximum_expectation_difference,
    }
    if not design.island_ids:
        return dict(base, status="no_eligible_islands", n_permutations=0,
                    observed_ratio=np.nan, p_two_sided=np.nan), np.empty(0)
    observed = median_profile_ratio(design.observed_profile_counts, design.expected_profiles)
    rng = np.random.default_rng(seed)
    draws = np.empty(permutations)
    for index in range(permutations):
        assignment = draw_profile_assignment(design, rng)
        draws[index] = median_profile_ratio(island_profile_counts(design, assignment),
                                           design.expected_profiles)
    inference = finite_two_sided_p(observed, draws)
    delta = observed - inference["null_median"]
    tolerance = 100 * np.finfo(float).eps * max(1.0, abs(observed))
    direction = "enriched" if delta > tolerance else "depleted" if delta < -tolerance else "equal"
    return dict(base, status="ok", n_permutations=permutations, observed_ratio=observed,
                **inference, observed_minus_null_median=delta, direction=direction), draws


def _file_record(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return {"path": str(path), "sha256": digest.hexdigest(), "bytes": path.stat().st_size}


def run_suite(results: Path, flora_path: Path, geography_path: Path, *,
              permutations: int = 999, min_species: int = 50, seed: int = 20261106) -> dict:
    """Run all sixteen comparisons and save aggregate-only outputs."""
    flora = pd.read_csv(flora_path, usecols=["island_id", "accepted_species"], dtype=str)
    geography = pd.read_csv(geography_path, usecols=["island_id", "analysis_regime"], dtype=str)
    sources = {"flora": _file_record(flora_path), "geography": _file_record(geography_path)}
    rows, draw_tables = [], []
    for scope_index, scope in enumerate(SCOPES):
        profile_path, metric_path = results / scope / "species_profiles.csv.gz", results / scope / "island_metrics.csv"
        profiles = pd.read_csv(profile_path, usecols=["accepted_species", "genus", "profile_id"])
        metrics = pd.read_csv(metric_path)
        sources[scope] = {"profiles": _file_record(profile_path), "island_metrics": _file_record(metric_path)}
        membership = build_membership(profiles, flora, geography)
        for region_index, region in enumerate(REGIMES):
            design = prepare_region(membership, metrics, region, min_species=min_species)
            stream = np.random.SeedSequence([seed, scope_index, region_index])
            row, draws = run_comparison(design, permutations=permutations, seed=stream)
            rows.append({"evidence_scope": scope, **row})
            draw_tables.append(pd.DataFrame({
                "evidence_scope": [scope] * len(draws), "analysis_regime": [region] * len(draws),
                "permutation": np.arange(1, len(draws) + 1), "median_profile_ratio": draws,
            }))
            print(f"{scope}/{region}: n={row['n_islands']} taxa={row['n_taxa']} "
                  f"T={row.get('observed_ratio', np.nan):.6f} "
                  f"null={row.get('null_median', np.nan):.6f} "
                  f"p={row['p_two_sided']:.4g}", flush=True)
    summary = pd.DataFrame(rows)
    assert len(summary) == len(SCOPES) * len(REGIMES) == 16
    summary["q_bh_16"] = bh_adjust_fixed(summary.p_two_sided)
    summary.to_csv(results / "permutation_summary.csv", index=False)
    pd.concat(draw_tables, ignore_index=True).to_csv(
        results / "permutation_null_draws.csv.gz", index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    manifest = {
        "schema_version": "floral_profile_joint_genus_permutation_v1", "status": "complete",
        "design_status": "exploratory conditional assemblage comparison", "fixed_comparison_family_size": 16,
        "evidence_scopes": list(SCOPES), "analysis_regimes": list(REGIMES),
        "permutations_requested_per_comparison": permutations, "minimum_complete_name_units_per_island": min_species,
        "seed": seed, "rng": "NumPy default_rng(SeedSequence([seed, scope_index, region_index]))",
        "statistic": "median across eligible islands of n_profiles / fixed analytic expected_profiles_genus_pool",
        "conditioning": [
            "observed island-by-taxon membership and hence island sample sizes and genus composition",
            "within-genus profile multisets in the union of all complete taxa recorded in each regime",
            "one permuted profile per taxon shared across all islands of the same regime",
            "pool includes complete taxa found only on islands below the eligibility threshold",
        ],
        "regional_interpretation": "independent conditional reference distributions; not a reassignment of global species traits",
        "tail_definition": "p_lower=(1+sum(T_perm<=T_obs+tol))/(B+1); p_upper=(1+sum(T_perm>=T_obs-tol))/(B+1); two-sided=min(1,2*min(p_lower,p_upper))",
        "tie_tolerance": "100 * machine_epsilon * max(1, abs(T_obs)); ties in both tails",
        "null_range": "2.5 and 97.5 percentiles of permutation statistics, not a confidence interval",
        "multiple_testing": "BH across fixed 16 comparisons; absent comparisons retain a family slot",
        "claim_limits": [
            "conditional descriptive assemblage association, not historical or causal selection",
            "observational coverage and source-dependent profile assignments remain part of the input",
            "overlapping evidence scopes are not independent replications",
            "BH adjustment is exploratory; arbitrary-dependence FDR control is not asserted",
        ],
        "sources": sources,
        "code": _file_record(Path(__file__)),
        "outputs": {
            name: _file_record(results / name)
            for name in ["permutation_summary.csv", "permutation_null_draws.csv.gz"]
        },
        "method_references": [
            "https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.permutation_test.html",
            "https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.false_discovery_control.html",
        ],
    }
    (results / "permutation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--flora", type=Path, required=True)
    parser.add_argument("--geography", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--min-species", type=int, default=50)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    run_suite(args.results, args.flora, args.geography, permutations=args.permutations,
              min_species=args.min_species, seed=args.seed)


if __name__ == "__main__":
    main()
