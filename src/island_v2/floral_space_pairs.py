"""Pair profile sharing after removing shared taxa, under a joint regional null.

Pair selection is supplied by the caller and never uses profile overlap here.
Every draw permutes the full regional within-genus pools once, then scores all
pair sides with that assignment. Recorded shared taxa are excluded before any
profile lookup; exclusive taxa carrying those same profile codes remain.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from island_v2.floral_space_null import (
    RegionPermutationDesign,
    draw_profile_assignment,
    finite_two_sided_p,
)


@dataclass(frozen=True)
class PairPermutationDesign:
    """Prepared exclusive membership using original regional taxon/profile IDs.

    ``pairs`` retains input row order, resets the row index, canonicalizes the
    endpoints and adds total, shared and nonshared taxon counts. Membership side
    ``2*i`` is endpoint a of row i; ``2*i+1`` is endpoint b. ``region`` retains the
    entire pool, including taxa on ineligible or unselected islands.
    """

    region: RegionPermutationDesign
    pairs: pd.DataFrame
    member_side_indices: np.ndarray
    member_taxon_indices: np.ndarray

    @property
    def n_pairs(self) -> int:
        return len(self.pairs)


def prepare_pair_design(
    region_design: RegionPermutationDesign, pairs: pd.DataFrame, *, min_unique: int = 20,
) -> PairPermutationDesign:
    """Validate fixed pairs and remove shared names before profile assignment.

    Empty, repeated unordered, self, missing and unknown/ineligible endpoints
    raise ValueError. Both exclusive sides must have at least ``min_unique``
    taxa. Input metadata columns are retained; derived counts are recomputed.
    """
    if not isinstance(min_unique, int) or isinstance(min_unique, bool) or min_unique < 1:
        raise ValueError("A positive integer minimum unique-taxon support is required")
    if not {"island_a", "island_b"}.issubset(pairs.columns):
        raise ValueError("Pair columns must include island_a and island_b")
    if pairs.empty:
        raise ValueError("Pair selection must not be empty")
    endpoints = pairs[["island_a", "island_b"]]
    if endpoints.isna().any().any():
        raise ValueError("Pair endpoints must not be missing")
    eligible = set(region_design.island_ids)
    if not set(endpoints.to_numpy().ravel()).issubset(eligible):
        raise ValueError("Pair endpoints contain unknown or ineligible islands")
    canonical = [tuple(sorted(pair)) for pair in endpoints.itertuples(index=False, name=None)]
    if any(a == b for a, b in canonical):
        raise ValueError("Pair selection must not contain self pairs")
    if len(set(canonical)) != len(canonical):
        raise ValueError("Pair selection contains repeated unordered pairs")
    table = pairs.copy().reset_index(drop=True)
    table[["island_a", "island_b"]] = pd.DataFrame(canonical, index=table.index)
    members = [set() for _ in region_design.island_ids]
    for island, taxon in zip(region_design.member_island_indices,
                             region_design.member_taxon_indices, strict=True):
        members[island].add(int(taxon))
    lookup = dict(zip(region_design.island_ids, members, strict=True))
    sides, taxa, counts = [], [], []
    for index, (a, b) in enumerate(canonical):
        shared = lookup[a] & lookup[b]
        exclusive_a, exclusive_b = sorted(lookup[a] - shared), sorted(lookup[b] - shared)
        if min(len(exclusive_a), len(exclusive_b)) < min_unique:
            raise ValueError(f"Pair {a!r}, {b!r} has insufficient nonshared taxon support")
        counts.append((len(lookup[a]), len(lookup[b]), len(shared),
                       len(exclusive_a), len(exclusive_b)))
        sides.extend([2 * index] * len(exclusive_a))
        sides.extend([2 * index + 1] * len(exclusive_b))
        taxa.extend(exclusive_a)
        taxa.extend(exclusive_b)
    count_columns = ["n_taxa_a", "n_taxa_b", "n_shared_taxa", "n_unique_a", "n_unique_b"]
    table[count_columns] = pd.DataFrame(counts, index=table.index)
    return PairPermutationDesign(region=region_design, pairs=table,
                                 member_side_indices=np.asarray(sides, dtype=np.int64),
                                 member_taxon_indices=np.asarray(taxa, dtype=np.int64))


def pair_overlap_metrics(design: PairPermutationDesign, assignment: np.ndarray) -> dict[str, np.ndarray]:
    """Score every fixed pair with a shared regional taxon-to-profile assignment.

    The temporary boolean occupancy matrix has 2*n_pairs rows and the original
    regional profile-code width; no pair-by-taxon-by-draw tensor is constructed.
    Positive nonshared support guarantees a positive union for every pair.
    """
    values = np.asarray(assignment)
    if values.shape != design.region.profile_ids.shape or not np.issubdtype(values.dtype, np.integer):
        raise ValueError("Assignment must contain one integer profile ID per regional taxon")
    if (values < 0).any() or (values >= design.region.n_profile_codes).any():
        raise ValueError("Assignment contains a profile outside the prepared reference codes")
    occupied = np.zeros((2 * design.n_pairs, design.region.n_profile_codes), dtype=bool)
    occupied[design.member_side_indices, values[design.member_taxon_indices]] = True
    a, b = occupied[0::2], occupied[1::2]
    intersection = np.count_nonzero(a & b, axis=1)
    union = np.count_nonzero(a | b, axis=1)
    return {
        "n_profiles_a": np.count_nonzero(a, axis=1),
        "n_profiles_b": np.count_nonzero(b, axis=1),
        "profile_intersection": intersection,
        "profile_union": union,
        "profile_jaccard": intersection / union,
    }


def run_pair_comparison(
    design: PairPermutationDesign, *, permutations: int = 999, seed=20261106,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    """Return finite-tail summary, per-draw medians and per-pair aggregates.

    Pair Jaccard mean and population SD (ddof=0) use streamed Welford moments.
    Only five comparison medians per draw are retained. Identities and repeated
    assignments are included. ``n_islands`` counts sampled pair endpoints.
    A zero null median is a valid run with status ``ok`` and an explicitly
    undefined observed/null ratio, while retaining the absolute difference.
    """
    if not isinstance(permutations, int) or isinstance(permutations, bool) or permutations < 1:
        raise ValueError("A positive integer permutation count is required")
    region = design.region
    observed = pair_overlap_metrics(design, region.profile_ids)
    metrics = ("profile_jaccard", "profile_intersection", "profile_union",
               "n_profiles_a", "n_profiles_b")
    draws = np.empty((permutations, len(metrics)))
    mean, m2 = np.zeros(design.n_pairs), np.zeros(design.n_pairs)
    rng = np.random.default_rng(seed)
    for index in range(permutations):
        result = pair_overlap_metrics(design, draw_profile_assignment(region, rng))
        draws[index] = [np.median(result[metric]) for metric in metrics]
        delta = result["profile_jaccard"] - mean
        mean += delta / (index + 1)
        m2 += delta * (result["profile_jaccard"] - mean)
    observed_median = float(np.median(observed["profile_jaccard"]))
    inference = finite_two_sided_p(observed_median, draws[:, 0])
    zero_null_median = inference["null_median"] == 0
    summary = {
        "analysis_regime": region.analysis_regime, "status": "ok",
        "n_pairs": design.n_pairs,
        "n_islands": len(set(design.pairs.island_a) | set(design.pairs.island_b)),
        "n_taxa_pool": len(region.taxon_names), "n_genera": len(region.genus_indices),
        "n_permutable_genera": len(region.permutable_genus_indices),
        "n_permutations": permutations, "observed_median_jaccard": observed_median,
        **inference,
        "observed_minus_null_median": observed_median - inference["null_median"],
        "observed_to_null_median": (np.nan if zero_null_median
                                    else observed_median / inference["null_median"]),
        "null_median_is_zero": zero_null_median,
        "ratio_status": "undefined_zero_null_median" if zero_null_median else "ok",
        "null_min": float(np.min(draws[:, 0])), "null_max": float(np.max(draws[:, 0])),
        "null_n_distinct_sampled_values": int(np.unique(draws[:, 0]).size),
    }
    for index, metric in enumerate(metrics[1:], start=1):
        summary[f"observed_median_{metric}"] = float(np.median(observed[metric]))
        summary[f"null_median_{metric}"] = float(np.median(draws[:, index]))
    draw_table = pd.DataFrame(draws, columns=[f"median_{metric}" for metric in metrics])
    draw_table.insert(0, "permutation", np.arange(1, permutations + 1))
    pair_table = design.pairs.copy()
    for metric, values in observed.items():
        pair_table[metric] = values
    pair_table["null_mean_jaccard"] = mean
    pair_table["null_sd_jaccard"] = np.sqrt(np.maximum(m2 / permutations, 0))
    return summary, draw_table, pair_table
