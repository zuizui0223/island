"""Private count-only functional turnover with sequential equal-carrier rarefaction.

Produces broad-region aggregates only, not independent spatial hypothesis tests.
Source islands/species/genera and all row-level counts stay in protected memory.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .functional_turnover import TRAITS

STAGES = (
    ("all", "species_exclusive"),
    ("species_exclusive", "genus_exclusive"),
)
KEYS = ["scope", "region", "island_a", "island_b"]


def _int_counts(frame: pd.DataFrame, key: str) -> np.ndarray:
    raw = pd.to_numeric(frame[key], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw).all() or (raw < 0).any() or (raw != np.floor(raw)).any():
        raise ValueError(f"Invalid integer carrier count: {key}")
    return raw.astype(np.int64)


def audited_rarefaction(
    private_counts: pd.DataFrame,
    private_blocks: pd.DataFrame,
    *,
    draws: int = 999,
    seed: int = 20261110,
    min_carriers: int = 20,
    min_pairs: int = 40,
    min_islands: int = 40,
    min_spatial_blocks: int = 10,
    traits: tuple[str, ...] = TRAITS,
) -> pd.DataFrame:
    """Evaluate strict binary trait prevalence differences after taxon removal.

    The private stage records have columns scope, region, island_a, island_b,
    trait_scope, n_TRAIT, n_TRAIT_b, k_TRAIT_a, k_TRAIT_b.  "n" is the
    observed binary-trait carrier count and "k" the number in state 1.
    All three stages must have identical, unique pair keys.  Only the
    broad-region grouped output is safe to share after rights review.

    For each pair, independently thin eligible baseline carriers of each
    island to the trait-known sample sizes surviving taxonomic exclusion.
    Hypergeometric expectation is a sample-size diagnostic, NOT a
    genus/source-conditioned evolutionary null.
    """
    if draws < 100 or min_carriers < 1 or min_pairs < 5 or min_islands < 5 or min_spatial_blocks < 5:
        raise ValueError("Invalid draws or report support gate")
    if not set(traits).issubset(TRAITS) or len(traits) != len(set(traits)):
        raise ValueError("Traits must be distinct predefined binary traits")
    required = set(KEYS + ["trait_scope"])
    for trait in traits:
        required.update((f"n_{trait}", f"n_{trait}_b", f"k_{trait}_a", f"k_{trait}_b"))
    if not required.issubset(private_counts.columns):
        raise ValueError(f"Missing private count fields: {sorted(required-set(private_counts.columns))}")
    if private_counts[KEYS + ["trait_scope"]].isna().any().any():
        raise ValueError("Missing cohort or stage identity")
    if (private_counts.island_a == private_counts.island_b).any():
        raise ValueError("Self-pair is invalid")
    if private_counts.duplicated(KEYS + ["trait_scope"]).any():
        raise ValueError("Duplicate island-pair-stage record")
    if not set(private_counts.trait_scope).issubset({"all", "species_exclusive", "genus_exclusive"}):
        raise ValueError("Unknown taxon exclusion stage")

    if not {"island_id", "spatial_block"}.issubset(private_blocks.columns):
        raise ValueError("Missing private spatial-block fields")
    blocks = private_blocks[["island_id", "spatial_block"]].copy()
    if blocks.isna().any().any() or blocks.duplicated("island_id").any():
        raise ValueError("Missing or duplicate spatial-block mapping")
    block_lookup = blocks.set_index("island_id")["spatial_block"]
    if not set(private_counts.island_a).issubset(block_lookup.index) or not set(private_counts.island_b).issubset(block_lookup.index):
        raise ValueError("Missing blocks for evaluated islands")
    if private_counts.empty:
        raise ValueError("Empty private stage table")
    rng = np.random.default_rng(seed)
    grouped = []
    for (scope, region), part in private_counts.groupby(["scope", "region"], sort=True):
        stages = {key: rows.set_index(["island_a", "island_b"]).sort_index()
                  for key, rows in part.groupby("trait_scope")}
        for source_stage, target_stage in STAGES:
            if source_stage not in stages or target_stage not in stages:
                raise ValueError(f"Missing {source_stage} or {target_stage} stage")
            source, target = stages[source_stage], stages[target_stage]
            if not source.index.equals(target.index):
                raise ValueError("Mismatched baseline and excluded island pairs")
            ids_a = np.asarray(target.index.get_level_values("island_a"))
            ids_b = np.asarray(target.index.get_level_values("island_b"))
            for trait in traits:
                n1, n2 = _int_counts(source, f"n_{trait}"), _int_counts(source, f"n_{trait}_b")
                k1, k2 = _int_counts(source, f"k_{trait}_a"), _int_counts(source, f"k_{trait}_b")
                m1, m2 = _int_counts(target, f"n_{trait}"), _int_counts(target, f"n_{trait}_b")
                j1, j2 = _int_counts(target, f"k_{trait}_a"), _int_counts(target, f"k_{trait}_b")
                if ((k1 > n1) | (k2 > n2) | (j1 > m1) | (j2 > m2)
                    | (m1 > n1) | (m2 > n2) | (j1 > k1) | (j2 > k2)
                    | (k1-j1 > n1-m1) | (k2-j2 > n2-m2)).any():
                    raise ValueError("Impossible baseline-to-subset carrier transition for " + trait)
                keep = (m1 >= min_carriers) & (m2 >= min_carriers)
                if not keep.any():
                    continue
                n1,n2,k1,k2,m1,m2,j1,j2 = (v[keep] for v in (n1,n2,k1,k2,m1,m2,j1,j2))
                actual = np.abs(j1/m1 - j2/m2)
                expected = np.empty(len(actual), dtype=float)
                for start in range(0, len(actual), 1200):
                    sl = slice(start, min(start + 1200, len(actual)))
                    xa = rng.hypergeometric(k1[sl], n1[sl]-k1[sl], m1[sl], size=(draws, len(actual[sl])))
                    xb = rng.hypergeometric(k2[sl], n2[sl]-k2[sl], m2[sl], size=(draws, len(actual[sl])))
                    expected[sl] = np.abs(xa/m1[sl]-xb/m2[sl]).mean(axis=0)
                a, b = ids_a[keep], ids_b[keep]
                n_islands = len(set(a) | set(b))
                n_blocks = len(set(block_lookup.loc[a]) | set(block_lookup.loc[b]))
                if len(actual) < min_pairs or n_islands < min_islands or n_blocks < min_spatial_blocks:
                    continue
                residual = actual-expected
                grouped.append({
                    "scope": str(scope), "region": str(region),
                    "contrast": source_stage+"_to_"+target_stage, "trait": trait,
                    "n_pair_comparisons": int(len(actual)),
                    "n_distinct_islands": int(n_islands), "n_spatial_blocks": int(n_blocks),
                    "median_observed_gap": float(np.median(actual)),
                    "median_expected_equal_carrier_gap": float(np.median(expected)),
                    "mean_paired_excess_gap": float(np.mean(residual)),
                    "median_paired_excess_gap": float(np.median(residual)),
                    "fraction_pairs_observed_below_reference": float(np.mean(residual < 0)),
                })
    columns = [
        "scope", "region", "contrast", "trait", "n_pair_comparisons",
        "n_distinct_islands", "n_spatial_blocks",
        "median_observed_gap", "median_expected_equal_carrier_gap",
        "mean_paired_excess_gap", "median_paired_excess_gap",
        "fraction_pairs_observed_below_reference",
    ]
    return pd.DataFrame.from_records(grouped, columns=columns)
