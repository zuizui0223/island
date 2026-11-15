"""Original-data, privacy-preserving geographical test of floral composition residuals.

Primary unit: island floras; repeated dyads are never independent replicates.
Only species-disjoint, different-10-degree-block, same-broad-region pairs are
eligible. Every response compares exact squared categorical frequency distance
with the expected distance under global within-genus whole-profile shuffling,
holding the original island-by-species mask unchanged. A genus-only global null
does NOT remove climate, geography, source-pool, phylogeny or observation bias.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from fd_island_functional_beta import (
    strict_traits, prepare_endpoint, compare_block, membership,
    geographical_distance,
)

MIN_DIRECT = 5
MIN_BLOCKS = 10
MIN_ISLANDS = 30
MIN_PAIRS = 100
BINS = (("under_500km", 0., 500.),
        ("500_to_2000km", 500., 2000.),
        ("at_least_2000km", 2000., np.inf),
        ("all_distances", 0., np.inf))


class Aggregate:
    """Sufficient statistics for pair-, island- and block-equal diagnostics."""
    def __init__(self, n_islands: int, island_block: np.ndarray, n_blocks: int):
        self.n_pairs = 0
        self.sum_observed = 0.
        self.sum_expected = 0.
        self.sum_delta = 0.
        self.node_count = np.zeros(n_islands, np.int64)
        self.node_sum = np.zeros(n_islands)
        self.block_count = np.zeros(n_blocks, np.int64)
        self.block_sum = np.zeros(n_blocks)
        self.island_block = island_block

    def add(self, left: np.ndarray, right: np.ndarray,
            observed: np.ndarray, expected: np.ndarray) -> None:
        assert len(left) == len(right) == len(observed) == len(expected)
        if not len(left):
            return
        assert np.all(left < right)
        assert np.all(self.island_block[left] != self.island_block[right])
        delta = observed - expected
        self.n_pairs += len(left)
        self.sum_observed += float(observed.sum())
        self.sum_expected += float(expected.sum())
        self.sum_delta += float(delta.sum())
        both = np.r_[left, right]
        weights = np.r_[delta, delta]
        self.node_count += np.bincount(both, minlength=len(self.node_count))
        self.node_sum += np.bincount(
            both, weights=weights, minlength=len(self.node_count))
        b = self.island_block[both]
        self.block_count += np.bincount(
            b, minlength=len(self.block_count))
        self.block_sum += np.bincount(
            b, weights=weights, minlength=len(self.block_count))

    def export(self, min_pairs=MIN_PAIRS, min_islands=MIN_ISLANDS,
               min_blocks=MIN_BLOCKS) -> dict | None:
        active = self.node_count > 0
        n_islands = int(active.sum())
        active_blocks = np.unique(self.island_block[active])
        n_blocks = len(active_blocks)
        if (self.n_pairs < min_pairs or n_islands < min_islands or
                n_blocks < min_blocks):
            return None
        island_mean = self.node_sum[active] / self.node_count[active]
        groups = self.island_block[active]
        block_means = [float(island_mean[groups == g].mean())
                       for g in active_blocks]
        loo = []
        for b in active_blocks:
            other_n = self.n_pairs - self.block_count[b]
            if other_n > 0:
                loo.append(float((self.sum_delta - self.block_sum[b]) /
                                 other_n))
        return dict(
            n_pairs=int(self.n_pairs),
            n_islands=n_islands,
            n_physical_blocks=n_blocks,
            pair_equal_mean_observed_D2=float(self.sum_observed / self.n_pairs),
            pair_equal_mean_genus_null_D2=float(self.sum_expected / self.n_pairs),
            pair_equal_mean_residual_D2=float(self.sum_delta / self.n_pairs),
            island_equal_mean_residual_D2=float(island_mean.mean()),
            physical_block_equal_mean_residual_D2=float(np.mean(block_means)),
            leave_one_physical_block_out_pair_mean_residual_D2_min=(
                float(min(loo)) if loo else None),
            leave_one_physical_block_out_pair_mean_residual_D2_max=(
                float(max(loo)) if loo else None),
            n_leave_one_block_checks=len(loo),
            n_leave_one_block_positive=int(sum(x > 0 for x in loo)),
        )


def analyze(states: dict, flora: pd.DataFrame, geo: pd.DataFrame, *,
            block_size: int = 96, min_pairs: int = MIN_PAIRS,
            min_islands: int = MIN_ISLANDS,
            min_blocks: int = MIN_BLOCKS) -> dict:
    """Compute identified descriptive same-region, species-disjoint strata.

    Coverage gates and distance bins are fixed without reading trait outcomes.
    The absolute-difference residual cannot prove selection or guild shifts.
    """
    if block_size < 1:
        raise ValueError("Invalid block size")
    if not {'island_id', 'accepted_species', 'origin_status'} <= set(flora):
        raise ValueError("Incomplete source botanical roster")
    if flora.duplicated(['island_id', 'accepted_species']).any():
        raise ValueError("Duplicate original accepted species per island")
    if not flora.origin_status.isin(('native', 'introduced', 'unresolved')).all():
        raise ValueError("Missing or invalid original botanical status")
    required = {'island_id','analysis_regime','spatial_block',
                'island_latitude','island_longitude'}
    if not required <= set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Verified original GSHHG geography required")
    islands = sorted(flora.island_id.unique())
    geography = geo.set_index('island_id').reindex(islands)
    lat = pd.to_numeric(geography.island_latitude, errors='coerce').to_numpy()
    lon = pd.to_numeric(geography.island_longitude, errors='coerce').to_numpy()
    blocks = geography.spatial_block.to_numpy()
    regions = geography.analysis_regime.fillna('geography_unresolved').astype(str).to_numpy()
    all_species = flora.groupby('island_id').accepted_species.agg(set).to_dict()
    native_species = flora.loc[flora.origin_status.eq('native')].groupby(
        'island_id').accepted_species.agg(set).to_dict()
    common = set(states['flower_primary_color']) & set(states['floral_form'])
    category = {sp: (states['flower_primary_color'][sp], states['floral_form'][sp])
                for sp in common}
    output = []
    scopes = (('all_recorded', all_species),
              ('original_native', native_species))
    eligible = {}
    for scope, mapping in scopes:
        originals = [mapping.get(i, set()) for i in islands]
        rosters = [r & common for r in originals]
        n = np.array([len(r) for r in rosters])
        use = np.flatnonzero(n >= MIN_DIRECT)
        eligible[scope] = int(len(use))
        if len(use) < 2:
            continue
        p = prepare_endpoint([rosters[i] for i in use], category)
        ids = [islands[i] for i in use]
        full_species = sorted(set().union(*originals))
        full_A = membership(originals, full_species)[use].tocsr()
        full_n = np.asarray(full_A.sum(axis=1)).ravel()
        region = regions[use]
        block = blocks[use]
        latitude = lat[use]
        longitude = lon[use]
        values = np.asarray(p['n'])
        fraction = values / full_n
        block_number, distinct_block = pd.factorize(block)
        n_blocks = len(distinct_block)
        sinks = {}
        for start in range(0, len(use), block_size):
            stop = min(start + block_size, len(use))
            mat = compare_block(p, start, stop)
            shared_full = (full_A[start:stop] @ full_A.T).toarray()
            row, col = np.nonzero(
                np.arange(len(use))[None, :] >
                np.arange(start, stop)[:, None])
            left, right = row + start, col
            if not len(left):
                continue
            distances = geographical_distance(latitude, longitude, start, stop)
            d = distances[row, col]
            valid = (
                (shared_full[row, col] == 0) &
                (region[left] == region[right]) &
                (region[left] != 'geography_unresolved') &
                (block_number[left] >= 0) &
                (block_number[right] >= 0) &
                (block_number[left] != block_number[right]) &
                np.isfinite(d) & np.isfinite(latitude[left]) &
                np.isfinite(latitude[right]) &
                np.isfinite(longitude[left]) &
                np.isfinite(longitude[right]) &
                (np.abs(latitude[left]) <= 90) &
                (np.abs(latitude[right]) <= 90) &
                (np.abs(longitude[left]) <= 180) &
                (np.abs(longitude[right]) <= 180)
            )
            if not valid.any():
                continue
            obs = mat['squared_frequency_distance'][row, col]
            ref = mat['genus_null_squared_distance'][row, col]
            gates = {
                'n5': valid,
                'both_n10': valid & (values[left] >= 10) & (values[right] >= 10),
                'both_direct_fraction_ge_0p10': (
                    valid & (fraction[left] >= .10) & (fraction[right] >= .10)),
            }
            for rg in sorted(set(region[left[valid]])):
                group = (region[left] == rg)
                for gate, criterion in gates.items():
                    if gate == 'n5':
                        dbins = BINS
                    else:
                        dbins = (BINS[-1],)
                    for name, lo, hi in dbins:
                        keep = group & criterion & (d >= lo) & (d < hi)
                        if not keep.any():
                            continue
                        key = (scope, rg, gate, name)
                        if key not in sinks:
                            sinks[key] = Aggregate(len(use), block_number, n_blocks)
                        sinks[key].add(left[keep], right[keep], obs[keep], ref[keep])
        for (s, rg, gate, dist), summary in sorted(sinks.items()):
            result = summary.export(min_pairs, min_islands, min_blocks)
            if result is not None:
                output.append(dict(scope=s, region=rg, support=gate,
                                   distance_class=dist, **result))
    return dict(
        status='ORIGINAL_GEOGRAPHIC_FLORAL_GENUS_NULL_DESCRIPTIVE',
        n_original_botanical_islands=len(islands),
        endpoint='observed_direct_colour__form',
        eligible_islands_by_scope=eligible,
        primary_pair_eligibility='zero_shared_accepted_species_within_current_floristic_scope__same_region__different_10degree_blocks',
        geographic_bins_km=[500,2000],
        comparisons=output,
        claims=dict(
            no_private_island_species_or_source_ids_exported=True,
            original_plant_morphology_state_fixed_globally_per_species=True,
            native_scope_species_disjoint_does_not_guarantee_full_flora_disjoint=True,
            no_pair_as_independent_replication=True,
            genus_null_not_region_climate_source_or_phylogeny_conditioned=True,
            no_adaptive_pollinator_syndrome_or_service_claim=True,
            missing_direct_traits_not_imputed_zero=True,
            leave_one_block_range_is_sensitivity_not_confidence_interval=True,
        ),
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as fp:
        for chunk in iter(lambda: fp.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    for arg in ('original-root','geology','output'):
        ap.add_argument('--' + arg, required=True, type=Path)
    ap.add_argument('--block-size', type=int, default=96)
    args = ap.parse_args()
    paths = {
        'traits': args.original_root / 'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz',
        'flora': args.original_root / 'fixed/canonical/input/chapter1_status_flora.csv.gz',
    }
    expected = {
        'traits': '9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745',
        'flora': '6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a',
    }
    for k in paths:
        if sha256(paths[k]) != expected[k]:
            raise ValueError('Original protected source checksum changed')
    trait = pd.read_csv(paths['traits'], dtype=str, usecols=[
        'accepted_species', 'trait_name', 'normalized_value',
        'evidence_scope', 'evidence_quality'])
    flora = pd.read_csv(paths['flora'], dtype=str, usecols=[
        'island_id', 'accepted_species', 'origin_status'])
    geo = pd.read_csv(args.geology, dtype=str)
    if flora.island_id.nunique() != 4505:
        raise ValueError('Source 4505-island botanical universe changed')
    if geo.geological_origin.isin(('oceanic','continental')).sum() != 434:
        raise ValueError('Independent geological source crosswalk changed')
    result = analyze(strict_traits(trait), flora, geo,
                     block_size=args.block_size)
    if result['eligible_islands_by_scope']['all_recorded'] != 2324:
        raise ValueError('Original floral Direct colour/form eligibility changed')
    body = json.dumps(result, sort_keys=True, indent=2, allow_nan=False)
    for forbidden in ('"island_id"', '"accepted_species"', '"source_signature"',
                      '"spatial_block"', '"latitude"', '"longitude"'):
        if forbidden in body:
            raise ValueError('Original protected botanical identity in output')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(body + '\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'],
                          n_original_islands=result['n_original_botanical_islands'],
                          n_exported_aggregate_strata=len(result['comparisons']),
                          note='No inferential p values from island pairs.')))


if __name__ == '__main__':
    main()
