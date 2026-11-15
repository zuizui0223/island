"""Identical-island, identical-pair source-native sensitivity for floral beta.

Paired outcomes come from the SAME original physical island pairs, selected
without consulting floral outcomes, and with ZERO FULL-FLORA shared species,
not just zero shared native-known species. Compare all recorded plants versus
independently original-source-labeled native plants. Independent equal-five
species thinnings address distinct Direct-carrier count differences; 29 draws
are Monte Carlo sensitivities, not 29 biological replicates.
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
from fd_island_beta_geographic_null import Aggregate, BINS


def _genus_residual(rosters: list[set], states: dict,
                    left: np.ndarray, right: np.ndarray):
    p = prepare_endpoint(rosters, states)
    m = compare_block(p, 0, len(rosters))
    observed = m['squared_frequency_distance'][left, right]
    expected = m['genus_null_squared_distance'][left, right]
    return np.asarray(observed), np.asarray(expected)


def _scores(agg_a: Aggregate, agg_n: Aggregate,
            agg_d: Aggregate, min_pairs=100,
            min_islands=30, min_blocks=10):
    q = agg_d.export(min_pairs, min_islands, min_blocks)
    if q is None:
        return None
    a = agg_a.export(min_pairs, min_islands, min_blocks)
    n = agg_n.export(min_pairs, min_islands, min_blocks)
    assert a['n_pairs'] == n['n_pairs'] == q['n_pairs']
    assert a['n_islands'] == n['n_islands'] == q['n_islands']
    return dict(
        n_pairs=q['n_pairs'], n_islands=q['n_islands'],
        n_physical_blocks=q['n_physical_blocks'],
        all_recorded_observed_D2=a['pair_equal_mean_observed_D2'],
        all_recorded_expected_D2=a['pair_equal_mean_genus_null_D2'],
        native_observed_D2=n['pair_equal_mean_observed_D2'],
        native_expected_D2=n['pair_equal_mean_genus_null_D2'],
        all_recorded_mean_residual_D2=a['pair_equal_mean_residual_D2'],
        native_mean_residual_D2=n['pair_equal_mean_residual_D2'],
        all_minus_native_paired_mean_residual_D2=q['pair_equal_mean_residual_D2'],
        island_equal_mean_paired_delta_D2=q['island_equal_mean_residual_D2'],
        physical_block_equal_mean_paired_delta_D2=(
            q['physical_block_equal_mean_residual_D2']),
        leave_one_block_out_paired_delta_min=(
            q['leave_one_physical_block_out_pair_mean_residual_D2_min']),
        leave_one_block_out_paired_delta_max=(
            q['leave_one_physical_block_out_pair_mean_residual_D2_max']),
        n_block_deletion_positive=q['n_leave_one_block_positive'],
        n_block_deletion_checks=q['n_leave_one_block_checks'])


def analyze(states: dict, flora: pd.DataFrame, geo: pd.DataFrame,
            *, draws: int = 29, seed: int = 20261115,
            min_pairs: int = 100, min_islands: int = 30,
            min_blocks: int = 10) -> dict:
    if draws < 1:
        raise ValueError('Need >=1 deterministic thinnings')
    need = {'island_id','accepted_species','origin_status'}
    if not need <= set(flora) or flora.duplicated(
            ['island_id','accepted_species']).any():
        raise ValueError('Missing/duplicate accepted original island species')
    if not flora.origin_status.isin(('native','introduced','unresolved')).all():
        raise ValueError('Original source nativeness not certified')
    required = {'island_id','analysis_regime','spatial_block',
                'island_latitude','island_longitude'}
    if not required <= set(geo) or geo.island_id.duplicated().any():
        raise ValueError('Missing independent geography')
    flora_all = flora.groupby('island_id').accepted_species.agg(set).to_dict()
    flora_native = flora.loc[flora.origin_status.eq('native')].groupby(
        'island_id').accepted_species.agg(set).to_dict()
    islands = sorted(flora_all)
    common = set(states['flower_primary_color']) & set(states['floral_form'])
    categories = {
        sp:(states['flower_primary_color'][sp],states['floral_form'][sp])
        for sp in common}
    qualified = [
        sp for sp in islands
        if len(flora_native.get(sp,set()) & common) >= 5
    ]
    n = len(qualified)
    if n < 2:
        raise ValueError('Not enough source-native morphological floras')
    roster_all = [flora_all[i] & common for i in qualified]
    roster_native = [flora_native.get(i,set()) & common for i in qualified]
    geos = geo.set_index('island_id').reindex(qualified)
    reg = geos.analysis_regime.fillna('geography_unresolved').astype(str).to_numpy()
    block = geos.spatial_block.to_numpy()
    latitude = pd.to_numeric(geos.island_latitude,errors='coerce').to_numpy()
    longitude = pd.to_numeric(geos.island_longitude,errors='coerce').to_numpy()
    block_number, unique_blocks = pd.factorize(block)
    full_names = sorted(set().union(*(flora_all[i] for i in qualified)))
    full_rosters = [flora_all[i] for i in qualified]
    A = membership(full_rosters, full_names)
    full_shared = (A @ A.T).toarray()
    distances = geographical_distance(latitude,longitude,0,n)
    left,right = np.triu_indices(n,1)
    distance = distances[left,right]
    valid = (
        (full_shared[left,right] == 0) &
        (reg[left] == reg[right]) &
        (reg[left] != 'geography_unresolved') &
        (block_number[left] >= 0) & (block_number[right] >= 0) &
        (block_number[left] != block_number[right]) &
        np.isfinite(distance) &
        np.isfinite(latitude[left]) & np.isfinite(latitude[right]) &
        np.isfinite(longitude[left]) & np.isfinite(longitude[right]) &
        (np.abs(latitude[left]) <= 90) &
        (np.abs(latitude[right]) <= 90) &
        (np.abs(longitude[left]) <= 180) &
        (np.abs(longitude[right]) <= 180)
    )
    left,right,distance = left[valid],right[valid],distance[valid]
    groups = {}
    for region in sorted(set(reg[left])):
        same = reg[left] == region
        for dist_class, low, high in BINS:
            mask = same & (distance >= low) & (distance < high)
            if mask.any():
                groups[(region,dist_class)] = np.flatnonzero(mask)
    if not groups:
        raise ValueError('No valid original same-region geographic pair support')
    # Public eligibility is metadata about islands, not pair-level identities.
    prelim = dict(n_original_botanical_islands=len(islands),
                  n_original_native_direct_pair_eligible_islands=n,
                  n_matched_eligible_same_region_other_block_full_species_disjoint_pairs=len(left))
    rng = np.random.default_rng(seed)
    protocols = [('full_direct',None),('equal_5_direct',draws)]
    results = []
    for protocol, repeats in protocols:
        ndraw = 1 if repeats is None else repeats
        local_scores = {}
        for k in range(ndraw):
            if protocol == 'full_direct':
                a_roster = roster_all
                n_roster = roster_native
            else:
                a_roster = [
                    set(rng.choice(sorted(r),size=5,replace=False))
                    for r in roster_all]
                n_roster = [
                    set(rng.choice(sorted(r),size=5,replace=False))
                    for r in roster_native]
            ao,ae = _genus_residual(a_roster,categories,left,right)
            no,ne = _genus_residual(n_roster,categories,left,right)
            for key,chosen in groups.items():
                if key not in local_scores:
                    local_scores[key] = dict(
                        a_ob=np.zeros(len(chosen)),a_ex=np.zeros(len(chosen)),
                        n_ob=np.zeros(len(chosen)),n_ex=np.zeros(len(chosen)),
                        montecarlo_paired_delta=[])
                bag = local_scores[key]
                bag['a_ob'] += ao[chosen]/ndraw
                bag['a_ex'] += ae[chosen]/ndraw
                bag['n_ob'] += no[chosen]/ndraw
                bag['n_ex'] += ne[chosen]/ndraw
                bag['montecarlo_paired_delta'].append(
                    float(np.mean((ao[chosen]-ae[chosen])-
                                  (no[chosen]-ne[chosen]))))
        for (rg, dist), bag in sorted(local_scores.items()):
            ix=groups[(rg,dist)]
            aa = Aggregate(n,block_number,len(unique_blocks))
            nn = Aggregate(n,block_number,len(unique_blocks))
            dd = Aggregate(n,block_number,len(unique_blocks))
            aa.add(left[ix],right[ix],bag['a_ob'],bag['a_ex'])
            nn.add(left[ix],right[ix],bag['n_ob'],bag['n_ex'])
            dd.add(left[ix],right[ix],
                   bag['a_ob']-bag['a_ex'],
                   bag['n_ob']-bag['n_ex'])
            row = _scores(aa,nn,dd,min_pairs,min_islands,min_blocks)
            if row is None:
                continue
            vals=bag['montecarlo_paired_delta']
            results.append(dict(
                comparison='all_recorded_minus_native_on_identical_full_species_disjoint_island_pairs',
                protocol=protocol, region=rg, distance_class=dist,
                n_montecarlo_thinnings=ndraw,
                montecarlo_paired_mean_delta_min=float(min(vals)),
                montecarlo_paired_mean_delta_max=float(max(vals)),
                **row))
    return dict(
        status='ORIGINAL_MATCHED_NATIVE_FLORAL_BETA_DIAGNOSTIC',
        **prelim,
        comparisons=results,
        claims=dict(
            identical_physical_island_pairs_within_each_stratum=True,
            no_species_shared_in_original_full_flora=True,
            globally_fixed_species_wide_phenotypes=True,
            source_native_subset_not_random_exposure_or_treatment=True,
            category_codes_and_global_genus_null_shared_across_scopes=True,
            equal_5_draws_not_independent_biological_replicates=True,
            no_pairwise_independence_p_value=True,
            no_reproductive_fitness_or_pollinator_shift_inference=True,
            no_private_original_island_or_plant_ids_exported=True,
        ),
    )


def _sha(path:Path)->str:
    h=hashlib.sha256()
    with path.open('rb') as fp:
        for chunk in iter(lambda:fp.read(1<<20),b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap=argparse.ArgumentParser()
    for k in ('original-root','geology','output'):
        ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--draws',type=int,default=29)
    args=ap.parse_args()
    paths={
        'traits':args.original_root/'fixed/source/current/evidence/direct_only_trait_ledger.csv.gz',
        'flora':args.original_root/'fixed/canonical/input/chapter1_status_flora.csv.gz',
    }
    expected={
        'traits':'9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745',
        'flora':'6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a',
    }
    for k in paths:
        if _sha(paths[k]) != expected[k]:
            raise ValueError('Protected original source checksum mismatch')
    traits=pd.read_csv(paths['traits'],dtype=str,usecols=[
        'accepted_species','trait_name','normalized_value',
        'evidence_scope','evidence_quality'])
    flora=pd.read_csv(paths['flora'],dtype=str,usecols=[
        'island_id','accepted_species','origin_status'])
    geo=pd.read_csv(args.geology,dtype=str)
    if flora.island_id.nunique()!=4505:
        raise ValueError('Original flora changed')
    out=analyze(strict_traits(traits),flora,geo,draws=args.draws)
    if out['n_original_native_direct_pair_eligible_islands'] != 327:
        raise ValueError('Original native direct support changed')
    body=json.dumps(out,allow_nan=False,sort_keys=True,indent=2)
    for forbidden in ('"island_id"','"accepted_species"','"spatial_block"',
                      '"source_signature"','"island_latitude"'):
        if forbidden in body:raise ValueError('Protected original identity exported')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+'\n',encoding='utf-8')
    print(json.dumps(dict(status=out['status'],
                          n_native_islands=out['n_original_native_direct_pair_eligible_islands'],
                          n_anonymous_strata=len(out['comparisons']),
                          no_pairwise_independence_claim=True)))


if __name__=='__main__':
    main()
