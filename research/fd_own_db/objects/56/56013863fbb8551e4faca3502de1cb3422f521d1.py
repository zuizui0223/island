"""Outcome-blind, source-group-aware functional retention contrasts.

This module deliberately makes only descriptive comparisons. It assumes input
pair matches were fixed from geography, not biological outcomes. Private keys
are never returned by `safe_matched_summary`.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .functional_turnover import TRAIT_GROUPS, audit_pairs

ROLES = ('local', 'distance_control')
PRIMARY_DOMAINS = ('reproductive_assurance', 'floral_accessibility')
INPUT_COLUMNS = ('match_id', 'island_a', 'island_b', 'source_group', 'donor_group', 'role', 'distance_km')


def validate_matches(pairs: pd.DataFrame, *, max_distance_ratio: float = 3.) -> pd.DataFrame:
    """Require two different donors per predefined focal match and group."""
    missing = set(INPUT_COLUMNS) - set(pairs.columns)
    if missing:
        raise ValueError(f'Match columns missing: {sorted(missing)}')
    if max_distance_ratio < 1:
        raise ValueError('max_distance_ratio must be >=1')
    x = pairs[list(INPUT_COLUMNS)].copy()
    for field in ('match_id', 'island_a', 'island_b', 'source_group', 'donor_group', 'role'):
        if x[field].isna().any() or x[field].astype(str).str.strip().eq('').any():
            raise ValueError(f'Missing {field}')
        x[field] = x[field].astype(str).str.strip()
    for field in ('source_group', 'donor_group'):
        if x[field].str.lower().isin(['na','nan','none','unknown']).any():
            raise ValueError(f'Unknown {field} label')
    if not x.role.isin(ROLES).all():
        raise ValueError('role must be local or distance_control')
    x['distance_km'] = pd.to_numeric(x['distance_km'], errors='coerce')
    if x.distance_km.isna().any() or (~np.isfinite(x.distance_km)).any() or (x.distance_km <= 0).any():
        raise ValueError('distance_km must be positive and finite')
    if ((x.role == 'local') & (x.donor_group != x.source_group)).any():
        raise ValueError('Local donor must share focal source_group')
    if ((x.role == 'distance_control') & (x.donor_group == x.source_group)).any():
        raise ValueError('Distance control must belong to another source_group')
    if (x.island_a == x.island_b).any():
        raise ValueError('Focal and donor must differ')
    for match, group in x.groupby('match_id'):
        if len(group) != 2 or set(group.role) != set(ROLES):
            raise ValueError(f'{match}: require one local and one distance_control')
        if group.island_a.nunique() != 1 or group.source_group.nunique() != 1:
            raise ValueError(f'{match}: mismatched focal or source group')
        if group.island_b.nunique() != 2:
            raise ValueError(f'{match}: donor reuse')
        ratio = group.distance_km.max() / group.distance_km.min()
        if ratio > max_distance_ratio:
            raise ValueError(f'{match}: geographic donor distance differs > {max_distance_ratio}x')
    if x.match_id.nunique() != x.island_a.nunique():
        raise ValueError('Each focal island may appear in only one match')
    if x.duplicated(['island_a', 'island_b']).any():
        raise ValueError('Repeated island donor pair')
    return x


def matched_retention(traits: pd.DataFrame, matches: pd.DataFrame,
                      *, min_carriers: int = 20, max_distance_ratio: float = 3.) -> pd.DataFrame:
    """Private match- and trait-specific local advantage; does not make p-values.

    Positive = distance-matched different-archipelago control has LARGER
    prevalence difference than the local donor.  Each domain includes only
    traits scorable on BOTH donor comparisons; >=2 components required.
    """
    m = validate_matches(matches, max_distance_ratio=max_distance_ratio)
    input_pairs = m.rename(columns={'role':'comparison'})[['island_a','island_b','comparison']]
    scores = audit_pairs(traits, input_pairs, min_carriers=min_carriers)
    meta = m[['match_id','island_a','island_b','source_group','role']]
    merged = scores.merge(meta, on=['island_a','island_b'], how='left', validate='many_to_one')
    if merged['match_id'].isna().any():
        raise AssertionError('Unassigned scored pair')
    key = ['match_id','source_group','trait_scope','domain','trait']
    wide = merged.pivot(index=key, columns='role', values='abs_prevalence_gap').reset_index()
    wide = wide.dropna(subset=list(ROLES)).copy()
    wide['local_advantage'] = wide['distance_control'] - wide['local']
    primary = wide[wide.domain.isin(PRIMARY_DOMAINS)]
    domain = primary.groupby(['match_id','source_group','trait_scope','domain'], as_index=False).agg(
        n_component_traits=('trait','nunique'),
        local_advantage=('local_advantage','mean'),
    )
    domain = domain[domain.n_component_traits >= 2].copy()
    domain['trait'] = 'DOMAIN_MEAN'
    return pd.concat([wide[['match_id','source_group','trait_scope','domain','trait',
                            'local_advantage']].assign(n_component_traits=1),
                      domain], ignore_index=True)


def safe_matched_summary(scores: pd.DataFrame, *, min_matches: int=10,
                         min_groups: int=5) -> pd.DataFrame:
    """Equal source-group mean; never exports focal or geographic identities.

    Block/group independence is not guaranteed. This is NOT a p-value.
    Group anonymization is intentionally strict; if counts fail return empty.
    """
    if min_matches < 5 or min_groups < 5:
        raise ValueError('minimum reporting support: >=5 matches and source groups')
    fields = ['trait_scope', 'domain', 'trait']
    records=[]
    for key, sub in scores.groupby(fields, dropna=False):
        if sub.match_id.duplicated().any():
            raise ValueError('Duplicated scored match')
        if sub.match_id.nunique()<min_matches or sub.source_group.nunique()<min_groups:
            continue
        groups = sub.groupby('source_group')['local_advantage'].mean()
        records.append({
            **dict(zip(fields,key)),
            'n_focal_matches':int(sub.match_id.nunique()),
            'n_source_groups':int(groups.size),
            'mean_equal_group_local_advantage':float(groups.mean()),
            'median_match_local_advantage':float(sub.local_advantage.median()),
            'share_groups_positive':float((groups>0).mean()),
        })
    return pd.DataFrame.from_records(records, columns=[*fields,'n_focal_matches',
        'n_source_groups','mean_equal_group_local_advantage',
        'median_match_local_advantage','share_groups_positive'])


def main(argv: list[str] | None = None) -> None:
    """Run locally with rights-authorized source data; release grouped output only."""
    import argparse
    import json
    from pathlib import Path

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--private-traits', type=Path, required=True)
    parser.add_argument('--private-matches', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--min-carriers', type=int, default=20)
    parser.add_argument('--minimum-matches', type=int, default=10)
    parser.add_argument('--minimum-groups', type=int, default=5)
    parser.add_argument('--max-distance-ratio', type=float, default=3.)
    args = parser.parse_args(argv)
    traits = pd.read_csv(args.private_traits, dtype={'island_id':str,'taxon_name':str,'genus':str})
    matches = pd.read_csv(args.private_matches, dtype={
        'match_id':str,'island_a':str,'island_b':str,
        'source_group':str,'donor_group':str,
    })
    scores = matched_retention(traits, matches, min_carriers=args.min_carriers,
                               max_distance_ratio=args.max_distance_ratio)
    summary = safe_matched_summary(scores, min_matches=args.minimum_matches,
                                   min_groups=args.minimum_groups)
    args.output.mkdir(parents=True,exist_ok=True)
    summary.to_csv(args.output/'functional_matched_grouped.csv',index=False)
    (args.output/'status.json').write_text(json.dumps({
        'status':'observational_grouped_only_not_a_null_test',
        'private_identifiers_released':False,
        'min_trait_carriers':args.min_carriers,
        'min_focal_matches':args.minimum_matches,
        'min_source_groups':args.minimum_groups,
        'max_distance_ratio':args.max_distance_ratio,
        'released_comparison_rows':len(summary),
    },indent=2)+'\n',encoding='utf-8')


if __name__ == '__main__':
    main()
