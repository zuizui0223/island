"""Materialize the pinned public checkpoint; never label private history reproduced."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import itertools
from functools import lru_cache
from pathlib import Path

import pandas as pd
import yaml

RUN_ID = 33605098044
ARTIFACT_ID = 9836760588
COVERAGE_SHA256 = 'c236ce2829ec21326ddc757321d34e9d2057ba7525e98c5db22c089cd64b5d9e'
AXES = {'flower_colour', 'floral_structural_complexity', 'reproductive_assurance'}
LOW_SIDECAR_SHA256 = 'b62d5ae133c7029ce7a67fb57afa230b18160ea6ea27c9fea5beb14a557b60a1'


def restore_low_values(coverage: pd.DataFrame, sidecar: Path, ontology: Path):
    if sha256(sidecar) != LOW_SIDECAR_SHA256:
        raise ValueError('Pinned historical Low sidecar checksum mismatch')
    allowed = {k: set(v['allowed_values']) for k, v in yaml.safe_load(ontology.read_text(encoding='utf-8'))['traits'].items()}
    @lru_cache(maxsize=None)
    def map_states(names_json, states_json):
        names = json.loads(names_json)
        states = [json.loads(s) for s in json.loads(states_json)]
        if not names or not states or len(states) > len(names) or len(names) > 4 or len(set(names)) != len(names):
            raise ValueError('Malformed historical trait/state lists')
        mappings = set()
        candidates = [[s for s in states if s and set(s).issubset(allowed.get(t, set()))] for t in names]
        supplied = {json.dumps(sorted(s)) for s in states}
        # Identical state sets were deduplicated in the historical sidecar;
        # multiple traits may share a set, but every supplied set must be used.
        for order in itertools.product(*candidates):
            if {json.dumps(sorted(s)) for s in order} == supplied:
                mappings.add(json.dumps(dict(zip(names, order)), sort_keys=True))
        # The old sidecar sorts traits and states independently. Never zip them
        # by position: SC must not become autonomous_selfing_capacity.
        if len(mappings) != 1:
            raise ValueError('Historical trait/state association is ambiguous')
        return json.loads(next(iter(mappings)))
    low = pd.read_csv(sidecar, dtype=str).fillna('')
    keys = ['accepted_species', 'axis']
    if low.duplicated(keys).any():
        raise ValueError('Duplicate historical Low sidecar key')
    lookup = low.set_index(keys)
    repaired = coverage.copy()
    affected = repaired.quality.eq('low') & repaired.trait_composition.eq('')
    audit = []
    for i, row in repaired.loc[affected].iterrows():
        key = (row.accepted_species, row.axis)
        if key not in lookup.index:
            raise ValueError(f'Missing Low sidecar for {key}')
        prior = lookup.loc[key]
        if prior.family_inference.lower() != 'false' or prior.global_fallback.lower() != 'false':
            raise ValueError('Forbidden fallback in historical sidecar')
        mapping = map_states(prior.trait_names, prior.predicted_state_sets)
        lineages = json.loads(prior.support_source_lineages)
        if not lineages:
            raise ValueError('Missing historical supporting lineage')
        repaired.loc[i, 'trait_composition'] = '|'.join(f'{trait}={json.dumps(states, separators=(",", ":"))}' for trait, states in sorted(mapping.items()))
        repaired.loc[i, 'trait_names'] = '|'.join(sorted(mapping))
        repaired.loc[i, 'source_groups'] = 'wave33_validated_probabilistic_genus_low'
        repaired.loc[i, 'source_lineages'] = '|'.join(lineages)
        audit.append({'accepted_species': key[0], 'axis': key[1], 'restoration': 'historical_secondary_low_sidecar', 'imputation_required': True, 'new_validation_claim': False})
    return repaired, pd.DataFrame(audit)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def recover(source: Path, output: Path, low_sidecar: Path | None = None, ontology: Path | None = None) -> dict:
    if output.exists():
        raise ValueError('Output must be new: existing checkpoints are never overwritten')
    coverage_path = source / 'wave53-trait-coverage/wave53_species_axis_coverage.csv.gz'
    if sha256(coverage_path) != COVERAGE_SHA256:
        raise ValueError('Pinned coverage checksum mismatch')
    manifests = list(source.glob('*/**/*summary.json'))
    checked = {}
    for manifest in manifests:
        data = json.loads(manifest.read_text(encoding='utf-8'))
        for name, expected in data.get('artifact_sha256', {}).items():
            path = manifest.parent / name
            if not path.is_file() or sha256(path) != expected:
                raise ValueError(f'Missing or altered source artifact: {name}')
            checked[path.relative_to(source).as_posix()] = expected
    required = [
        'wave53-all-evidence-audit/resolved_direct_species_trait.csv.gz',
        'wave53-all-evidence-audit/rebuilt_all_evidence_validated_low.csv.gz',
        'wave53-all-evidence-audit/trait_specific_genus_rule_audit.csv.gz',
    ]
    if not all(name in checked for name in required):
        raise ValueError('Required ledger files were not hash-verified')
    coverage = pd.read_csv(coverage_path, dtype=str).fillna('')
    if len(coverage) != 318885 or coverage.accepted_species.nunique() != 106295:
        raise ValueError('Wrong fixed denominator')
    if coverage.duplicated(['accepted_species', 'axis']).any() or set(coverage.axis) != AXES:
        raise ValueError('Duplicated or invalid species-axis keys')
    if not coverage.quality.isin(['', 'high', 'medium', 'low']).all():
        raise ValueError('Unknown quality tier')
    filled = coverage.quality.ne('')
    counts = coverage.loc[filled].groupby('axis').size().astype(int).to_dict()
    expected = {'flower_colour': 82502, 'floral_structural_complexity': 91613, 'reproductive_assurance': 48260}
    if counts != expected:
        raise ValueError('Pinned checkpoint counts mismatch')
    report = {
        'status': 'verified_public_checkpoint_recovered',
        'source_run_id': RUN_ID, 'source_artifact_id': ARTIFACT_ID,
        'denominator_species': 106295, 'denominator_species_axis': len(coverage),
        'filled': int(filled.sum()), 'unresolved': int((~filled).sum()),
        'coverage_percent': 100 * float(filled.mean()), 'by_axis': counts,
        'quality_counts': coverage.loc[filled, 'quality'].value_counts().to_dict(),
        'species_by_filled_axes': filled.groupby(coverage.accepted_species).sum().value_counts().sort_index().to_dict(),
        'historical_private_reported_filled': 222759,
        'historical_private_checkpoint_reproduced': False,
        'historical_difference_not_reproduced': 222759 - int(filled.sum()),
        'new_acquisition_gain': 0,
        'claim_limit': 'Restores the immutable public ledger; does not independently validate every original biological assertion or reproduce private TRY/Wave55 adjudication.',
        'verified_source_sha256': checked,
    }
    if low_sidecar is not None:
        if ontology is None:
            raise ValueError('Ontology required for trait-specific restoration')
        materialized, restoration = restore_low_values(coverage, low_sidecar, ontology)
        if not coverage[['accepted_species','axis','quality']].equals(materialized[['accepted_species','axis','quality']]):
            raise ValueError('Restoration changed coverage or quality')
        report['restored_missing_value_cells'] = len(restoration)
        report['remaining_labelled_cells_without_value'] = int((materialized.quality.ne('') & materialized.trait_composition.eq('')).sum())
        report['low_sidecar_sha256'] = sha256(low_sidecar)
        report['ontology_sha256'] = sha256(ontology)
        report['secondary_low_not_confirmatory'] = True
    shutil.copytree(source, output / 'public_source')
    if low_sidecar is not None:
        materialized.to_csv(output/'materialized_species_axis_coverage.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        restoration.to_csv(output/'low_value_restoration_audit.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        shutil.copy2(low_sidecar, output/'historical_secondary_low_sidecar.csv.gz')
    (output / 'recovery_manifest.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--low-sidecar', type=Path)
    parser.add_argument('--ontology', type=Path)
    args = parser.parse_args()
    print(json.dumps(recover(args.source, args.output, args.low_sidecar, args.ontology), indent=2))
