"""Materialize the pinned public checkpoint; never label private history reproduced."""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import pandas as pd

RUN_ID = 33605098044
ARTIFACT_ID = 9836760588
COVERAGE_SHA256 = 'c236ce2829ec21326ddc757321d34e9d2057ba7525e98c5db22c089cd64b5d9e'
AXES = {'flower_colour', 'floral_structural_complexity', 'reproductive_assurance'}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def recover(source: Path, output: Path) -> dict:
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
    shutil.copytree(source, output / 'public_source')
    (output / 'recovery_manifest.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(recover(args.source, args.output), indent=2))
