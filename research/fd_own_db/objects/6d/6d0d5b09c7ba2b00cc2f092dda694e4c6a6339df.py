import importlib.util
import json
from pathlib import Path

import pytest

spec=importlib.util.spec_from_file_location('taxon_audit',Path(__file__).resolve().parents[1]/'scripts'/'audit_public_island_taxon_frequency.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_source_freq_audit_matches_status_and_single_island_counts(tmp_path):
    source=tmp_path/'taxa.csv'; status=tmp_path/'status.json'
    source.write_text('accepted_species,genus,n_islands,n_records\n'
                      'A a,A,1,1\nB b,B,1,2\nC c,C,3,9\nD d,D,12,30\n')
    status.write_text(json.dumps({'n_accepted_species':4,
        'n_island_species_pairs':17,'n_islands_with_records':20,
        'n_collected_succeeded_blocks':2}))
    r=module.summarize(source,status)
    assert r['n_accepted_taxon_names']==4
    assert r['n_island_taxon_occurrences']==17
    assert r['taxa_recorded_on_one_island']==2
    assert r['fraction_recorded_on_one_island']==.5
    assert r['n_islands_with_any_record']==20


def test_coverage_and_duplicate_fail_closed(tmp_path):
    p=tmp_path/'taxa.csv'
    p.write_text('accepted_species,genus,n_islands,n_records\nA a,A,3,2\n')
    with pytest.raises(ValueError,match='Impossible'):
        module.summarize(p)
    p.write_text('accepted_species,genus,n_islands,n_records\nA a,A,1,2\nA a,A,1,3\n')
    with pytest.raises(ValueError,match='duplicated'):
        module.summarize(p)
