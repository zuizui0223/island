"""Add reviewed direct claims to empty cells; never overwrite or infer genus values."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import pandas as pd
import yaml

STRICT_AXES = {
    'flower_primary_color': 'flower_colour',
    **{t:'floral_structural_complexity' for t in ['floral_form','floral_symmetry','tube_depth_class','flower_size_class','inflorescence_display']},
    **{t:'reproductive_assurance' for t in ['self_incompatibility','autonomous_selfing_capacity','mating_system','cleistogamy']},
}


def integrate(c, d, records, ontology):
    c, d = c.copy(), d.copy()
    original = c.copy()
    if c.duplicated(['accepted_species', 'axis']).any():
        raise ValueError('Duplicate baseline cell')
    changes, new_direct = [], []
    seen = set()
    for r in records:
        key = (r['accepted_species'], r['trait_name'])
        if key in seen or ((d.accepted_species == key[0]) & (d.trait_name == key[1])).any():
            raise ValueError('Existing or duplicate trait: requires full conflict resolution')
        seen.add(key)
        trait = ontology['traits'][key[1]]
        states = r.get('state_set', [r['normalized_value']])
        if STRICT_AXES.get(key[1]) != r['axis'] or not states or not set(states).issubset(trait['allowed_values']) or 'unresolved' in states:
            raise ValueError('Invalid trait ontology')
        if r['quality'] not in ('high', 'medium') or r['review_status'] != 'accepted_direct_statement':
            raise ValueError('Unreviewed direct evidence')
        for field in ('source_url','source_excerpt','source_lineage','review_date','reviewer','name_match_method'):
            if not r.get(field):
                raise ValueError('Missing provenance')
        if r['cultivar_status'] not in ['natural_populations','species_database_record'] or r['genus_rule_training_allowed']:
            raise ValueError('Unsupported promotion scope')
        mask = c.accepted_species.eq(key[0]) & c.axis.eq(r['axis'])
        if mask.sum() != 1 or original.loc[mask, 'quality'].iloc[0] != '':
            raise ValueError('Not an empty baseline axis; requires full conflict resolution')
        states = sorted(set(states))
        state = json.dumps(states, separators=(',', ':'))
        old_quality = c.loc[mask, 'quality'].iloc[0]
        c.loc[mask, 'quality'] = 'high' if 'high' in [old_quality,r['quality']] else r['quality']
        for field, value in [('trait_names',r['trait_name']),('trait_composition',r['trait_name']+'='+state),('source_lineages',r['source_lineage'])]:
            existing = c.loc[mask,field].iloc[0]
            c.loc[mask,field] = '|'.join(sorted(set(filter(None,existing.split('|')+[value]))))
        c.loc[mask, 'source_groups'] = 'reviewed_restart_20260907'
        row = {col: '' for col in d.columns}
        row.update(accepted_species=key[0], axis=r['axis'], trait_name=key[1],
            classification='single_independent_lineage', resolution_status='resolved',
            selected_quality=r['quality'], quality=r['quality'], state_set=state,
            state_sets=json.dumps([states]), normalized_value='|'.join(states),
            n_independent_lineages='1', source_lineages=r['source_lineage'],
            source_groups='reviewed_restart_20260907', multistate=str(len(states)>1), genus=key[0].split()[0])
        new_direct.append(row)
        changes.append(dict(accepted_species=key[0], axis=r['axis'], before_quality='', after_quality=r['quality'], record_id=r['record_id']))
    return c, pd.concat([d, pd.DataFrame(new_direct)], ignore_index=True), pd.DataFrame(changes)


def correct_harrisia(c, d, receipt, ontology):
    """Pinned source correction; does not relax the generic add-only admission gate."""
    c, d = c.copy(), d.copy()
    species = 'Harrisia portoricensis'
    axis = 'reproductive_assurance'
    old_lineage = 'url:https://europepmc.org/article/MED/21622342'
    lineage = 'doi:10.3732/ajb.0900026'
    if (receipt['accepted_species'] != species or receipt['source_lineage'] != lineage
            or receipt['review_status'] != 'accepted_source_correction'
            or receipt['genus_rule_training_allowed'] is not False
            or receipt['compatibility_quote'] != 'H. portoricensis has a partially self-compatible breeding system'
            or receipt['autonomy_quote'] != 'it is not autogamous'):
        raise ValueError('Unapproved source correction')
    cm = c.accepted_species.eq(species) & c.axis.eq(axis)
    dm = d.accepted_species.eq(species) & d.trait_name.eq('self_incompatibility')
    am = d.accepted_species.eq(species) & d.trait_name.eq('autonomous_selfing_capacity')
    if (cm.sum() != 1 or dm.sum() != 1
            or c.loc[cm,'trait_composition'].iloc[0] != 'self_incompatibility=["SC"]'
            or c.loc[cm,'quality'].iloc[0] != 'high'
            or d.loc[dm,'normalized_value'].iloc[0] != 'SC'
            or d.loc[dm,'source_lineages'].iloc[0] != old_lineage
            or am.sum() != 1 or d.loc[am,'normalized_value'].iloc[0] != 'absent'
            or d.loc[am,'quality'].iloc[0] != 'medium'
            or d.loc[am,'source_lineages'].iloc[0] != 'doi:10.1038/ncomms13313'):
        raise ValueError('Correction precondition changed')
    for trait, value in [('self_incompatibility','mixed_or_variable'),('autonomous_selfing_capacity','absent')]:
        if value not in ontology['traits'][trait]['allowed_values']:
            raise ValueError('Correction ontology mismatch')
    d.loc[dm,'classification'] = 'partial_self_compatibility'
    d.loc[dm,'normalized_value'] = 'mixed_or_variable'
    d.loc[dm,'state_set'] = '["mixed_or_variable"]'
    d.loc[dm,'state_sets'] = '[["mixed_or_variable"]]'
    d.loc[dm,'source_lineages'] = lineage
    d.loc[am,'quality'] = 'high'
    d.loc[am,'selected_quality'] = 'high'
    d.loc[am,'source_lineages'] = lineage
    d.loc[am,'source_groups'] = 'reviewed_restart_source_correction'
    c.loc[cm,'trait_names'] = 'autonomous_selfing_capacity|self_incompatibility'
    c.loc[cm,'trait_composition'] = 'autonomous_selfing_capacity=["absent"]|self_incompatibility=["mixed_or_variable"]'
    c.loc[cm,'source_lineages'] = lineage
    return c,d


def correct_eulobus(c, d, receipt, ontology):
    """Replace the stranded Meyer SI polarity with Raven's tested SC state."""
    c, d = c.copy(), d.copy()
    species = 'Eulobus californicus'
    axis = 'reproductive_assurance'
    trait = 'self_incompatibility'
    old_lineage = 'dataset:dryad.cc2fqz6hr:v4'
    new_lineage = 'hdl:10088/27111'
    required = {
        'accepted_species': species,
        'trait_name': trait,
        'before_value': 'SI',
        'after_value': 'SC',
        'before_lineage': old_lineage,
        'source_lineage': new_lineage,
        'source_url': 'https://repository.si.edu/items/894bb186-9482-49f0-b0ba-d74736e330d5',
        'review_status': 'accepted_source_correction',
        'reviewer': 'source-scale-methodology-review-20260908',
        'review_date': '2026-09-08',
        'name_match_method': 'accepted_name_via_current_synonym_of_camissonia_californica',
        'support_basis': 'multiple diploid and tetraploid populations were tested by pollen-tube growth and found self-compatible',
    }
    if any(receipt.get(k) != v for k, v in required.items()) or receipt.get('genus_rule_training_allowed') is not False:
        raise ValueError('Unapproved Eulobus source correction')
    if 'SC' not in ontology['traits'][trait]['allowed_values']:
        raise ValueError('Eulobus correction ontology mismatch')
    cm = c.accepted_species.eq(species) & c.axis.eq(axis)
    dm = d.accepted_species.eq(species) & d.trait_name.eq(trait)
    if (cm.sum() != 1 or dm.sum() != 1
            or c.loc[cm,'trait_composition'].iloc[0] != 'self_incompatibility=["SI"]'
            or c.loc[cm,'quality'].iloc[0] != 'medium'
            or d.loc[dm,'normalized_value'].iloc[0] != 'SI'
            or d.loc[dm,'quality'].iloc[0] != 'medium'
            or d.loc[dm,'source_lineages'].iloc[0] != old_lineage):
        raise ValueError('Eulobus correction precondition changed')
    d.loc[dm,'classification'] = 'reviewed_source_correction'
    d.loc[dm,'normalized_value'] = 'SC'
    d.loc[dm,'state_set'] = '["SC"]'
    d.loc[dm,'state_sets'] = '[["SC"]]'
    d.loc[dm,'quality'] = 'high'
    d.loc[dm,'selected_quality'] = 'high'
    d.loc[dm,'source_lineages'] = new_lineage
    d.loc[dm,'source_groups'] = 'reviewed_restart_source_correction'
    c.loc[cm,'quality'] = 'high'
    c.loc[cm,'trait_names'] = trait
    c.loc[cm,'trait_composition'] = 'self_incompatibility=["SC"]'
    c.loc[cm,'source_lineages'] = new_lineage
    c.loc[cm,'source_groups'] = 'reviewed_restart_source_correction'
    return c,d


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--records', type=Path, required=True)
    p.add_argument('--additional-records', type=Path, action='append', default=[])
    p.add_argument('--ontology', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--harrisia-correction', type=Path)
    p.add_argument('--eulobus-correction', type=Path)
    a = p.parse_args()
    if a.output.exists():
        raise ValueError('Output already exists')
    cf = a.baseline/'materialized_species_axis_coverage.csv.gz'
    df = a.baseline/'public_source/wave53-all-evidence-audit/resolved_direct_species_trait.csv.gz'
    for f, expected in [(cf,'0dbe6632142536279bd40ec707f2f73e356fd1b213878f5a5e227c3a37d81677'),(df,'9a6b0a7bec119e68c8fbe05ceadddb84f97526ec3bbfb6f9d89e3b2fb761db79')]:
        if hashlib.sha256(f.read_bytes()).hexdigest() != expected:
            raise ValueError('Pinned baseline checksum mismatch')
    c, d = [pd.read_csv(f, dtype=str).fillna('') for f in (cf, df)]
    if len(c) != 318885 or c.accepted_species.nunique() != 106295 or int(c.quality.ne('').sum()) != 222375:
        raise ValueError('Wrong baseline')
    records = json.loads(a.records.read_text(encoding='utf-8'))
    for extra in a.additional_records:
        records.extend(json.loads(extra.read_text(encoding='utf-8')))
    ontology = yaml.safe_load(a.ontology.read_text(encoding='utf-8'))
    updated, ledger, changes = integrate(c, d, records, ontology)
    before, after = c.quality.ne(''), updated.quality.ne('')
    if not c.loc[before].equals(updated.loc[before]):
        raise ValueError('Modified existing evidence')

    corrections = []
    correction_audit = []
    if a.harrisia_correction:
        receipt=json.loads(a.harrisia_correction.read_text(encoding='utf-8'))
        prior=updated.copy(); prior_ledger=ledger.copy()
        updated,ledger=correct_harrisia(updated,ledger,receipt,ontology)
        allowed=prior.accepted_species.eq('Harrisia portoricensis') & prior.axis.eq('reproductive_assurance')
        if not prior.loc[~allowed].equals(updated.loc[~allowed]) or not prior.quality.equals(updated.quality):
            raise ValueError('Correction modified unrelated coverage')
        for trait in ('self_incompatibility','autonomous_selfing_capacity'):
            old=prior_ledger.loc[prior_ledger.accepted_species.eq('Harrisia portoricensis') & prior_ledger.trait_name.eq(trait)].iloc[0]
            new=ledger.loc[ledger.accepted_species.eq('Harrisia portoricensis') & ledger.trait_name.eq(trait)].iloc[0]
            correction_audit.append(dict(accepted_species=old.accepted_species,trait_name=trait,
                before_value=old.normalized_value,after_value=new.normalized_value,
                before_quality=old.quality,after_quality=new.quality,
                before_lineage=old.source_lineages,after_lineage=new.source_lineages))
        corrections.append(receipt)

    if a.eulobus_correction:
        receipt=json.loads(a.eulobus_correction.read_text(encoding='utf-8'))
        prior=updated.copy(); prior_ledger=ledger.copy()
        updated,ledger=correct_eulobus(updated,ledger,receipt,ontology)
        allowed=prior.accepted_species.eq('Eulobus californicus') & prior.axis.eq('reproductive_assurance')
        if not prior.loc[~allowed].equals(updated.loc[~allowed]):
            raise ValueError('Eulobus correction modified unrelated coverage')
        if int(updated.quality.ne('').sum()) != int(prior.quality.ne('').sum()):
            raise ValueError('Eulobus correction changed filled-cell count')
        old=prior_ledger.loc[prior_ledger.accepted_species.eq('Eulobus californicus') & prior_ledger.trait_name.eq('self_incompatibility')].iloc[0]
        new=ledger.loc[ledger.accepted_species.eq('Eulobus californicus') & ledger.trait_name.eq('self_incompatibility')].iloc[0]
        correction_audit.append(dict(accepted_species=old.accepted_species,trait_name='self_incompatibility',
            before_value=old.normalized_value,after_value=new.normalized_value,
            before_quality=old.quality,after_quality=new.quality,
            before_lineage=old.source_lineages,after_lineage=new.source_lineages))
        corrections.append(receipt)

    a.output.mkdir(parents=True)
    for name, table in [('species_axis_coverage',updated),('direct_species_trait_ledger',ledger),('change_audit',changes)]:
        table.to_csv(a.output/(name+'.csv.gz'),index=False,compression={'method':'gzip','mtime':0})
    (a.output/'reviewed_source_records.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    shutil.copy2(a.baseline/'recovery_manifest.json',a.output/'baseline_recovery_manifest.json')
    if a.harrisia_correction:
        shutil.copy2(a.harrisia_correction,a.output/'reviewed_source_correction.json')
    if corrections:
        (a.output/'reviewed_source_corrections.json').write_text(json.dumps(corrections,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        pd.DataFrame(correction_audit).to_csv(a.output/'source_correction_audit.csv.gz',index=False,compression={'method':'gzip','mtime':0})

    final_filled = updated.quality.ne('')
    correction_paths = [p for p in (a.harrisia_correction,a.eulobus_correction) if p]
    summary = dict(baseline_run=34093932899, baseline_artifact=10007872352,
        denominator_species=106295, denominator_species_axis=len(c),
        baseline_filled=int(before.sum()), filled=int(final_filled.sum()), unresolved=int((~final_filled).sum()),
        direct_species_trait_gain=len(ledger)-len(d), direct_species_axis_gain=int((~before & final_filled).sum()),
        gain=int((~before & final_filled).sum()), loss=int((before & ~final_filled).sum()),
        low_additions=0, low_invalidations=0, low_upgrades=0,
        genus_rebuild_status='deferred_no_genus_training_from_this_batch',
        corrected_existing_direct_values=len(corrections),
        correction_new_direct_traits=0,
        correction_medium_to_high_trait_upgrades=sum(1 for row in correction_audit if row['before_quality']=='medium' and row['after_quality']=='high'),
        unpromoted_harrisia_rule_candidates_held=5 if a.harrisia_correction else 0,
        claim_limit='Incremental direct promotion on historical secondary-Low baseline; not full all-evidence Low revalidation.',
        quality_counts=updated.loc[final_filled,'quality'].value_counts().to_dict(),
        axis_filled=updated.loc[final_filled].groupby('axis').size().to_dict(),
        input_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [cf,df,a.records,a.ontology,*a.additional_records,*correction_paths]},
        output_sha256={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in a.output.iterdir()})
    (a.output/'integration_summary.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    main()
