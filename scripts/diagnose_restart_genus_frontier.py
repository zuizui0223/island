"""Use the common trait-specific audit for a non-promoting affected-genus diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import pandas as pd

p=argparse.ArgumentParser()
p.add_argument('--common-root',type=Path,required=True)
p.add_argument('--recovered-raw',type=Path,required=True)
p.add_argument('--baseline',type=Path,required=True)
p.add_argument('--integrated',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--all-current-frontier',action='store_true',help='Diagnose all current min3/dominance candidates with unresolved axes; never promotes.')
a=p.parse_args()
if a.output.exists():
    raise ValueError('Existing output')
sys.path.insert(0,str(a.common_root/'src'))
from island_v2.all_evidence_trait_audit import build_rule_audit, dedupe_direct_lineages, load_ontology, apply_genus_rules, SETTINGS

records=json.loads((a.integrated/'reviewed_source_records.json').read_text(encoding='utf-8'))
groups={(r['accepted_species'].split()[0],r['trait_name']) for r in records}
if a.all_current_frontier:
    all_direct=pd.read_csv(a.integrated/'direct_species_trait_ledger.csv.gz',dtype=str).fillna('')
    all_direct['genus']=all_direct.accepted_species.str.split().str[0]
    all_coverage=pd.read_csv(a.integrated/'species_axis_coverage.csv.gz',dtype=str).fillna('')
    all_coverage['genus']=all_coverage.accepted_species.str.split().str[0]
    unresolved_groups=set(map(tuple,all_coverage.loc[all_coverage.quality.eq(''),['genus','axis']].to_numpy()))
    groups=set()
    thresholds={'flower_colour':.9,'reproductive_assurance':.95,'floral_structural_complexity':.8}
    for (genus,axis,trait),cells in all_direct.groupby(['genus','axis','trait_name']):
        if (genus,axis) not in unresolved_groups or cells.accepted_species.nunique()<3:
            continue
        if cells.state_set.value_counts().iloc[0]/len(cells)>=thresholds.get(axis,2):
            groups.add((genus,trait))
def affected(df):
    df=df.copy()
    df['genus']=df.accepted_species.str.split().str[0]
    return df.loc[[(g,t) in groups for g,t in zip(df.genus,df.trait_name)]]
before=affected(pd.read_csv(a.baseline/'public_source/wave53-all-evidence-audit/resolved_direct_species_trait.csv.gz',dtype=str).fillna(''))
after=affected(pd.read_csv(a.integrated/'direct_species_trait_ledger.csv.gz',dtype=str).fillna(''))
raw=affected(pd.read_csv(a.recovered_raw,dtype=str).fillna(''))
# Only raw lineage keys actually used by the baseline resolved cells may train.
expected={(r.accepted_species,r.trait_name,s) for r in before.itertuples() for s in r.source_lineages.split('|') if s}
raw=raw.loc[[(s,t,l) in expected for s,t,l in zip(raw.accepted_species,raw.trait_name,raw.source_lineage)]]
observed=set(zip(raw.accepted_species,raw.trait_name,raw.source_lineage))
missing=expected-observed
ontology=load_ontology(a.common_root/'config/trait_ontology.yml')
base_lines,_=dedupe_direct_lineages(raw,ontology)
base_lines['genus']=base_lines.accepted_species.str.split().str[0]
added=[]
for r in records:
    v={col:'' for col in raw.columns}
    v.update(r)
    v['source_group']='reviewed_restart_20260907'
    v['source_record_id']=r['record_id']
    v['source_excerpt']=r['source_excerpt']
    added.append(v)
new_lines,_=dedupe_direct_lineages(pd.DataFrame(added),ontology)
new_lines['genus']=new_lines.accepted_species.str.split().str[0]
old_low=pd.DataFrame(columns=['genus','trait_name','state_set'])
rb=build_rule_audit(before,base_lines,old_low,settings=SETTINGS[:2])
ra=build_rule_audit(after,pd.concat([base_lines,new_lines],ignore_index=True),old_low,settings=SETTINGS[:2])
coverage=pd.read_csv(a.integrated/'species_axis_coverage.csv.gz',dtype=str).fillna('')
unresolved=set(zip(coverage.loc[coverage.quality.eq(''),'accepted_species'],coverage.loc[coverage.quality.eq(''),'axis']))
summary={}
frontiers=[]
for setting in SETTINGS[:2]:
    x=apply_genus_rules(coverage[['accepted_species']].drop_duplicates(),after,ra,setting.name)
    if len(x):
        x=x.loc[[(s,axis) in unresolved for s,axis in zip(x.accepted_species,x.axis)]].copy()
        x['promotion_allowed']=False
        frontiers.append(x)
    oldkeys=set(map(tuple,rb.loc[(rb.setting==setting.name)&rb.eligible,['genus','trait_name']].to_numpy()))
    newkeys=set(map(tuple,ra.loc[(ra.setting==setting.name)&ra.eligible,['genus','trait_name']].to_numpy()))
    summary[setting.name]=dict(eligible_rules=len(newkeys),newly_eligible=len(newkeys-oldkeys),newly_ineligible=len(oldkeys-newkeys),unresolved_axis_candidates=len(x[['accepted_species','axis']].drop_duplicates()) if len(x) else 0)
a.output.mkdir(parents=True)
rb.to_csv(a.output/'rules_before.csv',index=False)
ra.to_csv(a.output/'rules_after.csv',index=False)
queue=ra.loc[ra.setting.eq('current_min3')].copy()
unresolved_counts=coverage.loc[coverage.quality.eq('')].copy()
unresolved_counts['genus']=unresolved_counts.accepted_species.str.split().str[0]
counts=unresolved_counts.groupby(['genus','axis']).size().rename('unresolved_axis_upper_bound').reset_index()
queue=queue.merge(counts,on=['genus','axis'],validate='many_to_one')
queue['priority_axis']=queue.axis.map({'reproductive_assurance':0,'flower_colour':1,'floral_structural_complexity':2})
queue['next_gate']=queue.apply(lambda r:'upstream_receipt_and_conflict_review' if r.eligible else ('independent_original_source_required' if r.lineage_loo_n==0 else 'masked_validation_or_conflict_review'),axis=1)
queue['promotion_allowed']=False
queue.sort_values(['priority_axis','unresolved_axis_upper_bound','genus','trait_name'],ascending=[True,False,True,True]).to_csv(a.output/'source_independence_acquisition_queue.csv',index=False)
pd.concat(frontiers,ignore_index=True).to_csv(a.output/'unpromoted_frontier.csv.gz',index=False) if frontiers else None
pd.DataFrame(sorted(missing),columns=['accepted_species','trait_name','source_lineage']).to_csv(a.output/'missing_lineage_receipts.csv',index=False)
report=dict(status='diagnostic_not_promotion',affected_genus_trait_pairs=len(groups),missing_source_lineage_keys=len(missing),settings=summary,
    common_module_sha256=hashlib.sha256((a.common_root/'src/island_v2/all_evidence_trait_audit.py').read_bytes()).hexdigest(),
    claim_limit='Genus-training flags in new source records remain false. This diagnostic does not authorize their use or revalidate all historical Low. Raw reconstruction and upstream source identities require audit before promotion.')
report['input_sha256']={str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in [
    a.recovered_raw,
    a.baseline/'public_source/wave53-all-evidence-audit/resolved_direct_species_trait.csv.gz',
    a.integrated/'reviewed_source_records.json',
    a.integrated/'direct_species_trait_ledger.csv.gz',
    a.integrated/'species_axis_coverage.csv.gz',
    a.common_root/'config/trait_ontology.yml',
]}
report['source_runs']={'baseline':34093932899,'integrated':34096408983}
report['source_artifacts']={'baseline':10007872352,'integrated':10008766395}
report['promoted_cells']=0
report['no_evaluable_lineage_holdout_rules']=int(queue.lineage_loo_n.eq(0).sum())
report['selection']='all_current_min3_dominance_frontier' if a.all_current_frontier else 'affected_added_records'
(a.output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
