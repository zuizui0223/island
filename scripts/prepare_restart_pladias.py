"""Select already accepted Pladias records with cached identity and excerpt receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from integrate_reviewed_restart import STRICT_AXES

p=argparse.ArgumentParser()
p.add_argument('--package',type=Path,required=True)
p.add_argument('--baseline',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.output.exists():
    raise ValueError('Existing output')
cf=a.package/'trait_candidates.csv.gz'
d=pd.read_csv(cf,dtype=str).fillna('')
c=pd.read_csv(a.baseline/'materialized_species_axis_coverage.csv.gz',dtype=str).fillna('')
d['axis']=d.trait_name.map(STRICT_AXES)
d=d.merge(c[['accepted_species','axis','quality']],on=['accepted_species','axis'],validate='many_to_one')
d=d[d.quality.eq('')]
cache={}
for f in sorted((a.package/'page_cache').glob('*.json')):
    obj=json.loads(f.read_text(encoding='utf-8'))
    cache[obj['accepted_species']]=(obj,f)
records=[]
for _,r in d.iterrows():
    if r.review_status!='accepted_source_package_contract' or r.needs_human_review!='false' or r.confidence!='high' or r.evidence_scope!='species_direct':
        raise ValueError('Pending or inferred candidate')
    page,f=cache[r.accepted_species]
    if page['identity_status']!='accepted_name_and_family_exact' or page['matched_page_name']!=r.accepted_species or page['page_family']!=r.family:
        raise ValueError('Identity mismatch')
    features=[x for x in page['features'] if x['supporting_excerpt']==r.evidence_excerpt and x['raw_value']==r.raw_value]
    if len(features)!=1:
        raise ValueError('Missing unique cached excerpt')
    states=json.loads(r.standardized_value) if r.standardized_value.startswith('[') else [r.standardized_value]
    records.append(dict(record_id=r.evidence_id,accepted_species=r.accepted_species,
        axis=r.axis,trait_name=r.trait_name,normalized_value='|'.join(states),state_set=states,
        quality='high',source_url=r.source_url,source_excerpt=r.evidence_excerpt,
        source_lineage=r.source_lineage,source_citation=r.source_citation,
        raw_value=r.raw_value,name_match_method='accepted_species_and_family_exact_cached',
        review_date='2026-09-07',retrieved_at_utc=page['retrieved_at_utc'],
        reviewer='Codex-cached-source-contract-review',review_status='accepted_direct_statement',
        cultivar_status='species_database_record',genus_rule_training_allowed=False,
        source_run_id=34086167061,source_artifact_id=10005354522,
        page_content_sha256=page['content_sha256'],cache_sha256=hashlib.sha256(f.read_bytes()).hexdigest(),
        package_candidates_sha256=hashlib.sha256(cf.read_bytes()).hexdigest(),
        audit_limit='Cached structured evidence and identity rechecked, not independent biological measurement. Pladias required citation retained; contact governing board before large-scale publication.'))
a.output.write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(records=len(records),cells=len({(r['accepted_species'],r['axis']) for r in records}))))
