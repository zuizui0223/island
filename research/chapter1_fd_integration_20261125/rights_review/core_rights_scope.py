#!/usr/bin/env python3
"""Overlap of frozen H1 Direct binary reproductive states and v1 axis-rights ledger.
Counts are audit proxies, not legal decisions or proof of original H1 lineage identity.
No trait values, taxa, or source URLs are exported in aggregate JSON. Optional private
scope CSVs contain only source references for custodian review; never auto-send them.
"""
import argparse, collections as C, hashlib, json, zipfile
from pathlib import Path
import pandas as pd
PINS={
 '10058653212':('10058653212-chapter1-progressive-analysis-34232450884.zip','b7d7f357d9d4d93062abc4671f27434725f6e0da07c48c1c574c2b601c64f999'),
 '10042377593':('10042377593-source-scale-batch-integration-34191508045.zip','5e9c887e093c42eb720f3abd2dc85701ce123fb7fe888c916dc4098ac898ccdd'),
 '10186293306':('10186293306-chapter1-database-v1-cell-rights-34566271397.zip','6c5f2577e85e95e841b8659427e27c72ce631f59a3da5a335f54f505de5ea32a'),
 '10154261753':('10154261753-chapter1-database-v1-rights-audit-34482428564.zip','5d3e5ea80e1cf9757cc7b3bbb64db2ca67f0b8ecce7f7f52668472511f07ad36')}
CORE={'self_incompatibility':{'SC','SI'},'mating_system':{'predominantly_outcrossing','predominantly_selfing','obligate_selfing'},'autonomous_selfing_capacity':{'autonomous','delayed','absent'}}
TARGETS={'biolflor':'database:floraweb_biolflor','pladias':'database:pladias','europe_pmc':'domain:europepmc.org','go_botany':'domain:gobotany.nativeplanttrust.org','ken_fern':'ken-fern-useful-plants:species'}
ALLOW=set(TARGETS.values())|{'unresolved:other_doi','unresolved:citation_hash','database:ecoflora'}
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(z,n):return pd.read_csv(z.open(n),compression='gzip' if n.endswith('.gz') else None,dtype=str,keep_default_na=False)
def tokens(s):return set(x for x in s.split('|') if x)
def label(s):return s if s in ALLOW else 'withheld_'+hashlib.sha256(s.encode()).hexdigest()[:16]
def run(root,out,private):
 zs={}
 for k,(n,h) in PINS.items():
  p=root/n;assert sha(p)==h,k;zs[k]=zipfile.ZipFile(p)
 state='snapshot/chapter1_trait_state_audit_direct.csv.gz'
 assert hashlib.sha256(zs['10058653212'].read(state)).hexdigest()=='4890da7df53b2e240c4b1493ef93c0814a68435a5e1670e3c0e7c3e6717b90e7'
 a=read(zs['10058653212'],state);v=read(zs['10042377593'],'direct_species_trait_ledger.csv.gz');r=read(zs['10186293306'],'CELL_RELEASE_RIGHTS.csv.gz');l=read(zs['10154261753'],'SOURCE_LINEAGE_DETAILS.csv');q=read(zs['10186293306'],'CITATION_PROVENANCE_AUDIT.csv')
 mp=dict(zip(l.source_lineage,l.source_family));mp.update(dict(zip(q.source_lineage,q.source_family)))
 key=['accepted_species','trait_name'];assert not a.duplicated(key).any() and not v.duplicated(key).any();assert not r.duplicated(['accepted_species','axis']).any()
 a=a[a.trait_name.isin(CORE) & a.resolved_for_primary.eq('True')].copy();a=a[a.apply(lambda s:s.canonical_signature in CORE[s.trait_name],axis=1)]
 v0=v[v.axis.eq('reproductive_assurance')].copy();v0['families']=v0.source_lineages.map(lambda s:{mp.get(x,'unknown') for x in tokens(s)})
 v=v0[v0.trait_name.isin(CORE)].copy();v['v1_signature']=v.state_set.map(lambda s:'|'.join(sorted(json.loads(s))))
 joined=a.merge(v,on=key,how='left',validate='one_to_one');m=joined[joined.canonical_signature.eq(joined.v1_signature)].copy();m=m.merge(r[['accepted_species','axis','release_status','blocked_source_families']],on=['accepted_species','axis'],how='left',validate='many_to_one');assert m.release_status.notna().all()
 def counts(f):return {str(k):int(n) for k,n in f.value_counts().items()}
 source={}
 for name,family in TARGETS.items():
  direct=v0[v0.families.map(lambda f:family in f)];same=m[m.families.map(lambda f:family in f)];affected=m[m.blocked_source_families.map(lambda s:family in tokens(s))];sole=affected[affected.blocked_source_families.eq(family)]
  source[name]={'family':family,'v1_direct_reproductive_trait_rows':len(direct),'v1_direct_by_trait':counts(direct.trait_name),'exact_state_overlap_core_rows':len(same),'exact_state_overlap_by_trait':counts(same.trait_name),'matched_core_rows_with_axis_blocker':len(affected),'matched_core_rows_sole_axis_blocker':len(sole),'historical_lineage_count':int(l.source_family.eq(family).sum())}
  if private is not None:
   items=[]
   for row in direct.itertuples():
    for t in sorted(tokens(row.source_lineages)):
     if mp.get(t)==family:items.append({'accepted_species':row.accepted_species,'trait_name':row.trait_name,'quality':row.quality,'own_source_lineage':t})
   pd.DataFrame(items,columns=['accepted_species','trait_name','quality','own_source_lineage']).drop_duplicates().to_csv(private/(name+'_PRIVATE_scope.csv'),index=False)
 sets=C.Counter(frozenset(tokens(s)) for s in m.loc[m.release_status.eq('rights_review_required'),'blocked_source_families']);selected=set();steps=[]
 for _ in range(8):
  gains=C.Counter()
  for b,n in sets.items():
   left=b-selected
   if len(left)==1:gains[next(iter(left))]+=n
  if not gains:break
  k,g=min(gains.items(),key=lambda kv:(-kv[1],kv[0]));selected.add(k);steps.append({'family':label(k),'marginal_matched_trait_rows':int(g)})
 quality=r[r.axis.eq('reproductive_assurance') & r.quality.isin(['high','medium'])]
 result={'pins':{k:h for k,(_,h) in PINS.items()},'frozen_H1_informative_single_state_rows':len(a),'H1_by_trait':counts(a.trait_name),'exact_v1_state_matches':len(m),'unmatched_rows':len(a)-len(m),'matched_unique_species':int(m.accepted_species.nunique()),'matched_status_counts':counts(m.release_status),'matched_by_trait':{t:counts(g.release_status) for t,g in m.groupby('trait_name')},'high_medium_reproductive_axis_cells':len(quality),'high_medium_authorized':int(quality.release_status.eq('redistributable').sum()),'sources':source,'greedy_matched_trait_rows':steps,'interpretation':'Same species, trait and exact single state only. v1 axis-cell rights status is a conservative proxy, not exact rights clearance for the historical H1/H2 records. The H1 audit lacks source lineages; unmatched rows are not adjudicated. Greedy steps are conditional hypothetical complete blocker removal, not approvals.'}
 assert len(a)==12401 and len(m)==10654 and result['matched_status_counts']=={'rights_review_required':9376,'redistributable':1278}
 out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--private-scope-dir',type=Path);a=p.parse_args();print(json.dumps(run(a.source_dir,a.out,a.private_scope_dir),sort_keys=True))
