"""Export diagnostic search leads, never validated rules or claimed gains."""
import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
if a.output.exists():
    raise ValueError('Output exists')
direct = a.baseline/'public_source/wave53-all-evidence-audit/resolved_direct_species_trait.csv.gz'
coverage = a.baseline/'materialized_species_axis_coverage.csv.gz'
d = pd.read_csv(direct, dtype=str).fillna('')
c = pd.read_csv(coverage, dtype=str).fillna('')
assert len(c) == 318885 and c.accepted_species.nunique() == 106295
assert not c.duplicated(['accepted_species','axis']).any()
d = d[d.axis.eq('reproductive_assurance') & d.quality.isin(['high','medium'])]
u = c[c.axis.eq('reproductive_assurance') & c.quality.eq('')].copy()
u['genus'] = u.accepted_species.str.split().str[0]
rows = []
for (genus, trait), g in d.groupby(['genus','trait_name'], sort=True):
    if g.accepted_species.nunique() != 2 or g.normalized_value.nunique() != 1:
        continue
    targets = sorted(u.loc[u.genus.eq(genus),'accepted_species'])
    if not targets:
        continue
    rows.append(dict(genus=genus, trait_name=trait, support_species=2,
        value=g.normalized_value.iloc[0], unresolved_axis_species=len(targets),
        supporting_species=json.dumps(sorted(g.accepted_species.unique().tolist())),
        source_lineages=json.dumps(sorted(g.source_lineages.unique().tolist())),
        unresolved_targets=json.dumps(targets), promotion_allowed=False,
        limitation='Diagnostic only; lineage independence, counterevidence and masked validation not established. Missing trait coverage is not measured by axis gaps.'))
q = pd.DataFrame(rows).sort_values(['unresolved_axis_species','genus','trait_name'],ascending=[False,True,True])
a.output.mkdir(parents=True)
q.to_csv(a.output/'support_two_reproductive_search_queue.csv',index=False)
(a.output/'manifest.json').write_text(json.dumps(dict(baseline_run=34093932899,
    baseline_artifact=10007872352, diagnostic_pairs=len(q), new_accepted_cells=0,
    hashes={str(f.relative_to(a.baseline)):hashlib.sha256(f.read_bytes()).hexdigest() for f in [direct,coverage]}),indent=2)+'\n')
print(q[['genus','trait_name','unresolved_axis_species']].head(10).to_string(index=False))
