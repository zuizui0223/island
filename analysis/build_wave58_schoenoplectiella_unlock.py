from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd


def main() -> None:
    p=argparse.ArgumentParser()
    p.add_argument('--coverage',type=Path,required=True); p.add_argument('--direct',type=Path,required=True)
    p.add_argument('--wave56',type=Path,required=True); p.add_argument('--out',type=Path,required=True)
    a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    cov=pd.read_csv(a.coverage,dtype=str).fillna(''); direct=pd.read_csv(a.direct,dtype=str).fillna(''); w56=pd.read_csv(a.wave56,dtype=str).fillna('')
    base=direct[(direct.axis=='reproductive_assurance')&direct.accepted_species.str.startswith('Schoenoplectiella ')&(direct.trait_name=='self_incompatibility')]
    base=base[['accepted_species','normalized_value','quality','source_lineages']].rename(columns={'source_lineages':'source_lineage'})
    add=w56[(w56.accepted_species=='Schoenoplectiella juncoides')&(w56.trait_name=='self_incompatibility')]
    assert len(add)==1 and add.iloc[0].normalized_value=='SC' and add.iloc[0].quality in {'high','medium'}
    add=pd.DataFrame([{'accepted_species':'Schoenoplectiella juncoides','normalized_value':'SC','quality':add.iloc[0].quality,'source_lineage':add.iloc[0].source_lineage}])
    support=pd.concat([base,add],ignore_index=True).drop_duplicates(['accepted_species','source_lineage'])
    assert support.accepted_species.nunique()==3 and set(support.normalized_value)=={'SC'} and support.source_lineage.nunique()==3
    dominance=1.0
    species_loo=sum(support.drop(index=i).normalized_value.mode().iat[0]==r.normalized_value for i,r in support.iterrows())/len(support)
    lineage_loo=sum(support[support.source_lineage!=r.source_lineage].normalized_value.mode().iat[0]==r.normalized_value for _,r in support.iterrows())/len(support)
    assert dominance>=.95 and species_loo>=.85 and lineage_loo>=.85
    unresolved=cov[(cov.axis=='reproductive_assurance')&(cov.quality=='')&cov.accepted_species.str.startswith('Schoenoplectiella ')]
    assert len(unresolved)==35 and 'Schoenoplectiella juncoides' in set(unresolved.accepted_species)
    low=unresolved[unresolved.accepted_species!='Schoenoplectiella juncoides'][['accepted_species']].copy()
    low['axis']='reproductive_assurance'; low['trait_name']='self_incompatibility'; low['normalized_value']='SC'; low['quality']='low'; low['evidence_scope']='genus_inference_validated_low'; low['support_n']='3'; low['dominance']='1.000'; low['species_loo_accuracy']=f'{species_loo:.3f}'; low['lineage_loo_accuracy']=f'{lineage_loo:.3f}'; low['promotion_allowed']='false'; low['canonical_collision_audit_complete']='false'
    support.to_csv(a.out/'schoenoplectiella_sc_support.csv',index=False); low.to_csv(a.out/'schoenoplectiella_public_validated_low_candidates.csv',index=False)
    summary={'contract':'wave58_schoenoplectiella_sc_unlock_v1','public_baseline':'Wave52','direct_species':3,'state':'SC','dominance':dominance,'species_loo_accuracy':species_loo,'lineage_loo_accuracy':lineage_loo,'wave52_unresolved_before':35,'direct_candidate_cells':1,'validated_low_candidate_cells':len(low),'public_candidate_cells_total':1+len(low),'canonical_net_gain_verified':None,'automatic_promotions':0}
    (a.out/'schoenoplectiella_unlock_summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
