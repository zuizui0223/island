from __future__ import annotations
import argparse, html, json, re, time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import pandas as pd

AXIS='reproductive_assurance'
TERMS=re.compile(r'self[ -]?(?:compatib|incompatib|pollinat|fertili)|selfing|autogam|cleistogam|breeding system|mating system|outcross|bagged|bagging',re.I)
BLOCKED={'Sideroxylon','Illicium','Portulaca','Eugenia','Liparis','Melaleuca','Callicarpa','Durio','Ourisia','Schoenoplectiella'}


def get_json(url: str) -> dict:
    req=Request(url,headers={'User-Agent':'island-wave58/1.0 (https://github.com/zuizui0223/island)'})
    with urlopen(req,timeout=25) as response: return json.load(response)


def query_epmc(species: str) -> list[tuple[str,str,str]]:
    query=f'TITLE_ABS:"{species}" AND ("self-compatible" OR "self-incompatible" OR selfing OR autogamy OR cleistogamy OR "breeding system" OR "mating system" OR outcrossing)'
    url='https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urlencode({'query':query,'format':'json','resultType':'core','pageSize':20})
    data=get_json(url); out=[]
    for item in data.get('resultList',{}).get('result',[]):
        text=re.sub(r'\s+',' ',html.unescape(re.sub(r'<[^>]+>',' ',item.get('title','')+' '+item.get('abstractText',''))))
        if species.casefold() not in text.casefold() or not TERMS.search(text): continue
        doi=str(item.get('doi','')).lower(); ident=f"{item.get('source','')}:{item.get('id','')}"
        source_url='https://doi.org/'+doi if doi else 'https://europepmc.org/article/'+ident.replace(':','/')
        out.append((('doi:'+doi) if doi else ('epmc:'+ident),source_url,text[:1400]))
    return out


def main() -> None:
    p=argparse.ArgumentParser(); p.add_argument('--coverage',type=Path,required=True); p.add_argument('--direct',type=Path,required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--max-species',type=int,default=160); p.add_argument('--acquire',action='store_true'); a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    cov=pd.read_csv(a.coverage,dtype=str).fillna(''); direct=pd.read_csv(a.direct,dtype=str).fillna(''); direct=direct[direct.axis.eq(AXIS)].copy(); direct['genus']=direct.accepted_species.str.split().str[0]
    groups=direct.groupby(['genus','trait_name']).agg(n_direct=('accepted_species','nunique'),values=('normalized_value',lambda s:'|'.join(sorted(set(s))))).reset_index()
    groups=groups[(groups.n_direct==2)&~groups.genus.isin(BLOCKED)&~groups['values'].str.contains(r'\|',regex=True)]
    unresolved=cov[(cov.axis.eq(AXIS))&(cov.quality.eq(''))].copy(); unresolved['genus']=unresolved.accepted_species.str.split().str[0]; counts=unresolved.groupby('genus').size(); groups['unresolved']=groups.genus.map(counts).fillna(0).astype(int); groups=groups[groups.unresolved>=10].sort_values(['unresolved','genus'],ascending=[False,True])
    groups[['genus','trait_name','values','unresolved']].to_csv(a.out/'coherent_support2_genera.csv',index=False)
    genera=list(dict.fromkeys(groups.genus)); targets=unresolved[unresolved.genus.isin(genera)][['accepted_species','genus']].drop_duplicates(); targets['leverage']=targets.genus.map(counts)
    targets=(targets.sort_values(['genus','accepted_species']).groupby('genus',as_index=False,group_keys=False).head(8).sort_values(['leverage','genus','accepted_species'],ascending=[False,True,True]).head(a.max_species)); targets.to_csv(a.out/'exact_species_targets.csv',index=False)
    rows=[]; logs=[]
    if a.acquire:
        for i,row in targets.reset_index(drop=True).iterrows():
            species=row.accepted_species
            try:
                papers=query_epmc(species); logs.append({'accepted_species':species,'status':'success','lead_documents':len(papers),'error':''})
                for lineage,url,excerpt in papers: rows.append({'accepted_species':species,'axis':AXIS,'source_lineage':lineage,'source_url':url,'source_excerpt':excerpt,'evidence_scope':'unreviewed_species_document','normalized_value':'','evidence_quality':'unreviewed','promotion_allowed':'false','genus_rule_training_allowed':'false','review_status':'needs_fulltext_trait_attribution'})
            except Exception as exc: logs.append({'accepted_species':species,'status':'error','lead_documents':0,'error':type(exc).__name__+': '+str(exc)[:200]})
            if i%20==19: time.sleep(1)
    cols=['accepted_species','axis','source_lineage','source_url','source_excerpt','evidence_scope','normalized_value','evidence_quality','promotion_allowed','genus_rule_training_allowed','review_status']
    pd.DataFrame(rows,columns=cols).drop_duplicates().to_csv(a.out/'species_literature_review_queue.csv',index=False); pd.DataFrame(logs).to_csv(a.out/'query_audit.csv',index=False)
    summary={'contract':'wave58_exact_species_reproductive_discovery_v1','coherent_support2_genus_trait_pairs':len(groups),'target_species':len(targets),'acquisition_executed':a.acquire,'query_success':sum(x['status']=='success' for x in logs),'query_errors':sum(x['status']=='error' for x in logs),'unreviewed_document_leads':len(rows),'lead_species':len({x['accepted_species'] for x in rows}),'automatic_promotions':0,'canonical_net_gain_verified':None}
    (a.out/'wave58_search_summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
