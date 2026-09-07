from __future__ import annotations
import argparse, html, json, re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import pandas as pd

AXIS='reproductive_assurance'
TERMS=re.compile(r'self[ -]?(?:compatib|incompatib|pollinat|fertili)|selfing|autogam|cleistogam|breeding system|mating system|outcross|bagged|bagging',re.I)
# Conflict blocks plus genera whose Wave52 support-two gaps were already recovered in Waves53-55.
BLOCKED={'Sideroxylon','Illicium','Portulaca','Eugenia','Liparis','Melaleuca','Callicarpa','Durio','Ourisia','Schoenoplectiella','Spermacoce','Dicliptera','Gomphrena','Marcgravia','Melicytus','Onopordum'}


def get_json(url: str) -> dict:
    req=Request(url,headers={'User-Agent':'island-wave58/1.2 (https://github.com/zuizui0223/island)'})
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
    p=argparse.ArgumentParser(); p.add_argument('--coverage',type=Path,required=True); p.add_argument('--rules',type=Path,required=True); p.add_argument('--out',type=Path,required=True); p.add_argument('--max-species',type=int,default=160); p.add_argument('--workers',type=int,default=8); p.add_argument('--acquire',action='store_true'); a=p.parse_args(); a.out.mkdir(parents=True,exist_ok=True)
    cov=pd.read_csv(a.coverage,dtype=str).fillna(''); rules=pd.read_csv(a.rules,dtype=str).fillna('')
    unresolved=cov[(cov.axis.eq(AXIS))&(cov.quality.eq(''))].copy(); unresolved['genus']=unresolved.accepted_species.str.split().str[0]; counts=unresolved.groupby('genus').size()
    pairs=rules[(rules.setting.eq('current_min2_diagnostic'))&(rules.axis.eq(AXIS))&(rules.eligible.str.casefold().eq('true'))].copy()
    pairs=pairs[~pairs.genus.isin(BLOCKED)]; pairs['unresolved']=pairs.genus.map(counts).fillna(0).astype(int); pairs=pairs[pairs.unresolved>=10].sort_values(['unresolved','genus','trait_name'],ascending=[False,True,True])
    pairs[['genus','trait_name','inferred_value','n_direct_species','dominance','species_loo_accuracy','lineage_loo_accuracy','unresolved']].to_csv(a.out/'validated_support2_genera.csv',index=False)
    genera=list(dict.fromkeys(pairs.genus)); targets=unresolved[unresolved.genus.isin(genera)][['accepted_species','genus']].drop_duplicates(); targets['leverage']=targets.genus.map(counts)
    targets=(targets.sort_values(['genus','accepted_species']).groupby('genus',as_index=False,group_keys=False).head(8).sort_values(['leverage','genus','accepted_species'],ascending=[False,True,True]).head(a.max_species)); targets.to_csv(a.out/'exact_species_targets.csv',index=False)
    rows=[]; logs=[]
    def one(species: str):
        try: return species,query_epmc(species),''
        except Exception as exc: return species,[],type(exc).__name__+': '+str(exc)[:200]
    if a.acquire:
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futures=[pool.submit(one,s) for s in targets.accepted_species]
            for future in as_completed(futures):
                species,papers,error=future.result(); logs.append({'accepted_species':species,'status':'error' if error else 'success','lead_documents':len(papers),'error':error})
                for lineage,url,excerpt in papers: rows.append({'accepted_species':species,'axis':AXIS,'source_lineage':lineage,'source_url':url,'source_excerpt':excerpt,'evidence_scope':'unreviewed_species_document','normalized_value':'','evidence_quality':'unreviewed','promotion_allowed':'false','genus_rule_training_allowed':'false','review_status':'needs_fulltext_trait_attribution'})
    cols=['accepted_species','axis','source_lineage','source_url','source_excerpt','evidence_scope','normalized_value','evidence_quality','promotion_allowed','genus_rule_training_allowed','review_status']
    pd.DataFrame(rows,columns=cols).drop_duplicates().sort_values(['accepted_species','source_lineage']).to_csv(a.out/'species_literature_review_queue.csv',index=False); pd.DataFrame(logs).sort_values('accepted_species').to_csv(a.out/'query_audit.csv',index=False)
    summary={'contract':'wave58_exact_species_reproductive_discovery_v3','validated_support2_genus_trait_pairs':len(pairs),'target_species':len(targets),'workers':a.workers,'blocked_or_already_recovered_genera':sorted(BLOCKED),'acquisition_executed':a.acquire,'query_success':sum(x['status']=='success' for x in logs),'query_errors':sum(x['status']=='error' for x in logs),'unreviewed_document_leads':len(rows),'lead_species':len({x['accepted_species'] for x in rows}),'automatic_promotions':0,'canonical_net_gain_verified':None}
    (a.out/'wave58_search_summary.json').write_text(json.dumps(summary,indent=2)+'\n'); print(json.dumps(summary,indent=2))


if __name__=='__main__': main()
