#!/usr/bin/env python3
"""Export public independent bird species incidence and exact diet/family/genus labels.
No plant data, no restricted botanical identifiers."""
import argparse,hashlib,json,pathlib
from tropical_bird_diet_label_null import prepare_incidence
from tropical_three_guild_candidate_audit import fetch
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    b=fetch("bee_bird");b=b[b.group!="Ryukyu sensitivity"].copy()
    assert len(b)==38 and b.bee_id.is_unique and b.group.nunique()==7
    presence,sp,labels,fam,gen=prepare_incidence(b)
    assert len(sp)==377 and len(presence)==38
    index={s:i for i,s in enumerate(sp)}
    # Use anonymous integer positions instead of transmitting bird names; retain
    # public family/genus labels only as integer block IDs. Shared species remain shared.
    f={v:i for i,v in enumerate(sorted(set(fam)))}
    g={v:i for i,v in enumerate(sorted(set(gen)))}
    out={"source_artifact":14901918629,"source_bird_count":377,"source_island_count":38,
         "diet_codes":labels.tolist(),"family_blocks":[f[x] for x in fam],"genus_blocks":[g[x] for x in gen],
         "positive_incidence":[[str(i),sorted(index[s] for s in v)] for i,v in sorted(presence.items())]}
    text=json.dumps(out,separators=(",",":"),allow_nan=False)
    (a.out/"public_bird_incidence.json").write_text(text)
    print("PUBLIC_BIRD_INCIDENCE_SHA256",hashlib.sha256(text.encode()).hexdigest(),flush=True)
    print("PUBLIC_BIRD_INCIDENCE_JSON",text,flush=True)
if __name__=="__main__":main()
