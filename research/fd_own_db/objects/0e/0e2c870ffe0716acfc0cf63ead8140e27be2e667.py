#!/usr/bin/env python3
"""Exploratory matched specialist-versus-partial bird floral contrast.
All animal count columns are disjoint. Public historical state-set traits only."""
import argparse,json,pathlib
import numpy as np
from tropical_public_flora_guild_association import read_artifact,get_public_flora,fit
from tropical_three_guild_candidate_audit import fetch
from tropical_public_flora_guild_pair_only import prepare
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=pathlib.Path);a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    fl,sha=get_public_flora(read_artifact());bird=fetch("bee_bird").rename(columns={"bee_area_km2":"area","dist_km":"distance"})
    assert len(bird)==52 and bird.bee_id.is_unique
    out={}
    for scope,f in fl.items():
        for n in (5,50):
            b=prepare(f[(f.analysis_regime=="tropical")&(f.n_names>=n)].merge(bird,on="bee_id",validate="one_to_one"))
            b["specialist"]=b.nectar_ge50;b["partial"]=b.nectar_gt0-b.nectar_ge50;b["other"]=b.all_birds-b.nectar_gt0
            assert (b[["specialist","partial","other"]]>=0).all().all()
            assert (b[["specialist","partial","other"]].sum(axis=1)==b.all_birds).all()
            for c in ("specialist","partial","other"):b["log_"+c]=np.log1p(b[c])
            b["log_feeder_component_sum"]=b.log_specialist+b.log_partial
            for cohort,z in [("all",b),("tight",b[b.area_ratio.between(.9,1.1)&b.match_type.eq("normalized_name")]),("exclude_shared_boundary",b[~b.bird_name.isin(["St.Martin","Guadeloupe"])])]:
                cv=["log_area","log_distance","log_n_names"]
                key=f"{scope}_n{n}_{cohort}"
                out[key+"_disjoint"]=fit(z,["log_bees","log_specialist","log_partial","log_other"],cv,draws=499)
                out[key+"_contrast"]=fit(z,["log_bees","log_specialist","log_feeder_component_sum","log_other"],cv,draws=499)
    result={"source_profile_sha256":sha,"panels":out,"contrast":"In contrast model coefficient of log_specialist equals disjoint beta_specialist minus beta_partial; not a fixed-total-count substitution or causal estimand.","warning":"Exploratory sparse groups, bird geography provisional. Category state-set membership not actual measured depth or pollination."}
    (a.out/"specialist_partial_flora.json").write_text(json.dumps(result,indent=2,allow_nan=False))
    for k,v in out.items():
        print("SPECIALIST_PANEL",k,json.dumps({x:y for x,y in v.items() if x!="endpoints"}|{"endpoints":{x:y for x,y in v.get("endpoints",{}).items() if x in ("tubular_only","extended_access_shape","red_pink_and_extended_shape","rao_q_without_replacement")}}),flush=True)
if __name__=="__main__":main()
