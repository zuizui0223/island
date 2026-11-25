#!/usr/bin/env python3
"""Fast source-pinned pair-only flower analysis; no re-download of the 262 MB bee archive."""
import argparse,json,pathlib
import numpy as np
from tropical_public_flora_guild_association import read_artifact,get_public_flora,fit
from tropical_three_guild_candidate_audit import fetch
def prepare(d):
    d=d.copy()
    d=d[(d.area>0)&(d.distance>=0)&(d.native_bees>=0)]
    d["log_area"]=np.log2(d.area)
    d["log_distance"]=np.log2(1+d.distance)
    d["log_n_names"]=np.log2(d.n_names)
    d["log_bees"]=np.log1p(d.native_bees)
    return d
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=pathlib.Path,required=True);a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    br=read_artifact();fl,profile_sha=get_public_flora(br)
    bird=fetch("bee_bird");but=fetch("butterfly")
    assert len(bird)==52 and len(but)==98
    bird=bird.rename(columns={"bee_area_km2":"area","dist_km":"distance"})
    but=but.rename(columns={"archipelago":"group","area_km2":"area","distance_km":"distance"})
    results={}
    cv=["log_area","log_distance","log_n_names"]
    for scope,v in fl.items():
        for n in (5,50):
            f=v[(v.analysis_regime=="tropical")&(v.n_names>=n)]
            z=prepare(f.merge(bird,on="bee_id",validate="one_to_one"))
            z["log_bird"]=np.log1p(z.nectar_gt0);z["log_nonnectar"]=np.log1p(z.all_birds-z.nectar_gt0)
            z["log_bird50"]=np.log1p(z.nectar_ge50)
            assert (z.all_birds>=z.nectar_gt0).all()
            for label,fs in [("nectar",["log_bees","log_bird"]),("nonnectar",["log_bees","log_nonnectar"]),("paired_control",["log_bees","log_bird","log_nonnectar"]),("specialist",["log_bees","log_bird50"])]:
                results[f"{scope}_n{n}_bird_{label}"]=fit(z,fs,cv,draws=299)
            x=prepare(f.merge(but,on="bee_id",validate="one_to_one"))
            for mask in ("all","precision_le1km","no_gbifid_with_reference"):
                rr=mask+"_records";ss=mask+"_recorded_species"
                q=x[(x[rr]>=10)&(x[ss]>0)].copy()
                q["log_butterfly"]=np.log1p(q[ss]);q["log_butterfly_records"]=np.log1p(q[rr])
                for label,c in [("raw",cv),("effort",cv+["log_butterfly_records"])]:
                    results[f"{scope}_n{n}_butterfly_{mask}_{label}"]=fit(q,["log_bees","log_butterfly"],c,draws=299)
    out={"profile_sha256":profile_sha,"panels":results,"warning":"Historical public flower state sets, provisional bird matches, positive-only butterfly records; no effective pollination. Exploratory conditional regression; no confirmatory intervals, no multiple-testing correction."}
    (a.out/"public_pair_results.json").write_text(json.dumps(out,indent=2,allow_nan=False))
    for k,v in results.items():
        compact={x:y for x,y in v.items() if x!="endpoints"}
        compact["endpoints"]={x:y for x,y in v.get("endpoints",{}).items() if x in ("tubular_only","extended_access_shape","red_pink_and_extended_shape","rao_q_without_replacement")}
        print("PAIR_FLORA_PANEL",k,json.dumps(compact),flush=True)
if __name__=="__main__":main()
