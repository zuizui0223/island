"""Support-only intersection of the existing signed bee/bird/butterfly sources
and newly computed bat incidence. No fitting on underpowered common cohorts.
Call only in protected temporary analysis; emit anonymous counts, not IDs.
"""
import numpy as np
import pandas as pd
from tropical_public_flora_guild_association import read_artifact
from tropical_three_guild_candidate_audit import fetch

MIN_BLOCKS=30
MIN_INFORMATIVE_ISLANDS=200

def block_support(d,features):
    d=d.dropna(subset=["block",*features])
    groups=[v for _,v in d.groupby("block") if len(v)>=2 and all(v[x].nunique()>1 for x in features)]
    k=len(groups);n=sum(len(v) for v in groups)
    return {"islands":len(d),"informative_blocks":k,"informative_islands":n,
            "gate_pass":k>=MIN_BLOCKS and n>=MIN_INFORMATIVE_ISLANDS}

def same_cohort_support(tab,geo):
    flora=read_artifact()
    bird=fetch("bee_bird")
    butterfly=fetch("butterfly")
    assert len(bird)==52 and bird.bee_id.is_unique
    assert len(butterfly)==98 and butterfly.bee_id.is_unique
    bird=bird[bird.group.ne("Ryukyu sensitivity")&bird.native_bees.gt(0)].copy()
    bird=bird[["bee_id","native_bees","nectar_ge50","nectar_gt0","all_birds","area_ratio","match_type"]]
    bird=bird.rename(columns={"native_bees":"bee_bird_count"})
    assert (bird.all_birds>=bird.nectar_gt0).all() and (bird.nectar_gt0>=bird.nectar_ge50).all()
    butterfly=butterfly[butterfly.native_bees.gt(0)][["bee_id","native_bees","all_records","all_recorded_species","precision_le1km_records","precision_le1km_recorded_species"]].copy()
    butterfly=butterfly.rename(columns={"native_bees":"bee_butterfly_count"})
    out={}
    for scope,f in flora.items():
        assert f.bee_id.is_unique and f.island_id.is_unique
        meta=geo[["island_id","island_latitude","island_longitude"]]
        m=f[["bee_id","island_id","analysis_regime","n_names"]].merge(meta,on="island_id",validate="one_to_one")
        m=m[m.analysis_regime.eq("tropical")].copy()
        m["block"]=(np.floor(m.island_latitude/10).astype("Int64").astype(str)+"_"+np.floor(m.island_longitude/10).astype("Int64").astype(str))
        for nmin in (5,50):
            q=m[m.n_names.ge(nmin)]
            b=q.merge(bird,on="bee_id",validate="one_to_one")
            t=q.merge(butterfly,on="bee_id",validate="one_to_one")
            j=b.merge(butterfly,on="bee_id",validate="one_to_one")
            assert (j.bee_bird_count==j.bee_butterfly_count).all()
            j=j[j.all_records.gt(0)&j.all_recorded_species.gt(0)]
            # These are pre-bat upper bounds. Adding an eligible bat island can
            # never enlarge an existing source intersection.
            pre={"bee_geometry_flora":len(q),"bird_candidate":len(b),
                 "butterfly_positive":int((t.all_records.gt(0)&t.all_recorded_species.gt(0)).sum()),
                 "three_guild_positive_record_cohort":len(j)}
            for cover,c in tab.items():
                for amin in (20,100):
                    bats=c[c.complete&c.range_area.ge(amin)][["island_id","nectar","other","specialist","partial"]]
                    a=q.merge(bats,on="island_id",validate="one_to_one")
                    ab=b.merge(bats,on="island_id",validate="one_to_one")
                    at=t.merge(bats,on="island_id",validate="one_to_one")
                    aj=j.merge(bats,on="island_id",validate="one_to_one")
                    at=at[at.all_records.gt(0)&at.all_recorded_species.gt(0)]
                    assert len(aj)<=len(j) and len(ab)<=len(b)
                    key=f"{scope}_n{nmin}_cover{cover}_area{amin}"
                    out[key]={"before_bat_filter":pre,
                              "bat_bee_geometry_overlap":block_support(a,["nectar","other"]),
                              "bat_bird_bee_candidate":block_support(ab,["nectar","other","nectar_ge50","bee_bird_count"]),
                              "bat_butterfly_bee_positive":block_support(at,["nectar","other","all_recorded_species","bee_butterfly_count"]),
                              "all_four_source_intersection":block_support(aj,["nectar","other","nectar_ge50","bee_bird_count","all_recorded_species"]),
                              "strict_bird_name_area":len(ab[ab.match_type.eq("normalized_name")&ab.area_ratio.between(.9,1.1)]),
                              "butterfly_10_records":int((at.all_records>=10).sum()),
                              "four_source_10_records":int((aj.all_records>=10).sum())}
    return {"panels":out,"thresholds":{"informative_blocks":MIN_BLOCKS,"informative_islands":MIN_INFORMATIVE_ISLANDS},
            "warning":"Eligibility and variation census only. Bird identity provisional; butterfly counts positive records, not complete richness; bat range opportunity not confirmed current occurrence. Bee geometry match alone does not certify a positive native list. No zeros from missing sources, no model fit or causal substitution claim."}

def self_test():
    x=pd.DataFrame({"block":np.repeat(np.arange(30),8),"a":np.tile(np.arange(8),30),"b":np.tile(np.arange(8)[::-1],30)})
    assert block_support(x,["a","b"])=={"islands":240,"informative_blocks":30,"informative_islands":240,"gate_pass":True}
    x["b"]=0
    assert block_support(x,["a","b"])["informative_blocks"]==0
    print("MULTIGUILD_SUPPORT_TEST_PASS",flush=True)

if __name__=="__main__":self_test()
