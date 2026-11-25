"""Exact-category atlas and finite morphology/reproductive capacity partition.

Original GBIF occurrence keys, own integrated Direct trait ledger, corrected
geography only. No GIFT, status, candidate mainland or external network input.
All outputs are anonymous summaries of recorded species, not effective service.
"""
from __future__ import annotations
import argparse
from collections import Counter
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import hypergeom
from scripts.fd_self_db_capacity_space import normalize_roster, A, SC
from scripts.fd_self_db_isolation import C,F,S,D,DIST,COV,LOCK,REGIONS,digest,fit_cluster
from scripts.fd_self_db_genus_null import genus_key

SIZE="flower_size_class"
DISPLAY="inflorescence_display"
EP={"colour":(C,), "form":(F,), "symmetry":(S,), "depth":(D,),
    "size":(SIZE,), "display":(DISPLAY,),
    "colour_form":(C,F), "colour_symmetry":(C,S), "form_symmetry":(F,S),
    "symmetry_depth":(S,D), "form_depth":(F,D), "colour_depth":(C,D),
    "colour_form_symmetry":(C,F,S)}
N=10
REGION_ALL=("global",*REGIONS)
V_METRICS=("total_mass","positive_only_mass","negative_only_mass","shared_mass",
           "positive_share","excess_random_share","rare_common_excess_gap")
PRED=[DIST,*COV,"log_joint_n","joint_coverage"]


def states(ledger):
    cols=set(t for ts in EP.values() for t in ts)|{A,SC}
    x=ledger.loc[ledger.trait_name.isin(cols)&ledger.evidence_scope.eq("species_direct")
        &ledger.evidence_quality.isin(["high","medium"])&ledger.normalized_value.notna()
        &~ledger.normalized_value.astype(str).str.contains("|",regex=False),
        ["accepted_species","trait_name","normalized_value"]]
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("duplicate trait states")
    d=x.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    d[A]=d[A].map({"autonomous":1.,"delayed":1.,"absent":0.})
    d[SC]=d[SC].map({"SC":1.,"SI":0.})
    return d


def no_hit(pop, count, sample):
    """Exact finite draw zero-success probability, vectorized; sample<=pop."""
    k=np.asarray(count,dtype=float)
    if sample>pop or np.any(k<0) or np.any(k>pop):
        raise ValueError("invalid counts")
    if sample==0:
        return np.ones_like(k)
    # Small sample: product stable, including k=0 or k=pop.
    if sample<=30:
        out=np.ones_like(k)
        for j in range(sample):
            out*=np.maximum(pop-k-j,0)/(pop-j)
        return out
    return hypergeom.pmf(0,pop,k,sample)


def partition(n, near, pos, neg, weights, rarity):
    """Exact expected set masses; ratio of expectations, not expected ratio."""
    t=near.sum(1)
    p=near[:,pos].sum(1)
    m=near[:,neg].sum(1)
    if not np.array_equal(t,p+m):
        raise ValueError("unknown label in complete-carrier partition")
    q=1-no_hit(n,t,N)
    u=1-no_hit(n,p,N)
    v=1-no_hit(n,m,N)
    shared=u+v-q
    # Expected positive coverage under random label placement preserving count.
    j=np.arange(N+1)
    sample_near=hypergeom.pmf(j[None,:],n,t[:,None],N)
    # no_hit has scalar sample size, so precompute each possible draw count.
    expected=1-sample_near@np.array([float(no_hit(n,int(pos.sum()),int(z))) for z in j])
    def mass(x):
        return float(weights@x)
    den=mass(q)
    if den<=0:
        raise ValueError("empty profile space")
    def excess(mask):
        d=mass(q*mask)
        return mass((u-expected)*mask)/d if d>1e-12 else None
    rare=excess(rarity==0); common=excess(rarity==2)
    return {"total_mass":den,"positive_only_mass":mass(q-v),
            "negative_only_mass":mass(q-u),"shared_mass":mass(shared),
            "positive_share":mass(u)/den,"excess_random_share":mass(u-expected)/den,
            "rare_common_excess_gap":rare-common if rare is not None and common is not None else np.nan}


def exact_genus_count(labels, profile, genera):
    """Conditional exact within-genus positive-label shuffle for a category."""
    obs=int(labels[profile].sum())
    expected=0.; offset=0; pmf=np.ones(1); movable=0
    for g in np.unique(genera[profile]):
        idx=genera==g
        n=int(idx.sum()); m=int(labels[idx].sum()); k=int((idx&profile).sum())
        expected+=k*m/n
        low=max(0,m+k-n); high=min(m,k)
        offset+=low
        if high>low:
            movable+=k
            d=hypergeom.pmf(np.arange(low,high+1),n,m,k)
            pmf=np.convolve(pmf,d)
    pmf=pmf/pmf.sum()
    at=obs-offset
    if not 0<=at<len(pmf):
        raise ValueError("observed outside null support")
    p=float(min(1.,2*min(pmf[:at+1].sum(),pmf[at:].sum())))
    return expected,p,movable


def adjust(rows, keys, pkey="p"):
    groups={}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in keys),[]).append(r)
    for group in groups.values():
        m=len(group); harm=sum(1/i for i in range(1,m+1)); prev=1.
        for rank,r in reversed(list(enumerate(sorted(group,key=lambda z:z.get(pkey,1.)),1))):
            prev=min(prev,m*r.get(pkey,1.)/rank)
            r["family_slots"]=m
            if pkey in r:
                r["q_bh"]=float(prev); r["q_by"]=float(min(1,prev*harm))


def quant(x):
    v=np.asarray(x,dtype=float)
    v=v[np.isfinite(v)]
    return [float(np.quantile(v,q)) for q in (.25,.5,.75)] if len(v) else None


def run(roster,traits,geo):
    if roster.island_id.nunique()!=4505 or len(roster)!=1039757 or len(geo)!=8264:
        raise ValueError("source denominator drift")
    if geo.island_id.duplicated().any():
        raise ValueError("duplicate geography")
    geo=geo[["island_id","analysis_regime","spatial_block",DIST,*COV]].copy()
    for c in (DIST,*COV):
        geo[c]=pd.to_numeric(geo[c],errors="coerce")
    totals=roster.groupby("island_id").size()
    joined=roster.merge(traits,left_on="accepted_species",right_index=True,validate="many_to_one")
    joined=joined.merge(geo[["island_id","analysis_regime"]],on="island_id",validate="many_to_one")
    cat=[]; slopes=[]; volumes=[]; volmodels=[]; geom=[]; geommodels=[]
    for ep,cols in EP.items():
        fam="single" if len(cols)==1 else "pair" if len(cols)==2 else "triple"
        morph=traits.dropna(subset=list(cols))
        ref, cnt=np.unique(morph[list(cols)].to_numpy(str),axis=0,return_counts=True)
        names=[" × ".join(r) for r in ref]
        idx={tuple(r):j for j,r in enumerate(ref)}
        rarity=np.where(cnt<=2,0,np.where(cnt>=26,2,1))
        all_morph=joined.dropna(subset=list(cols)).copy()
        all_morph["profile"]=[idx[tuple(r)] for r in all_morph[list(cols)].to_numpy(str)]
        # All true observed profiles, including those not observed in a region;
        # support failures explicitly emitted; no missing data coded as absence.
        for region in REGION_ALL:
            sub=all_morph if region=="global" else all_morph.loc[all_morph.analysis_regime.eq(region)]
            unique=sub.drop_duplicates("accepted_species")
            for flag,mode in ((A,"autonomy"),(SC,"SC")):
                known=unique.loc[unique[flag].notna()]
                labs=known[flag].to_numpy(int)
                g=np.array([genus_key(x) for x in known.accepted_species],dtype=str)
                cp=known.profile.to_numpy(int)
                for k,label in enumerate(names):
                    mask=cp==k
                    total=int(unique.profile.eq(k).sum())
                    n=int(mask.sum()); p=int(labs[mask].sum())
                    row={"endpoint":ep,"family":fam,"region":region,"category":label,"capacity":mode,
                         "global_reference_species":int(cnt[k]),"recorded_species":total,
                         "known_species":n,"positive_species":p,"negative_species":n-p,
                         "positive_fraction":p/n if n else None,
                         "lower_unknown_bound":p/total if total else None,
                         "upper_unknown_bound":(p+total-n)/total if total else None}
                    if n>=10:
                        e,pv,movable=exact_genus_count(labs,mask,g)
                        row.update(genus_expected_positive=e,genus_excess_fraction=(p-e)/n,
                                   movable_species=movable)
                        if movable>=5:
                            row["p"]=pv
                    cat.append(row)
            if region=="global":
                continue
            counts=pd.crosstab(sub.island_id,sub.profile).reindex(columns=range(len(ref)),fill_value=0)
            counts=counts.loc[counts.sum(1)>=N]
            n=counts.sum(1)
            base=pd.DataFrame({"island_id":counts.index,"log_joint_n":np.log(n.to_numpy()),
                "joint_coverage":n.to_numpy()/totals.loc[counts.index].to_numpy()})
            base=base.merge(geo,on="island_id",validate="one_to_one")
            # Full source-positive phenotype share, not GIFT or native comparison.
            for k,label in enumerate(names):
                v=counts[k].to_numpy()/n.to_numpy()
                part=base.assign(y=v)
                positive=v>0
                pos_n=int(positive.sum())
                pos_blocks=int(part.loc[positive,"spatial_block"].nunique())
                fit=fit_cluster(part,"y",PRED) if pos_n>=10 and pos_blocks>=5 else {
                    "status":"rare_category_not_testable","n":len(part),"blocks":int(part.spatial_block.nunique())}
                slopes.append({"endpoint":ep,"family":fam,"region":region,"category":label,
                    "positive_islands":pos_n,"positive_blocks":pos_blocks,
                    "median_share":float(np.median(v)) if len(v) else None,**fit})
        # Finite metric mass, all morphology independent of capacity. Exact
        # category radius and one mismatch radius. Reuse same global reference.
        for radius in ([0] if len(cols)==1 else [0,1]):
            geomrec=[]; rec=[]
            # profile x profile lookup fixed once
            adj=(ref[:,None,:]!=ref[None,:,:]).sum(2)<=radius
            degree=adj.sum(1)
            for island,p in all_morph.groupby("island_id",sort=True):
                if len(p)<N:
                    continue
                ids=p.profile.to_numpy(int)
                near=adj[:,ids]
                k=len(set(ids))
                unif=np.ones(len(ref))/len(ref)
                volume=float(unif@(1-no_hit(len(p),near.sum(1),N)))
                full=float(unif@near.any(1))
                exp=float(unif@(1-no_hit(len(ref),degree,k)))
                geomrec.append({"island_id":island,"log_joint_n":np.log(len(p)),
                    "joint_coverage":len(p)/totals[island],"mass_n10":volume,
                    "geometry_excess_full":full-exp})
                # Same known capacity carriers; retain 0/1 all-positive/negative.
                for flag,mode in ((A,"autonomy"),(SC,"SC")):
                    arr=p[flag].to_numpy(float)
                    known=np.isfinite(arr); n=int(known.sum())
                    if n<N:
                        continue
                    pos=arr[known]==1; neg=~pos
                    for weight,w in (("uniform",unif),("species_frequency",cnt/cnt.sum())):
                        d=partition(n,near[:,known],pos,neg,w,rarity)
                        rec.append({"island_id":island,"capacity":mode,"weight":weight,
                            "log_joint_n":np.log(n),"joint_coverage":n/totals[island],**d})
            g=pd.DataFrame(geomrec,columns=["island_id","log_joint_n","joint_coverage","mass_n10","geometry_excess_full"]).merge(geo,on="island_id",validate="one_to_one")
            df=pd.DataFrame(rec,columns=["island_id","capacity","weight","log_joint_n","joint_coverage",*V_METRICS]).merge(geo,on="island_id",validate="many_to_one")
            for region in REGIONS:
                a=g.loc[g.analysis_regime.eq(region)]
                geom.append({"endpoint":ep,"family":fam,"radius":radius,"region":region,
                             "n_islands":len(a),"reference_profiles":len(ref),
                             "mass_quartiles":quant(a.mass_n10),
                             "geometry_excess_quartiles":quant(a.geometry_excess_full)})
                for metric in ("mass_n10","geometry_excess_full"):
                    geommodels.append({"endpoint":ep,"family":fam,"radius":radius,"region":region,
                        "metric":metric,**fit_cluster(a,metric,PRED)})
                for mode in ("autonomy","SC"):
                    for weight in ("uniform","species_frequency"):
                        a=df.loc[df.analysis_regime.eq(region)&df.capacity.eq(mode)&df.weight.eq(weight)]
                        volumes.append({"endpoint":ep,"family":fam,"radius":radius,"region":region,
                            "capacity":mode,"weight":weight,"n_islands":len(a),
                            "reference_profiles":len(ref),"quartiles":{m:quant(a[m]) for m in V_METRICS}})
                        for metric in V_METRICS:
                            volmodels.append({"endpoint":ep,"family":fam,"radius":radius,"region":region,
                                "capacity":mode,"weight":weight,"metric":metric,**fit_cluster(a,metric,PRED)})
        print("ATLAS_ENDPOINT "+ep,flush=True)
    adjust(cat,["region","capacity","family"])
    adjust(slopes,["family"])
    adjust(geommodels,["family","radius"])
    adjust(volmodels,["family","radius","weight","capacity"])
    return {"status":"OWN_GBIF_CATEGORY_CAPACITY_ATLAS_NOT_REALIZED_FUNCTION",
        "source":"raw GBIF occurrence keys, integrated species-Direct evidence, corrected geography",
        "roster_islands":4505,"roster_rows":len(roster),"geo_islands":8264,
        "category_atlas":cat,"category_slopes":slopes,"geometry":geom,"geometry_models":geommodels,
        "capacity_partitions":volumes,"capacity_models":volmodels,
        "source_sha256":{"ledger":LOCK["ledger"],"geography":LOCK["geography"]},
        "raw_gbif_blob":"d311557334579a7efbef838cc01f927339d16f77",
        "no_GIFT_tables":True,"no_row_level_export":True,"not_causal":True}


def main():
    p=argparse.ArgumentParser()
    for k in ("roster","ledger","geography","output"):
        p.add_argument("--"+k,required=True,type=Path)
    a=p.parse_args()
    if digest(a.ledger)!=LOCK["ledger"] or digest(a.geography)!=LOCK["geography"]:
        raise ValueError("source hashes")
    data=run(normalize_roster(pd.read_csv(a.roster,dtype=str)),
             states(pd.read_csv(a.ledger,dtype=str)),
             pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str}))
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(data,allow_nan=False,indent=2)+"\n")
    for key in ("category_atlas","category_slopes","geometry","geometry_models","capacity_partitions","capacity_models"):
        rows=data[key]
        print("ATLAS_COUNT "+json.dumps({"key":key,"rows":len(rows),
            "BY_supported":sum(r.get("q_by",1)<.05 for r in rows)}))
        for r in rows:
            if r.get("q_by",1)<.05:
                print("ATLAS_BY "+json.dumps({"key":key,**r}))
    print("ATLAS_COMPLETE "+json.dumps({"roster_islands":data["roster_islands"],"roster_rows":data["roster_rows"]}))

if __name__=="__main__":
    main()
