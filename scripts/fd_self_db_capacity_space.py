"""Own GBIF roster × Direct floral/reproductive labels: potential autonomy coverage.

NO GIFT island checklist, mainland donor, or floristic status read. The metric
measures recorded species-wide capacity and floral category support, NOT
pollinator-free realized seed set or observed ecosystem function.
"""
from __future__ import annotations
import argparse
from collections import Counter
from functools import lru_cache
import json
from math import log
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import hypergeom

from scripts.fd_self_db_isolation import C,F,S,D,DIST,COV,LOCK,REGIONS,digest,fit_cluster
from scripts.fd_self_db_genus_null import genus_key

A="autonomous_selfing_capacity"
SC="self_incompatibility"
ENDPOINTS={"colour":(C,), "form":(F,), "symmetry":(S,), "depth":(D,),
           "colour_form":(C,F), "form_symmetry":(F,S), "symmetry_depth":(S,D),
           "colour_form_symmetry":(C,F,S)}
N=10
MODES=("autonomy_known","SC_known","joint_known_autonomy","joint_known_SC")
OUTCOMES=("share_n10","excess_random_n10","excess_genus_full")
PRED=[DIST,*COV,"log_joint_n","joint_coverage"]


def load_states(ledger):
    x=ledger.loc[
        ledger.trait_name.isin([C,F,S,D,A,SC])
        & ledger.evidence_scope.eq("species_direct")
        & ledger.evidence_quality.isin(["high","medium"])
        & ledger.normalized_value.notna()
        & ~ledger.normalized_value.astype(str).str.contains("|",regex=False),
        ["accepted_species","trait_name","normalized_value"]].copy()
    if x.duplicated(["accepted_species","trait_name"]).any():
        raise ValueError("duplicate species-trait labels")
    out=x.pivot(index="accepted_species",columns="trait_name",values="normalized_value")
    out[A]=out[A].map({"autonomous":1.,"delayed":1.,"absent":0.})
    out[SC]=out[SC].map({"SC":1.,"SI":0.})
    return out


def normalize_roster(raw):
    if "accepted_species" not in raw and "species" in raw:
        raw=raw.rename(columns={"species":"accepted_species"})
    x=raw[["island_id","accepted_species"]].copy()
    for c in x.columns:
        x[c]=x[c].fillna("").astype(str).str.strip()
    return x.loc[x.island_id.ne("")&x.accepted_species.ne("")].drop_duplicates()


@lru_cache(maxsize=200000)
def p_no(N,K,n):
    """Probability zero of K marked individuals in n without replacement."""
    if not (0<=K<=N and 0<=n<=N):
        raise ValueError("bad hypergeometric counts")
    if K==0 or n==0:
        return 1.
    if N-K<n:
        return 0.
    def lc(a,b):
        return gammaln(a+1)-gammaln(b+1)-gammaln(a-b+1)
    return float(np.clip(np.exp(lc(N-K,n)-lc(N,n)),0,1))


def expected_n10(N,counts):
    return np.array([1-p_no(N,int(k),10) for k in counts])


def calc_metrics(near,labels,genera):
    """near: reference profile × scored species boolean incidence.
    labels: 0 or 1, no missingness. Uniform reference profile weighting.
    Two nulls: exact same total positives randomized over whole roster, and
    exact within-genus positives held on the full roster. No causal null.
    """
    labels=np.asarray(labels,dtype=int)
    n=len(labels)
    if n<N or not np.isin(labels,[0,1]).all():
        raise ValueError("insufficient or nonbinary")
    total=near.sum(1)
    k=near[:,labels==1].sum(1)
    m=int(labels.sum())
    den=expected_n10(n,total).sum()
    if den<=0:
        raise ValueError("empty morphology support")
    observed=expected_n10(n,k).sum()/den
    # Number of positive-labelled species near each reference profile under
    # globally uniform fixed-m assignment is hypergeometric.
    rand=0.
    for t in total:
        t=int(t)
        if t==0:
            continue
        if m in (0,n):
            rand+=1-p_no(n,t if m==n else 0,N)
            continue
        z=np.arange(max(0,m+t-n),min(m,t)+1)
        pmf=hypergeom.pmf(z,n,t,m)
        rand+=float(pmf@expected_n10(n,z))
    random_share=rand/den
    # Full roster genus-conditioned coverage: retain positive count of every
    # genus. Product of no-positive probabilities across genera is exact.
    occupied=total>0
    no=np.ones(near.shape[0])
    g=np.asarray(genera)
    for key in np.unique(g):
        mask=g==key
        ng=int(mask.sum())
        mg=int(labels[mask].sum())
        if mg==0:
            continue
        for j,t in enumerate(near[:,mask].sum(1)):
            no[j]*=p_no(ng,int(t),mg)
    full_obs=float((k>0).sum()/occupied.sum())
    full_exp=float((1-no).sum()/occupied.sum())
    return {"share_n10":float(observed),"random_share_n10":float(random_share),
            "excess_random_n10":float(observed-random_share),
            "share_full":full_obs,"genus_expected_share_full":full_exp,
            "excess_genus_full":full_obs-full_exp,"positive_fraction":float(m/n),
            "volume_all_n10":float(den/len(total))}


def coverage_bounds(near,labels):
    """Bounds on reported capacity-supported coverage with unknown labels.
    Not population-level bounds or a missingness model."""
    occupied=near.any(1)
    if not occupied.any():
        raise ValueError("empty support")
    pos=np.asarray(labels)==1
    unknown=pd.isna(labels)
    lo=near[:,pos].any(1).sum()/occupied.sum()
    hi=near[:,pos|unknown].any(1).sum()/occupied.sum()
    return float(lo),float(hi)


def fixed_adjust(rows):
    # All 8 endpoints x 4 regions x 3 response slots per mode/radius family;
    # retain untestable slots as p=1.
    for key in {(r["mode"],r["radius"]) for r in rows}:
        group=[r for r in rows if (r["mode"],r["radius"])==key]
        ranked=sorted(group,key=lambda r:r.get("p",1.))
        m=len(group); harmonic=sum(1/i for i in range(1,m+1)); prev=1.
        for rank,r in reversed(list(enumerate(ranked,1))):
            prev=min(prev,m*r.get("p",1.)/rank)
            r["family_slots"]=m
            if r["status"]=="fit":
                r["q_bh"]=float(prev)
                r["q_by"]=float(min(1,prev*harmonic))


def interval_summary(a):
    if not len(a):
        return None
    return [float(np.quantile(a,q)) for q in (.25,.5,.75)]


def run(roster,states,geo):
    if roster.island_id.nunique()!=4505 or len(geo)!=8264 or geo.island_id.duplicated().any():
        raise ValueError("frozen source or geography count changed")
    geo=geo[["island_id","analysis_regime","spatial_block",DIST,*COV]].copy()
    for c in (DIST,*COV):
        geo[c]=pd.to_numeric(geo[c],errors="coerce")
    totals=roster.groupby("island_id").size()
    joined=roster.merge(states,left_on="accepted_species",right_index=True,how="inner",validate="many_to_one")
    support=[]; rows=[]; bounds=[]
    for endpoint,cols in ENDPOINTS.items():
        complete=states.dropna(subset=list(cols))
        ref=np.unique(complete[list(cols)].to_numpy(str),axis=0)
        mor=joined.dropna(subset=list(cols))
        regions=geo.set_index("island_id").analysis_regime
        radius_list=[0] if len(cols)==1 else [0,1]
        parts=[(i,p.sort_values("accepted_species")) for i,p in mor.groupby("island_id",sort=True) if len(p)>=N]
        for radius in radius_list:
            private=[]
            bound_records=[]
            for island,p in parts:
                if island not in regions.index:
                    continue
                # Distinct actual global Direct morphology combinations; no
                # Cartesian profiles invented. Radius 1 = one category mismatch.
                x=p[list(cols)].to_numpy(str)
                near=(ref[:,None,:]!=x[None,:,:]).sum(2)<=radius
                av=p[A].to_numpy(float); sc=p[SC].to_numpy(float)
                lo,hi=coverage_bounds(near,av)
                bound_records.append({"region":regions[island],"lower":lo,"upper":hi,
                    "unknown_frac":float(np.isnan(av).mean()),"n":len(p)})
                for mode in MODES:
                    flag=av if mode.endswith("autonomy") or mode=="autonomy_known" else sc
                    mask=np.isfinite(flag)
                    if mode.startswith("joint_known"):
                        mask &= np.isfinite(av)&np.isfinite(sc)
                    n=int(mask.sum())
                    if n<N:
                        continue
                    g=[genus_key(s) for s in p.accepted_species.to_numpy()[mask]]
                    metrics=calc_metrics(near[:,mask],flag[mask],g)
                    private.append({"island_id":island,"mode":mode,"log_joint_n":log(n),
                        "joint_coverage":n/totals[island],"joint_n":n,**metrics})
            pf=pd.DataFrame(private,columns=["island_id","mode","log_joint_n","joint_coverage","joint_n","share_n10","random_share_n10","excess_random_n10","share_full","genus_expected_share_full","excess_genus_full","positive_fraction","volume_all_n10"])
            pf=pf.merge(geo,on="island_id",validate="many_to_one")
            bf=pd.DataFrame(bound_records)
            for region in REGIONS:
                br=bf.loc[bf.region.eq(region)]
                bounds.append({"endpoint":endpoint,"radius":radius,"region":region,"n_morphology_islands":len(br),
                    "lower_quartiles":interval_summary(br.lower),"upper_quartiles":interval_summary(br.upper),
                    "unknown_fraction_quartiles":interval_summary(br.unknown_frac)})
                for mode in MODES:
                    p=pf.loc[pf.analysis_regime.eq(region)&pf["mode"].eq(mode)]
                    known=complete[A].notna() if mode=="autonomy_known" else complete[SC].notna()
                    if mode.startswith("joint_known"):
                        known=complete[A].notna()&complete[SC].notna()
                    support.append({"endpoint":endpoint,"radius":radius,"region":region,"mode":mode,
                        "global_morph_species":len(complete),"global_known_joint_species":int(known.sum()),
                        "reference_profiles":len(ref),"n_islands":len(p),"n_blocks":int(p.spatial_block.nunique()),
                        "joint_n_quartiles":interval_summary(p.joint_n),
                        "positive_fraction_quartiles":interval_summary(p.positive_fraction),
                        "share_n10_quartiles":interval_summary(p.share_n10),
                        "excess_random_n10_quartiles":interval_summary(p.excess_random_n10),
                        "excess_genus_full_quartiles":interval_summary(p.excess_genus_full),
                        "saturated_n10_share":float((p.volume_all_n10>=.98).mean()) if len(p) else None})
                    for outcome in OUTCOMES:
                        rows.append({"endpoint":endpoint,"radius":radius,"region":region,"mode":mode,"outcome":outcome,**fit_cluster(p,outcome,PRED)})
    fixed_adjust(rows)
    return {"status":"OWN_GBIF_CAPACITY_MORPHOLOGY_NOT_REALIZED_FUNCTION","source":"GBIF original roster + own Direct ledger; no GIFT or floristic status",
            "n_roster_islands":roster.island_id.nunique(),"n_roster_rows":len(roster),
            "n_geo":len(geo),"minimum_joint_species":N,"support":support,"bounds":bounds,"models":rows,
            "no_row_level_export":True,"unknown_not_absent":True,"not_causal":True}


def main():
    p=argparse.ArgumentParser()
    for key in ("roster","ledger","geography","output"):
        p.add_argument("--"+key,type=Path,required=True)
    a=p.parse_args()
    if digest(a.ledger)!=LOCK["ledger"] or digest(a.geography)!=LOCK["geography"]:
        raise ValueError("frozen source hash mismatch")
    roster=normalize_roster(pd.read_csv(a.roster,dtype=str))
    states=load_states(pd.read_csv(a.ledger,dtype=str))
    geo=pd.read_csv(a.geography,dtype={"island_id":str,"spatial_block":str})
    d=run(roster,states,geo)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(d,indent=2,allow_nan=False)+"\n")
    for r in d["support"]:
        print("CAP_SUPPORT "+json.dumps(r))
    for r in d["bounds"]:
        print("CAP_BOUNDS "+json.dumps(r))
    for r in d["models"]:
        print("CAP_MODEL "+json.dumps(r))
    print("CAP_COMPLETE "+json.dumps({k:d[k] for k in ("n_roster_islands","n_roster_rows","n_geo","status")}))


if __name__=="__main__":
    main()
