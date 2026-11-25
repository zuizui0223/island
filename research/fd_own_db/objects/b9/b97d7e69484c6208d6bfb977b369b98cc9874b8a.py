"""Exact Direct floral reference-mass coverage: original island-native WCVP sensitivity.

Published regional TDWG-L3 WCVP native distribution is NOT proof of island
native establishment. Original GIFT introduced names are never upgraded.
Source-rights restricted accepted plant/island/source IDs stay in CI only.

Expected one-trait-radius reference coverage under a RANDOM n=30 finite
population species sample is computed ANALYTICALLY, not by finite MC:
  1 - choose(N-m_r,30)/choose(N,30)
where m_r counts plants whose original Direct phenotypes are within one
recorded categorical mismatch of original global profile r.

The k-profile-count reference excess is a SEPARATE MC approximation
using multivariate hypergeometric n30 draws; it is NOT analytically exact.
"""
from __future__ import annotations
from collections import Counter
import hashlib
import numpy as np
import pandas as pd
from scipy.special import gammaln

TRAITS=("flower_primary_color","floral_form","floral_symmetry")


def frozen_seed(*labels:object)->int:
    return int.from_bytes(hashlib.sha256(
        ":".join(map(str,labels)).encode()).digest()[:8],"big")


def regional_native_candidates(original_flora:pd.DataFrame,
                               island_tdwg:pd.DataFrame,
                               source_wcvp_ranges:pd.DataFrame)->pd.DataFrame:
    for frame,cols in (
        (original_flora,{"island_id","accepted_species","origin_status"}),
        (island_tdwg,{"island_id","tdwg_l3_code","tdwg_match_status"}),
        (source_wcvp_ranges,{"accepted_species","native_l3_codes"})):
        if not cols.issubset(frame):
            raise ValueError("Source-pinned original species/status geography missing")
    if original_flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original focal island accepted species duplicated")
    if island_tdwg.island_id.duplicated().any():
        raise ValueError("Original island TDWG mapping duplicated")
    ranges=source_wcvp_ranges[["accepted_species","native_l3_codes"]].drop_duplicates()
    if ranges.accepted_species.duplicated().any():
        raise ValueError("Original published WCVP candidate native codes contradict")
    if not original_flora.origin_status.isin((
        "native","introduced","unresolved","endemic")).all():
        raise ValueError("Original island GIFT source nativeness states changed")
    x=original_flora.merge(
        island_tdwg,on="island_id",how="left",validate="many_to_one").merge(
        ranges,on="accepted_species",how="left",validate="many_to_one")
    confirmed=np.fromiter((
        isinstance(code,str) and bool(code) and stat=="accepted" and
        isinstance(regions,str) and code in regions.split("|")
        for code,stat,regions in zip(
            x.tdwg_l3_code,x.tdwg_match_status,x.native_l3_codes)
    ),dtype=bool,count=len(x))
    x["strict_source_native"]=x.origin_status.eq("native")
    x["regional_WCVP_upgrade_only"]=x.origin_status.eq("unresolved")&confirmed
    x["regional_native_compatible"]=(
        x.strict_source_native|x.regional_WCVP_upgrade_only)
    if (x.regional_native_compatible&x.origin_status.eq("introduced")).any():
        raise AssertionError("GIFT source-introduced taxa silently treated as native")
    return x


class ExactFiniteProfileSupport:
    def __init__(self,source_direct:pd.DataFrame, *,
                 reference_draws:int=1999,reference_seed:int=20261113):
        if list(source_direct.columns)!=list(TRAITS) or source_direct.isna().any().any():
            raise ValueError("Original Direct three-trait complete 100-profile schema changed")
        tuples=list(map(tuple,source_direct.astype(str).to_numpy()))
        c=Counter(tuples)
        self.profiles=tuple(sorted(c))
        if len(tuples)!=2234 or len(self.profiles)!=100:
            raise ValueError("Original 2234 Direct plant/100 triple reference drift")
        self.counts=np.asarray([c[p] for p in self.profiles],dtype=int)
        self.species_profile=dict(zip(source_direct.index,tuples))
        self.profile_idx={p:i for i,p in enumerate(self.profiles)}
        p=np.asarray(self.profiles,dtype=object)
        self.within_one=((p[:,None,:]!=p[None,:,:]).sum(axis=2)<=1)
        self.rare=np.flatnonzero(self.counts<=2)
        self.common=np.flatnonzero(self.counts>=26)
        if len(self.rare)!=29 or len(self.common)!=25:
            raise ValueError("Original source global rarity reference changed")
        if reference_draws<1999:
            raise ValueError("Cannot relax category reference null support")
        masks=[]
        for j in range(100):
            bits=0
            for i in np.flatnonzero(self.within_one[:,j]):
                bits|=1<<int(i)
            masks.append(bits)
        rm=sum(1<<int(i) for i in self.rare)
        cm=sum(1<<int(i) for i in self.common)
        rng=np.random.default_rng(reference_seed)
        self.null={}
        for k in range(1,31):
            values=[]
            for _ in range(reference_draws):
                ids=rng.choice(100,size=k,replace=False)
                covered=0
                for j in ids:covered|=masks[int(j)]
                values.append((covered&rm).bit_count()/29-
                              (covered&cm).bit_count()/25)
            self.null[k]=float(np.mean(values))
        self.reference_seed=reference_seed
        self.reference_draws=reference_draws

    def source_scores(self,names:list[str],*,seed:int,
                      k_distribution_draws:int=2000)->dict[str,float]:
        if len(names)<30 or k_distribution_draws<2000:
            raise ValueError("Original equal-n30 Direct source opportunity missing")
        frequencies=np.zeros(100,dtype=int)
        for n in names:
            frequencies[self.profile_idx[self.species_profile[n]]]+=1
        N=int(frequencies.sum())
        covered=self.within_one@frequencies
        p_none=np.zeros(100,dtype=float)
        possible=N-covered>=30
        remaining=N-covered[possible]
        p_none[possible]=np.exp(
            gammaln(remaining+1)-gammaln(remaining-30+1)+
            gammaln(N-30+1)-gammaln(N+1))
        mass=1-p_none
        rare=float(np.mean(mass[self.rare]))
        common=float(np.mean(mass[self.common]))
        rng=np.random.default_rng(seed)
        sampled=rng.multivariate_hypergeometric(
            frequencies[frequencies>0],30,size=k_distribution_draws)
        unique=np.count_nonzero(sampled>0,axis=1)
        mean_null=float(np.mean([self.null[int(k)] for k in unique]))
        gap=rare-common
        return dict(
            rare_coverage=rare,
            common_coverage=common,
            rare_minus_common_gap=gap,
            k_uniform_residual_gap=gap-mean_null,
            estimated_mean_distinct_k=float(np.mean(unique)),
            original_source_pool_species=N,
        )


def cluster_equal_bootstrap(rows:pd.DataFrame,column:str,cluster:str,*,
                            draws:int=1999,seed:int=20261114)->dict|None:
    values=rows.groupby(cluster)[column].mean().dropna().to_numpy(float)
    if len(values)<5 or draws<999:return None
    rng=np.random.default_rng(seed)
    sampled=values[rng.integers(len(values),size=(draws,len(values)))].mean(axis=1)
    return dict(mean=float(np.mean(values)),
                ci95=[float(np.quantile(sampled,.025)),
                      float(np.quantile(sampled,.975))],
                n_independent_clusters=int(len(values)))
