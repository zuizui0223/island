"""Rarity-resolved finite three-trait *category* support coverage.

Different from profile *number*: reference consists of 100 distinct
actually recorded global Direct colour×form×symmetry triples, and each
island/sample's one-trait Hamming neighborhood covers different subsets.

Rare/common bins are defined **before island results**, by number of global
Direct species carrying each profile, using quartile thresholds with ties.
Outcomes are fractions WITHIN a frozen reference-bin, so bin size and
category frequency are not silently confounded.

Null 1: exactly k distinct reference profile centers selected uniformly.
Null 2: k distinct reference profiles selected without replacement with
probability proportional to ORIGINAL GLOBAL Direct species frequencies.
Neither null preserves genus/species phylogeny or original observation source.
Pair controls: exactly n=30 sampled species on each island and five-nearest
candidate mainland native flora; exclude exact island species from mainland.
The additional strict identical-genus source support gate is set by
taxonomic names only, BEFORE trait state/volume comparisons.

Every output is regional, spatial-block/source-signature aggregated only.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import hashlib
import numpy as np
import pandas as pd

from .discrete_functional_volume import (
    TRAITS, FiniteProfileReference, _resample_group, _seed
)

REGIONS=("northern_high_latitude","northern_midlatitude",
         "southern_extratropical","tropical")
BIN_NAMES=("rare","middle","common")
OUTCOMES=("covered_rare","covered_middle","covered_common",
          "rare_minus_common","k_residual_rare","k_residual_common",
          "k_residual_rare_minus_common",
          "speciesfreq_k_residual_rare_minus_common")


def _seedstr(*label:object)->int:
    return int.from_bytes(hashlib.sha256(
       ":".join(map(str,label)).encode()).digest()[:8],"big")


def genus(name:str)->str:
    words=str(name).strip().split()
    if len(words)<2 or not words[0][0:1].isalpha():
        return ""
    return words[0]


@dataclass
class RarityReference:
    volume: FiniteProfileReference
    bins: dict[str,int]
    n_bins: dict[str,int]
    low_quartile_species_count:int
    high_quartile_species_count:int
    n_global_species: int
    profile_n_species:np.ndarray
    expected_uniform:dict[int,np.ndarray]
    expected_speciesfreq:dict[int,np.ndarray]
    boundary_tie_count:dict[str,int]

    @classmethod
    def from_original(cls,triple:pd.DataFrame,*,
                      null_draws:int=399,seed:int=20261113):
        if null_draws<199:
            raise ValueError("Insufficient predeclared k-specific reference null draws")
        ref=FiniteProfileReference.from_direct_species(
           triple,null_draws=199,seed=seed)
        species=Counter(map(tuple,triple[list(TRAITS)].astype(str).to_numpy()))
        counts=np.asarray([species[x] for x in ref.profiles],dtype=int)
        if len(counts)!=100 or len(triple)!=2234:
            raise ValueError("Original full global 100-profile/2234-Direct source drift")
        # Tie-inclusive quartile split, not arbitrary genome-name rank.
        low=int(np.quantile(counts,.25,method="lower"))
        high=int(np.quantile(counts,.75,method="higher"))
        if not (low<high):
            raise ValueError("Rare and common original Direct profile bins overlap")
        masks={}
        labels=np.where(counts<=low,0,np.where(counts>=high,2,1))
        for label,j in zip(BIN_NAMES,range(3)):
            mask=0
            for k in np.flatnonzero(labels==j):
                mask|=1<<int(k)
            masks[label]=mask
        sizes={k:v.bit_count() for k,v in masks.items()}
        if any(size<10 for size in sizes.values()):
            raise ValueError("Insufficient ties-inclusive observed reference rarity group")
        if sum(sizes.values())!=100:
            raise ValueError("Rarity bins do not partition the actual original profiles")
        probabilities=counts.astype(float)/counts.sum()
        rng=np.random.default_rng(seed+51)
        null_u={}
        null_w={}
        for k in range(1,31):
            uniform=[]
            weighted=[]
            for _ in range(null_draws):
                for destination,ids in (
                   (uniform,rng.choice(len(counts),size=k,replace=False)),
                   (weighted,rng.choice(len(counts),size=k,replace=False,
                          p=probabilities))
                ):
                    covered=0
                    for index in ids:
                        covered|=ref.masks[ref.profiles[int(index)]][1]
                    destination.append([
                      (covered & masks[band]).bit_count()/sizes[band]
                      for band in BIN_NAMES
                    ])
            null_u[k]=np.mean(np.asarray(uniform,float),axis=0)
            null_w[k]=np.mean(np.asarray(weighted,float),axis=0)
        return cls(ref,masks,sizes,low,high,len(triple),counts,
                   null_u,null_w,
                   {"rare_low_threshold_ties":int(np.sum(counts==low)),
                    "common_high_threshold_ties":int(np.sum(counts==high))})

    def measure(self,triples:np.ndarray):
        selected=set(map(tuple,np.asarray(triples,dtype=str)))
        covered=0
        for p in selected:
            if p not in self.volume.masks:
                raise ValueError("Original observed flower state not in Direct global reference")
            covered |= self.volume.masks[p][1]
        fractions=np.asarray([
            (covered & self.bins[band]).bit_count()/self.n_bins[band]
            for band in BIN_NAMES
        ],dtype=float)
        return fractions,len(selected)

    def sample(self,triples:np.ndarray, *,
               repeats:int=40, n:int=30,seed:int=20261113)->dict[str,float]:
        records=np.asarray(triples,dtype=str)
        if n!=30 or repeats<20 or len(records)<n:
            raise ValueError("Insufficient fixed 30 Direct plant carriers")
        rng=np.random.default_rng(seed)
        results=[]
        for _ in range(repeats):
            group=records[rng.choice(len(records),size=n,replace=False)]
            seen,k=self.measure(group)
            null=self.expected_uniform[k]
            wnull=self.expected_speciesfreq[k]
            results.append([
                *seen,
                seen[0]-seen[2],
                seen[0]-null[0],
                seen[2]-null[2],
                (seen[0]-seen[2])-(null[0]-null[2]),
                (seen[0]-seen[2])-(wnull[0]-wnull[2]),
            ])
        return dict(zip(OUTCOMES,np.asarray(results).mean(axis=0).astype(float)))


def _taxonomy_eligibility(island_names:set[str],source_names:set[str])->int:
    isl=Counter(map(genus,island_names))
    pool=Counter(map(genus,source_names))
    return int(sum(min(n,pool[k]) for k,n in isl.items() if k and pool[k]>0))


def original_rare_category_coverage(
    triple:pd.DataFrame,original_flora:pd.DataFrame,
    geo:pd.DataFrame,assignments:pd.DataFrame,mainland:pd.DataFrame,*,
    n:int=30,repeats:int=40,null_draws:int=399,
    bootstraps:int=1999,seed:int=20261113,
):
    if n!=30 or repeats<20 or null_draws<199 or bootstraps<999:
        raise ValueError("Prespecified rare phenotype original-source protocol cannot be relaxed")
    for label,frame,cols in (
      ("flora",original_flora,("island_id","accepted_species","origin_status")),
      ("geography",geo,("island_id","analysis_regime","spatial_block")),
      ("source",assignments,("island_id","source_mode","source_rank","entity_ID")),
      ("mainland",mainland,("entity_ID","work_species"))):
        if not set(cols).issubset(frame):
            raise ValueError("Missing original "+label+" source contract")
    ref=RarityReference.from_original(triple,null_draws=null_draws,seed=seed)
    if original_flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original source island×plant duplicate")
    flora=original_flora.loc[original_flora.accepted_species.isin(triple.index)].copy()
    if not flora.origin_status.isin(("native","endemic","introduced","unresolved")).all():
        raise ValueError("Original source nativeness labels changed")
    island_rich=flora.groupby("island_id").accepted_species.nunique()
    eligible=set(island_rich.loc[island_rich>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original 928-island Direct cohort changed")
    flora=flora.loc[flora.island_id.isin(eligible)]
    geo=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates()
    if geo.island_id.duplicated().any() or not eligible.issubset(set(geo.island_id)):
        raise ValueError("Original geographic island block not unique")
    loc=geo.set_index("island_id")
    src=assignments.loc[
         assignments.island_id.isin(eligible)&
         assignments.source_mode.eq("geo_k5")&
         assignments.source_rank.between(1,5),
         ["island_id","source_rank","entity_ID"]].copy()
    if src.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original rank within five-nearest source")
    nearest=src.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(nearest.get(i,()))!=5 or len(set(nearest[i]))!=5 for i in eligible):
        raise ValueError("Original five GIFT donor slots incomplete")
    g=mainland.loc[mainland.work_species.isin(triple.index),
       ["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    byisland={key:frame for key,frame in flora.groupby("island_id")}
    rows=[]
    for island in sorted(eligible):
        own=byisland[island]
        all_names=set(own.accepted_species)
        nearby=set().union(*(g.get(id_,set()) for id_ in nearest[island]))
        nearby &= set(triple.index)
        disjoint=nearby-all_names
        source_scores={}
        for label,names in (("all",nearby),("disjoint",disjoint)):
            if len(names)>=n:
                source_scores[label]=ref.sample(
                    triple.loc[sorted(names),list(TRAITS)].to_numpy(dtype=str),
                    repeats=repeats,n=n,seed=_seedstr(seed,island,label))
        for scope in ("all_recorded","source_native"):
            selected=(own if scope=="all_recorded"
                      else own.loc[own.origin_status.eq("native")])
            names=set(selected.accepted_species)
            if len(names)<n:
                continue
            scores=ref.sample(triple.loc[sorted(names),list(TRAITS)].to_numpy(dtype=str),
                repeats=repeats,n=n,seed=_seedstr(seed,island,scope))
            entry=dict(island_id=island,scope=scope,
                regime=str(loc.at[island,"analysis_regime"]),
                spatial_block=str(loc.at[island,"spatial_block"]),
                source_signature=nearest[island][0],
                n_source_direct=len(nearby),
                n_source_disjoint=len(disjoint),
                n_local_direct=len(names),
                n_genus_donors_available_for_30_sample_source_disjoint=(
                    _taxonomy_eligibility(names,disjoint)),
                source_all_eligible=("all" in source_scores),
                source_disjoint_eligible=("disjoint" in source_scores))
            for key in OUTCOMES:
                entry[key]=scores[key]
                for which,s in source_scores.items():
                    entry[which+"_delta_"+key]=scores[key]-s[key]
            rows.append(entry)
    priv=pd.DataFrame(rows)
    region=[]
    for (scope,regime),q in priv.groupby(["scope","regime"]):
        if len(q)<10 or q.spatial_block.nunique()<5:
            continue
        entry=dict(scope=scope,regime=regime,n_islands=len(q),
            n_spatial_blocks=q.spatial_block.nunique(),
            n_first_nearest_GIFT_source_signatures=q.source_signature.nunique(),
            n_source_all_eligible=int(q.source_all_eligible.sum()),
            n_source_disjoint_eligible=int(q.source_disjoint_eligible.sum()),
            n_source_disjoint_genus_matching_max_30plus=int(
                q.n_genus_donors_available_for_30_sample_source_disjoint.ge(30).sum()),
            n_source_disjoint_genus_matching_max_15plus=int(
                q.n_genus_donors_available_for_30_sample_source_disjoint.ge(15).sum()),
            median_number_source_disjoint_same_genus_donors=float(
                q.n_genus_donors_available_for_30_sample_source_disjoint.median()),
            median_island_rare_coverage=float(q.covered_rare.median()),
            median_island_common_coverage=float(q.covered_common.median()),
            median_island_rare_minus_common=float(q.rare_minus_common.median()))
        for k,srcname in enumerate(("all","disjoint")):
            keys=("covered_rare","covered_common","rare_minus_common",
                  "k_residual_rare_minus_common",
                  "speciesfreq_k_residual_rare_minus_common")
            for j,key in enumerate(keys):
                col=srcname+"_delta_"+key
                if col not in q:
                    continue
                sub=q.dropna(subset=[col])
                if len(sub)<10 or sub.spatial_block.nunique()<5:
                    continue
                entry["n_"+srcname+"_supported_for_"+key]=len(sub)
                entry["median_"+col]=float(sub[col].median())
                for level,weight in (("block","spatial_block"),("source","source_signature")):
                    mean,lo,hi=_resample_group(sub,col,weight,
                        seed+100*k+20*j+(level=="source"),draws=bootstraps)
                    entry[level+"_mean_"+col]=mean
                    entry[level+"_ci025_"+col]=lo
                    entry[level+"_ci975_"+col]=hi
        region.append(entry)
    # Biological phylogenetic conditioning requires actual genetic history.
    # Genus donor opportunity is an independent pre-outcome *identifiability*
    # gate, not successful genus-adjusted effect estimation.
    return dict(
      status="ORIGINAL_RARITY_RESOLVED_CATEGORICAL_VOLUME_DESCRIPTIVE",
      n_original_reference_complete_Direct_species=len(triple),
      n_original_reference_triple_profiles=len(ref.volume.profiles),
      n_original_928_island_flora=len(eligible),
      n_original_island_scope_rows_scored=len(priv),
      n_strict_native_source_flora_scored=int(priv.scope.eq("source_native").sum()),
      rarity_bins=dict(
          n_original_profile_categories_by_tier=ref.n_bins,
          rare_profile_global_direct_species_frequency_max=ref.low_quartile_species_count,
          common_profile_global_direct_species_frequency_min=ref.high_quartile_species_count,
          exact_frequency_boundary_ties=ref.boundary_tie_count,
          reference_tier_labels_predeclared=True),
      rare_vs_common_island_source_aggregate=region,
      sampling=dict(equal_direct_species_island_and_GIFT_source=30,
          original_direct_random_rarefactions=repeats,
          equal_k_uniform_profile_null_draws=null_draws,
          equal_k_global_speciesfrequency_profile_null_draws=null_draws,
          block_and_mainland_source_bootstrap_draws=bootstraps),
      genus_matched_null_effect_estimated=False,
      measured_actual_pollinator_function_or_adaptation=False,
      source_is_verified_ancestral_pool=False,
      main_sensitivity_removes_all_island_species_from_nearby_GIFT_source=True,
      nonrecorded_flower_profiles_not_assumed_biologically_absent=True,
      source_species_island_taxa_genera_coordinates_not_released=True,
    )
