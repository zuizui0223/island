"""Finite-reference categorical floral FUNCTIONAL coverage volume.

A discrete category trait space has no intrinsic Euclidean Lebesgue volume.
Here volume is the REFERENCE MEASURE of a union of fixed Hamming balls:
 V_r(A) = [sum_{x in ORIGINAL_GLOBAL_REFERENCE} 1[min_{a in A}d(x,a)<=r]] / K
for r=0,1/3,2/3 and d=(1/3)sum of three categorical mismatches.

r=0 is observed unique profile richness / K, NOT a new hypervolume finding.
r=1/3 is an observable *one-trait-neighbourhood* volume: unlike r=0 it
depends on how profiles are placed relative to one another in shared
three-trait geometry. Crucially both the n=30 original-species carrier count
AND the number k of distinct recorded profiles are controlled.

All reference profiles and weights come from original Direct single-state
species (no hypothetical colour×form×symmetry profiles are called observed).
Permutation combinations that were not globally observed may be evaluated as
counterfactual null CENTERS, not added to the reference measure.

This tests source-relative flora assembly, not evolutionary trait plasticity,
colonization history, pollination services or actual island visitor turnover.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import hashlib

import numpy as np
import pandas as pd

TRAITS=("flower_primary_color","floral_form","floral_symmetry")
GEO=("island_id","analysis_regime","spatial_block")
SOURCE=("island_id","source_mode","source_rank","entity_ID")
MEASURES=("occupancy_r0","metric_volume_r1","metric_volume_r2",
          "metric_volume_r1_weighted","metric_volume_r1_richness_residual",
          "colour_architecture_pairing_excess_r1")
ALLOWED_ORIGINS={"native","endemic","introduced","unresolved"}


def _seed(*labels:object)->int:
    h=hashlib.sha256(":".join(map(str,labels)).encode()).digest()
    return int.from_bytes(h[:8],"big")


@dataclass
class FiniteProfileReference:
    """Original Direct population reference, with three binary mismatch balls."""
    profiles:tuple[tuple[str,str,str],...]
    state_probabilities:np.ndarray
    masks:dict[tuple[str,str,str],tuple[int,int,int]]
    expected_volume_r1_given_k:dict[int,float]

    @classmethod
    def from_direct_species(cls,triples:pd.DataFrame, *,
                            null_draws:int=399,seed:int=20261113):
        if list(triples.columns)!=list(TRAITS):
            raise ValueError("Direct original 3-trait reference column contract changed")
        if len(triples)<100 or triples.isna().any().any():
            raise ValueError("Original Direct three-trait reference incomplete")
        if any(triples[c].astype(str).str.contains("|",regex=False).any() for c in TRAITS):
            raise ValueError("Ambiguous Direct phenotype state cannot be expanded")
        tuples=list(map(tuple,triples.astype(str).to_numpy()))
        counts=pd.Series(tuples).value_counts()
        profiles=tuple(sorted(counts.index))
        n=len(profiles)
        if n<8 or n>500:
            raise ValueError("Expected finite original global 3-trait combination support")
        weights=np.asarray([float(counts[p])/len(tuples) for p in profiles])
        # All combinations are allowed solely as hypothetical permutation
        # centers. The reference mass is NEVER expanded into those categories.
        states=[sorted(triples[t].astype(str).unique()) for t in TRAITS]
        candidates=tuple(product(*states))
        if len(candidates)>10000:
            raise ValueError("Global original categorical product has implausible cardinality")
        masks={}
        for candidate in candidates:
            mismatch=[
                int(candidate[0]!=p[0])+int(candidate[1]!=p[1])+int(candidate[2]!=p[2])
                for p in profiles
            ]
            masks[candidate]=tuple(
                sum(1<<j for j,x in enumerate(mismatch) if x<=radius)
                for radius in range(3)
            )
        expected={}
        rng=np.random.default_rng(seed)
        # Independent UNIQUE-profile geometric null, conditional on the
        # observed number k of occupied profiles after n=30 rarefaction.
        # This null does not hold genus or trait marginal prevalence fixed.
        for k in range(1,min(30,n)+1):
            vals=[]
            for _ in range(null_draws):
                selected=rng.choice(n,size=k,replace=False)
                mask=0
                for j in selected:
                    mask |= masks[profiles[int(j)]][1]
                vals.append(mask.bit_count()/n)
            expected[k]=float(np.mean(vals))
        return cls(profiles,weights,masks,expected)

    def coverage(self,arr:np.ndarray)->tuple[float,float,float,float,int]:
        arr=np.asarray(arr,dtype=str)
        if arr.ndim!=2 or arr.shape[1]!=3 or len(arr)==0:
            raise ValueError("One recorded phenotype triple per source species required")
        selected=set(map(tuple,arr))
        bits=[0,0,0]
        for profile in selected:
            if profile not in self.masks:
                raise ValueError("Source category absent from original Direct universe")
            cover=self.masks[profile]
            for j in range(3):
                bits[j] |= cover[j]
        n=len(self.profiles)
        r0,r1,r2=(v.bit_count()/n for v in bits)
        indices=[j for j in range(n) if (bits[1]>>j)&1]
        weighted=float(self.state_probabilities[indices].sum())
        return float(r0),float(r1),float(r2),weighted,len(selected)

    def sample_scores(self,species_profiles:np.ndarray,*,n:int=30,
                      repeats:int=60,null_shuffles:int=29,
                      rng:np.random.Generator)->dict[str,float]:
        a=np.asarray(species_profiles,dtype=str)
        if n!=30 or repeats<10 or null_shuffles<9 or len(a)<n:
            raise ValueError("Insufficient fixed-30 Direct richness / null support")
        vals=[]
        for _ in range(repeats):
            sample=a[rng.choice(len(a),size=n,replace=False)].copy()
            r0,r1,r2,r1w,k=self.coverage(sample)
            expected=self.expected_volume_r1_given_k.get(k)
            if expected is None:
                raise ValueError("Profile richness k outside frozen null grid")
            shuffle=[]
            for _ in range(null_shuffles):
                changed=sample.copy()
                changed[:,0]=rng.permutation(sample[:,0])
                shuffle.append(self.coverage(changed)[1])
            vals.append((r0,r1,r2,r1w,r1-expected,r1-float(np.mean(shuffle))))
        mean=np.mean(np.asarray(vals,float),axis=0)
        return dict(zip(MEASURES,map(float,mean)))


def _resample_group(group:pd.DataFrame,field:str,
                    grouping:str,seed:int,draws:int=1999)->tuple[float,float,float]:
    levels=group.groupby(grouping)[field].mean().dropna().to_numpy(float)
    if len(levels)<5:
        return (float("nan"),float("nan"),float("nan"))
    rand=np.random.default_rng(seed)
    boot=levels[rand.integers(len(levels),size=(draws,len(levels)))].mean(axis=1)
    return tuple(map(float,(levels.mean(),np.quantile(boot,.025),np.quantile(boot,.975))))


def original_categorical_flora_volume(
    triple:pd.DataFrame,original_flora:pd.DataFrame,geo:pd.DataFrame,
    source_assignment:pd.DataFrame,source_mainland:pd.DataFrame, *,
    n:int=30,repeats:int=60,null_shuffles:int=29,random_profiles:int=399,
    bootstrap:int=1999,seed:int=20261113,
)->dict:
    """Run protected real source within CI; return ONLY safe region aggregates."""
    if n!=30 or repeats<20 or null_shuffles<19 or random_profiles<199 or bootstrap<999:
        raise ValueError("Cannot relax original registered categorical volume audit")
    for label,table,cols in (
        ("flora",original_flora,{"island_id","accepted_species","origin_status"}),
        ("geography",geo,set(GEO)),
        ("sources",source_assignment,set(SOURCE)),
        ("mainland",source_mainland,{"entity_ID","work_species"})):
        if not cols.issubset(table):
            raise ValueError("Original source "+label+" missing mandatory columns")
    reference=FiniteProfileReference.from_direct_species(
        triple,null_draws=random_profiles,seed=seed)
    if original_flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate original source island × accepted species")
    if not original_flora.origin_status.isin(ALLOWED_ORIGINS).all():
        raise ValueError("Unknown original flora origin state")
    flora=original_flora.loc[original_flora.accepted_species.isin(triple.index)].copy()
    counts=flora.groupby("island_id").accepted_species.nunique()
    eligible=set(counts.loc[counts.ge(30)].index)
    if len(eligible)!=928:
        raise ValueError("Original source 928 Direct morphological island cohort drifted")
    flora=flora.loc[flora.island_id.isin(eligible)].copy()
    geo=geo[list(GEO)].drop_duplicates()
    if geo[list(GEO)].isna().any().any() or geo.island_id.duplicated().any():
        raise ValueError("Missing original geographic block or duplicated island ID")
    if not eligible.issubset(set(geo.island_id)):
        raise ValueError("One or more original 928 island IDs not geographically qualified")
    geography=geo.set_index("island_id")
    source=source_assignment.loc[
        source_assignment.island_id.isin(eligible)&
        source_assignment.source_mode.eq("geo_k5")&
        source_assignment.source_rank.between(1,5),
        list(SOURCE)
    ].copy()
    if source.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Duplicate original five-nearest GIFT source rank")
    sources=source.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(sources.get(i,()))!=5 or len(set(sources[i]))!=5 for i in eligible):
        raise ValueError("Original five native mainland region source slots incomplete")
    mainland=source_mainland.loc[
       source_mainland.work_species.isin(triple.index),
       ["entity_ID","work_species"]].dropna().drop_duplicates()
    gift=mainland.groupby("entity_ID").work_species.agg(set).to_dict()
    byisland={i:x for i,x in flora.groupby("island_id")}
    rows=[]
    for island in sorted(eligible):
        f=byisland[island]
        all_names=set(f.accepted_species)
        donor=set().union(*(gift.get(j,set()) for j in sources[island]))
        donor=donor.intersection(triple.index)
        disjoint=donor-all_names
        src={}
        for src_scope,names in (("source_all",donor),
                                ("source_excluding_same_island_species",disjoint)):
            if len(names)>=n:
                src[src_scope]=reference.sample_scores(
                    triple.loc[sorted(names),list(TRAITS)].to_numpy(dtype=str),
                    n=n,repeats=repeats,null_shuffles=null_shuffles,
                    rng=np.random.default_rng(_seed(seed,island,src_scope)))
        for scope in ("all_recorded","source_native"):
            own=f if scope=="all_recorded" else f.loc[f.origin_status.eq("native")]
            if len(own)<n:
                continue
            values=reference.sample_scores(
                triple.loc[sorted(set(own.accepted_species)),list(TRAITS)].to_numpy(dtype=str),
                n=n,repeats=repeats,null_shuffles=null_shuffles,
                rng=np.random.default_rng(_seed(seed,island,scope))
            )
            record=dict(
                island_id=island,scope=scope,
                analysis_regime=str(geography.at[island,"analysis_regime"]),
                spatial_block=str(geography.at[island,"spatial_block"]),
                source_signature=sources[island][0],
                n_eligible_direct_species=len(own),
                n_global_reference_profiles=len(reference.profiles),
                n_candidate_source_species=len(donor),
                n_exclusive_candidate_source_species=len(disjoint),
                source_all_available=("source_all" in src),
                source_disjoint_available=("source_excluding_same_island_species" in src))
            for measure in MEASURES:
                record[measure]=values[measure]
                for donor_name,score in src.items():
                    record[donor_name+"_delta_"+measure]=values[measure]-score[measure]
            rows.append(record)
    private=pd.DataFrame.from_records(rows)
    output=[]
    for (scope,regime),p in private.groupby(["scope","analysis_regime"]):
        if len(p)<10 or p.spatial_block.nunique()<5:
            continue
        rec=dict(scope=scope,regime=regime,n_islands=int(len(p)),
           n_spatial_blocks=int(p.spatial_block.nunique()),
           n_nearest_mainland_source_signatures=int(p.source_signature.nunique()),
           n_source_all_supported=int(p.source_all_available.sum()),
           n_source_disjoint_supported=int(p.source_disjoint_available.sum()))
        for key in MEASURES:
            rec["median_"+key]=float(p[key].median())
        for k,donor in enumerate(("source_all","source_excluding_same_island_species")):
            for j,metric in enumerate(("occupancy_r0","metric_volume_r1",
                   "metric_volume_r1_weighted",
                   "metric_volume_r1_richness_residual",
                   "colour_architecture_pairing_excess_r1")):
                col=donor+"_delta_"+metric
                if col not in p:
                    continue
                ok=p.dropna(subset=[col])
                if len(ok)<10 or ok.spatial_block.nunique()<5:
                    continue
                rec["n_"+donor+"_"+metric+"_eligible"]=int(len(ok))
                rec["median_"+col]=float(ok[col].median())
                for grouping,tier in (("spatial_block","block"),
                                       ("source_signature","source")):
                    mean,lo,hi=_resample_group(ok,col,grouping,
                                    seed+100*k+20*j+(tier=="source"),draws=bootstrap)
                    rec[tier+"_equal_mean_"+col]=mean
                    rec[tier+"_ci025_"+col]=lo
                    rec[tier+"_ci975_"+col]=hi
        rec["fraction_islands_r1_near_saturation"]=float(
            p.metric_volume_r1.ge(.98).mean())
        rec["fraction_islands_r2_near_saturation"]=float(
            p.metric_volume_r2.ge(.98).mean())
        output.append(rec)
    # Compare ALL and strict native on the same island, not two unrelated
    # geographic source samples; unknown origin is never renamed nonnative.
    paired=[]
    raw_all=private.loc[private.scope.eq("all_recorded")].set_index("island_id")
    raw_native=private.loc[private.scope.eq("source_native")].set_index("island_id")
    for regime,g in raw_native.groupby("analysis_regime"):
        same=g.index.intersection(raw_all.index)
        if len(same)<10 or g.spatial_block.nunique()<5:
            continue
        a=raw_all.loc[same]
        b=raw_native.loc[same]
        x=b[["spatial_block"]].copy()
        entry=dict(regime=regime,n_native_eligible_paired_islands=len(same),
                   n_spatial_blocks=g.spatial_block.nunique(),
                   all_recorded_and_native_are_nested_not_independent=True)
        for j,metric in enumerate(("occupancy_r0","metric_volume_r1",
                     "metric_volume_r1_richness_residual")):
            col="native_minus_all_"+metric
            x[col]=(b[metric]-a[metric]).to_numpy(float)
            mean,lo,hi=_resample_group(x,col,"spatial_block",
                                      seed+500+j,draws=bootstrap)
            entry["block_equal_mean_"+col]=mean
            entry["block_ci025_"+col]=lo
            entry["block_ci975_"+col]=hi
        paired.append(entry)
    return {
      "status":"CATEGORICAL_FINITE_REFERENCE_METRIC_VOLUME_NOT_17D_KDE",
      "reference":"original_global_direct_high_medium_single_state_triple",
      "n_original_reference_direct_plant_species":int(len(triple)),
      "n_original_distinct_triple_profiles":len(reference.profiles),
      "n_928_direct_islands":len(eligible),
      "n_all_recorded_islands_scored":int(private.scope.eq("all_recorded").sum()),
      "n_original_strict_native_islands_scored":int(private.scope.eq("source_native").sum()),
      "floral_metric_volume_by_geographic_regime":output,
      "native_same_island_sensitivity":paired,
      "measurement":{
        "metric":"three_category_slot_Hamming",
        "mismatch_weights_per_trait":[1/3,1/3,1/3],
        "radius_values_in_mismatch_count":[0,1,2],
        "reference_measure":"uniform_over_original_observed_global_Direct_profiles",
        "secondary_measure":"original_global_Direct_species_prevalence_weighted",
        "radius0_is_only_unique_observed_profile_richness_over_K":True,
        "radius1_interpreted_as_one_trait_neighbourhood_not_continuous_Euclidean_volume":True,
        "radius2_only_saturation_diagnostic":True,
        "fixed_species_per_island_or_source":n,
        "species_rarefaction_repeats":repeats,
        "colour_shuffle_per_rarefaction":null_shuffles,
        "null_preserves_observed_colour_margins_and_form_symmetry_pairs":True,
        "richness_conditional_geometry_null_unique_profiles":random_profiles,
        "richness_conditional_geometry_null_preserves_taxonomic_identity":False,
        "spatial_and_GIFT_source_bootstrap_draws":bootstrap
      },
      "no_implicit_source_species_or_island_taxa_export":True,
      "observed_source_is_not_verified_ancestor":True,
      "demonstrated_causal_pollinator_loss_or_flower_evolution":False,
      "not_a_species_independent_or_phylogenetically_adjusted_null":True,
    }
