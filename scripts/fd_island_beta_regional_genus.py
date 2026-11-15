"""Region-local versus global within-genus floral beta references on EXACT dyads.

Main source universe: original 4505 GIFT island floras; endpoint: original
species-direct high/medium colour + form on the SAME accepted plant species.
The region-local null reuses botanical observations from its target region;
it is an in-sample sensitivity, NOT an external mainland source prediction.
All pair masks exclude shared accepted species in the FULL botanical flora,
regardless of native-specific endpoint status. No independent dyad p-values.
"""
from __future__ import annotations
import argparse
import json
import hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd
from fd_island_functional_beta import (
    strict_traits, prepare_endpoint, compare_block, membership,
    geographical_distance,
)
from fd_island_beta_geographic_null import Aggregate, BINS

GATES = ("n5", "both_n10", "both_recorded_fraction_ge_0p10")
MIN_PAIRS, MIN_ISLANDS, MIN_BLOCKS = 100, 30, 10


def _subset_prepared(ref:dict, rows:np.ndarray) -> dict:
    """Slice target rows but retain the EXACT original global genus donors."""
    out=dict(ref)
    for k in ("A","B","n","counts","P","genus_mean","genus_trace"):
        out[k]=ref[k][rows]
    return out


def _genus_name(name:str)->str:
    terms=name.split()
    return terms[0] if len(terms)>=2 else "single_species:"+name


def _regional_donor_audit(states:dict, used_species:set) -> dict:
    groups=defaultdict(set)
    for sp in used_species:
        groups[_genus_name(sp)].add(states[sp])
    var={k for k,v in groups.items() if len(v)>1}
    return dict(
        n_region_reference_species=len(used_species),
        n_region_reference_genera=len(groups),
        n_region_genera_with_multiple_direct_profiles=len(var),
        fraction_distinct_region_species_in_variable_genus=(
            sum(_genus_name(sp) in var for sp in used_species) /
            len(used_species)),
    )


def analyze(states:dict, flora:pd.DataFrame, geo:pd.DataFrame, *,
            batch:int=96, min_pairs:int=MIN_PAIRS,
            min_islands:int=MIN_ISLANDS, min_blocks:int=MIN_BLOCKS)->dict:
    """Both reference expectations evaluated on identical eligible island dyads."""
    if batch<1: raise ValueError("batch must be positive")
    if not {"island_id","accepted_species","origin_status"}<=set(flora):
        raise ValueError("Original accepted flora fields required")
    if flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate island×accepted species")
    if flora.origin_status.isna().any() or not flora.origin_status.isin(
            ("native","introduced","unresolved")).all():
        raise ValueError("Unknown botanical origin token")
    required={"island_id","analysis_regime","spatial_block",
              "island_latitude","island_longitude"}
    if not required<=set(geo) or geo.island_id.duplicated().any():
        raise ValueError("Original physical geography required")
    if not {"flower_primary_color","floral_form"}<=set(states):
        raise ValueError("Both measured flower traits required")
    both=set(states["flower_primary_color"])&set(states["floral_form"])
    joint={sp:(states["flower_primary_color"][sp],states["floral_form"][sp])
           for sp in both}
    all_f=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native_f=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    all_islands=sorted(all_f)
    geo=geo.set_index("island_id").reindex(all_islands)
    reg_all=geo.analysis_regime.fillna("geography_unresolved").astype(str).to_numpy()
    block_all=geo.spatial_block.to_numpy()
    lat_all=pd.to_numeric(geo.island_latitude,errors="coerce").to_numpy()
    lon_all=pd.to_numeric(geo.island_longitude,errors="coerce").to_numpy()
    out=[]
    coverage={}
    for scope, source in (("all_recorded",all_f),("original_native",native_f)):
        source_rosters=[source.get(i,set()) for i in all_islands]
        flower_rosters=[r&both for r in source_rosters]
        observed_n=np.array([len(x) for x in flower_rosters])
        qualified=np.flatnonzero(observed_n>=5)
        coverage[scope]=dict(n_eligible_islands=int(len(qualified)),
                             n_source_islands=len(all_islands))
        if len(qualified)<2:continue
        targets=[flower_rosters[i] for i in qualified]
        full_rosters=[all_f[all_islands[i]] for i in qualified]
        n_qual=len(qualified)
        # Source-wide reference retains the ORIGINAL global Direct trait donor pool.
        global_p=prepare_endpoint(targets,joint)
        target_region=reg_all[qualified]
        target_block=block_all[qualified]
        target_lat=lat_all[qualified]
        target_lon=lon_all[qualified]
        target_frac=observed_n[qualified]/np.array(
            [len(source_rosters[i]) for i in qualified],float)
        for region in sorted(set(target_region)-{"geography_unresolved"}):
            pick=np.flatnonzero(target_region==region)
            if len(pick)<2:continue
            regional_targets=[targets[i] for i in pick]
            regional_used=set().union(*regional_targets)
            # Critical local null: ONLY recorded source-region target names are
            # considered exchangeable; this is in-sample and not causal.
            regional_p=prepare_endpoint(
                regional_targets,{sp:joint[sp] for sp in regional_used})
            global_subset=_subset_prepared(global_p,pick)
            audit=_regional_donor_audit(joint,regional_used)
            orig_f=[full_rosters[i] for i in pick]
            full_sp=sorted(set().union(*orig_f))
            full_A=membership(orig_f,full_sp)
            blocks=target_block[pick]
            block_id,block_labels=pd.factorize(blocks)
            longitude=target_lon[pick]
            latitude=target_lat[pick]
            carrier_n=observed_n[qualified[pick]]
            carrier_fraction=target_frac[pick]
            accum={}
            for start in range(0,len(pick),batch):
                stop=min(start+batch,len(pick))
                original_ref=compare_block(global_subset,start,stop)
                local_ref=compare_block(regional_p,start,stop)
                # The two reference universes must not change observed D².
                if not np.allclose(
                        original_ref["squared_frequency_distance"],
                        local_ref["squared_frequency_distance"],atol=1e-11):
                    raise AssertionError("Null donor selection changed observed phenotypes")
                shared_full=(full_A[start:stop]@full_A.T).toarray()
                ii,jj=np.nonzero(
                    np.arange(len(pick))[None,:] >
                    np.arange(start,stop)[:,None])
                if not len(ii):continue
                I,J=ii+start,jj
                km=geographical_distance(latitude,longitude,start,stop)[ii,jj]
                valid=(
                    (shared_full[ii,jj]==0) &
                    (block_id[I]>=0) & (block_id[J]>=0) &
                    (block_id[I]!=block_id[J]) &
                    np.isfinite(km) &
                    np.isfinite(latitude[I]) & np.isfinite(latitude[J]) &
                    np.isfinite(longitude[I]) & np.isfinite(longitude[J]) &
                    (np.abs(latitude[I])<=90)&(np.abs(latitude[J])<=90) &
                    (np.abs(longitude[I])<=180)&(np.abs(longitude[J])<=180))
                if not valid.any():continue
                obs=original_ref["squared_frequency_distance"][ii,jj]
                expected_global=original_ref["genus_null_squared_distance"][ii,jj]
                expected_region=local_ref["genus_null_squared_distance"][ii,jj]
                gates={
                    "n5":valid,
                    "both_n10":valid&(carrier_n[I]>=10)&(carrier_n[J]>=10),
                    "both_recorded_fraction_ge_0p10":(
                        valid & (carrier_fraction[I]>=.10) &
                        (carrier_fraction[J]>=.10)),
                }
                for gate,selection in gates.items():
                    bins=BINS if gate=="n5" else (BINS[-1],)
                    for label,lo,hi in bins:
                        keep=selection&(km>=lo)&(km<hi)
                        if not keep.any():continue
                        key=(gate,label)
                        if key not in accum:
                            accum[key]=[Aggregate(len(pick),block_id,len(block_labels))
                                        for _ in range(3)]
                        g,r,d=accum[key]
                        g.add(I[keep],J[keep],obs[keep],expected_global[keep])
                        r.add(I[keep],J[keep],obs[keep],expected_region[keep])
                        # delta=(observed - regional expected) -
                        #       (observed - global expected), same dyads.
                        d.add(I[keep],J[keep],
                              obs[keep]-expected_region[keep],
                              obs[keep]-expected_global[keep])
            for (gate,dist),(g,r,d) in sorted(accum.items()):
                G=g.export(min_pairs,min_islands,min_blocks)
                R=r.export(min_pairs,min_islands,min_blocks)
                D=d.export(min_pairs,min_islands,min_blocks)
                if G is None or R is None or D is None:continue
                if G["n_pairs"]!=R["n_pairs"] or G["n_pairs"]!=D["n_pairs"]:
                    raise AssertionError("Reference comparisons use different pairs")
                row=dict(scope=scope,region=region,support=gate,
                         distance_class=dist,
                         n_pairs=G["n_pairs"],n_active_islands=G["n_islands"],
                         n_physical_blocks=G["n_physical_blocks"],
                         mean_observed_D2=G["pair_equal_mean_observed_D2"],
                         global_genus_expected_D2=G["pair_equal_mean_genus_null_D2"],
                         region_genus_expected_D2=R["pair_equal_mean_genus_null_D2"],
                         global_genus_residual_D2=G["pair_equal_mean_residual_D2"],
                         region_genus_residual_D2=R["pair_equal_mean_residual_D2"],
                         region_minus_global_residual_D2=D["pair_equal_mean_residual_D2"],
                         island_equal_region_minus_global_residual_D2=D["island_equal_mean_residual_D2"],
                         block_equal_region_minus_global_residual_D2=D["physical_block_equal_mean_residual_D2"],
                         region_residual_leave_one_block_min=R["leave_one_physical_block_out_pair_mean_residual_D2_min"],
                         region_residual_leave_one_block_max=R["leave_one_physical_block_out_pair_mean_residual_D2_max"],
                         reference_shift_leave_one_block_min=D["leave_one_physical_block_out_pair_mean_residual_D2_min"],
                         reference_shift_leave_one_block_max=D["leave_one_physical_block_out_pair_mean_residual_D2_max"],
                         **audit)
                out.append(row)
    return dict(
        status="FLORAL_BETA_SAME_DYAD_REGION_GENUS_REFERENCE_SENSITIVITY",
        n_original_botanical_flora=len(all_islands),
        primary_endpoint="direct_flower_colour_x_form",
        eligible_by_scope=coverage,
        comparisons=out,
        claims=dict(
            all_pairs_reuse_identical_islands_across_two_nulls=True,
            exclude_shared_species_on_full_original_botanical_rosters=True,
            same_broad_region_and_different_physical_10degree_blocks=True,
            global_genus_reference_uses_all_direct_species=True,
            regional_genus_reference_in_sample_not_external=True,
            recorded_source_native_not_verified_evolutionary_origin=True,
            no_dyad_independent_inference=True,
            no_plant_adaptation_or_pollination_service_claim=True,
            no_private_taxon_island_source_geography_identifiers_exported=True,
        ))


def main():
    p=argparse.ArgumentParser()
    for name in ("original-root","geology","output"):
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args()
    paths={
        "traits":args.original_root/"fixed/source/current/evidence/direct_only_trait_ledger.csv.gz",
        "flora":args.original_root/"fixed/canonical/input/chapter1_status_flora.csv.gz",
    }
    sha={
        "traits":"9ddac1cedf2bf6235f3b426f37688550391ee8d96586e8d6a095ea3150bbe745",
        "flora":"6b5006322d15dd87aa3aad3b4dd050d4c6420485521be314ae511f39e05b859a",
    }
    for k,path in paths.items():
        h=hashlib.sha256()
        with path.open("rb") as fp:
            for chunk in iter(lambda:fp.read(1<<20),b""):h.update(chunk)
        if h.hexdigest()!=sha[k]:raise ValueError("Protected source checksum mismatch")
    trait=pd.read_csv(paths["traits"],dtype=str,usecols=[
        "accepted_species","trait_name","normalized_value",
        "evidence_scope","evidence_quality"])
    flora=pd.read_csv(paths["flora"],dtype=str,usecols=[
        "island_id","accepted_species","origin_status"])
    geo=pd.read_csv(args.geology,dtype=str)
    if flora.island_id.nunique()!=4505:
        raise ValueError("Botanical source universe changed")
    outcome=analyze(strict_traits(trait),flora,geo)
    if outcome["eligible_by_scope"]["all_recorded"]["n_eligible_islands"]!=2324:
        raise ValueError("Source pair-trait opportunity changed")
    if outcome["eligible_by_scope"]["original_native"]["n_eligible_islands"]!=327:
        raise ValueError("Source original-native opportunity changed")
    if not outcome["comparisons"]:
        raise ValueError("Original analytical result empty")
    body=json.dumps(outcome,sort_keys=True,indent=2,allow_nan=False)
    for private in ('"island_id"','"accepted_species"','"spatial_block"',
                    '"island_latitude"','"island_longitude"',
                    '"source_signature"'):
        if private in body:raise ValueError("Private identifiers in output")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(body+"\n",encoding="utf-8")
    print(json.dumps(dict(
        status=outcome["status"],n_original_flora=len(all_islands) if False else 4505,
        n_anonymous_strata=len(outcome["comparisons"]),
        not_a_causal_comparison=True)))


if __name__=="__main__":
    main()
