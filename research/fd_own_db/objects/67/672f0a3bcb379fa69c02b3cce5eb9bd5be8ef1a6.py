"""Potential mainland-source floral repertoire sensitivity (private inputs only).

This analysis removes a *descriptive* GIFT mainland source-repertoire vector
from each observed island repertoire vector, then tests paired origin residuals.
It is not a fully source-/genus-conditioned null, nor evidence of colonist
ancestry or geological selection.

Never write private species, island identities, mainland source IDs or taxon
profile assignments. Results contain only aggregate contrasts and input hashes.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .core import ProfileSpace, load_profiles
from .genus_source_null import _sha, load_private_maps
from .g1_fixed_margins import _bh
from .repertoire import (
    REGIONS, match_origin_islands, _paired_swap_statistic, _upper_tail,
)

SOURCE_MODES = ("geo_k5", "geo_k10", "geo_k20", "geo50_climate10")
SCORES = ("exact_profile_L2", "categorical_similarity_kernel")
CONTRASTS = ("original_observed", "source_prediction", "source_subtracted_residual")
COLUMNS = ("scope","source_mode","region","endpoint","contrast",
           "n_pairs","observed","paired_null_median","paired_null_q025",
           "paired_null_q975","p_pair_swap_upper","q_bh_primary4",
           "q_by_primary4","permutations")


def load_candidate_sources(assignments_path: Path, source_flora: Path,
                           island_ids: set[str]) -> tuple[dict, pd.DataFrame, dict]:
    """Read only needed source assignments and GIFT species for their entities."""
    assignments = pd.read_csv(
        assignments_path, dtype={"island_id":str,"entity_ID":str},
        usecols=["island_id","source_mode","entity_ID","source_rank"],
    )
    assignments = assignments.loc[assignments.island_id.isin(island_ids)].copy()
    if len(assignments) == 0:
        raise ValueError("No source assignment for matching islands")
    if not set(SOURCE_MODES).issubset(set(assignments.source_mode)):
        raise ValueError("Required fixed source-matching modes missing")
    seen = assignments.groupby(["source_mode","island_id"]).size()
    if any(set(assignments.loc[assignments.source_mode==m,"island_id"]) != island_ids
           for m in SOURCE_MODES):
        raise ValueError("One or more matched islands lack mainland source assignments")
    entities = set(assignments.entity_ID)
    chunks = []
    for chunk in pd.read_csv(
        source_flora, dtype={"entity_ID":str,"work_species":str},
        usecols=["entity_ID","work_species"], chunksize=200000,
    ):
        part = chunk.loc[chunk.entity_ID.isin(entities)]
        if not part.empty:
            chunks.append(part.copy())
    if not chunks:
        raise ValueError("No mainland source flora records")
    native = pd.concat(chunks,ignore_index=True).drop_duplicates()
    if native.work_species.isna().any():
        raise ValueError("Blank source native taxon names")
    return (
        {key: group.entity_ID.tolist()
         for key,group in assignments.groupby(["island_id","source_mode"])},
        native.rename(columns={"work_species":"accepted_species"}),
        {"assigned_source_entities":len(entities),
         "source_taxon_entity_rows":len(native),
         "n_islands":len(island_ids)},
    )


def _normalize_profile_set(profiles: set[int], n: int) -> np.ndarray:
    if not profiles:
        raise ValueError("Empty floral profile support cannot be interpreted as zero")
    if not set(profiles).issubset(set(range(n))):
        raise ValueError("Unknown reference profile ID")
    vector = np.zeros(n,dtype=float)
    vector[list(profiles)] = 1/len(profiles)
    return vector


def _source_and_observed(
    original_flora: Path, private_profiles: pd.DataFrame,
    source_native: pd.DataFrame,
    assignments: dict[tuple[str,str],list[str]],
    islands: set[str], n_reference: int, *, min_source_profiles: int = 10,
) -> tuple[dict[str,np.ndarray],dict[tuple[str,str],np.ndarray],dict]:
    # All source joins stay in memory; raw names/source IDs are never written.
    observed = pd.read_csv(
        original_flora,usecols=["island_id","accepted_species"],
        dtype=str,
    )
    observed = observed.loc[observed.island_id.isin(islands)].merge(
        private_profiles[["accepted_species","profile_id"]],
        on="accepted_species",how="inner",validate="many_to_one",
    )
    unique_obs = observed.groupby("island_id").profile_id.agg(lambda s:set(s))
    if len(unique_obs) != len(islands):
        raise ValueError("One matched island has no complete floral profiles")
    all_profiles = {
        island:_normalize_profile_set(ps,n_reference)
        for island,ps in unique_obs.items()
    }
    sourced = source_native.merge(
        private_profiles[["accepted_species","profile_id"]],
        on="accepted_species",how="inner",validate="many_to_one",
    )
    by_entity = sourced.groupby("entity_ID").profile_id.agg(lambda s:set(s)).to_dict()
    vectors,coverage = {}, {}
    for key,entities in assignments.items():
        island,mode=key
        if island not in islands:
            continue
        pool=set()
        for entity in entities:
            pool.update(by_entity.get(entity,set()))
        if len(pool)<min_source_profiles:
            raise ValueError("Insufficient source profiles for a matched island/mode")
        vectors[key]=_normalize_profile_set(pool,n_reference)
        coverage[key]={"observed_profiles":len(unique_obs[island]),
                       "source_profiles":len(pool),
                       "source_overlap":len(pool&unique_obs[island])/len(unique_obs[island])}
    return all_profiles,vectors,coverage


def run(private_all: Path, private_direct: Path, flora: Path,
        source_assignments: Path, mainland_native: Path,
        crosswalk: Path, g0_root: Path, reference: Path, output: Path,
        *, permutations: int=9999,seed: int=20261108,strict:bool=True) -> dict:
    if permutations < 999:
        raise ValueError("At least 999 label swaps required")
    paths = [Path(p) for p in (
        private_all,private_direct,flora,source_assignments,
        mainland_native,crosswalk,reference,
    )]
    private_all,private_direct,flora,source_assignments,mainland_native,crosswalk,reference = paths
    g0_root,output = Path(g0_root),Path(output)
    output.mkdir(parents=True,exist_ok=True)
    space=ProfileSpace(load_profiles(reference))
    profile_maps={
        "all":load_private_maps(private_all,"all",set(space.ids),strict=strict),
        "direct":load_private_maps(private_direct,"direct",set(space.ids),strict=strict),
    }
    geology=pd.read_csv(crosswalk,dtype={"island_id":str}).set_index("island_id")
    g0_all=pd.read_csv(g0_root/"all_islands.csv",dtype={"island_id":str})
    g0_direct=pd.read_csv(g0_root/"direct_islands.csv",dtype={"island_id":str})
    matched,_=match_origin_islands(
        g0_all.to_dict("records"),g0_direct.to_dict("records"),
        geology.to_dict("index"),
    )
    if strict and len(matched)!=37:
        raise ValueError("Matched cohort differs from frozen PR4 37 pairs")
    island_ids={p[key] for p in matched
                for key in ("continental_id","oceanic_id")}
    sources,source_native,audit=load_candidate_sources(
        source_assignments,mainland_native,island_ids,
    )
    D=space.matrix
    reference_mats={"exact_profile_L2":np.eye(len(space.ids)),
                    "categorical_similarity_kernel":1-D}
    rng=np.random.default_rng(seed)
    global_signs=rng.choice([-1.0,1.0],size=(permutations,len(matched)))
    results, diagnostic = [], []
    for scope,mapping in profile_maps.items():
        observed,source,coverage=_source_and_observed(
            flora,mapping,source_native,sources,island_ids,len(space.ids),
        )
        for mode in SOURCE_MODES:
            items=[coverage[island,mode] for island in island_ids]
            diagnostic.append({
                "scope":scope,"source_mode":mode,"n_islands":len(items),
                "median_source_profiles":float(np.median([x["source_profiles"] for x in items])),
                "median_observed_profiles":float(np.median([x["observed_profiles"] for x in items])),
                "median_observed_profile_overlap_with_source":float(
                    np.median([x["source_overlap"] for x in items])),
            })
            observed_delta=np.vstack([
                observed[p["oceanic_id"]]-observed[p["continental_id"]]
                for p in matched
            ])
            source_delta=np.vstack([
                source[p["oceanic_id"],mode]-source[p["continental_id"],mode]
                for p in matched
            ])
            residual=observed_delta-source_delta
            for region in ("global",*REGIONS):
                ix=np.arange(len(matched)) if region=="global" else np.asarray([
                    i for i,p in enumerate(matched) if p["region"]==region
                ],dtype=int)
                if not len(ix):
                    continue
                signs=global_signs if region=="global" else rng.choice(
                    [-1.0,1.0],size=(permutations,len(ix))
                )
                for endpoint,kernel in reference_mats.items():
                    for kind,arr in (
                        ("original_observed",observed_delta[ix]),
                        ("source_prediction",source_delta[ix]),
                        ("source_subtracted_residual",residual[ix]),
                    ):
                        stat,distribution=_paired_swap_statistic(arr,kernel,signs)
                        results.append({
                            "scope":scope,"source_mode":mode,"region":region,
                            "endpoint":endpoint,"contrast":kind,"n_pairs":len(ix),
                            "observed":float(stat),
                            "paired_null_median":float(np.median(distribution)),
                            "paired_null_q025":float(np.quantile(distribution,0.025)),
                            "paired_null_q975":float(np.quantile(distribution,0.975)),
                            "p_pair_swap_upper":float(_upper_tail(stat,distribution)),
                            "permutations":permutations,
                        })
    # Predeclared exploratory primary: geo_k5 source residual, 2 scopes x 2 endpoints.
    primary=[r for r in results if
             r["source_mode"]=="geo_k5" and r["region"]=="global"
             and r["contrast"]=="source_subtracted_residual"]
    if len(primary)!=4:
        raise AssertionError("Fixed source-remainder comparison family changed")
    q=_bh([r["p_pair_swap_upper"] for r in primary])
    for r,adjusted in zip(primary,q):
        r["q_bh_primary4"]=float(adjusted)
        r["q_by_primary4"]=float(min(1,adjusted*(1+1/2+1/3+1/4)))
    pd.DataFrame(results).reindex(columns=COLUMNS).to_csv(
        output/"source_repertoire_residual_summary.csv",index=False)
    pd.DataFrame(diagnostic).to_csv(
        output/"source_coverage_summary.csv",index=False)
    manifest={
        "schema":"fd_gift_mainland_source_repertoire_sensitivity_v1",
        "status":"exploratory_source_subtracted_no_causal_inference",
        "source_modes":list(SOURCE_MODES),
        "n_origin_matched_pairs":len(matched),
        "n_islands":len(island_ids),
        "source_audit":audit,
        "source_type":"GIFT 3.2 native mainland, archived in original island PR138 checkpoint",
        "species_source_matching":"exact work_species-to-accepted_species string match, no inferred synonyms",
        "comparison":"observed island profile vector minus union of nearby mainland entity profile vectors, then paired geological-origin sign swaps",
        "not_a_genus_null":True,
        "not_conditioned_on":["observed island genus richness","true colonist ancestor populations",
                              "phylogenetic branch lengths","nonrandom habitat/trait observation",
                              "floristic native status on islands"],
        "inference":"geo_k5 source residual BH/BY 4-family; other modes and regions are unadjusted sensitivities",
        "privacy":"no mainland entity IDs, private taxa, island IDs or restricted source rows exported",
        "input_sha256":{
            "private_all_profiles":_sha(private_all),
            "private_direct_profiles":_sha(private_direct),
            "private_flora":_sha(flora),
            "mainland_source_assignments":_sha(source_assignments),
            "mainland_native_flora":_sha(mainland_native),
            "geological_crosswalk":_sha(crosswalk),
            "public_profile_reference":_sha(reference),
        },
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"n_pairs":len(matched),"primary":[
        {"scope":r["scope"],"endpoint":r["endpoint"],
         "p":r["p_pair_swap_upper"],"BH":r["q_bh_primary4"],
         "BY":r["q_by_primary4"]} for r in primary]}))
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for label in ("private-all","private-direct","flora","source-assignments",
                  "mainland-native","crosswalk","g0","reference","output"):
        parser.add_argument("--"+label,type=Path,required=True)
    parser.add_argument("--permutations",type=int,default=9999)
    parser.add_argument("--seed",type=int,default=20261108)
    parser.add_argument("--allow-synthetic",action="store_true")
    a=parser.parse_args()
    run(a.private_all,a.private_direct,a.flora,a.source_assignments,
        a.mainland_native,a.crosswalk,a.g0,a.reference,a.output,
        permutations=a.permutations,seed=a.seed,strict=not a.allow_synthetic)


if __name__=="__main__":
    main()
