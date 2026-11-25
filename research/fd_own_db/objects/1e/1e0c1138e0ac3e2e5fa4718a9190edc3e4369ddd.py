"""Independent Curveball versus edge-swap null on IDENTICAL original DoPI sites.

Only original source publication×published site/grid×year×month positive
accepted plant × true visitor species edges. Curveball trades incident
visitor sets of two plants while preserving BOTH botanical and animal
positive-degrees exactly. Same original sites must pass both algorithms'
chain/swap-support gates before a paired publication estimate is reported.

Biological target is documented European insect-guild breadth, not local
GIFT island pollination, true nondetection, effective pollen delivery, or
adaptive evolution. Every original study/plant/animal/site ID stays private.
"""
from __future__ import annotations

from collections import Counter,defaultdict
from hashlib import sha256
import numpy as np

from .dopi_degree_fixed_visitor_null import (
    FOCAL_MOTIFS, N_CHAINS, N_SNAPSHOTS_PER_CHAIN,
    BURN_EDGE_MULTIPLIER, STEP_EDGE_MULTIPLIER, MIN_STUDIES,
    MAX_CHAIN_DISAGREEMENT_GUILDS, MIN_ACCEPTED_SWAPS_PER_CHAIN,
    MIN_MOTIF_PLANTS, MIN_DISTINCT_VISITORS, MIN_SITE_GUILDS,
    MAX_EDGES_PER_SITE, site_degree_preserving_null,
    _degrees, _positive_guild_breadth,
)

METHODS=("edge_swap","curveball")
FIELDS=("observed_guild_breadth_difference","edge_swap_expected",
        "curveball_expected","edge_swap_excess","curveball_excess",
        "curveball_minus_edge_swap_excess")
MAX_PAIRED_METHOD_EXPECTATION_DIFFERENCE=.25


def curveball_trades(rows:dict[str,set[str]],rng:np.random.Generator,
                     attempts:int)->tuple[dict[str,set[str]],int]:
    """True bipartite plant-neighbor Curveball trade (not label permutation).

    The *unique* visitor species of both chosen plants are allocated uniformly
    to the same original row lengths; shared visitor species remain in both
    rows. Column degrees and row degrees are unchanged. A no-change draw
    is not counted as an accepted state change.
    """
    if attempts<0:raise ValueError("negative Curveball trade proposals")
    current={k:set(v) for k,v in rows.items()}
    plants=sorted(current)
    if len(plants)<2:return current,0
    accepted=0
    for _ in range(attempts):
        i=int(rng.integers(len(plants)))
        j=int(rng.integers(len(plants)))
        if i==j:continue
        p,q=plants[i],plants[j]
        shared=current[p]&current[q]
        p_only=current[p]-shared
        q_only=current[q]-shared
        if not p_only or not q_only:continue
        pool=sorted(p_only|q_only)
        rng.shuffle(pool)
        new_p=set(pool[:len(p_only)])
        if new_p==p_only:continue
        new_q=set(pool[len(p_only):])
        current[p]=shared|new_p
        current[q]=shared|new_q
        accepted+=1
    return current,accepted


def _rows_to_edges(rows):
    return {(p,v) for p,vs in rows.items() for v in vs}


def site_curveball_null(edges:set[tuple[str,str]],visitor_guilds:dict,
                        motifs:dict,*,seed:str,draws:int=N_SNAPSHOTS_PER_CHAIN):
    if draws<20 or not motifs:
        raise ValueError("source Curveball requires at least 20 null draws")
    if len(edges)>MAX_EDGES_PER_SITE:
        return dict(status="NOT_IDENTIFIABLE_SOURCE_NETWORK_TOO_LARGE")
    if any(v not in visitor_guilds for _,v in edges):
        raise ValueError("real positive source visitor guild identity missing")
    plant_degrees,visitor_degrees=_degrees(edges)
    if (len(plant_degrees)<4 or len(visitor_degrees)<3 or
        len({visitor_guilds[v] for v in visitor_degrees})<MIN_SITE_GUILDS):
        return dict(status="NOT_IDENTIFIABLE_SOURCE_BIPARTITE_OPPORTUNITY")
    eligible={}
    for name,profile in motifs.items():
        known=set(profile["known"])
        joint=set(profile["joint"])
        positive=[p for p,n in sorted(plant_degrees.items())
                  if n>=MIN_DISTINCT_VISITORS and p in known and p in joint]
        negative=[p for p,n in sorted(plant_degrees.items())
                  if n>=MIN_DISTINCT_VISITORS and p in known and p not in joint]
        if min(len(positive),len(negative))>=MIN_MOTIF_PLANTS:
            eligible[name]=(positive,negative)
    if not eligible:
        return dict(status="NOT_IDENTIFIABLE_REAL_DIRECT_MOTIF_STATES")
    original_rows=defaultdict(set)
    for p,v in edges:original_rows[p].add(v)
    observed_breadth=_positive_guild_breadth(edges,visitor_guilds)
    observed={name:float(
        np.mean([observed_breadth[p] for p in pos])-
        np.mean([observed_breadth[p] for p in neg])
    ) for name,(pos,neg) in eligible.items()}
    n_edges=len(edges)
    accepted_trades=[]
    by_name={name:[] for name in eligible}
    for chain in range(N_CHAINS):
        seed_value=int.from_bytes(sha256(
            f"curveball|{seed}|chain{chain}".encode()).digest()[:8],"big")
        rng=np.random.default_rng(seed_value)
        state,naccept=curveball_trades(original_rows,rng,
                                     BURN_EDGE_MULTIPLIER*n_edges)
        chains={m:[] for m in eligible}
        for _ in range(draws):
            state,acc=curveball_trades(state,rng,STEP_EDGE_MULTIPLIER*n_edges)
            naccept+=acc
            real_edges=_rows_to_edges(state)
            if _degrees(real_edges)!=(plant_degrees,visitor_degrees):
                raise AssertionError("Curveball did not preserve observed degrees")
            breadth=_positive_guild_breadth(real_edges,visitor_guilds)
            for m,(pos,neg) in eligible.items():
                chains[m].append(float(
                    np.mean([breadth[p] for p in pos])-
                    np.mean([breadth[p] for p in neg])))
        accepted_trades.append(naccept)
        for m in eligible:by_name[m].append(float(np.mean(chains[m])))
    result={}
    for name,(pos,neg) in eligible.items():
        chain_values=by_name[name]
        d=max(chain_values)-min(chain_values)
        row=dict(
            original_observed_guild_richness_difference=observed[name],
            n_positive_motif_plants=len(pos),n_negative_motif_plants=len(neg),
            n_accepted_curveball_trades_each_chain=list(accepted_trades),
            chain_mean_difference_abs=float(d))
        if min(accepted_trades)<max(MIN_ACCEPTED_SWAPS_PER_CHAIN,n_edges//2):
            row["status"]="NOT_IDENTIFIABLE_CURVEBALL_NO_TRADABLE_OPPORTUNITY"
        elif d>MAX_CHAIN_DISAGREEMENT_GUILDS:
            row["status"]="NOT_IDENTIFIABLE_CURVEBALL_CHAIN_DISAGREEMENT"
        else:
            expectation=float(np.mean(chain_values))
            row.update(status="DESCRIPTIVE_CURVEBALL_BOTH_DEGREES_FIXED",
                       expected_curveball_guild_richness_difference=expectation,
                       excess_over_curveball_null=observed[name]-expectation)
        result[name]=row
    return dict(
        status="ORIGINAL_SOURCE_CURVEBALL_SITE_PROFILES",
        n_original_edges=n_edges,n_original_plants=len(plant_degrees),
        n_original_visitors=len(visitor_degrees),motif_profiles=result)


def _same_publication_bootstrap(by_study:dict,*,seed:str,draws:int)->dict:
    pubs=sorted(by_study)
    n=len(pubs)
    out=dict(n_original_source_publications=n,
             n_original_study_site_months=sum(map(len,by_study.values())))
    if n<MIN_STUDIES:
        out["status"]="NOT_IDENTIFIABLE_PAIRED_SOURCE_PUBLICATIONS"
        return out
    arr=np.asarray([np.mean([list(z[k] for k in FIELDS)
                             for z in by_study[pub]],axis=0)
                    for pub in pubs],float)
    if arr.shape!=(n,len(FIELDS)) or not np.isfinite(arr).all():
        raise ValueError("invalid exact original paired source article cohort")
    rng=np.random.default_rng(int.from_bytes(
        sha256(seed.encode()).digest()[:8],"big"))
    resample=arr[rng.integers(n,size=(draws,n))].mean(axis=1)
    mean=arr.mean(axis=0)
    for j,name in enumerate(FIELDS):
        out[name]=float(mean[j])
        out[name+"_ci95"]=[float(z) for z in
                           np.quantile(resample[:,j],[.025,.975])]
    # The pair difference in expected networks controls algorithm sensitivity,
    # not a second biological effect estimate.
    diffs=arr[:,4]-arr[:,3]
    out["n_source_papers_curveball_excess_negative"]=int((arr[:,4]<0).sum())
    out["n_source_papers_edge_swap_excess_negative"]=int((arr[:,3]<0).sum())
    out["n_source_papers_methods_have_different_effect_sign"]=int(
        (np.sign(arr[:,3])!=np.sign(arr[:,4])).sum())
    out["source_pub_leave_one_out_curveball_min"]=float(min(
        (arr[:,4].sum()-z)/(n-1) for z in arr[:,4]))
    out["source_pub_leave_one_out_curveball_max"]=float(max(
        (arr[:,4].sum()-z)/(n-1) for z in arr[:,4]))
    out["difference_between_null_algorithms_small_vs_0p25"]=bool(
        abs(mean[5])<=MAX_PAIRED_METHOD_EXPECTATION_DIFFERENCE)
    out["status"]="DESCRIPTIVE_IDENTICAL_SOURCE_SITE_TWO_DEGREE_NULLS"
    return out


def audit_curveball_vs_swap(index:dict,*,motifs:dict,draws:int=999,
                            site_draws:int=N_SNAPSHOTS_PER_CHAIN)->dict:
    """Same source sites for both null methods, no publication/site IDs exported."""
    if draws<199 or site_draws<20:
        raise ValueError("too few source or null resamples")
    if not index.get("documented_guild_absence_not_inferred"):
        raise ValueError("GIFT original positive observer source missing")
    edges_by_site=index.get("private_positive_edges",{})
    if not edges_by_site or not set(FOCAL_MOTIFS)<=set(motifs):
        raise ValueError("original plant×true insect source and registered motifs required")
    selected={m:motifs[m] for m in FOCAL_MOTIFS}
    by_motif={m:defaultdict(list) for m in FOCAL_MOTIFS}
    support=Counter()
    source_motif_support=Counter()
    for key,raw in edges_by_site.items():
        perpair={}
        visitor_guilds={}
        for plant,visitor,guild in raw:
            if visitor in visitor_guilds and visitor_guilds[visitor]!=guild:
                raise ValueError("source original visitor species changed guild identity")
            visitor_guilds[visitor]=guild
            pair=(plant,visitor)
            if pair in perpair and perpair[pair]!=guild:
                raise ValueError("same original plant insect edge has conflicting guild")
            perpair[pair]=guild
        edges=set(perpair)
        # Both independent algorithms start from EXACT same source-positive
        # incidence matrix, with its original GIFT Direct phenotype masks.
        original=site_degree_preserving_null(
            edges,visitor_guilds,selected,seed=repr(key),draws=site_draws)
        alternative=site_curveball_null(
            edges,visitor_guilds,selected,seed=repr(key),draws=site_draws)
        if (original["status"]!="ORIGINAL_SOURCE_DEGREE_FIXED_SITE_PROFILES"
            or alternative["status"]!="ORIGINAL_SOURCE_CURVEBALL_SITE_PROFILES"):
            support["source_sites_preflight_ineligible"]+=1
            continue
        publication=key[0]
        for name in FOCAL_MOTIFS:
            a=original["motif_profiles"].get(name)
            b=alternative["motif_profiles"].get(name)
            if a is None or b is None:continue
            source_motif_support[name]+=1
            if (a["status"]!="DESCRIPTIVE_FIXED_BOTH_SIDES_DEGREES_NETWORK_NULL"
                or b["status"]!="DESCRIPTIVE_CURVEBALL_BOTH_DEGREES_FIXED"):
                support["site_motif_not_both_chain_qualified"]+=1
                continue
            obs=float(a["original_observed_guild_richness_difference"])
            if abs(obs-b["original_observed_guild_richness_difference"])>1e-12:
                raise AssertionError("different original source phenotype or visitor cohort")
            swap=float(a["expected_from_fixed_species_degree_null"])
            curve=float(b["expected_curveball_guild_richness_difference"])
            row=dict(observed_guild_breadth_difference=obs,
                     edge_swap_expected=swap,curveball_expected=curve,
                     edge_swap_excess=obs-swap,curveball_excess=obs-curve,
                     curveball_minus_edge_swap_excess=swap-curve)
            by_motif[name][publication].append(row)
    panels=[]
    for name in FOCAL_MOTIFS:
        info=dict(motif=name,
                  n_source_sites_with_both_motif_states=source_motif_support[name])
        info.update(_same_publication_bootstrap(by_motif[name],
                seed=f"source|same-site|curveball|swap|{name}",draws=draws))
        panels.append(info)
    return dict(
        status="GIFT4505_DOPI2026_IDENTICAL_SITE_CURVEBALL_VS_EDGE_SWAP",
        n_original_GIFT_botanical_islands=4505,
        n_original_DoPI_v2_raw_rows=index["source_row_counts"].get(
            "original_DoPI_v2_source_rows",0),
        n_original_site_month_positive_networks=len(edges_by_site),
        n_source_study_keys=index["n_original_source_studies_qualified"],
        source_site_gate_diagnostics=dict(support),
        source_null_methods=list(METHODS),
        chains_per_method=N_CHAINS,
        snapshots_per_chain=site_draws,
        tested_morphologies=list(FOCAL_MOTIFS),
        panels=panels,
        limitations=dict(
            actual_accepted_botanical_and_visitor_species_degrees_both_fixed=True,
            same_exact_original_source_site_papers_both_null_algorithms=True,
            no_unobserved_visitors_or_guild_absences_inferred=True,
            study_month_location_originally_documented_not_same_flower_session=True,
            curveball_and_swap_mixing_not_proven_uniform=True,
            original_GIFT_island_flora_4505_not_local_visitor_observations=True,
            source_publication_bootstrap_not_independent_evolution=True,
            botanical_lineage_and_observation_effort_uncontrolled=True,
            four_focal_morphologies_posthoc_selected_from_27=True,
            no_island_pollinator_switch_or_effective_pollen_service_identified=True,
            original_taxon_site_study_ids_private_not_released=True))
