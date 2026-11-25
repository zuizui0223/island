"""Outcome-blind within-focal matching of strict whole-flora genus-disjoint
versus accepted-species-disjoint, genus-sharing island floral ensembles.

The floral response is never passed into matching. Source botanical island
parent N=4505, per-endpoint high/medium Direct SAME-species colour×form, with
separate all-recorded and original-native scopes. All private island/source/
species identifiers remain within original-source CI. No pollinator causality.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
import numpy as np
import pandas as pd
from fd.island_functional_beta import (
    prepare_endpoint, membership, geographical_distance,
    certified_genus_sets,
)

PRIMARY_FACTOR = 2.0
LOOSE_FACTOR = 4.0
RANDOM_SEED = 20261118


def _safe_logs(x):
    v = np.asarray(x, dtype=float)
    out = np.full(v.shape, np.nan)
    ok = np.isfinite(v) & (v > 0)
    out[ok] = np.log(v[ok])
    return out


def select_focal_triplets(
    regions, blocks, genus_certified, shared_species, shared_genera, km,
    n_direct, n_full, n_genera, area_km2, direct_fraction, *,
    factor=PRIMARY_FACTOR, donor_blocks_distinct=False,
    climate_pc1=None, require_climate=False,
):
    """Select at most one (focal, G0 donor, G1 donor) per focal ISLAND.

    Both donors have zero full-flora shared accepted species with focal,
    are independently genus-certified and in a different physical block.
    G0 has no shared genus and G1 shares >=1 genus with focal. Control
    donor is nearest by OUTCOME-BLIND covariate distance to G0 donor.
    Both donors and focal are in same broad region. Pairwise calipers
    apply to the two donors, NOT to the focal morphology.
    """
    r=np.asarray(regions,dtype=object)
    b=np.asarray(blocks,dtype=object)
    certified=np.asarray(genus_certified,bool)
    n=len(r)
    a=np.asarray(shared_species)
    g=np.asarray(shared_genera)
    km=np.asarray(km,float)
    if any(mat.shape!=(n,n) for mat in (a,g,km)):
        raise ValueError("Original island-by-island source matrices misaligned")
    if len(b)!=n or len(certified)!=n:
        raise ValueError("Unmatched certified island identities")
    if factor<=1:raise ValueError("Matching factor must exceed 1")
    logs=np.stack([_safe_logs(n_direct),_safe_logs(n_full),
        _safe_logs(n_genera),_safe_logs(area_km2),
        _safe_logs(direct_fraction)],axis=1)
    lc=np.log(factor)
    pc=(np.full(n,np.nan) if climate_pc1 is None
        else np.asarray(climate_pc1,float))
    if len(pc)!=n:raise ValueError("Climate source alignment lost")
    chosen=[]
    opportunity=dict(n_focal_genus_disjoint_candidate=0,
        n_focal_genus_shared_control_candidate=0,
        n_focal_having_both_unmatched=0)
    for i in range(n):
        if not certified[i] or pd.isna(r[i]) or pd.isna(b[i]) or not (
                np.isfinite(logs[i]).all()):
            continue
        valid=(certified & (r==r[i]) & (b!=b[i]) &
            pd.notna(b) & (a[i]==0) & np.isfinite(km[i]) &
            (km[i]>0) & np.isfinite(logs).all(axis=1))
        valid[i]=False
        g0=np.flatnonzero(valid & (g[i]==0))
        g1=np.flatnonzero(valid & (g[i]>0))
        if len(g0):opportunity["n_focal_genus_disjoint_candidate"]+=1
        if len(g1):opportunity["n_focal_genus_shared_control_candidate"]+=1
        if not len(g0) or not len(g1):continue
        opportunity["n_focal_having_both_unmatched"]+=1
        best=None
        for j in g0:
            # Five donor-level matching covariates plus focal-donor distance.
            # Island area is given a 4x factor, other factors use 'factor'.
            dlog=np.abs(logs[g1]-logs[j])
            lim=np.array([lc,lc,lc,np.log(max(4.,factor)),lc])
            dist=np.abs(np.log(km[i,g1]/km[i,j]))
            ok=np.isfinite(dist) & (dist<=lc+1e-12) & (
                dlog<=lim+1e-12).all(axis=1)
            if donor_blocks_distinct:
                ok &= b[g1]!=b[j]
            if require_climate:
                ok &= np.isfinite(pc[j]) & np.isfinite(pc[g1]) & (
                    np.abs(pc[g1]-pc[j])<=1.5)
            if not ok.any():continue
            eligible=g1[ok]
            # Equal-weight normalized (outcome-blind) covariate discrepancy.
            cost=(dlog[ok]/lim).sum(axis=1)+dist[ok]/lc
            order=np.lexsort((eligible,cost))
            k=int(order[0])
            test=(float(cost[k]),int(j),int(eligible[k]))
            if best is None or test<best:best=test
        if best is not None:chosen.append((int(i),best[1],best[2],
                                           float(best[0])))
    return chosen, opportunity


def _connected_components(edges):
    # Physical donor and focal block graph: a stress diagnostic only.
    parent={}
    def find(x):
        parent.setdefault(x,x)
        while parent[x]!=x:
            parent[x]=parent[parent[x]]
            x=parent[x]
        return x
    def join(a,b):
        a=find(a);b=find(b)
        if a!=b:parent[a]=b
    for a,b,c in edges:
        join(a,b);join(a,c)
    return len({find(x) for x in parent})


def _greedy_block_disjoint_stress(triples, blocks, observed_contrast,
                                  *, draws=199):
    """Outcome-blind three-block-disjoint hyperedge reselections.

    Every selected focal, G0 and G1 donor occupies a unique spatial 10deg
    block globally within one reselection. Distinct matches never reuse a
    block. Greedy choices use independent seeded priorities only, NEVER
    observed floral states. The reselection interval is NOT a sampling CI.
    """
    b=np.asarray(blocks,dtype=object)
    contrast=np.asarray(observed_contrast,float)
    if len(triples)!=len(contrast):
        raise ValueError("Unmatched source outcome and original candidate triplets")
    admissible=[]
    for k,(i,j0,j1,score) in enumerate(triples):
        ids=(str(b[i]),str(b[j0]),str(b[j1]))
        if all(pd.notna(b[x]) for x in (i,j0,j1)) and len(set(ids))==3:
            admissible.append((k,ids))
    if not admissible:
        return dict(n_three_block_eligible_triplets=0,
            n_greedy_reselections=draws,
            three_block_disjoint_support_unavailable=True)
    rng=np.random.default_rng(RANDOM_SEED+17)
    n=[];effect=[]
    for iteration in range(draws):
        selected=[]
        occupied=set()
        for pos in rng.permutation(len(admissible)):
            k,ids=admissible[int(pos)]
            if not occupied.intersection(ids):
                selected.append(k)
                occupied.update(ids)
        if not selected:raise AssertionError("Original geographic disjoint match empty")
        n.append(len(selected))
        effect.append(float(contrast[selected].mean()))
    return dict(n_three_block_eligible_triplets=len(admissible),
        n_greedy_reselections=draws,
        independent_block_triplets_min=int(min(n)),
        independent_block_triplets_median=float(np.median(n)),
        independent_block_triplets_max=int(max(n)),
        independent_block_triplet_effect_median=float(np.median(effect)),
        independent_block_triplet_effect_p05=float(np.quantile(effect,.05)),
        independent_block_triplet_effect_p95=float(np.quantile(effect,.95)),
        proportion_outcome_blind_reselections_positive=float(
            np.mean(np.asarray(effect)>0)),
        reselection_interval_not_statistical_confidence_interval=True,
        spatial_blocks_disjoint_not_independent_genus_or_source_histories=True)


def _summarize(triples,prepared,blocks,archipelago,*,boot=499):
    if not triples:return None
    I=np.array([x[0] for x in triples],int)
    G0=np.array([x[1] for x in triples],int)
    G1=np.array([x[2] for x in triples],int)
    P=prepared["P"]
    o0=np.minimum(P[I],P[G0]).sum(axis=1)
    o1=np.minimum(P[I],P[G1]).sum(axis=1)
    d=o0-o1
    blocks=np.asarray(blocks,dtype=object)
    arch=np.asarray(archipelago,dtype=object)
    groups=np.array([str(x) for x in blocks[I]],dtype=object)
    unique=sorted(set(groups))
    means=np.array([d[groups==x].mean() for x in unique],float)
    mean0=float(o0.mean());mean1=float(o1.mean())
    res=dict(
        n_matched_focals=int(len(I)),
        n_matched_distinct_focal_blocks=int(len(unique)),
        n_distinct_treatment_donors=int(len(set(G0))),
        n_distinct_control_donors=int(len(set(G1))),
        n_distinct_donor_physical_blocks=int(
            len(set(blocks[np.r_[G0,G1]]))),
        n_named_focal_archipelagos=int(len({
            str(x) for x in arch[I] if pd.notna(x) and str(x).strip()})),
        n_spatial_block_network_components=_connected_components(
            zip(blocks[I],blocks[G0],blocks[G1])),
        mean_G0_overlap=float(mean0),mean_G1_overlap=float(mean1),
        focal_equal_G0_minus_G1=float(d.mean()),
        focal_block_equal_G0_minus_G1=float(means.mean()),
        mean_outcome_blind_matching_score=float(
            np.mean([x[3] for x in triples])),
        source_gift_first_mainland_region_unavailable=True,
        **_greedy_block_disjoint_stress(triples,blocks,d),
    )
    if len(unique)>=10 and boot>=99:
        rng=np.random.default_rng(RANDOM_SEED)
        take=rng.integers(len(means),size=(boot,len(means)))
        m=means[take].mean(axis=1)
        res["focal_block_bootstrap_95"]=list(
            map(float,np.quantile(m,[.025,.975])))
        res["focal_block_bootstrap_ignores_reused_donor_blocks"]=True
        res["bootstrap_is_exploratory_not_causal"]=True
    else:
        res["focal_block_bootstrap_95"]=None
        res["insufficient_focal_blocks_for_bootstrap"]=True
    return res


def analyze_matched(states:dict, flora:pd.DataFrame, geo:pd.DataFrame,
                    *,calipers=(2.,4.),boot=499)->dict:
    required={"island_id","accepted_species","origin_status"}
    cols={"island_id","analysis_regime","spatial_block",
        "island_latitude","island_longitude","area_km2","climate_pc1",
        "weigelt_archipelago"}
    if not required<=set(flora) or not cols<=set(geo):
        raise ValueError("Original source flora and geography fields required")
    if (flora.duplicated(["island_id","accepted_species"]).any()
        or geo.island_id.duplicated().any()):
        raise ValueError("Original accepted species/geography IDs duplicate")
    if not flora.origin_status.isin(("native","introduced","unresolved")).all():
        raise ValueError("Original nativeness token unqualified")
    if not {"flower_primary_color","floral_form"}<=set(states):
        raise ValueError("Original direct traits missing")
    two=set(states["flower_primary_color"]) & set(states["floral_form"])
    joint={s:(states["flower_primary_color"][s],states["floral_form"][s])
           for s in two}
    full_all=flora.groupby("island_id").accepted_species.agg(set).to_dict()
    native=flora.loc[flora.origin_status.eq("native")].groupby(
        "island_id").accepted_species.agg(set).to_dict()
    population=sorted(full_all)
    joined=geo.set_index("island_id").reindex(population)
    output=[]
    scope_n={}
    for scope,sources in (("all_recorded",full_all),("original_native",native)):
        elig=[island for island in population if len(
            sources.get(island,set())&two)>=5]
        scope_n[scope]=len(elig)
        if len(elig)<3:continue
        full=[sources.get(x,set()) for x in elig]
        observed=[x&two for x in full]
        pg=prepare_endpoint(observed,joint)
        genus,cert=certified_genus_sets(full)
        # Only source-accepted whole-flora species and genera define exclusions.
        fids=sorted(set().union(*full))
        Gids=sorted(set().union(*genus))
        FS=membership(full,fids)
        GS=membership(genus,Gids)
        shared_species=(FS@FS.T).toarray()
        shared_genera=(GS@GS.T).toarray()
        fr=joined.reindex(elig)
        regions=fr.analysis_regime.to_numpy()
        blocks=fr.spatial_block.to_numpy()
        lat=pd.to_numeric(fr.island_latitude,errors="coerce").to_numpy()
        lon=pd.to_numeric(fr.island_longitude,errors="coerce").to_numpy()
        km=geographical_distance(lat,lon,0,len(elig))
        n_direct=np.array([len(x) for x in observed],float)
        n_full=np.array([len(x) for x in full],float)
        n_genus=np.array([len(x) for x in genus],float)
        area=pd.to_numeric(fr.area_km2,errors="coerce").to_numpy()
        pc=pd.to_numeric(fr.climate_pc1,errors="coerce").to_numpy()
        frac=n_direct/n_full
        arch=fr.weigelt_archipelago.to_numpy()
        for factor in calipers:
            for mode,block_disjoint,climate in (
                    ("primary",False,False),
                    ("donor_blocks_distinct",True,False),
                    ("climate_pc1_caliper",False,True)):
                triples,opportunity=select_focal_triplets(
                    regions,blocks,cert,shared_species,shared_genera,km,
                    n_direct,n_full,n_genus,area,frac,
                    factor=factor,donor_blocks_distinct=block_disjoint,
                    climate_pc1=pc,require_climate=climate)
                for region in ("GLOBAL_ALL","northern_midlatitude",
                               "tropical","southern_extratropical",
                               "northern_highlatitude"):
                    group=[t for t in triples if (region=="GLOBAL_ALL" or
                           regions[t[0]]==region)]
                    report=_summarize(group,pg,blocks,arch,boot=boot)
                    if report is None:continue
                    output.append(dict(scope=scope,region=region,mode=mode,
                        covariate_factor=float(factor),
                        **(opportunity if region=="GLOBAL_ALL" else {}),
                        **report))
    return dict(status="GENUS_DISJOINT_MATCHED_FLORAL_BETA_DESCRIPTIVE",
        n_source_accepted_island_floras=len(population),
        eligible_direct_islands=scope_n,
        results=output,
        claims=dict(
            matched_outcome_blind=True,
            all_flora_roster_genus_not_just_direct_species=True,
            same_focal_G0_vs_genus_shared_G1=True,
            same_region_full_species_disjoint_and_diff_block=True,
            donor_km_direct_n_full_n_genus_fraction_area_calipers=True,
            double_donor_block_network_nonindependence_audited=True,
            first_gift_source_donor_unavailable_in_geographic_crosswalk=True,
            no_matched_island_pairs_considered_independent_replicates=True,
            no_reproductive_fitness_effective_visitation_or_adaptation=True,
            no_private_island_taxon_source_or_geo_identifiers_exported=True,
        ))
