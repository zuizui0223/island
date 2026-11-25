"""Degree-preserving positive DoPI flower-visitor network null.

Randomize only true accepted source positive plant×visitor-species EDGES
inside the SAME published study×site×year×month, preserving:
  * EACH accepted plant's actual distinct visitor species degree;
  * EACH real named visitor species' actual positive plant degree;
  * the guild identity of every visitor species; and
  * the exact set of positive source botanical and insect names.

This isolates a network configuration reference, not plant visitation
choice conditional on local floral availability, true animal absence,
within-island evolution, effective pollen deposition or source effort.
All source identities are private and never serialized into output.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from math import comb
import numpy as np

from .dopi_gift_guild_motif_bridge import botanical_genus

FOCAL_MOTIFS=("open_radial","tubular_bell","zygomorphic","deep_tube")
MIN_STUDIES=10
MIN_MOTIF_PLANTS=2
MIN_DISTINCT_VISITORS=2
MIN_SITE_GUILDS=2
MAX_EDGES_PER_SITE=400
N_CHAINS=4
N_SNAPSHOTS_PER_CHAIN=80
BURN_EDGE_MULTIPLIER=25
STEP_EDGE_MULTIPLIER=4
MIN_ACCEPTED_SWAPS_PER_CHAIN=10
MAX_CHAIN_DISAGREEMENT_GUILDS=.20


def bipartite_degree_swap(edges:set[tuple[str,str]],rng:np.random.Generator,
                          attempts:int)->tuple[set[tuple[str,str]],int]:
    """Swap (plant1,visitor1)+(plant2,visitor2) into crossed unobserved edges.

    Edges retain unique plant and visitor degrees EXACTLY (no duplicates).
    Visitor guild is an immutable species attribute held outside this model.
    """
    if attempts<0:raise ValueError("invalid swap attempts")
    work=set(edges)
    if not work:return work,0
    if any(len(x)!=2 for x in work):raise ValueError("malformed original positive edge")
    rows=sorted(work)
    accepted=0
    n=len(rows)
    for _ in range(attempts):
        i=int(rng.integers(n));j=int(rng.integers(n))
        if i==j:continue
        p1,v1=rows[i];p2,v2=rows[j]
        if p1==p2 or v1==v2:continue
        x=(p1,v2);y=(p2,v1)
        if x in work or y in work:continue
        work.remove(rows[i]);work.remove(rows[j])
        work.add(x);work.add(y)
        rows[i]=x;rows[j]=y
        accepted+=1
    if len(work)!=len(edges):raise AssertionError("positive source edge count drift")
    return work,accepted


def _degrees(edges):
    plants=defaultdict(int);visitors=defaultdict(int)
    for p,v in edges:
        plants[p]+=1;visitors[v]+=1
    return dict(plants),dict(visitors)


def _positive_guild_breadth(edges,visitor_guilds):
    per_plant=defaultdict(set)
    for p,v in edges:
        per_plant[p].add(visitor_guilds[v])
    return {p:len(g) for p,g in per_plant.items()}


def site_degree_preserving_null(
    edges:set[tuple[str,str]],visitor_guilds:dict[str,str],motifs:dict,
    *,seed:str,draws:int=N_SNAPSHOTS_PER_CHAIN,
)->dict:
    """Source-private site data with real fauna-degree-preserving null estimates."""
    if draws<20 or not motifs:
        raise ValueError("null requires >=20 snapshots and registered Direct flower states")
    if any(v not in visitor_guilds for _,v in edges):
        raise ValueError("original positive source visitor species missing guild identity")
    if len(edges)>MAX_EDGES_PER_SITE:
        return dict(status="NOT_IDENTIFIABLE_SOURCE_NETWORK_TOO_LARGE",
                    n_edges=len(edges))
    dp,dv=_degrees(edges)
    guilds=set(visitor_guilds[v] for v in dv)
    if len(dp)<4 or len(dv)<3 or len(guilds)<MIN_SITE_GUILDS:
        return dict(status="NOT_IDENTIFIABLE_SOURCE_BIPARTITE_OPPORTUNITY",
                    n_edges=len(edges),n_visitor_species=len(dv),
                    n_positive_plant_species=len(dp),n_source_guilds=len(guilds))
    eligible={}
    for name,profile in motifs.items():
        known=set(profile["known"])
        yes=set(profile["joint"])
        positive=sorted(p for p in dp if dp[p]>=MIN_DISTINCT_VISITORS
                        and p in known and p in yes)
        negative=sorted(p for p in dp if dp[p]>=MIN_DISTINCT_VISITORS
                        and p in known and p not in yes)
        if min(len(positive),len(negative))>=MIN_MOTIF_PLANTS:
            eligible[name]=(positive,negative)
    if not eligible:
        return dict(status="NOT_IDENTIFIABLE_REAL_DIRECT_MOTIF_STATES",
                    n_edges=len(edges),n_visitor_species=len(dv),
                    n_positive_plant_species=len(dp),n_source_guilds=len(guilds))

    actual=_positive_guild_breadth(edges,visitor_guilds)
    obs={k:float(np.mean([actual[p] for p in a])-
                 np.mean([actual[p] for p in b]))
         for k,(a,b) in eligible.items()}
    # Run two independently seeded chains from the true observed source
    # network, using different initial random generators. Never scramble
    # botanical trait labels as substitutes for network rewiring.
    chain_means={k:[] for k in eligible}
    acceptance=[]
    E=len(edges)
    for chain in range(N_CHAINS):
        digest=sha256(f"{seed}|chain{chain}".encode()).digest()
        rng=np.random.default_rng(int.from_bytes(digest[:8],"big"))
        current,swaps=bipartite_degree_swap(
            edges,rng,BURN_EDGE_MULTIPLIER*E)
        acc=swaps
        observed_samples={k:[] for k in eligible}
        for _ in range(draws):
            current,s=bipartite_degree_swap(
                current,rng,STEP_EDGE_MULTIPLIER*E)
            acc+=s
            if _degrees(current)!=(dp,dv):
                raise AssertionError("source visitor or plant degree broken")
            guild_counts=_positive_guild_breadth(current,visitor_guilds)
            for k,(a,b) in eligible.items():
                observed_samples[k].append(
                    float(np.mean([guild_counts[p] for p in a])-
                          np.mean([guild_counts[p] for p in b])))
        acceptance.append(acc)
        for k in eligible:
            chain_means[k].append(float(np.mean(observed_samples[k])))
    result={}
    for name,(a,b) in eligible.items():
        vals=chain_means[name]
        disagreement=float(max(vals)-min(vals))
        rec=dict(n_positive_motif_plants=len(a),
                 n_negative_motif_plants=len(b),
                 n_original_visitor_species=len(dv),
                 n_original_positive_source_edges=E,
                 original_observed_guild_richness_difference=obs[name],
                 n_accepted_degree_swaps_each_chain=list(acceptance),
                 chain_mean_difference_abs=float(disagreement))
        if min(acceptance)<max(MIN_ACCEPTED_SWAPS_PER_CHAIN,E//2):
            rec["status"]="NOT_IDENTIFIABLE_SOURCE_NETWORK_NO_SWAPPABLE_OPPORTUNITY"
        elif disagreement>MAX_CHAIN_DISAGREEMENT_GUILDS:
            rec["status"]="NOT_IDENTIFIABLE_SOURCE_NETWORK_CHAIN_DISAGREEMENT"
        else:
            expected=float(np.mean(vals))
            rec.update(
                status="DESCRIPTIVE_FIXED_BOTH_SIDES_DEGREES_NETWORK_NULL",
                expected_from_fixed_species_degree_null=expected,
                excess_over_network_degree_null=obs[name]-expected,
            )
        result[name]=rec
    return dict(status="ORIGINAL_SOURCE_DEGREE_FIXED_SITE_PROFILES",
                n_original_edges=E,
                n_plant_species=len(dp),
                n_visitor_species=len(dv),
                n_guilds=len(guilds),
                motif_profiles=result)


def _publication_sign_pvalue(values:list[float])->dict:
    """Exact signed *publication*-median test; not a causal or prespecified p."""
    npos=sum(x>1e-12 for x in values)
    nneg=sum(x<-1e-12 for x in values)
    ntie=len(values)-npos-nneg
    n=npos+nneg
    if n==0:
        return dict(n_positive=0,n_negative=0,n_exact_zero_ties=ntie,
                    two_sided_sign_p=None)
    tail=sum(comb(n,j) for j in range(min(npos,nneg)+1))/2**n
    return dict(n_positive=npos,n_negative=nneg,n_exact_zero_ties=ntie,
                two_sided_sign_p=float(min(1,2*tail)))


def _familywise_study_bootstrap(by_motif:dict,draws:int)->dict:
    """Correlated max-t bootstrap over shared original study labels.

    Align support by original publication, not by individual network edges.
    Each motif's statistic remains conditional on its own study availability.
    Family is these FOUR explicitly selected exploratory morphology profiles,
    not the 27 earlier searched, so DO NOT claim globally corrected evidence.
    """
    good=[name for name in FOCAL_MOTIFS if len(by_motif[name])>=MIN_STUDIES]
    if not good:return dict(status="NOT_IDENTIFIABLE_FOUR_PROFILE_FAMILY")
    keys=sorted(set().union(*(set(by_motif[m]) for m in good)))
    if len(keys)<MIN_STUDIES:return dict(status="NOT_IDENTIFIABLE_SOURCE_STUDIES")
    # Critical: missing original paper/motif observations are *missing*,
    # never manufactured source-negative visitor observations.
    rows=np.full((len(keys),len(good)),np.nan)
    for i,pub in enumerate(keys):
        for j,motif in enumerate(good):
            if pub in by_motif[motif]:
                rows[i,j]=float(np.mean([z["excess_over_network_degree_null"]
                                          for z in by_motif[motif][pub]]))
    original=np.nanmean(rows,axis=0)
    seed=int.from_bytes(sha256(b"fd|DoPI|degree-fixed|source-pub|maxT").digest()[:8],"big")
    rng=np.random.default_rng(seed)
    sampled=rng.integers(len(keys),size=(draws,len(keys)))
    # All resampled studies must include >=1 actual original source paper
    # with available data for each eligible morphology; otherwise resample
    # cannot estimate the same target and must be excluded.
    bs=np.empty((draws,len(good)),float)
    bs[:]=np.nan
    for d,ix in enumerate(sampled):
        x=rows[ix,:]
        good_cols=np.isfinite(x).any(axis=0)
        if good_cols.all():
            bs[d,:]=np.nanmean(x,axis=0)
    valid=np.isfinite(bs).all(axis=1)
    if int(valid.sum())<int(draws*.99):
        return dict(status="NOT_IDENTIFIABLE_STUDY_FAMILY_RESAMPLING",
                    n_resamples_valid=int(valid.sum()))
    bs=bs[valid]
    std=bs.std(axis=0,ddof=1)
    if np.any(std<=0):
        return dict(status="NOT_IDENTIFIABLE_ZERO_STUDY_VARIANCE",
                    n_resamples_valid=int(valid.sum()))
    max_t=np.max(np.abs((bs-original)/std),axis=1)
    threshold=float(np.quantile(max_t,.95))
    out={}
    for j,m in enumerate(good):
        out[m]=dict(
            estimate=float(original[j]),
            simultaneous_selected_four_profiles_ci95=[
                float(original[j]-threshold*std[j]),
                float(original[j]+threshold*std[j])],
            bootstrap_sd=float(std[j]),
            n_studies=int(np.isfinite(rows[:,j]).sum()))
    return dict(status="EXPLORATORY_FOUR_PROFILE_STUDY_MAXT_INTERVALS",
                n_unique_source_publications=len(keys),
                n_valid_publication_resamples=int(valid.sum()),
                max_studentized_threshold95=threshold,
                profiles=out,
                four_selected_profiles_not_global_multiplicity=True,
                nonexchangeable_source_studies_and_lineages=True)


def audit_original_degree_fixed_networks(index:dict,*,motifs:dict,
                                         draws:int=999,site_draws:int=40)->dict:
    """Only source-evidenced study/site/month cohorts; never actual GIFT island visits."""
    if draws<199 or site_draws<20:
        raise ValueError("insufficient source site/republication re-sampling")
    if not index.get("documented_guild_absence_not_inferred"):
        raise ValueError("original positive-only DoPI index missing")
    private=index.get("private_positive_edges")
    if not private or not motifs or not set(FOCAL_MOTIFS)<=set(motifs):
        raise ValueError("original degree-fixed source support or focal motif missing")

    selected={k:motifs[k] for k in FOCAL_MOTIFS}
    by_motif={name:defaultdict(list) for name in FOCAL_MOTIFS}
    sites_with_candidates=defaultdict(set)
    n_sites=defaultdict(int)
    n_rejected=defaultdict(int)
    n_source_sites=len(private)
    for site_key,raw in private.items():
        per_pair={}
        guild_by_visitor={}
        for plant,visitor,guild in raw:
            if visitor in guild_by_visitor and guild_by_visitor[visitor]!=guild:
                raise ValueError("one original visitor accepted binomial has multiple mutually conflicting guilds")
            guild_by_visitor[visitor]=guild
            pair=(plant,visitor)
            if pair in per_pair and per_pair[pair]!=guild:
                raise ValueError("source plant visitor pair has conflicting taxon family")
            per_pair[pair]=guild
        edge_set=set(per_pair)
        diag=site_degree_preserving_null(
            edge_set,guild_by_visitor,selected,seed=repr(site_key),
            draws=site_draws,
        )
        if diag["status"]!="ORIGINAL_SOURCE_DEGREE_FIXED_SITE_PROFILES":
            n_rejected[diag["status"]]+=1
            continue
        publication=site_key[0]
        for name,rec in diag["motif_profiles"].items():
            n_sites[name]+=1
            if rec["status"]!="DESCRIPTIVE_FIXED_BOTH_SIDES_DEGREES_NETWORK_NULL":
                n_rejected[rec["status"]]+=1
                continue
            by_motif[name][publication].append(rec)
            sites_with_candidates[name].add(publication)
    panels=[]
    for name in FOCAL_MOTIFS:
        by_pub=by_motif[name]
        n=len(by_pub)
        item=dict(motif=name,
            n_real_original_study_site_month_null_support=n_sites[name],
            n_source_publications_with_degree_fixed_null=n,
            n_source_qualifying_site_months=sum(map(len,by_pub.values())),
            n_accepted_positive_plant_edges_in_qualifying_sites=sum(
                r["n_original_positive_source_edges"]
                for rows in by_pub.values() for r in rows),
            n_distinct_positive_source_visitor_species_site_units=sum(
                r["n_original_visitor_species"]
                for rows in by_pub.values() for r in rows))
        if n<MIN_STUDIES:
            item["status"]="NOT_IDENTIFIABLE_INDEPENDENT_SOURCE_PUBLICATIONS"
            panels.append(item)
            continue
        arr=np.asarray([[np.mean([r[k] for r in by_pub[pub]])
            for k in ("original_observed_guild_richness_difference",
                      "expected_from_fixed_species_degree_null",
                      "excess_over_network_degree_null")]
            for pub in sorted(by_pub)],dtype=float)
        if arr.shape!=(n,3) or not np.isfinite(arr).all():
            raise ValueError("original source null cannot aggregate finite degree-fixed values")
        digest=sha256(f"{name}|publication".encode()).digest()
        rng=np.random.default_rng(int.from_bytes(digest[:8],"big"))
        boot=arr[rng.integers(n,size=(draws,n))].mean(axis=1)
        item["status"]="DESCRIPTIVE_SOURCE_VISITOR_DEGREE_FIXED_FLORAL_GUILD_ASSOCIATION"
        for j,k in enumerate(("observed_guild_difference",
                              "expected_degree_fixed_network_difference",
                              "residual_over_degree_fixed_network_null")):
            item[k]=float(arr[:,j].mean())
            item[k+"_bootstrap95"]=[
                float(x) for x in np.quantile(boot[:,j],[.025,.975])]
        pub_residuals=arr[:,2].tolist()
        item["publication_sign_test"]=_publication_sign_pvalue(pub_residuals)
        item["source_publication_loo_residual_min"]=float(
            min((sum(pub_residuals)-v)/(n-1) for v in pub_residuals))
        item["source_publication_loo_residual_max"]=float(
            max((sum(pub_residuals)-v)/(n-1) for v in pub_residuals))
        site_counts=[len(by_pub[pub]) for pub in sorted(by_pub)]
        item["largest_one_publication_share_of_source_sites"]=float(
            max(site_counts)/sum(site_counts))
        item["n_source_papers_with_opposite_direction_than_full_mean"]=int(
            sum(np.sign(v)!=np.sign(item["residual_over_degree_fixed_network_null"])
                for v in pub_residuals))
        panels.append(item)
    # A transparent 4-profile family check, preserving overlapping
    # publication identities across *different* available motif cohorts.
    family=_familywise_study_bootstrap(by_motif,draws)
    # Holm over exact source-paper sign tests is another descriptive control.
    available=[(p["motif"],p["publication_sign_test"]["two_sided_sign_p"])
        for p in panels if p.get("status")==
        "DESCRIPTIVE_SOURCE_VISITOR_DEGREE_FIXED_FLORAL_GUILD_ASSOCIATION"
        and p["publication_sign_test"]["two_sided_sign_p"] is not None]
    ordered=sorted(available,key=lambda pair:pair[1])
    adjusted={}
    running=0.
    for i,(name,p) in enumerate(ordered):
        running=max(running,min(1.,float(p)*(len(ordered)-i)))
        adjusted[name]=float(running)
    for p in panels:
        if p["motif"] in adjusted:
            p["holm_four_focal_publication_sign_p"]=adjusted[p["motif"]]
    return dict(
        status="DOPI_GIFT_4505_SOURCE_VISITOR_SPECIES_DEGREE_FIXED_NETWORK_NULL",
        selected_four_profile_familywise_source_study_bootstrap=family,
        original_GIFT_botanical_island_population=4505,
        original_DoPI_raw_source_rows=index["source_row_counts"].get(
            "original_DoPI_v2_source_rows",0),
        n_original_source_studies=index["n_original_source_studies_qualified"],
        n_original_published_site_months=n_source_sites,
        tested_flower_profiles=list(FOCAL_MOTIFS),
        source_site_null_chains=N_CHAINS,
        source_site_null_snapshots_per_chain=site_draws,
        n_source_sites_without_degree_fixed_comparability=dict(n_rejected),
        panels=panels,
        caveats=dict(
            plant_positive_visitor_species_degrees_preserved=True,
            each_insect_species_positive_host_degree_preserved=True,
            species_to_original_guild_fixed=True,
            only_positive_observed_source_edges_are_randomized=True,
            publications_bootstrapped_not_network_edges=True,
            site_study_year_month_required=True,
            no_source_absences_or_effective_pollination_inferred=True,
            network_null_not_actual_pollinator_preference=True,
            original_plant_abundance_and_visit_intensity_unknown=True,
            original_4505_islands_parent_not_visitor_study_sites=True,
            European_visitors_not_tropical_island_guilds=True,
            no_original_site_study_accepted_plant_or_insect_names_exported=True,
            four_profile_study_cluster_multiplicity_only_not_full_prior_search=True,
            publication_leave_one_out_is_influence_not_cross_dataset_replication=True))
