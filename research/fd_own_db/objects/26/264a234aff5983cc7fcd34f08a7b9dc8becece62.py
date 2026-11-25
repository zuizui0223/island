"""Exact botanical genus partition of TRUE SAME-SPECIES floral trait coupling.

For each observed island and independently constructed potential mainland
source species list, never multiply marginals from distinct plants. On the
identical direct two-axis accepted-species cohort:

J = P(A and B);
C = J - P(A)P(B);
E_genus[J] = SUM_g (n_g/N)(A_g/n_g)(B_g/n_g);
C_within = J - E_genus[J]; C_between = E_genus[J] - P(A)P(B).

Thus C = C_within + C_between EXACTLY.
The variance of J under independent finite fixed-margin permutation of
B state WITHIN each original botanical genus is analytically available.

These are descriptive partition/counterfactual margins conditioned on
recorded genera and species-level phenotypes. NO phylogenetic controls
within genera, causal pollinator effects, historical lineage filtering,
or effective visitor-syndrome labels are inferred.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import numpy as np
import pandas as pd

from .real_joint_syndrome_access import (
    JOINT_SPECS, REGIMES, SCOPES, MIN_ISLAND, MIN_SOURCE,
)

MODES=("geo5","climate10")
MIN_REGION_ISLANDS=30
MIN_BLOCKS=10
MIN_FIRST_SOURCES=10


def _genus(name:str)->str:
    parts=str(name).split()
    return parts[0] if len(parts)>=2 else "invalid_no_binomial_"+str(name)


def partition_genus_coupling(roster:set[str],profile:dict)->dict:
    known=roster & set(profile["known"])
    if (not set(profile["joint"]) <= set(profile["A"]) &
        set(profile["B"]) or
        not (set(profile["A"])|set(profile["B"]))<=set(profile["known"])):
        raise ValueError("original same-accepted-plant phenotype ontology inconsistent")
    N=len(known)
    out=dict(n_direct=N,n_genera=0,n_variable_genera=0,
             n_plants_in_variable_genera=0)
    if not N:
        out.update({k:np.nan for k in (
            "joint","pA","pB","coupling","genus_independent_joint",
            "within_genus_coupling","between_genus_coupling",
            "finite_within_genus_joint_variance")})
        return out
    posA=known & set(profile["A"])
    posB=known & set(profile["B"])
    joint=known & set(profile["joint"])
    by_g=defaultdict(lambda:[0,0,0,0])
    for sp in known:
        row=by_g[_genus(sp)]
        row[0]+=1
        row[1]+=int(sp in posA)
        row[2]+=int(sp in posB)
        row[3]+=int(sp in joint)
    expected_count=0.
    variance_count=0.
    var_g=0
    nvar_sp=0
    for n,a,b,j in by_g.values():
        expected_count+=a*b/n
        if 0<a<n and 0<b<n:
            var_g+=1
            nvar_sp+=n
            # Hypergeometric overlap of two recorded binary traits with
            # EXACT genus-stratified finite accepted-plant margins.
            variance_count+=a*b*(n-a)*(n-b)/(n*n*(n-1))
    pa=len(posA)/N;pb=len(posB)/N
    j=len(joint)/N
    expected=expected_count/N
    between=expected-pa*pb
    within=j-expected
    coupling=j-pa*pb
    if abs(coupling-between-within)>1e-11:
        raise AssertionError("exact genus within/between flower coupling partition failed")
    if variance_count<0 or variance_count/N**2<0:
        raise AssertionError("invalid fixed-flower-genus finite variance")
    return dict(n_direct=N,n_genera=len(by_g),
                n_variable_genera=var_g,
                n_plants_in_variable_genera=nvar_sp,
                joint=j,pA=pa,pB=pb,coupling=coupling,
                genus_independent_joint=expected,
                between_genus_coupling=between,
                within_genus_coupling=within,
                finite_within_genus_joint_variance=variance_count/N**2)


def _genera_with_real_within_variation(roster:set[str],profile:dict)->frozenset:
    """Private genus identifiers; source-safe output contains COUNTS ONLY."""
    known=set(roster)&set(profile["known"])
    A=set(profile["A"]);B=set(profile["B"])
    by_g=defaultdict(lambda:[0,0,0])
    for sp in known:
        v=by_g[_genus(sp)]
        v[0]+=1;v[1]+=int(sp in A);v[2]+=int(sp in B)
    return frozenset(g for g,(n,a,b) in by_g.items()
                     if n>1 and 0<a<n and 0<b<n)


def partition_island_source(*,full_flora:set[str],
                            scope_rosters:dict[str,set[str]],
                            source_rosters:dict[str,set[str]],
                            pairs:dict)->list[dict]:
    """Return 33 rows for one original island, no missing-phenotype imputation."""
    if (not set(SCOPES)<=set(scope_rosters) or
        not set(MODES)<=set(source_rosters) or
        not set(JOINT_SPECS)<=set(pairs)):
        raise ValueError("original flora/source/morphological axes missing")
    if any(not set(scope_rosters[s])<=full_flora for s in SCOPES):
        raise ValueError("GIFT origin status source is outside focal flora")
    if any(set(source_rosters[m])&full_flora for m in MODES):
        raise ValueError("mainland donor accepted plant leaked into focal GIFT flora")
    src={}
    source_var={}
    for p in JOINT_SPECS:
        profile=pairs[p]
        src[p]={m:partition_genus_coupling(source_rosters[m],profile)
                for m in MODES}
        source_var[p]={m:_genera_with_real_within_variation(source_rosters[m],profile)
                       for m in MODES}
    rows=[]
    for scope in SCOPES:
        for name,p in JOINT_SPECS.items():
            island=partition_genus_coupling(scope_rosters[scope],pairs[name])
            island_var=_genera_with_real_within_variation(scope_rosters[scope],
                                                              pairs[name])
            common_var=(island_var & source_var[name]["geo5"] &
                        source_var[name]["climate10"])
            z=dict(scope=scope,sector=name,
                   island_total=len(scope_rosters[scope]),
                   n_shared_informative_genera_both_source_modes=len(common_var))
            z.update({f"island__{k}":v for k,v in island.items()})
            for mode in MODES:
                m=src[name][mode]
                z[f"{mode}__total"]=len(source_rosters[mode])
                z.update({f"{mode}__{k}":v for k,v in m.items()})
                for metric in ("joint","coupling","between_genus_coupling",
                               "within_genus_coupling"):
                    x,y=island[metric],m[metric]
                    z[f"{mode}__delta_{metric}"]=x-y if np.isfinite(x+y) else np.nan
                if np.isfinite(z[f"{mode}__delta_coupling"]):
                    diff=z[f"{mode}__delta_coupling"]
                    parts=(z[f"{mode}__delta_between_genus_coupling"]+
                           z[f"{mode}__delta_within_genus_coupling"])
                    if abs(diff-parts)>1e-11:
                        raise AssertionError("original source-relative flower coupling failed genus accounting")
            rows.append(z)
    if len(rows)!=len(SCOPES)*len(JOINT_SPECS):
        raise AssertionError("global GIFT status/structural panel missing")
    return rows


def aggregate_pair_genus_partition(frame:pd.DataFrame, *,
                                   population:int=4505,draws:int=999)->dict:
    """Source-block / first-GIFT-source bootstrap, fail-closed infovar cohort."""
    needed={"island_id","regime","spatial_block","first_source",
            "scope","sector","island__n_direct","island__n_genera",
            "island__n_variable_genera","island_total",
            "n_shared_informative_genera_both_source_modes"}
    for mode in MODES:
        needed|={f"{mode}__total",f"{mode}__n_direct",
                 f"{mode}__n_variable_genera"}
        needed|={f"{mode}__delta_{m}" for m in (
            "joint","coupling","between_genus_coupling","within_genus_coupling")}
    if missing:=needed-set(frame):
        raise ValueError("missing source genus within/between data "+repr(sorted(missing)))
    if (frame.island_id.nunique()!=population or
        len(frame)!=population*len(SCOPES)*len(JOINT_SPECS) or
        frame.duplicated(["island_id","scope","sector"]).any()):
        raise ValueError("original complete 4505 island population lost")
    if (not frame.regime.isin(REGIMES).all() or
        not frame.scope.isin(SCOPES).all() or
        not frame.sector.isin(JOINT_SPECS).all()):
        raise ValueError("unknown original GIFT region/trait/status")
    if draws<199:
        raise ValueError("not enough source-block draws")
    panels=[]
    for region in REGIMES:
        for scope in SCOPES:
            for sector in JOINT_SPECS:
                pool=frame.loc[
                    frame.regime.eq(region)&frame.scope.eq(scope)&
                    frame.sector.eq(sector)]
                valid=pool.island__n_direct.ge(MIN_ISLAND)
                for m in MODES:
                    valid &=pool[f"{m}__n_direct"].ge(MIN_SOURCE)
                    valid &=np.isfinite(pool[f"{m}__delta_joint"])
                q=pool.loc[valid].copy()
                # Strictly require actual within-genus species trait variation
                # on BOTH sides under BOTH potential mainland source definitions.
                # Otherwise within-genus = 0 is forced by the taxon composition.
                informative=(q.island__n_variable_genera.gt(0)&
                    q.geo5__n_variable_genera.gt(0)&
                    q.climate10__n_variable_genera.gt(0))
                common_informative=q.n_shared_informative_genera_both_source_modes.gt(0)
                # This is a TRUE source/genus-matched follow-up; the previous
                # broad "variation on either side" cohort may have completely
                # different botanical genera on island and mainland.
                # No missing phenotype value is replaced by a source zero.
                # Additional documentation-imbalance check, scored separately.
                balanced=(q.island__n_direct.div(q.island_total.replace(0,np.nan))
                    .sub(q.geo5__n_direct.div(q.geo5__total.replace(0,np.nan)))
                    .abs().le(.10)&
                    q.island__n_direct.div(q.island_total.replace(0,np.nan))
                    .sub(q.climate10__n_direct.div(q.climate10__total.replace(0,np.nan)))
                    .abs().le(.10))
                z=dict(region=region,scope=scope,sector=sector,
                       n_original_region_islands=int(len(pool)),
                       n_both_mainland_source_supported=int(len(q)),
                       n_islands_with_true_variation_on_all_three_rosters=int(informative.sum()),
                       n_islands_with_same_variable_genus_on_island_and_both_sources=int(common_informative.sum()),
                       n_annotation_rate_le10pp=int(balanced.sum()),
                       n_informative_and_annotation_rate_le10pp=int((informative&balanced).sum()))
                cols=[f"{m}__delta_{metric}" for m in MODES for metric in (
                    "joint","coupling","between_genus_coupling",
                    "within_genus_coupling")]
                for mask,sub in (("all",q),("genuine_within_genus_var",q.loc[informative]),
                                 ("genus_var_and_doc_le10pp",q.loc[informative&balanced]),
                                 ("same_variable_genus_across_island_and_sources",
                                   q.loc[common_informative])):
                    for group,threshold in (("spatial_block",MIN_BLOCKS),
                                            ("first_source",MIN_FIRST_SOURCES)):
                        grouped=sub.dropna(subset=[group]).groupby(group,sort=True)
                        n_groups=grouped.ngroups
                        item=dict(n_islands=int(len(sub)),n_groups=int(n_groups))
                        if len(sub)<MIN_REGION_ISLANDS or n_groups<threshold:
                            item["status"]=(
                                "NOT_IDENTIFIABLE_TRUE_WITHIN_GENUS_SOURCE_SUPPORT"
                                if mask!="all" else "NOT_IDENTIFIABLE_DIRECT_PAIR_SOURCE_SUPPORT")
                        else:
                            a=grouped[cols].mean().to_numpy(float)
                            mean=a.mean(axis=0)
                            seed=int.from_bytes(sha256(
                                f"{region}|{scope}|{sector}|{mask}|{group}".encode()
                            ).digest()[:8],"big")
                            rng=np.random.default_rng(seed)
                            boot=a[rng.integers(len(a),size=(draws,len(a)))].mean(axis=1)
                            for j,key in enumerate(cols):
                                item[key+"__equal_group_mean"]=float(mean[j])
                                item[key+"__bootstrap95"]=[
                                    float(v) for v in np.quantile(boot[:,j],[.025,.975])]
                            for m in MODES:
                                # exact partition must also survive equal cluster averaging
                                lhs=item[f"{m}__delta_coupling__equal_group_mean"]
                                rhs=(item[f"{m}__delta_between_genus_coupling__equal_group_mean"]+
                                     item[f"{m}__delta_within_genus_coupling__equal_group_mean"])
                                if abs(lhs-rhs)>1e-11:
                                    raise AssertionError("original mean genus partition drifted")
                            item["status"]="DESCRIPTIVE_EXACT_GENUS_WITHIN_BETWEEN_PAIR_PARTITION"
                        z[mask+"__"+group]=item
                panels.append(z)
    return dict(
        status="ORIGINAL_GIFT_4505_EXACT_WITHIN_VS_BETWEEN_GENUS_FLOWER_PAIRING",
        original_botanical_island_population=population,
        n_island_scope_pair_records=len(frame),
        n_distinct_true_same_accepted_plant_pairs=len(JOINT_SPECS),
        n_region_scope_pair_panels=len(panels),
        panels=panels,
        limits=dict(
            exact_two_flower_traits_on_same_accepted_plant=True,
            within_genus_expectation_conditions_on_real_species_fixed_margins=True,
            between_genus_component_is_flora_compositional_not_evolution=True,
            strict_variation_requires_original_trait_variation_on_island_and_both_sources=True,
            matched_same_variable_genus_required_for_strongest_contrast=True,
            previous_simple_full_flora_10pp_documentation_gate_is_not_proven_missing_at_random=True,
            high_missingness_weakens_identification_not_zero_effect=True,
            potential_mainland_sources_not_actual_historical_propagules=True,
            genus_not_phylogeny_and_original_source_provider_not_independent=True,
            no_visitor_guild_or_effective_pollination_inferred=True,
            all_4505_island_parent_floras_preserved=True))
