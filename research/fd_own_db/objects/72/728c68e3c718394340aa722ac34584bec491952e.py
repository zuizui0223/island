"""Analytic null for actual same-accepted-species two-trait flower coupling.

For each fixed plant roster, independently exchange trait B labels WITHIN its
documented botanical genera, preserving each genus' number of source-recorded
trait-A and trait-B states, carrier n and globally observed margins. No plants,
flowers, or morphological values are imputed.

This is a compositional explanation / negative control, NOT pollinator
effectiveness, phylogenetic sister/ancestor reconstruction, or evolution.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
import numpy as np
import pandas as pd

from .real_joint_syndrome_access import JOINT_SPECS

PRIORITY=("actinomorphic_deep","zygomorphic_deep","tubular_bell_deep")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
NAMES=("island","geo5","climate10")
MIN_ISLAND=5
MIN_SOURCE=10
MIN_ISLANDS=30
MIN_BLOCKS=10


def genus_of(sp: str) -> str:
    t=str(sp).split()
    return t[0] if len(t)>=2 else "single_original_accepted_name:"+str(sp)


def one_roster_analytic_genus_pair(roster: set[str], profile: dict)->dict:
    """Exact conditional E(C), SD(C), and residual J-E(J), no random shuffles.

    C = P(A and B) - P(A)P(B). The source reference preserves within each
    genus species count plus marginal A and B frequencies, and computes the
    expectation under all equiprobable within-genus B-label assignments.
    """
    known=roster & profile["known"]
    n=len(known)
    base=dict(
        n_direct=n,n_A=0,n_B=0,n_joint=0,
        joint_prevalence=np.nan,coupling=np.nan,
        genus_null_expected_joint_prevalence=np.nan,
        genus_null_expected_coupling=np.nan,
        excess_coupling_beyond_genus_null=np.nan,
        genus_null_conditional_sd=np.nan,
        n_recorded_genera=0,n_genus_with_both_marginals_variable=0,
        n_carriers_in_variable_genera=0,
    )
    if n==0:return base
    bygenus=defaultdict(lambda:[0,0,0,0])
    for sp in known:
        a=int(sp in profile["A"]);b=int(sp in profile["B"])
        g=bygenus[genus_of(sp)]
        g[0]+=1
        g[1]+=a
        g[2]+=b
        g[3]+=a*b
    A=sum(q[1] for q in bygenus.values())
    B=sum(q[2] for q in bygenus.values())
    J=sum(q[3] for q in bygenus.values())
    null_j=0.
    null_var=0.
    variable=0
    carriers=0
    for k,a,b,j in bygenus.values():
        null_j += a*b/k
        if k>1:
            null_var += a*(b/k)*(1-b/k)*(k-a)/(k-1)
            if 0<a<k and 0<b<k:
                variable+=1
                carriers+=k
    p_joint=J/n
    fixed_marginal=A*B/(n*n)
    expected_j=null_j/n
    base.update(
        n_A=A,n_B=B,n_joint=J,
        joint_prevalence=float(p_joint),
        coupling=float(p_joint-fixed_marginal),
        genus_null_expected_joint_prevalence=float(expected_j),
        genus_null_expected_coupling=float(expected_j-fixed_marginal),
        excess_coupling_beyond_genus_null=float((J-null_j)/n),
        genus_null_conditional_sd=float(np.sqrt(max(0.,null_var))/n),
        n_recorded_genera=len(bygenus),
        n_genus_with_both_marginals_variable=variable,
        n_carriers_in_variable_genera=carriers,
    )
    if abs(
        base["coupling"]-base["genus_null_expected_coupling"]-
        base["excess_coupling_beyond_genus_null"]
    )>1e-11:
        raise AssertionError("analytical within-genus coupling partition error")
    return base


def _valid_source_rows(frame:pd.DataFrame)->pd.Series:
    return (frame.island_n_direct.ge(MIN_ISLAND)&
            frame.geo5_n_direct.ge(MIN_SOURCE)&
            frame.climate10_n_direct.ge(MIN_SOURCE)&
            frame.island_coupling.notna()&
            frame.geo5_coupling.notna()&
            frame.climate10_coupling.notna())


def _group_summary(q:pd.DataFrame,*,group:str,seed:str,boots:int)->dict:
    clusters=q.dropna(subset=[group]).groupby(group,sort=True)
    out=dict(n_islands=int(len(q)),n_groups=int(clusters.ngroups))
    if len(q)<MIN_ISLANDS or clusters.ngroups<MIN_BLOCKS:
        out["status"]="NOT_IDENTIFIABLE_GEO_OR_SOURCE_REPLICATION"
        return out
    columns=[
        f"{mode}__{metric}"
        for mode in ("geo5","climate10")
        for metric in (
            "source_relative_original_coupling",
            "source_relative_genus_expected_coupling",
            "source_relative_beyond_genus_null",
        )
    ]
    X=clusters[columns].mean().to_numpy(float)
    avg=X.mean(axis=0)
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    draws=X[rng.integers(len(X),size=(boots,len(X)))].mean(1)
    for i,col in enumerate(columns):
        out[col+"_group_equal_mean"]=float(avg[i])
        out[col+"_group_bootstrap95"]=[
            float(z) for z in np.quantile(draws[:,i],[.025,.975])
        ]
    out["status"]="DESCRIPTIVE_GENUS_CONDITIONAL_COUPLING_NOT_CAUSAL"
    return out


def summarize_original_genus_joint_coupling(
    frame:pd.DataFrame,*,population:int=4505,boots:int=999
)->dict:
    necessary={"island_id","regime","scope","sector",
               "spatial_block","source_signature",
               "island_total","geo5_total","climate10_total"}
    metrics=(
        "n_direct","coupling","genus_null_expected_coupling",
        "excess_coupling_beyond_genus_null",
        "genus_null_conditional_sd",
        "n_genus_with_both_marginals_variable",
        "n_carriers_in_variable_genera",
    )
    for mode in NAMES:
        necessary.update(f"{mode}_{z}" for z in metrics)
    if missing:=necessary-set(frame):
        raise ValueError(f"missing real-pair genus source frame {sorted(missing)}")
    if frame.duplicated(["island_id","scope","sector"]).any():
        raise ValueError("duplicate original island×scope×real pair")
    if frame.island_id.nunique()!=population or (
        len(frame)!=population*len(SCOPES)*len(PRIORITY)):
        raise ValueError("original GIFT study frame narrowed")
    if not frame.sector.isin(PRIORITY).all() or (
        not frame.scope.isin(SCOPES).all()) or (
        not frame.regime.isin(REGIMES).all()):
        raise ValueError("unexpected floral pair or original status/regime")
    if boots<199:
        raise ValueError("insufficient fixed original geographic block draws")
    for mode in NAMES:
        if (frame[f"{mode}_n_direct"]>frame[f"{mode}_total"]).any():
            raise ValueError("observed Direct pair species exceed original source roster")
        if (
            frame[f"{mode}_coupling"].notna() !=
            frame[f"{mode}_n_direct"].gt(0)
        ).any():
            raise ValueError("missing trait or source category improperly treated as zero")
        if (
            frame[f"{mode}_genus_null_expected_coupling"].notna() !=
            frame[f"{mode}_n_direct"].gt(0)
        ).any():
            raise ValueError("genus null missingness inconsistent")
        defined=frame[f"{mode}_n_direct"].gt(0)
        if defined.any() and not np.allclose(
            frame.loc[defined,f"{mode}_coupling"]-
            frame.loc[defined,f"{mode}_genus_null_expected_coupling"],
            frame.loc[defined,f"{mode}_excess_coupling_beyond_genus_null"],
            rtol=0,atol=1e-11,
        ):
            raise AssertionError("per-island genus conditional coupling not additive")

    result=[]
    for regime in REGIMES:
        for scope in SCOPES:
            for sector in PRIORITY:
                p=frame.loc[
                    frame.regime.eq(regime)&frame.scope.eq(scope)&
                    frame.sector.eq(sector)
                ].copy()
                q=p.loc[_valid_source_rows(p)].copy()
                for mode in ("geo5","climate10"):
                    q[f"{mode}__source_relative_original_coupling"]=(
                        q.island_coupling-q[f"{mode}_coupling"])
                    q[f"{mode}__source_relative_genus_expected_coupling"]=(
                        q.island_genus_null_expected_coupling-
                        q[f"{mode}_genus_null_expected_coupling"])
                    q[f"{mode}__source_relative_beyond_genus_null"]=(
                        q.island_excess_coupling_beyond_genus_null-
                        q[f"{mode}_excess_coupling_beyond_genus_null"])
                    if len(q) and not np.allclose(
                        q[f"{mode}__source_relative_original_coupling"],
                        q[f"{mode}__source_relative_genus_expected_coupling"]+
                        q[f"{mode}__source_relative_beyond_genus_null"],
                        rtol=0,atol=1e-11):
                        raise AssertionError("source-relative coupling sum differs")
                if len(q):
                    fi=q.island_n_direct/q.island_total.replace(0,np.nan)
                    fg=q.geo5_n_direct/q.geo5_total.replace(0,np.nan)
                    fc=q.climate10_n_direct/q.climate10_total.replace(0,np.nan)
                    balanced=fi.sub(fg).abs().le(.10)&fi.sub(fc).abs().le(.10)
                    informative=(
                        q.island_n_genus_with_both_marginals_variable.ge(1)&
                        q.geo5_n_genus_with_both_marginals_variable.ge(1)&
                        q.climate10_n_genus_with_both_marginals_variable.ge(1)&
                        q.island_n_carriers_in_variable_genera.ge(3)&
                        q.geo5_n_carriers_in_variable_genera.ge(3)&
                        q.climate10_n_carriers_in_variable_genera.ge(3)
                    )
                else:
                    balanced=pd.Series(False,index=q.index)
                    informative=pd.Series(False,index=q.index)
                # The genuinely informative part of the genus null requires
                # variation in BOTH focal traits inside a genus. Report
                # opportunity separately on island and both source lists:
                # a near-zero residual can otherwise be mechanically fixed.
                inform_count={
                    mode: int((
                        q[f"{mode}_n_genus_with_both_marginals_variable"].ge(1)&
                        q[f"{mode}_n_carriers_in_variable_genera"].ge(3)
                    ).sum()) for mode in NAMES
                }
                panel=dict(
                    region=regime,scope=scope,
                    genuine_joint_sector=sector,
                    n_original_region_flora=len(p),
                    n_both_real_source_mode_qualified_islands=int(len(q)),
                    n_with_both_source_and_island_informative_genera=int(informative.sum()),
                    n_islands_with_informative_genus_on_island=inform_count["island"],
                    n_islands_with_informative_genus_in_geo5_source=inform_count["geo5"],
                    n_islands_with_informative_genus_in_climate10_source=inform_count["climate10"],
                    n_direct_documentation_gap_both_le10pp=int(balanced.sum()),
                    n_both_informative_and_documentation_balanced=int((balanced&informative).sum()),
                    # Null cannot identify a within-genus effect in monomorphic
                    # or singleton genera, so report that opportunity explicitly.
                    original_genus_conditional_null_not_plant_adaptation=True,
                )
                for name,sub in (
                    ("all_direct_qualified",q),
                    ("trait_documentation_balanced",q.loc[balanced]),
                    ("informative_genera_both_sources",q.loc[informative]),
                    ("both_balanced_and_informative",q.loc[balanced&informative]),
                ):
                    for g in ("spatial_block","source_signature"):
                        panel[name+"__"+g]=_group_summary(
                            sub,group=g,
                            seed=f"{regime}|{scope}|{sector}|{name}|{g}",
                            boots=boots
                        )
                result.append(panel)
    return dict(
        status="ORIGINAL_4505_REAL_JOINT_GENUS_CONDITIONAL_COUPLING",
        n_original_GIFT_botanical_islands=population,
        n_exact_same_species_morphological_pairs=len(PRIORITY),
        n_region_scope_structural_pair_panels=len(result),
        n_original_private_island_scope_pair_rows=len(frame),
        panels=result,
        epistemic_limits=dict(
            true_same_species_two_trait_pairs_only=True,
            analytic_within_genus_fixed_margin_null_no_random_seeds=True,
            genus_null_different_from_genus_composition_effect_on_total_pair_share=True,
            no_ancestral_phylogeny_or_local_population_morphological_change=True,
            no_effective_bird_butterfly_moth_guild_observed_in_this_analysis=True,
            documentation_and_congener_opportunity_not_imputed=True,
            all_4505_islands_universe_kept=True,
            posthoc_after_tropical_symmetry_depth_result=True,
        ),
    )
