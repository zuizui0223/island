"""Out-of-focal-island accepted-name congener prediction of REAL floral pairs.

For each botanical island, do NOT train on any accepted species in that
island's entire flora, even if it has a source GIFT regional record. Predict
original documented same-accepted-species joint flower states from different
source-only congeneric species. Compare three models on the SAME supported
island plants: global source frequency, shared-genus donor-pooled frequency,
and island-genus-composition weighted congener means.

No causal pollinator guild, local phenotypic adaptation, actual ancestral
population, arrival filter, or genuine unseen-genus prediction is identified.
"""
from __future__ import annotations

from collections import defaultdict
from hashlib import sha256

import numpy as np
import pandas as pd

from .real_joint_syndrome_access import JOINT_SPECS

SECTORS=("actinomorphic_deep","zygomorphic_deep","tubular_bell_deep")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
MODES=("geo5","climate10")
DONORS=("any_congener","at_least_two_congeners")
MIN_TARGET=5
MIN_GENERA=2
MIN_DONORS=10
MIN_ISLANDS=30
MIN_BLOCKS=10


def _genus(sp:str)->str:
    bits=str(sp).split()
    return bits[0] if len(bits)>=2 else "single_source_name:"+str(sp)


def _carrier(sp:str,profile:dict)->tuple[int,int,int]:
    a=int(sp in profile["A"])
    b=int(sp in profile["B"])
    return (a,b,a*b)


def _build_source_donors(source:set[str],profile:dict)->dict:
    q=defaultdict(list)
    for sp in source & profile["known"]:
        q[_genus(sp)].append(_carrier(sp,profile))
    return q


def _prediction_on_fixed_target(
    observed:set[str],source_by_genus:dict,all_source_states:list,
    profile:dict, *,min_donor_per_genus:int,
)->dict:
    """Same target species for all three baseline models within each mode."""
    eligible={g for g,x in source_by_genus.items()
              if len(x)>=min_donor_per_genus}
    target_by_genus=defaultdict(list)
    for sp in observed:
        genus=_genus(sp)
        if genus in eligible:
            target_by_genus[genus].append(sp)
    n=sum(len(q) for q in target_by_genus.values())
    base=dict(n_target_direct=n,n_target_genera=len(target_by_genus),
              n_source_matched_direct=sum(len(source_by_genus[g])
                                          for g in target_by_genus),
              n_source_matched_total=None,
              joint_observed=np.nan,coupling_observed=np.nan,
              joint_genus_predicted=np.nan,coupling_genus_predicted=np.nan,
              joint_factorized_predicted=np.nan,coupling_factorized_predicted=np.nan,
               n_source_genera_with_both_traits_variable=0,
               n_target_carriers_on_variable_source_genera=0,
               joint_pooled_predicted=np.nan,coupling_pooled_predicted=np.nan,
              joint_global_predicted=np.nan,coupling_global_predicted=np.nan,
              brier_genus=np.nan,brier_factorized=np.nan,
               brier_pooled=np.nan,brier_global=np.nan,
              n_target_positive=0,
              target_species_per_source_genus_support=None)
    if not n:
        return base
    if not all_source_states:
        raise ValueError("cannot predict from an empty source-only donor pool")
    source_by_target={g:np.asarray(source_by_genus[g],float)
                      for g in target_by_genus}
    merged=np.concatenate(list(source_by_target.values()),axis=0)
    global_state=np.asarray(all_source_states,float)
    pglobal=np.mean(global_state,axis=0)
    ppool=np.mean(merged,axis=0)
    observed_state=np.asarray([
        _carrier(sp,profile) for g in target_by_genus
        for sp in target_by_genus[g]
    ],float)
    target_means=np.mean(observed_state,axis=0)
    pred=np.zeros(3)
    brier=0.
    joint_factorized=0.
    brier_factorized=0.
    informative_genera=0
    informative_carriers=0
    for g,sp in target_by_genus.items():
        n_g=len(sp)
        p=np.mean(source_by_target[g],axis=0)
        pred += p*n_g/n
        factorized=float(p[0]*p[1])
        joint_factorized+=factorized*n_g/n
        if 0<p[0]<1 and 0<p[1]<1:
            informative_genera+=1
            informative_carriers+=n_g
        positive=sum(int(s in profile["joint"]) for s in sp)
        brier += (positive*(1.-p[2])**2+
                  (n_g-positive)*p[2]**2)/n
        brier_factorized += (
            positive*(1.-factorized)**2+
            (n_g-positive)*factorized**2
        )/n
    obsJ=target_means[2]
    def c(x):
        return float(x[2]-x[0]*x[1])
    base.update(
        joint_observed=float(obsJ),
        coupling_observed=c(target_means),
        joint_genus_predicted=float(pred[2]),
        coupling_genus_predicted=c(pred),
        joint_factorized_predicted=float(joint_factorized),
        coupling_factorized_predicted=float(
            joint_factorized-pred[0]*pred[1]),
        n_source_genera_with_both_traits_variable=informative_genera,
        n_target_carriers_on_variable_source_genera=informative_carriers,
        joint_pooled_predicted=float(ppool[2]),
        coupling_pooled_predicted=c(ppool),
        joint_global_predicted=float(pglobal[2]),
        coupling_global_predicted=c(pglobal),
        brier_genus=float(brier),
        brier_factorized=float(brier_factorized),
        brier_pooled=float(obsJ*(1-ppool[2])**2+(1-obsJ)*ppool[2]**2),
        brier_global=float(obsJ*(1-pglobal[2])**2+(1-obsJ)*pglobal[2]**2),
        n_target_positive=int(target_means[2]*n+1e-9),
    )
    if not (0<=base["brier_genus"]<=1 and
            0<=base["brier_factorized"]<=1 and
            0<=base["brier_pooled"]<=1 and
            0<=base["brier_global"]<=1):
        raise AssertionError("invalid species name held-out source-only predictions")
    return base


def predict_joint_by_source_scopes(
    *,full_island:set[str],scope_rosters:dict[str,set[str]],
    geo5_source:set[str],climate10_source:set[str],profile:dict,
)->dict:
    """Reuse original name-exclusive donors for the three nested status scopes."""
    if any(not s<=full_island for s in scope_rosters.values()):
        raise ValueError("focal plant status roster exceeds original island flora")
    if (geo5_source|climate10_source)&full_island:
        raise ValueError("focal accepted island species leaked into source model")
    if not scope_rosters or not profile["known"]:
        raise ValueError("empty original status/flower pair source specification")
    source={
        "geo5":_build_source_donors(geo5_source,profile),
        "climate10":_build_source_donors(climate10_source,profile),
    }
    global_pool={
        m:[z for x in source[m].values() for z in x]
        for m in MODES
    }
    scores={}
    for scope,roster in scope_rosters.items():
        observed=roster & profile["known"]
        # To avoid a moving same-island target for geo5 vs climate10,
        # only include plant names with congeneric Direct source donors in
        # BOTH modes, at the required minimum donor class.
        scores[scope]={}
        for donor_key,threshold in (("any_congener",1),
                                    ("at_least_two_congeners",2)):
            both_g={
                g for g,x in source["geo5"].items() if len(x)>=threshold
            }&{
                g for g,x in source["climate10"].items() if len(x)>=threshold
            }
            fixed={sp for sp in observed if _genus(sp) in both_g}
            item={}
            for mode in MODES:
                result=_prediction_on_fixed_target(
                    fixed,source[mode],global_pool[mode],profile,
                    min_donor_per_genus=threshold,
                )
                # Per-source trait documentation fraction in exactly the
                # focal GENUS set, including unscored distinct species, is
                # handled by separate source roster counts in the runner.
                item[mode]=result
            if item["geo5"]["n_target_direct"]!=item["climate10"]["n_target_direct"]:
                raise AssertionError("source mode changed original focal species")
            scores[scope][donor_key]=item
    return scores


def _columns(prefix:str):
    return [
        f"{prefix}__bias_genus",f"{prefix}__bias_pooled",
        f"{prefix}__bias_global",f"{prefix}__brier_gain_vs_pooled",
        f"{prefix}__brier_gain_joint_vs_factorized",
        f"{prefix}__bias_factorized",
        f"{prefix}__brier_gain_vs_global",
        f"{prefix}__coupling_residual_congener",
    ]


def aggregate_source_congener_prediction(
    frame:pd.DataFrame,*,population:int=4505,draws:int=999
)->dict:
    keys={"island_id","regime","scope","sector","donor_gate",
          "spatial_block","source_signature","n_target_direct",
          "n_target_genera", "island_direct_coverage_supported",
          "geo5_direct_coverage_supported",
          "climate10_direct_coverage_supported"}
    for mode in MODES:
        keys |= {
            f"{mode}__n_source_matched_direct",
            f"{mode}__n_source_matched_total",
            f"{mode}__n_target_positive",
            f"{mode}__joint_observed",
            f"{mode}__joint_genus_predicted",
            f"{mode}__joint_factorized_predicted",
            f"{mode}__n_source_genera_with_both_traits_variable",
            f"{mode}__n_target_carriers_on_variable_source_genera",
            f"{mode}__joint_pooled_predicted",
            f"{mode}__joint_global_predicted",
            f"{mode}__coupling_observed",
            f"{mode}__coupling_genus_predicted",
            f"{mode}__brier_genus",f"{mode}__brier_factorized",
            f"{mode}__brier_pooled",
            f"{mode}__brier_global",
        }
    if missing:=keys-set(frame):
        raise ValueError(f"missing original source congener result columns {sorted(missing)}")
    if frame.duplicated(["island_id","scope","sector","donor_gate"]).any():
        raise ValueError("duplicate original island/status/real pair/donor gate")
    if frame.island_id.nunique()!=population or (
        len(frame)!=population*len(SCOPES)*len(SECTORS)*len(DONORS)
    ):
        raise ValueError("full 4505-island source/flower frame contracted")
    if not frame.scope.isin(SCOPES).all() or not frame.sector.isin(SECTORS).all() or (
        not frame.donor_gate.isin(DONORS).all()
    ) or not frame.regime.isin(REGIMES).all():
        raise ValueError("unexpected native statuses, flower pair or source donor mode")
    if draws<199:
        raise ValueError("too few physical-source cluster draws")
    for mode in MODES:
        if not (
            (frame[f"{mode}__joint_observed"].notna() ==
             frame.n_target_direct.gt(0)).all()
        ):
            raise ValueError("missing actual Direct joint pair not preserved")
        if (frame[f"{mode}__n_source_matched_direct"]<0).any():
            raise ValueError("negative original congener donor count")
        defined=frame.n_target_direct.gt(0)
        for kind in ("genus","factorized","pooled","global"):
            b=frame.loc[defined,f"{mode}__brier_{kind}"]
            if b.isna().any() or not b.between(0,1).all():
                raise ValueError("invalid source-only Brier metric")
        if (frame.loc[defined,f"{mode}__n_target_positive"]>
            frame.loc[defined,"n_target_direct"]).any():
            raise ValueError("target marked positive beyond known Direct species")

    panels=[]
    for reg in REGIMES:
        for scope in SCOPES:
            for sector in SECTORS:
                for donor_gate in DONORS:
                    x=frame.loc[
                        frame.regime.eq(reg)&frame.scope.eq(scope)&
                        frame.sector.eq(sector)&frame.donor_gate.eq(donor_gate)
                    ].copy()
                    valid=x.loc[
                        x.n_target_direct.ge(MIN_TARGET)&
                        x.n_target_genera.ge(MIN_GENERA)&
                        x.geo5__n_source_matched_direct.ge(MIN_DONORS)&
                        x.climate10__n_source_matched_direct.ge(MIN_DONORS)
                    ].copy()
                    if len(valid):
                        for mode in MODES:
                            valid[f"{mode}__bias_genus"]=(
                                valid[f"{mode}__joint_observed"]-
                                valid[f"{mode}__joint_genus_predicted"])
                            valid[f"{mode}__bias_factorized"]=(
                                valid[f"{mode}__joint_observed"]-
                                valid[f"{mode}__joint_factorized_predicted"])
                            valid[f"{mode}__brier_gain_joint_vs_factorized"]=(
                                valid[f"{mode}__brier_factorized"]-
                                valid[f"{mode}__brier_genus"])
                            valid[f"{mode}__bias_pooled"]=(
                                valid[f"{mode}__joint_observed"]-
                                valid[f"{mode}__joint_pooled_predicted"])
                            valid[f"{mode}__bias_global"]=(
                                valid[f"{mode}__joint_observed"]-
                                valid[f"{mode}__joint_global_predicted"])
                            valid[f"{mode}__brier_gain_vs_pooled"]=(
                                valid[f"{mode}__brier_pooled"]-
                                valid[f"{mode}__brier_genus"])
                            valid[f"{mode}__brier_gain_vs_global"]=(
                                valid[f"{mode}__brier_global"]-
                                valid[f"{mode}__brier_genus"])
                            valid[f"{mode}__coupling_residual_congener"]=(
                                valid[f"{mode}__coupling_observed"]-
                                valid[f"{mode}__coupling_genus_predicted"])
                        # Balanced observed Direct known trait proportion on
                        # island among the supported genera vs counterpart
                        # source-only congener references is separately audited.
                        coverage=(
                            (valid.island_direct_coverage_supported-
                             valid.geo5_direct_coverage_supported).abs().le(.10)&
                            (valid.island_direct_coverage_supported-
                             valid.climate10_direct_coverage_supported).abs().le(.10)
                        )
                    else:
                        coverage=pd.Series(False,index=valid.index)
                    item=dict(
                        region=reg,scope=scope,
                        genuine_same_species_joint_sector=sector,
                        congener_donor_requirement=donor_gate,
                        n_original_region_flora=int(len(x)),
                        n_same_islands_with_both_source_congeners=int(len(valid)),
                        n_coverage_balanced_10pp=int(coverage.sum()),
                        median_supported_original_island_Direct_carriers=(
                            float(valid.n_target_direct.median()) if len(valid) else None
                        ),
                        median_supported_original_island_genera=(
                            float(valid.n_target_genera.median()) if len(valid) else None
                        ),
                        n_islands_with_source_genus_variable_in_both_axes_and_modes=int((
                            valid.geo5__n_source_genera_with_both_traits_variable.ge(1)&
                            valid.climate10__n_source_genera_with_both_traits_variable.ge(1)
                        ).sum()),
                    )
                    informative=(
                        valid.geo5__n_source_genera_with_both_traits_variable.ge(1)&
                        valid.climate10__n_source_genera_with_both_traits_variable.ge(1)
                    )
                    item["n_islands_both_informative_and_coverage10"]=int(
                        (informative&coverage).sum()
                    )
                    for name,sub in (
                        ("all",valid),("coverage10",valid.loc[coverage]),
                        ("source_trait_variation_both",valid.loc[informative]),
                        ("source_trait_variation_and_coverage10",
                         valid.loc[informative&coverage]),
                    ):
                        for cluster in ("spatial_block","source_signature"):
                            group=sub.dropna(subset=[cluster]).groupby(cluster,sort=True)
                            a=dict(n_original_islands=int(len(sub)),
                                   n_independent_source_or_spatial_groups=int(group.ngroups))
                            if len(sub)<MIN_ISLANDS or group.ngroups<MIN_BLOCKS:
                                a["status"]="NOT_IDENTIFIABLE_GEOGRAPHICAL_SOURCE_SUPPORT"
                            else:
                                names=_columns("geo5")+_columns("climate10")
                                vals=group[names].mean().to_numpy(float)
                                m=vals.mean(0)
                                seed=int.from_bytes(sha256(
                                    f"{reg}|{scope}|{sector}|{donor_gate}|{name}|{cluster}".encode()
                                ).digest()[:8],"big")
                                rng=np.random.default_rng(seed)
                                boots=vals[rng.integers(len(vals),size=(draws,len(vals)))].mean(1)
                                for k,key in enumerate(names):
                                    a[key+"__mean"]=float(m[k])
                                    a[key+"__bootstrap95"]=[
                                        float(q) for q in np.quantile(boots[:,k],[.025,.975])
                                    ]
                                a["status"]="SOURCE_ONLY_NAME_HELDOUT_CONGENERIC_PREDICTION"
                            item[f"{name}__{cluster}"]=a
                    panels.append(item)
    return dict(
        status="ORIGINAL_4505_SOURCE_ONLY_REAL_FLORAL_PAIR_CONGENER_PREDICTION",
        n_original_GIFT_botanical_flora_islands=population,
        n_original_island_scope_sector_donor_cells=len(frame),
        n_region_scope_sector_donor_panels=len(panels),
        panels=panels,
        limits=dict(
            focal_accepted_species_names_excluded_from_both_source_pools=True,
            source_congeners_are_different_accepted_species=True,
            original_island_genus_composition_known_not_predicted=True,
            all_methods_evaluated_on_exact_same_supported_island_plants=True,
            no_true_out_of_lineage_or_independent_trait_provider=True,
            no_true_effective_pollination_guild_measured=True,
            source_coverage_and_nativeness_nonrandomly_missing=True,
            prediction_increment_of_true_joint_over_independent_marginal_products_tested=True,
            within_source_genus_both_traits_variable_support_explicit=True,
            global_original_4505_denominator_retained=True,
        )
    )
