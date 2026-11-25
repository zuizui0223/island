"""External NHM tube-length calibration and source-only congener prediction.

All training species are GIFT native mainland names that occur on NONE of
4,505 original island checklists. Predictions for each island use different
accepted species with NHM corolla-tube length recorded in its two fixed GIFT
candidate source pools. Original ordinal island tube-depth values NEVER enter
source model fitting. This tests predictability, NOT local adaptation or a
bird/butterfly pollination syndrome.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache
from hashlib import sha256

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit

REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
TIERS=("any_valid_NHM_numeric_tube_length",
       "two_source_corpora_consistent_tube_length")
MODES=("geo5","climate10")
MIN_TARGET=5
MIN_GENERA=2
MIN_SOURCE_DONORS=10
MIN_ISLANDS=30
MIN_BLOCKS=10


def botanical_genus(name:str)->str:
    t=str(name).split()
    return t[0] if len(t)>=2 else "single_source_name:"+str(name)


def fit_source_only_length_calibrator(
    external_cm:dict[str,float],depth:dict[str,str],
    mainland_names:set[str],all_island_names:set[str],
)->dict:
    """Simple pre-fixed one-axis logistic with source-only accepted names.

    The model has log1p(corolla length cm) as its only informative predictor.
    Mild L2 slope regularization prevents complete separation. Strict source
    opportunity guard: no focal or any original island accepted names in fit.
    """
    keys=(set(external_cm)&set(depth)&mainland_names)-all_island_names
    if keys&all_island_names:
        raise ValueError("focal Direct depth information entered calibration")
    names=sorted(s for s in keys if depth[s] in ("shallow","intermediate","deep"))
    if not names:raise ValueError("independent external source calibration missing")
    x=np.array([np.log1p(float(external_cm[s])) for s in names],float)
    y=np.array([int(depth[s]=="deep") for s in names],float)
    if not np.isfinite(x).all() or (x<=0).any():
        raise ValueError("invalid original numeric NHM tube length")
    if len(names)<100 or min(sum(y),len(names)-sum(y))<30:
        raise ValueError("insufficient source-only deep/nondeep calibration cases")
    def objective(beta):
        eta=beta[0]+beta[1]*x
        loss=np.logaddexp(0.,eta).sum()-y@eta+.5*.1*beta[1]**2
        g=np.array([
            (expit(eta)-y).sum(),
            ((expit(eta)-y)*x).sum()+.1*beta[1]
        ])
        return float(loss),g
    optim=minimize(objective,
                   x0=np.array([np.log(np.mean(y)/(1-np.mean(y))),0.]),
                   jac=True,method="L-BFGS-B",
                   options={"maxiter":1000,"ftol":1e-12})
    if not optim.success or not np.isfinite(optim.x).all():
        raise ValueError("source-only NHM length calibration failed")
    if optim.x[1]<=0:
        raise ValueError("longer source tube does not predict deeper Direct class")
    return dict(
        n_calibration_distinct_original_accepted_species=len(names),
        n_original_source_only_deep=int(sum(y)),
        n_original_source_only_not_deep=int(len(y)-sum(y)),
        intercept=float(optim.x[0]),log1p_cm_slope=float(optim.x[1]),
        original_island_accepted_species_fully_excluded=True,
        source_only_train_loglik=float(-objective(optim.x)[0]),
    )


def make_genus_heldout_length_calibrator(
    external_cm:dict[str,float],depth:dict[str,str],
    mainland_names:set[str],
):
    """Build a *cached* exact leave-one-BOTANICAL-GENUS-out calibration.

    Source-only species entirely absent from ALL GIFT islands had zero
    numeric×ordinal overlap in original data. Do not silently relax that
    failed preflight and claim an island-species-free global fit. Instead,
    omit the ENTIRE target genus (including any species found on other
    islands) from the logistic calibration on each prediction. All input
    observations to this calibration remain independent released NHM numeric
    records paired to the original Direct ordinal species labels, possibly
    provider-related; this is NOT an independent study.
    """
    names=sorted(
        (set(external_cm)&set(depth)&mainland_names)
    )
    names=[x for x in names if depth[x] in ("shallow","intermediate","deep")]
    X=np.array([np.log1p(float(external_cm[x])) for x in names],float)
    Y=np.array([int(depth[x]=="deep") for x in names],float)
    GEN=np.array([botanical_genus(x) for x in names],dtype=str)
    if len(names)<100 or min(Y.sum(),len(Y)-Y.sum())<30 or (
            not np.isfinite(X).all()) or (X<=0).any():
        raise ValueError("source-wide genus-heldout numeric length calibration insufficient")
    all_counts=(int(len(names)),int(Y.sum()),int(len(Y)-Y.sum()))
    @lru_cache(maxsize=None)
    def calibrate(heldout_genus:str)->dict:
        mask=GEN!=heldout_genus
        x=X[mask];y=Y[mask]
        if len(x)<100 or min(y.sum(),len(y)-y.sum())<30:
            raise ValueError("no independent botanical genera for heldout tube calibration")
        def objective(beta):
            eta=beta[0]+beta[1]*x
            resid=expit(eta)-y
            return (
                float(np.logaddexp(0.,eta).sum()-y@eta+.05*beta[1]**2),
                np.array([resid.sum(),resid@x+.1*beta[1]]),
            )
        p=y.mean()
        fit=minimize(
            objective,x0=np.array([np.log(p/(1-p)),0.]),
            jac=True,method="L-BFGS-B",
            options={"maxiter":1000,"ftol":1e-12}
        )
        if not fit.success or fit.x[1]<=0 or not np.isfinite(fit.x).all():
            raise ValueError("leave-one-genus-out length calibration not monotone")
        return dict(
            intercept=float(fit.x[0]),
            log1p_cm_slope=float(fit.x[1]),
            n_calibration_distinct_original_accepted_species=int(mask.sum()),
            heldout_target_genus_never_in_calibration=True,
        )
    metadata=dict(
        calibration_strategy="leave_entire_target_genus_out",
        n_source_accepted_species_with_NHM_numeric_and_Direct_depth=all_counts[0],
        n_source_deep_direct_categories=all_counts[1],
        n_source_non_deep_direct_categories=all_counts[2],
        complete_source_only_calibration_absent_any_original_island_is_not_available=True,
        every_prediction_excludes_its_own_genus_from_all_original_direct_training=True,
        possible_shared_botanical_study_documentation_not_excluded=True,
    )
    return calibrate,metadata


def calibrated_prob(cm:float,cal:dict)->float:
    if not (np.isfinite(cm) and 0<float(cm)<=100):
        raise ValueError("unknown/unreliable original NHM numeric tube cm")
    return float(expit(cal["intercept"]+cal["log1p_cm_slope"]*np.log1p(cm)))


def source_predict_same_island(
    *, full_island:set[str],plant_scopes:dict[str,set[str]],
    source_geo5:set[str],source_climate10:set[str],
    target_depth:dict[str,str],length_cm:dict[str,float],cal
)->dict:
    """Estimate original Direct deep for same observed island species.

    Source-based genus-specific and same-donor pooled models share exactly
    the same outcome plant names, Direct trait axis and donor species set.
    """
    if source_geo5&full_island or source_climate10&full_island:
        raise ValueError("same island accepted name leaked into source donors")
    if any(not z<=full_island for z in plant_scopes.values()):
        raise ValueError("target plant provenance outside original checklist")
    donor_sets={"geo5":source_geo5,"climate10":source_climate10}
    # Never calibrate unrelated source genera that have no original target
    # botanical accepted species in this island.
    needed_genera={botanical_genus(x) for x in full_island
                   if x in target_depth}
    source_genus_numeric={}
    source_genus_all={}
    for mode,source in donor_sets.items():
        any_genus=defaultdict(int)
        donor=defaultdict(list)
        for sp in source:
            g=botanical_genus(sp)
            if g not in needed_genera:continue
            any_genus[g]+=1
            if sp in length_cm:
                genus_model=(cal(g) if callable(cal) else cal)
                if genus_model is None:
                    # Failed single-genus provider-disjoint model support:
                    # explicit source-donor omission, NEVER a zero phenotype.
                    continue
                if callable(cal) and not genus_model[
                    "heldout_target_genus_never_in_calibration"
                ]:
                    raise AssertionError("heldout genus leakage")
                donor[g].append(calibrated_prob(
                    float(length_cm[sp]),genus_model
                ))
        source_genus_numeric[mode]=donor
        source_genus_all[mode]=any_genus
    genera=set(source_genus_numeric["geo5"])&set(source_genus_numeric["climate10"])
    out={}
    for scope,plants in plant_scopes.items():
        names=sorted(
            x for x in plants if x in target_depth and botanical_genus(x) in genera
            and target_depth[x] in ("shallow","intermediate","deep")
        )
        target_genera={botanical_genus(x) for x in names}
        answer=dict(
            n_same_island_direct_target_species=len(names),
            n_original_target_genera=len(target_genera),
            n_target_deep=int(sum(target_depth[s]=="deep" for s in names)),
            island_documentation_fraction=(
                len(names)/sum(botanical_genus(s) in target_genera for s in plants)
                if target_genera else None
            ),
        )
        for mode in MODES:
            bygenus=source_genus_numeric[mode]
            donors=[p for g in target_genera for p in bygenus[g]]
            nsource=len(donors)
            nfull=sum(source_genus_all[mode][g] for g in target_genera)
            baseline=float(np.mean(donors)) if donors else None
            score_genus=score_pool=0.
            for sp in names:
                actual=int(target_depth[sp]=="deep")
                p_g=float(np.mean(bygenus[botanical_genus(sp)]))
                score_genus+=(actual-p_g)**2
                score_pool+=(actual-baseline)**2
            answer[mode]=dict(
                n_distinct_original_source_numeric_congeners=nsource,
                n_original_source_species_in_focal_genera=nfull,
                source_direct_numeric_coverage=nsource/nfull if nfull else None,
                brier_length_genus=score_genus/len(names) if names else None,
                brier_length_pooled=score_pool/len(names) if names else None,
                observed_island_deep_fraction=(
                    answer["n_target_deep"]/len(names) if names else None),
                predicted_pooled_source_deep_probability=baseline,
                predicted_genus_weighted_deep_probability=(
                    np.mean([np.mean(bygenus[botanical_genus(sp)])
                             for sp in names]) if names else None
                ),
                n_original_source_genera_with_at_least_two_distinct_numeric_donors=sum(
                    len(bygenus[g])>=2 for g in target_genera
                ),
            )
        if answer["geo5"]["n_distinct_original_source_numeric_congeners"]<0:
            raise AssertionError("negative accepted source donor opportunity")
        out[scope]=answer
    return out


def aggregate_external_source_numeric_prediction(
    table:pd.DataFrame,*,population=4505,bootstrap=999
)->dict:
    necessary={"island_id","region","scope","tier","spatial_block",
               "first_source","n_target_species","n_target_genera",
               "island_documentation_fraction"}
    for m in MODES:
        necessary|={
            m+"__n_numeric_donors",m+"__source_numeric_coverage",
            m+"__brier_genus",m+"__brier_pooled",
        }
    if x:=necessary-set(table):
        raise ValueError(f"missing source numeric prediction columns {sorted(x)}")
    if table.duplicated(["island_id","scope","tier"]).any():
        raise ValueError("duplicate original botanical source numeric prediction")
    if table.island_id.nunique()!=population or len(table)!=population*len(SCOPES)*len(TIERS):
        raise ValueError("full original 4505 island roster reduced")
    if (not table.region.isin(REGIMES).all() or
        not table.scope.isin(SCOPES).all() or not table.tier.isin(TIERS).all()):
        raise ValueError("source region/status/provenance mismatch")
    if bootstrap<199:raise ValueError("bootstrap draws too small")
    panels=[]
    for r in REGIMES:
        for scope in SCOPES:
            for tier in TIERS:
                allq=table.loc[table.region.eq(r)&table.scope.eq(scope)&table.tier.eq(tier)]
                good=allq.loc[
                    allq.n_target_species.ge(MIN_TARGET)&
                    allq.n_target_genera.ge(MIN_GENERA)&
                    allq.geo5__n_numeric_donors.ge(MIN_SOURCE_DONORS)&
                    allq.climate10__n_numeric_donors.ge(MIN_SOURCE_DONORS)
                ].copy()
                for m in MODES:
                    good[m+"__brier_gain"]=good[m+"__brier_pooled"]-good[m+"__brier_genus"]
                if len(good):
                    balanced=(
                        (good.island_documentation_fraction-good.geo5__source_numeric_coverage).abs().le(.1)&
                        (good.island_documentation_fraction-good.climate10__source_numeric_coverage).abs().le(.1)
                    )
                else:
                    balanced=pd.Series(False,index=good.index)
                out=dict(
                    region=r,scope=scope,tier=tier,
                    n_original_regional_floras=int(len(allq)),
                    n_same_island_both_source_mode_qualified=int(len(good)),
                    n_original_documentation_balanced_both_10pp=int(balanced.sum()),
                )
                for label,p in (("all",good),("coverage10",good.loc[balanced])):
                    for cluster in ("spatial_block","first_source"):
                        grp=p.dropna(subset=[cluster]).groupby(cluster,sort=True)
                        cell=dict(n_islands=int(len(p)),n_groups=int(grp.ngroups))
                        if len(p)<MIN_ISLANDS or grp.ngroups<MIN_BLOCKS:
                            cell["status"]="NOT_IDENTIFIABLE_ORIGINAL_SPATIAL_SOURCE_SUPPORT"
                        else:
                            names=["geo5__brier_gain","climate10__brier_gain"]
                            value=grp[names].mean().to_numpy(float)
                            avg=value.mean(0)
                            seed=int.from_bytes(sha256(
                                f"{r}|{scope}|{tier}|{label}|{cluster}".encode()
                            ).digest()[:8],"big")
                            rng=np.random.default_rng(seed)
                            boots=value[rng.integers(len(value),size=(bootstrap,len(value)))].mean(1)
                            for idx,name in enumerate(names):
                                cell[name+"__equal_group_mean"]=float(avg[idx])
                                cell[name+"__bootstrap95"]=[float(y) for y in np.quantile(boots[:,idx],[.025,.975])]
                            cell["status"]="EXPLORATORY_EXTERNAL_NUMERIC_SOURCE_PREDICTION"
                        out[label+"__"+cluster]=cell
                panels.append(out)
    return dict(
        status="ORIGINAL_4505_NHM_NUMERIC_LENGTH_SOURCE_CONGENER_DEEP_PREDICTION",
        n_original_botanical_flora_islands=population,
        n_original_island_status_source_tier_rows=len(table),
        n_region_status_tier_panels=len(panels),
        panels=panels,
        limitations=dict(
            real_original_NHM_numeric_length_not_local_island_population=True,
            initial_all_island_name_excluded_calibration_failed_due_to_zero_overlap=True,
            each_predicted_genus_direct_ordinal_labels_excluded_from_its_calibration=True,
            other_genera_from_island_floras_can_be_in_calibration=True,
            original_island_ordinal_target_never_used_in_calibration=True,
            same_island_focal_species_across_two_source_modes=True,
            observed_island_genus_roster_given_not_predicted=True,
            original_Direct_annotation_provider_bias_uncontrolled=True,
            no_directly_measured_bird_butterfly_moth_visitors=True,
            numeric_length_not_demonstrated_floral_evolution=True,
            all_4505_island_population_retained=True,
        )
    )
