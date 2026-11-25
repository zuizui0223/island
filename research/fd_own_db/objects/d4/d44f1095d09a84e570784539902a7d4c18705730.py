"""Negative-control design: does floral genus predictability occur on MAINLAND too?

For each GIFT island compare (A) the actual island's SAME triply-direct
plant species to (B) other accepted mainland source donor species from the
SAME genera. All models leave ONE particular source donor out, so neither
actual island forecasts nor mainland pseudo-target forecasts can memorize
their held-out donor's categorical outcome. Plant species in the actual
focal island's full flora are never source donors.

This is NOT independent data-provider validation, historical colonization,
island local population phenotype, or a pollinator effect.
"""
from __future__ import annotations

from hashlib import sha256
import numpy as np
import pandas as pd

MODES=("geo5","climate10")
STRATA=("triple_direct","triple_direct_deep_only")
SECTORS=("radial","tubular")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")


def heldout_source_placebo(*, targets: list[str], donors_by_genus: dict[str,set[str]],
                           radial: dict, tubular: dict) -> dict:
    """Weight source leave-one-out outcomes by real target GENUS frequencies.

    Each target and each donor has three measured traits in both profiles.
    Each target genus needs >=2 distinct mainland source donor species.
    The island prediction is also averaged over the SAME donor omissions
    as the source pseudo-target control to equalize model training sizes.
    """
    if not targets or len(set(targets))!=len(targets):
        raise ValueError("nonempty unique source-independent target names required")
    if not donors_by_genus or any(len(q)<2 for q in donors_by_genus.values()):
        raise ValueError("source placebo needs >=2 donor species per target genus")
    if sum(map(len,donors_by_genus.values()))<3:
        raise ValueError("insufficient source donor species")
    all_donors=set().union(*donors_by_genus.values())
    if set(targets)&all_donors:
        raise ValueError("focal target species leaked into mainland donor pool")
    if sum(map(len,donors_by_genus.values()))!=len(all_donors):
        raise ValueError("same mainland donor assigned to multiple plant genera")
    from .triple_matched_congener import genus
    target_groups={}
    for sp in targets:
        g=genus(sp)
        if g not in donors_by_genus:
            raise ValueError("focal target genus lacks matched source donor")
        target_groups.setdefault(g,[]).append(sp)
    if set(target_groups)!=set(donors_by_genus):
        raise ValueError("source donors outside focal botanical genera")
    for prof in (radial,tubular):
        common=set(prof["known"])
        if not (set(targets)|all_donors)<=common:
            raise ValueError("three Direct morphological traits not known on all plants")
    out={}
    for name,prof in (("radial",radial),("tubular",tubular)):
        present=set(prof["joint"])
        k_total=len(all_donors&present)
        n_total=len(all_donors)
        lg_island=lp_island=lg_source=lp_source=0.0
        for g,targets_g in target_groups.items():
            donors_g=sorted(donors_by_genus[g])
            ys=np.fromiter((sp in present for sp in donors_g),dtype=float)
            ky=ys.sum(); ng=len(ys)
            p_genus=(ky-ys)/(ng-1)
            p_pool=(k_total-ys)/(n_total-1)
            y_target=float(np.mean([sp in present for sp in targets_g]))
            weight=len(targets_g)/len(targets)
            # Binomial sufficient statistics for each omitted donor.
            loss_ig=(y_target*(1-p_genus)**2 +
                     (1-y_target)*p_genus**2).mean()
            loss_ip=(y_target*(1-p_pool)**2 +
                     (1-y_target)*p_pool**2).mean()
            loss_sg=np.square(ys-p_genus).mean()
            loss_sp=np.square(ys-p_pool).mean()
            lg_island+=weight*loss_ig
            lp_island+=weight*loss_ip
            lg_source+=weight*loss_sg
            lp_source+=weight*loss_sp
        out[f"island_gain_{name}"]=float(lp_island-lg_island)
        out[f"source_loo_gain_{name}"]=float(lp_source-lg_source)
        out[f"island_pooled_brier_{name}"]=float(lp_island)
        out[f"source_pooled_brier_{name}"]=float(lp_source)
    out["island_gain_difference"]=(
        out["island_gain_radial"]-out["island_gain_tubular"])
    out["source_loo_gain_difference"]=(
        out["source_loo_gain_radial"]-out["source_loo_gain_tubular"])
    out["island_excess_over_source_loo"]=(
        out["island_gain_difference"]-out["source_loo_gain_difference"])
    return out


def summarize_source_loo_placebo(frame:pd.DataFrame, *,
                                 population:int=4505,draws:int=999) -> dict:
    required={"island_id","regime","spatial_block","source_signature",
              "scope","donor_gate","stratum","n_target","n_target_genera",
              "island_direct_coverage"}
    for m in MODES:
        required|={f"{m}__n_donor",f"{m}__donor_direct_coverage"}
        required|={
            f"{m}__{label}" for label in
            ("island_gain_difference","source_loo_gain_difference",
             "island_excess_over_source_loo","island_gain_radial",
             "island_gain_tubular","source_loo_gain_radial",
             "source_loo_gain_tubular")
        }
    if missing:=required-set(frame):
        raise ValueError(f"missing real source negative-control fields {sorted(missing)}")
    if frame.duplicated(["island_id","scope","donor_gate","stratum"]).any() or (
        frame.island_id.nunique()!=population or len(frame)!=population*12
    ):
        raise ValueError("all 4505 original GIFT island panels must be retained")
    if draws<199:
        raise ValueError("source-LOO group bootstrap too short")
    out=[]
    for region in REGIMES:
        for scope in SCOPES:
            for stratum in STRATA:
                p=frame.loc[
                    frame.regime.eq(region)&
                    frame.scope.eq(scope)&
                    frame.stratum.eq(stratum)&
                    frame.donor_gate.eq("at_least_two_congeners")
                ].copy()
                good=p.n_target.ge(5)&p.n_target_genera.ge(2)&p.geo5__n_donor.ge(10)&p.climate10__n_donor.ge(10)
                for m in MODES:
                    good&=np.isfinite(p[f"{m}__island_excess_over_source_loo"])
                q=p.loc[good].copy()
                balanced=(
                    q.island_direct_coverage.sub(q.geo5__donor_direct_coverage).abs().le(.10)
                    &q.island_direct_coverage.sub(q.climate10__donor_direct_coverage).abs().le(.10)
                )
                for support,data in (("all",q),("doc_gap_le10pp",q.loc[balanced])):
                    for mode in MODES:
                        x=data.dropna(subset=["spatial_block","source_signature"])
                        fields=[f"{mode}__{f}" for f in
                                ("island_gain_radial","island_gain_tubular",
                                 "source_loo_gain_radial","source_loo_gain_tubular",
                                 "island_gain_difference",
                                 "source_loo_gain_difference",
                                 "island_excess_over_source_loo")]
                        item=dict(region=region,scope=scope,stratum=stratum,
                                  source_mode=mode,support=support,
                                  n_original_regional_islands=len(p),
                                  n_supported_islands=len(x),
                                  n_geographic_blocks=x.spatial_block.nunique(),
                                  n_first_sources=x.source_signature.nunique(),
                                  n_documentation_balanced=len(q.loc[balanced]))
                        if (len(x)<30 or item["n_geographic_blocks"]<10 or
                                item["n_first_sources"]<10):
                            item["status"]="NOT_IDENTIFIABLE_SOURCE_LOO_SUPPORT"
                        else:
                            item["status"]="DESCRIPTIVE_REAL_ISLAND_VS_SOURCE_LOO"
                            for grouping,col in (("block","spatial_block"),
                                                 ("first_source","source_signature")):
                                X=x.groupby(col,sort=True)[fields].mean().to_numpy(float)
                                mean=X.mean(axis=0)
                                seed=int.from_bytes(sha256(
                                    f"{region}|{scope}|{stratum}|{mode}|{support}|{grouping}".encode()
                                ).digest()[:8],"big")
                                rng=np.random.default_rng(seed)
                                boot=X[rng.integers(len(X),size=(draws,len(X)))].mean(1)
                                for j,k in enumerate(fields):
                                    key=k.split("__",1)[1]
                                    item[f"{grouping}__{key}"]=float(mean[j])
                                    item[f"{grouping}__{key}_ci95"]=[
                                        float(v) for v in np.quantile(boot[:,j],[.025,.975])]
                        out.append(item)
    return dict(status="ORIGINAL_4505_TRUE_FLORAL_SOURCE_LOO_NEGATIVE_CONTROL",
                n_original_botanical_islands=population,
                n_original_island_scope_gate_stratum_rows=len(frame),
                n_reporting_panels=len(out),
                panels=out,
                limits=dict(
                    same_accepted_target_plant_species=True,
                    same_distinct_mainland_donor_species_for_both_flower_pairs=True,
                    same_donor_leave_one_out_for_island_and_source_control=True,
                    island_genus_composition_known_not_predicted=True,
                    mainland_donors_not_observed_arrivals=True,
                    source_provider_disjoint=False,
                    island_phenotype_provider_disjoint=False,
                    documentation_bias_strong=True,
                    not_causal_guild_or_adaptation=True))
