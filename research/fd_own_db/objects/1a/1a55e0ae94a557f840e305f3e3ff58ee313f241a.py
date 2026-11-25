"""Full documented island-source floral difference: exact signed genus accounting.

Each mode's observed difference is partitioned algebraically into four signed
terms using actual Direct-labelled accepted species. These terms may cancel;
they are NOT causal variance fractions or a pollinator mechanism.
"""
from __future__ import annotations

from hashlib import sha256
import numpy as np
import pandas as pd

MORPHOLOGY=("zygomorphic","tube_or_bell","open_radial","deep_tube",
            "bilateral_tube_or_bell_joint")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic",
        "native_nonendemic_source_listed_both")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
TERMS=("shared_genus_composition","shared_genus_within_species_states",
       "native_exclusive_genus","source_exclusive_genus_negative")
MIN_ISLAND_DIRECT=5
MIN_SOURCE_DIRECT=10
MIN_ISLANDS=30
MIN_BLOCKS=10


def _validate(frame, population):
    needed={"island_id","regime","scope","sector","spatial_block",
            "source_signature","n_scope_total","n_scope_direct"}
    for mode in ("geo5","climate10"):
        needed|={
            f"{mode}__n_source_total",f"{mode}__n_source_direct",
            f"{mode}__observed_difference",
            f"{mode}__n_shared_genera",
            f"{mode}__n_native_carriers_from_shared_genera",
            *(f"{mode}__{term}" for term in TERMS),
        }
    missing=needed-set(frame)
    if missing:
        raise ValueError(f"missing original exact genus parts: {sorted(missing)}")
    if frame.island_id.nunique()!=population:
        raise ValueError("all original botanical island IDs must remain")
    if frame.duplicated(["island_id","scope","sector"]).any():
        raise ValueError("duplicate island scope sector")
    if len(frame)!=population*len(SCOPES)*len(MORPHOLOGY):
        raise ValueError("original plant statuses or direct floral categories collapsed")
    if not frame.regime.isin(REGIMES).all() or not frame.scope.isin(SCOPES).all() or not frame.sector.isin(MORPHOLOGY).all():
        raise ValueError("unknown source regime/provenance/real floral category")
    for mode in ("geo5","climate10"):
        part=frame[[f"{mode}__{term}" for term in TERMS]]
        y=frame[f"{mode}__observed_difference"]
        defined=part.notna().all(axis=1)
        if ((part.notna().any(axis=1))!=defined).any():
            raise ValueError("partial genus-accounting terms are not an observed zero")
        if ((y.notna()) != defined).any():
            raise ValueError("observed difference and accounting terms disagree on missingness")
        if defined.any() and not np.allclose(
            part.loc[defined].sum(axis=1).to_numpy(float),
            y.loc[defined].to_numpy(float),rtol=0,atol=1e-10
        ):
            raise ValueError("four signed genus components fail exact reconciliation")


def _group_summary(q,group,seed,bootstrap):
    g=q.dropna(subset=[group]).groupby(group,sort=True)
    n=g.ngroups
    o=dict(n_original_islands=len(q),n_original_source_clusters=n)
    if len(q)<MIN_ISLANDS or n<MIN_BLOCKS:
        o["status"]="NOT_IDENTIFIABLE_INDEPENDENT_CLUSTER_SUPPORT"
        return o
    # Explicitly preserve every signed (potentially cancelling) term.
    cols=[f"{mode}__{key}" for mode in ("geo5","climate10")
          for key in ("observed_difference",*TERMS)]
    x=g[cols].mean().to_numpy(float)
    mean=x.mean(axis=0)
    rng=np.random.default_rng(int.from_bytes(sha256(seed.encode()).digest()[:8],"big"))
    simulated=x[rng.integers(n,size=(bootstrap,n))].mean(axis=1)
    for k,col in enumerate(cols):
        o[col+"__block_equal"]=float(mean[k])
        o[col+"__group_bootstrap95"]=[float(v) for v in np.quantile(simulated[:,k],[.025,.975])]
    o["status"]="EXACT_SIGNED_GENUS_ACCOUNTING_DESCRIPTIVE_NOT_CAUSAL"
    return o


def summarize_exact_genus_partition(frame,population=4505,bootstrap=999):
    _validate(frame,population)
    if bootstrap<199:
        raise ValueError("bootstrap draws insufficient")
    result=[]
    for reg in REGIMES:
        for scope in SCOPES:
            for sector in MORPHOLOGY:
                t=frame.loc[
                    frame.regime.eq(reg)&frame.scope.eq(scope)&frame.sector.eq(sector)
                ]
                support=(
                    t.n_scope_direct.ge(MIN_ISLAND_DIRECT)
                    & t.geo5__n_source_direct.ge(MIN_SOURCE_DIRECT)
                    & t.climate10__n_source_direct.ge(MIN_SOURCE_DIRECT)
                    & t.geo5__observed_difference.notna()
                    & t.climate10__observed_difference.notna()
                )
                q=t.loc[support].copy()
                if len(q):
                    p_island=q.n_scope_direct/q.n_scope_total.replace(0,np.nan)
                    p_geo=q.geo5__n_source_direct/q.geo5__n_source_total.replace(0,np.nan)
                    p_clim=q.climate10__n_source_direct/q.climate10__n_source_total.replace(0,np.nan)
                    balanced=(p_island.sub(p_geo).abs().le(.1)&
                              p_island.sub(p_clim).abs().le(.1))
                else:
                    balanced=pd.Series(False,index=q.index)
                row=dict(
                    region=reg,scope=scope,structural_sector=sector,
                    n_region_original_botanical_islands=int(len(t)),
                    n_same_original_islands_with_both_sources=int(len(q)),
                    n_documentation_gap_both_le10pp=int(balanced.sum()),
                    median_number_shared_source_genera_geo5=(
                        float(q.geo5__n_shared_genera.median()) if len(q) else None),
                    median_number_shared_source_genera_climate10=(
                        float(q.climate10__n_shared_genera.median()) if len(q) else None),
                    island_direct_species_fraction_from_shared_genera_geo5=(
                        float((q.geo5__n_native_carriers_from_shared_genera/
                               q.n_scope_direct.clip(lower=1)).median())
                        if len(q) else None),
                    island_direct_species_fraction_from_shared_genera_climate10=(
                        float((q.climate10__n_native_carriers_from_shared_genera/
                               q.n_scope_direct.clip(lower=1)).median())
                        if len(q) else None),
                )
                for subset_name,z in (("all",q),("documentation_gap_le10pp",q.loc[balanced])):
                    for cluster in ("spatial_block","source_signature"):
                        row[f"{subset_name}__{cluster}"]=_group_summary(
                            z,cluster,
                            seed=f"{reg}|{scope}|{sector}|{subset_name}|{cluster}",
                            bootstrap=bootstrap,
                        )
                result.append(row)
    return dict(
        status="ORIGINAL_GIFT_4505_GEO5_CLIMATE10_EXACT_SIGNED_GENUS_PARTITION",
        n_original_botanical_floras=population,
        n_region_scope_sector_panels=len(result),
        n_original_private_island_scope_sector_cells=int(len(frame)),
        panels=result,
        scientific_limits=dict(
            signed_terms_can_cancel_and_are_not_causal_shares=True,
            same_focal_plants_in_two_source_mode_comparisons=True,
            source_candidate_species_are_not_known_arrivals=True,
            no_shared_genus_does_not_equal_zero_evolution=True,
            no_historical_pollinator_community_or_effective_visits_measured=True,
            four_status_scopes_nested_and_nonrandom=True,
            original_4505_island_population_unrestricted=True,
        ),
    )
