"""Source-only, fixed-island, exact signed botanical genus accounting for
REAL same-accepted-species floral trait pair configurations.

Not a causal decomposition of pollinator-associated evolution. In particular,
paired flower morphology can turn over because different plant lineages occur
on islands; the source-only congeneric record is a documentation contrast.
"""
from __future__ import annotations
from hashlib import sha256
import numpy as np
import pandas as pd

TARGETS=("actinomorphic_deep","zygomorphic_deep","tubular_bell_deep")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIMES=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
PARTS=("shared_genus_composition","shared_genus_within_species_states",
       "native_exclusive_genus","source_exclusive_genus_negative")
MIN_ISLAND=5
MIN_MAINLAND=10
MIN_REGION_ISLANDS=30
MIN_BLOCKS=10


def summarize_real_joint_genus(frame:pd.DataFrame, *, population:int=4505,
                               draws:int=999) -> dict:
    required={"island_id","region","scope","sector","spatial_block",
              "first_source","n_island_total","n_island_direct"}
    for mode in ("geo5","climate10"):
        required |= {
            f"{mode}__n_mainland_total", f"{mode}__n_mainland_direct",
            f"{mode}__observed_difference",
            f"{mode}__n_shared_genera",
            *(f"{mode}__{z}" for z in PARTS),
        }
    missing=required-set(frame)
    if missing:
        raise ValueError(f"missing source-JOINT genus accounting fields {sorted(missing)}")
    if frame.duplicated(["island_id","scope","sector"]).any():
        raise ValueError("duplicate original island scope joint phenotype")
    if frame.island_id.nunique()!=population or (
            len(frame)!=population*len(SCOPES)*len(TARGETS)):
        raise ValueError("full original botanical population or paired trait set reduced")
    if not frame.scope.isin(SCOPES).all() or not frame.sector.isin(TARGETS).all():
        raise ValueError("invalid provenance or confirmed genuine pair target")
    if not frame.region.isin(REGIMES).all():
        raise ValueError("unrecognized GSHHG geography")
    if draws<199:
        raise ValueError("insufficient geographical block resamples")
    for mode in ("geo5","climate10"):
        matrix=frame[[f"{mode}__{part}" for part in PARTS]]
        defined=matrix.notna().all(axis=1)
        if (matrix.notna().any(axis=1)!=defined).any():
            raise ValueError("partially missing genus decomposition")
        raw=frame[f"{mode}__observed_difference"]
        if (raw.notna()!=defined).any():
            raise ValueError("observed flower pair and exact terms missingness mismatch")
        if defined.any() and not np.allclose(
                matrix.loc[defined].sum(axis=1),raw.loc[defined],
                atol=1e-10,rtol=0):
            raise AssertionError("four signed botanical genus contributions do not sum")

    summary=[]
    for reg in REGIMES:
        for scope in SCOPES:
            for sector in TARGETS:
                records=frame.loc[
                    frame.region.eq(reg)&frame.scope.eq(scope)
                    &frame.sector.eq(sector)
                ].copy()
                both=records.loc[
                    records.n_island_direct.ge(MIN_ISLAND)
                    &records.geo5__n_mainland_direct.ge(MIN_MAINLAND)
                    &records.climate10__n_mainland_direct.ge(MIN_MAINLAND)
                    &records.geo5__observed_difference.notna()
                    &records.climate10__observed_difference.notna()
                ].copy()
                if len(both):
                    i=both.n_island_direct/both.n_island_total.replace(0,np.nan)
                    g=both.geo5__n_mainland_direct/both.geo5__n_mainland_total.replace(0,np.nan)
                    c=both.climate10__n_mainland_direct/both.climate10__n_mainland_total.replace(0,np.nan)
                    balanced=i.sub(g).abs().le(.10)&i.sub(c).abs().le(.10)
                else:
                    balanced=pd.Series(False,index=both.index)
                cell=dict(
                    region=reg,scope=scope,sector=sector,
                    n_original_botanical_region_floras=len(records),
                    n_same_original_islands_both_mainland_sources=int(len(both)),
                    n_original_islands_documentation_gap_both_le10pp=int(balanced.sum()),
                    median_shared_genus_geo5=float(both.geo5__n_shared_genera.median())
                    if len(both) else None,
                    median_shared_genus_climate10=float(both.climate10__n_shared_genera.median())
                    if len(both) else None,
                )
                for k,p in (("all",both),("doc_gap_le10pp",both.loc[balanced])):
                    for axis in ("spatial_block","first_source"):
                        group=p.dropna(subset=[axis]).groupby(axis,sort=True)
                        ngrp=group.ngroups
                        d=dict(n_islands=int(len(p)),n_geographical_or_source_groups=int(ngrp))
                        if len(p)<MIN_REGION_ISLANDS or ngrp<MIN_BLOCKS:
                            d["status"]="NOT_IDENTIFIABLE_BLOCK_OR_SOURCE_SUPPORT"
                        else:
                            keys=[f"{mode}__{part}"
                                  for mode in ("geo5","climate10")
                                  for part in ("observed_difference",*PARTS)]
                            values=group[keys].mean().to_numpy(float)
                            ave=values.mean(0)
                            seed=int.from_bytes(sha256(
                                f"{reg}|{scope}|{sector}|{k}|{axis}".encode()
                            ).digest()[:8],"big")
                            rng=np.random.default_rng(seed)
                            bs=values[rng.integers(ngrp,size=(draws,ngrp))].mean(1)
                            for j,key in enumerate(keys):
                                d[key+"__equal_group_mean"]=float(ave[j])
                                d[key+"__bootstrap95"]=[
                                    float(z) for z in np.quantile(bs[:,j],[.025,.975])
                                ]
                            d["status"]="EXACT_REAL_PAIRED_FLORA_GENUS_ACCOUNTING_NOT_CAUSAL"
                        cell[f"{k}__{axis}"]=d
                summary.append(cell)
    return dict(
        status="ORIGINAL_4505_REAL_SYNDROME_PAIR_SOURCE_GENUS_PARTITION",
        n_original_GIFT_botanical_islands=population,
        n_original_island_scope_true_pair_rows=len(frame),
        n_regional_scope_true_pair_panels=len(summary),
        panels=summary,
        limits={
            "same_species_binary_joint_state_only":True,
            "within_shared_genus_different_accepted_species_not_ancestor":True,
            "signed_genus_terms_can_cancel_and_cannot_be_called_causal_shares":True,
            "source_candidate_not_observed_seed_arrival_or_establishment_failure":True,
            "no_effective_bird_butterfly_hawkmoth_pollinator_identity":True,
            "all_4505_island_floras_keep_denominator":True,
            "posthoc_joint_morphology_source_validation_not_preregistered":True,
        },
    )
