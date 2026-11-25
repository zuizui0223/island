"""Exact SAME-island SAME-plant contrast between two original NHM source corpora.

Target botanical genera must have numerical congener donors in BOTH GIFT
mainland source definitions AND BOTH heldout named NHM corpora. Thus all
scores below concern the identical focal Direct-phenotyped island species.
A corpus name is not a distinct underlying original botanical study.
"""
from __future__ import annotations
from hashlib import sha256
import numpy as np
import pandas as pd

PAIR=("BHL_text","FloraChina_text")
SCOPES=("all_recorded","strict_native","strict_native_nonendemic")
REGIONS=("northern_midlatitude","tropical","southern_extratropical",
         "northern_high_latitude","geography_unresolved")
MODES=("geo5","climate10")
MIN_ISLAND=5
MIN_GENUS=2
MIN_DONORS=10
MIN_REGIONAL_ISLANDS=30
MIN_BLOCKS=10


def summarize_paired_nhm_providers(frame:pd.DataFrame,*,population=4505,draws=999):
    cols={"island_id","region","scope","spatial_block","first_source",
          "n_target","n_target_genera","island_documentation_fraction"}
    for provider in PAIR:
        for mode in MODES:
            cols |= {
                f"{provider}__{mode}__n_numeric",
                f"{provider}__{mode}__source_coverage",
                f"{provider}__{mode}__brier_genus",
                f"{provider}__{mode}__brier_pooled",
                f"{provider}__{mode}__n_accepted_donor_species_shared_with_other_provider",
            }
    if missing:=cols-set(frame):
        raise ValueError(f"missing same-plant provider pair fields: {sorted(missing)}")
    if frame.duplicated(["island_id","scope"]).any():
        raise ValueError("original island status provider pair duplicated")
    if len(frame)!=population*len(SCOPES) or frame.island_id.nunique()!=population:
        raise ValueError("original 4505 island floral universe contracted")
    if not frame.scope.isin(SCOPES).all() or not frame.region.isin(REGIONS).all():
        raise ValueError("unknown botanical nativeness or GSHHG region")
    if draws<199:raise ValueError("cluster bootstrap draws insufficient")
    for provider in PAIR:
        for mode in MODES:
            loss=frame[f"{provider}__{mode}__brier_genus"]
            if (loss.notna()!=frame.n_target.gt(0)).any():
                raise ValueError("missing original provider pair score treated as biological zero")
            if not loss.dropna().between(0,1).all():
                raise ValueError("invalid named-provider Brier score")
            if ((frame[f"{provider}__{mode}__n_numeric"]<0)|
                (frame[f"{provider}__{mode}__n_accepted_donor_species_shared_with_other_provider"]<0)).any():
                raise ValueError("negative original provider source opportunity")
    summary=[]
    for region in REGIONS:
        for scope in SCOPES:
            panel=frame.loc[frame.region.eq(region)&frame.scope.eq(scope)].copy()
            enough=(panel.n_target.ge(MIN_ISLAND)&
                    panel.n_target_genera.ge(MIN_GENUS))
            for provider in PAIR:
                for mode in MODES:
                    enough &= panel[f"{provider}__{mode}__n_numeric"].ge(MIN_DONORS)
            paired=panel.loc[enough].copy()
            for provider in PAIR:
                for mode in MODES:
                    paired[f"{provider}__{mode}__gain"]=(
                        paired[f"{provider}__{mode}__brier_pooled"]-
                        paired[f"{provider}__{mode}__brier_genus"]
                    )
            for mode in MODES:
                paired[mode+"__FloraChina_minus_BHL_gain"]=(
                    paired[f"FloraChina_text__{mode}__gain"]-
                    paired[f"BHL_text__{mode}__gain"]
                )
            if len(paired):
                balanced=pd.Series(True,index=paired.index)
                for provider in PAIR:
                    for mode in MODES:
                        balanced &= (
                            paired.island_documentation_fraction-
                            paired[f"{provider}__{mode}__source_coverage"]
                        ).abs().le(.10)
            else:balanced=pd.Series(False,index=paired.index)
            one=dict(
                region=region,scope=scope,
                n_original_region_islands=int(len(panel)),
                n_common_islands_same_direct_plants_with_both_NHM_donor_corpora=int(len(paired)),
                n_original_islands_balanced_across_all_four_NHM_sources=int(balanced.sum()),
                median_target_direct_species=float(paired.n_target.median()) if len(paired) else None,
                median_target_genera=float(paired.n_target_genera.median()) if len(paired) else None,
                n_provider_donor_overlaps_median_by_geo5=float(
                    paired[f"BHL_text__geo5__n_accepted_donor_species_shared_with_other_provider"].median()
                ) if len(paired) else None,
                same_exact_focal_island_accepted_names_in_both_models=True,
                same_exact_focal_island_accepted_names_in_both_source_modes=True,
            )
            for label,sub in (("all",paired),("documentation_balance10",paired.loc[balanced])):
                for groupby in ("spatial_block","first_source"):
                    groups=sub.dropna(subset=[groupby]).groupby(groupby,sort=True)
                    res=dict(n_original_islands=int(len(sub)),n_groups=int(groups.ngroups))
                    if len(sub)<MIN_REGIONAL_ISLANDS or groups.ngroups<MIN_BLOCKS:
                        res["status"]="NOT_IDENTIFIABLE_MATCHED_PROVIDER_GEO_SUPPORT"
                    else:
                        columns=[
                            f"{provider}__{mode}__gain"
                            for provider in PAIR for mode in MODES
                        ]+[mode+"__FloraChina_minus_BHL_gain" for mode in MODES]
                        X=groups[columns].mean().to_numpy(float)
                        mean=X.mean(axis=0)
                        seed=int.from_bytes(sha256(
                            f"{region}|{scope}|{label}|{groupby}".encode()
                        ).digest()[:8],"big")
                        rng=np.random.default_rng(seed)
                        b=X[rng.integers(len(X),size=(draws,len(X)))].mean(1)
                        for i,k in enumerate(columns):
                            res[k+"__equal_group_mean"]=float(mean[i])
                            res[k+"__bootstrap95"]=[
                                float(z) for z in np.quantile(b[:,i],[.025,.975])
                            ]
                        res["status"]="EXPLORATORY_SAME_PLANT_TWO_NHM_CORPUS_SOURCE_PREDICTION"
                    one[label+"__"+groupby]=res
            summary.append(one)
    return dict(
        status="ORIGINAL_4505_SAME_PLANT_BHL_VS_FLORACHINA_TUBE_SOURCE",
        n_original_botanical_GIFT_islands=population,
        n_original_private_island_scope_matched_provider_rows=len(frame),
        n_original_region_provenance_panels=len(summary),
        compared_named_NHM_corpora=list(PAIR),
        panels=summary,
        limitations=dict(
            original_same_target_species_within_island_across_both_corpora=True,
            original_same_target_species_across_geo5_climate10_source_modes=True,
            distinct_heldout_NHM_corpora_not_independent_field_studies=True,
            both_source_genus_model_cohorts_strongly_selected=True,
            source_heldout_calibration_excludes_predicted_genus=True,
            no_plant_specific_bird_butterfly_moth_effective_visit=True,
            no_directed_ecological_or_evolutionary_causality=True,
            all_original_4505_island_floras_in_denominator=True,
        ),
    )
