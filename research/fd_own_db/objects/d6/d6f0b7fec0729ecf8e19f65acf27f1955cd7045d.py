"""NHM original five-corpus *held-out provider* numerical tube prediction.

A tested NHM corpus provides only source-donor centimetre lengths. Every
calibration uses the OTHER FOUR original NHM source corpora and excludes
the ENTIRE predicted botanical genus. Outputs are only original-stratum
aggregate predictions. This removes numerical-provider reuse, NOT overlapping
underlying botanical publications or differential Direct annotation.
"""
from __future__ import annotations

from hashlib import sha256
import numpy as np
import pandas as pd

from .nhm_numeric_source_congener_predict import (
    REGIMES,SCOPES,MODES
)

MIN_ISLANDS=30
MIN_BLOCKS=10
PROVIDER_NAMES=("BHL_text","Ecoflora_text","FloraChina_text",
                "FloraNorthAmerica_text","FloraPakistan_text")


def summarize_provider_disjoint(rows:pd.DataFrame,*,population=4505,draws=999):
    required={"island_id","region","scope","provider","spatial_block",
              "first_source","calibration_supported","n_target",
              "n_target_genera","island_documentation_fraction"}
    for mode in MODES:
        required|={mode+"__n_numeric_donors",
                   mode+"__brier_genus",mode+"__brier_pooled",
                   mode+"__source_coverage"}
    if missing:=required-set(rows):
        raise ValueError(f"missing provider disjoint numerical source fields: {sorted(missing)}")
    if rows.duplicated(["island_id","scope","provider"]).any():
        raise ValueError("duplicate original island scope NHM provider combination")
    if rows.island_id.nunique()!=population or (
        len(rows)!=population*len(SCOPES)*len(PROVIDER_NAMES)
    ):
        raise ValueError("original full GIFT island×plant scope×NHM provider population lost")
    if not rows.scope.isin(SCOPES).all() or (
        not rows.provider.isin(PROVIDER_NAMES).all()
    ) or not rows.region.isin(REGIMES).all():
        raise ValueError("unexpected NHM original provider, geography or plant scope")
    if draws<199:raise ValueError("too few original independent spatial/source group draws")
    if not rows.calibration_supported.isin([True,False]).all():
        raise ValueError("NHM train-support state must be explicitly resolved")
    if (rows.loc[~rows.calibration_supported,"n_target"]>0).any():
        raise ValueError("failed NHM disjoint calibration must not produce predictions")
    for mode in MODES:
        scores=rows[mode+"__brier_genus"]
        if (scores.notna() != rows.n_target.gt(0)).any():
            raise ValueError("source donor prediction undefined is not a biological zero")
        if not scores.dropna().between(0,1).all():
            raise ValueError("non-probabilistic source-provider Brier metric")
        if ((rows[mode+"__n_numeric_donors"]<0)|
            (rows.n_target<0)).any():
            raise ValueError("negative original GIFT provider support count")
    out=[]
    for region in REGIMES:
        for scope in SCOPES:
            for provider in PROVIDER_NAMES:
                p=rows.loc[
                    rows.region.eq(region)&rows.scope.eq(scope)&
                    rows.provider.eq(provider)
                ]
                supported=p.loc[
                    p.calibration_supported&
                    p.n_target.ge(5)&p.n_target_genera.ge(2)&
                    p.geo5__n_numeric_donors.ge(10)&
                    p.climate10__n_numeric_donors.ge(10)
                ].copy()
                for mode in MODES:
                    supported[mode+"__gain"]=(
                        supported[mode+"__brier_pooled"]-
                        supported[mode+"__brier_genus"]
                    )
                if len(supported):
                    balanced=(
                        (supported.island_documentation_fraction-
                         supported.geo5__source_coverage).abs().le(.10)&
                        (supported.island_documentation_fraction-
                         supported.climate10__source_coverage).abs().le(.10)
                    )
                else:balanced=pd.Series(False,index=supported.index)
                record=dict(
                    region=region,scope=scope,source_numeric_heldout_provider=provider,
                    n_original_region_floras=int(len(p)),
                    n_qualified_source_provider_disjoint_islands=int(len(supported)),
                    n_direct_NHM_provider_coverage_balanced_both_le10pp=int(balanced.sum()),
                    n_original_islands_without_numeric_model=int((~p.calibration_supported).sum()),
                )
                for name,data in (("all",supported),("documentation_balance10",supported.loc[balanced])):
                    for grouping in ("spatial_block","first_source"):
                        groups=data.dropna(subset=[grouping]).groupby(grouping,sort=True)
                        r=dict(n_original_island_floras=int(len(data)),
                               n_groups=int(groups.ngroups))
                        if len(data)<MIN_ISLANDS or groups.ngroups<MIN_BLOCKS:
                            r["status"]="NOT_IDENTIFIABLE_SOURCE_PROVIDER_SUPPORT"
                        else:
                            cols=[mode+"__gain" for mode in MODES]
                            X=groups[cols].mean().to_numpy(float)
                            means=X.mean(axis=0)
                            seed=int.from_bytes(sha256(
                                f"{region}|{scope}|{provider}|{name}|{grouping}".encode()
                            ).digest()[:8],"big")
                            rng=np.random.default_rng(seed)
                            boot=X[rng.integers(len(X),size=(draws,len(X)))].mean(1)
                            for k,col in enumerate(cols):
                                r[col+"__block_equal"]=float(means[k])
                                r[col+"__bootstrap95"]=[
                                    float(v) for v in np.quantile(boot[:,k],[.025,.975])
                                ]
                            r["status"]="EXPLORATORY_ORIGINAL_PROVIDER_DISJOINT_NUMERIC_SOURCE_PREDICTION"
                        record[name+"__"+grouping]=r
                out.append(record)
    return dict(
        status="ORIGINAL_4505_NHM_PROVIDER_DISJOINT_TUBE_SOURCE_PREDICTION",
        n_original_GIFT_island_floras=population,
        n_original_island_scope_provider_rows=int(len(rows)),
        n_regional_status_source_provider_panels=len(out),
        panels=out,
        limitations=dict(
            each_NHM_source_donor_length_from_one_heldout_named_corpus=True,
            calibration_numeric_lengths_from_other_four_pinned_corpora=True,
            each_target_genus_depth_categories_excluded_from_its_own_fit=True,
            source_provider_disjoint_not_underlying_publication_independent=True,
            no_island_local_morphology_or_effective_visitation_measured=True,
            no_undocumented_source_or_island_trait_filled=True,
            strict_native_provenance_is_selective=True,
            original_all_4505_botanical_islands_retained=True,
            no_named_bird_butterfly_moth_pollinator_inference=True,
        ),
    )
