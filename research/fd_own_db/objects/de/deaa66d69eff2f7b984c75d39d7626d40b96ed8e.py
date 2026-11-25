"""Outcome-blind island-vs-mainland same-plant VISITOR NETWORK eligibility.

Uses corrected original EuPPollNet network geography and the original fd
plant species × island flora. All original taxon, network, study, polygon and
coordinate identities remain private. This is a cross-sectional *comparability
gate*, not historical pollinator loss, adaptive floral evolution or pollen
delivery.

A source study sampling the SAME plant in the SAME calendar year on an fd
island and on a verified continental landmass is a candidate comparison. The
nearest mainland site must be <=500 km by great-circle distance; an additional
<=250 km sensitivity is reported. Neither matching is ancestry.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .euppollnet_island_overlap import validate_original_visits

SITE_KEY=("Study_id","Network_id","Latitude","Longitude")
GROUP=("Study_id","Year","Plant_accepted_name")
LOCATION=("habitat_class","island_id","analysis_regime","spatial_block")
HABITATS=frozenset(("fd_island","mainland","unresolved"))
MIN_STUDIES=3
MIN_BLOCKS=5
MIN_ISLANDS=10
MIN_NATIVE_PAIR_GROUPS=20


def haversine_km(lat:float,lon:float,lats:np.ndarray,lons:np.ndarray)->np.ndarray:
    if not (-90<=lat<=90 and -180<=lon<=180):
        raise ValueError("Invalid source visit latitude/longitude")
    rad=np.pi/180
    a=np.sin((lats-lat)*rad/2.)**2+np.cos(lat*rad)*np.cos(lats*rad)*(
        np.sin((lons-lon)*rad/2.)**2)
    return 6371.0088*(2*np.arcsin(np.sqrt(np.maximum(0,np.minimum(a,1.)))))


def _location_join(visits:pd.DataFrame,site_class:pd.DataFrame)->pd.DataFrame:
    required=set((*SITE_KEY,*LOCATION))
    if not required.issubset(site_class):
        raise ValueError("Exact original GSHHG site classification absent")
    loc=site_class[[*SITE_KEY,*LOCATION]].copy()
    if loc[list(SITE_KEY)].isna().any().any() or loc.duplicated(list(SITE_KEY)).any():
        raise ValueError("Original study/network/geosite is ambiguous")
    if not loc.habitat_class.isin(HABITATS).all():
        raise ValueError("Unresolved island/mainland class")
    island=loc.habitat_class.eq("fd_island")
    if (loc.loc[island,["island_id","spatial_block","analysis_regime"]]
            .isna().any().any()):
        raise ValueError("Original fd island lacks exact island ID and geographic block")
    if loc.loc[~island,["island_id","spatial_block"]].notna().any().any():
        raise ValueError("Mainland or unresolved cannot carry an fd island ID")
    return visits.merge(loc,on=list(SITE_KEY),how="inner",validate="many_to_one")


def original_same_plant_island_mainland_opportunity(
    interaction:pd.DataFrame, site_class:pd.DataFrame, original_flora:pd.DataFrame,
    source_metadata:pd.DataFrame,
) -> dict:
    """Original source-qualified *coverage* only, without visitor outcomes.

    No differences in insect composition are calculated or tested here; both
    original plant-name and the study's actual sampling scope must match.
    """
    if not {"island_id","accepted_species","origin_status"}.issubset(original_flora):
        raise ValueError("Original fd plant island checklist missing")
    if original_flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Duplicate original island × accepted plant name")
    required={"Study_id","Year","Sampling_method","Taxa_recorded"}
    if not required.issubset(source_metadata):
        raise ValueError("Original EuPPollNet study metadata schema drift")
    x=validate_original_visits(interaction)
    x=x.loc[x.Plant_rank.eq("SPECIES")].copy()
    m=source_metadata[list(required)].copy()
    m.Year=pd.to_numeric(m.Year,errors="coerce")
    m=m.loc[m.Year.notna()].copy()
    m.Year=m.Year.astype(int)
    # If original metadata has two contradictory reported protocols for the
    # same study/year, do NOT call those sampling frames comparable.
    conflicts=(m.groupby(["Study_id","Year"])
               [["Sampling_method","Taxa_recorded"]].nunique(dropna=False).gt(1).any(axis=1))
    bad=set(conflicts[conflicts].index)
    m=m.loc[~pd.MultiIndex.from_frame(m[["Study_id","Year"]]).isin(bad)]
    m=m.drop_duplicates(["Study_id","Year"])
    x=_location_join(x,site_class)
    x=x.loc[x.habitat_class.isin(("fd_island","mainland"))].copy()
    if x.empty:
        return dict(status="NO_VERIFIED_ORIGINAL_CONTINENT_ISLAND_VISITS",
            n_mapped_source_plant_visit_rows=0,n_matching_same_study_year_species_island_mainland=0,
            actual_pollinator_replacement_measured=False,
            original_island_study_species_names_or_coordinates_exported=False)
    x=x.merge(m,on=["Study_id","Year"],how="inner",validate="many_to_one")
    x=x.loc[x.Sampling_method.notna() & x.Taxa_recorded.notna()].copy()
    flora=original_flora[["island_id","accepted_species","origin_status"]].rename(
        columns={"accepted_species":"Plant_accepted_name"})
    island=x.loc[x.habitat_class.eq("fd_island")].copy()
    island=island.merge(flora,on=["island_id","Plant_accepted_name"],
        how="inner",validate="many_to_one")
    mainland=x.loc[x.habitat_class.eq("mainland")].copy()
    if island.empty or mainland.empty:
        return dict(
            status="NO_SAME_STUDY_ISLAND_MAINLAND_ACTUAL_VISIT_SUPPORT",
            n_original_fd_islands_with_visit_plant_flora=int(island.island_id.nunique()),
            n_continental_source_visit_studies=int(mainland.Study_id.nunique()),
            n_matching_same_study_year_species_island_mainland=0,
            actual_pollinator_replacement_measured=False,
            original_island_study_species_names_or_coordinates_exported=False,
        )
    # Distinct site×plant×study×year only. NEVER count multiple flower
    # interaction records as independent sampling sites.
    i=island.drop_duplicates([*GROUP,*SITE_KEY,"island_id"])
    msites=mainland.drop_duplicates([*GROUP,*SITE_KEY])
    mdict={key:frame[["Latitude","Longitude"]].drop_duplicates().to_numpy(float)
           for key,frame in msites.groupby(list(GROUP),sort=False)}
    rows=[]
    for row in i.itertuples(index=False):
        key=(row.Study_id,row.Year,row.Plant_accepted_name)
        targets=mdict.get(key)
        if targets is None or len(targets)==0:
            continue
        distances=haversine_km(float(row.Latitude),float(row.Longitude),
                 targets[:,0],targets[:,1])
        smallest=float(np.min(distances))
        rows.append(dict(
            island_id=row.island_id,
            spatial_block=row.spatial_block,
            analysis_regime=row.analysis_regime,
            study_id=row.Study_id,year=row.Year,
            plant=row.Plant_accepted_name,
            origin_status=row.origin_status,
            same_study_sampling_method=row.Sampling_method,
            same_study_taxa_scope=row.Taxa_recorded,
            minimum_original_mainland_site_distance_km=smallest))
    if not rows:
        return dict(
            status="NO_EXACT_SAME_PLANT_SAME_STUDY_YEAR_ISLAND_MAINLAND",
            n_original_fd_islands_with_visit_plant_flora=int(island.island_id.nunique()),
            n_continental_source_visit_studies=int(mainland.Study_id.nunique()),
            n_matching_same_study_year_species_island_mainland=0,
            actual_pollinator_replacement_measured=False,
            original_island_study_species_names_or_coordinates_exported=False,
        )
    pairs=pd.DataFrame(rows)
    # Deduplicate networks/sites before calculating independent pair support;
    # multiple network visits of the same original plant are not new species.
    matched=pairs.drop_duplicates(["island_id","study_id","year","plant"])
    near=matched.loc[matched.minimum_original_mainland_site_distance_km.le(500)]
    native=near.loc[near.origin_status.eq("native")]
    near250=matched.loc[matched.minimum_original_mainland_site_distance_km.le(250)]
    # A study with explicit bee-and-nonbee target scope is at least a
    # candidate for cross-guild comparison; cannot equalize visit effort.
    broad=near.loc[near.same_study_taxa_scope.fillna("").astype(str).str.lower().map(
        lambda s:("hymenoptera" in s or "bees" in s) and
              ("diptera" in s or "flies" in s or "syrphid" in s))]
    support=(
       near.island_id.nunique()>=MIN_ISLANDS
       and near.spatial_block.nunique()>=MIN_BLOCKS
       and near.study_id.nunique()>=MIN_STUDIES
       and len(native)>=MIN_NATIVE_PAIR_GROUPS
       and broad.island_id.nunique()>=MIN_ISLANDS
    )
    output={
        "status":("SOURCE_REPEATED_OBSERVATIONAL_GUILD_COMPARISON_ELIGIBLE"
                  if support else "SAME_PLANT_ISLAND_MAINLAND_COMPARISON_SUPPORT_INSUFFICIENT"),
        "source":"EuPPollNet_v1.3.0_corrected_original",
        "units":"island_x_original_study_x_year_x_exact_accepted_plant_species",
        "n_original_fd_islands_with_any_visited_plant_in_flora":int(island.island_id.nunique()),
        "n_source_continental_visit_studies":int(mainland.Study_id.nunique()),
        "n_same_study_year_plant_candidate_groups_any_distance":int(len(matched)),
        "n_pair_groups_nearest_mainland_at_most_500_km":int(len(near)),
        "n_pair_groups_nearest_mainland_at_most_250_km":int(len(near250)),
        "n_eligible_strict_native_pair_groups_500km":int(len(native)),
        "n_eligible_broad_sampling_scope_pair_groups_500km":int(len(broad)),
        "n_eligible_islands_500km":int(near.island_id.nunique()),
        "n_eligible_spatial_blocks_500km":int(near.spatial_block.nunique()),
        "n_eligible_original_studies_500km":int(near.study_id.nunique()),
        "n_eligible_unique_plant_names_500km":int(near.plant.nunique()),
        "n_broad_study_scope_islands_500km":int(broad.island_id.nunique()),
        "n_native_pairs_by_regime_500km":[{
            "regime":str(name),
            "n_native_comparisons":int(len(g)),
            "n_islands":int(g.island_id.nunique()),
            "n_spatial_blocks":int(g.spatial_block.nunique())
        } for name,g in native.groupby("analysis_regime")],
        "n_independent_study_plant_year_pairs_not_visitor_records":True,
        "same_study_same_year_does_not_prove_equivalent_visit_effort":True,
        "no_verified_ancestor_or_counterfactual_colonist_source":True,
        "actual_pollinator_replacement_measured":False,
        "historical_decrease_or_pollen_delivery_not_estimated":True,
        "source_taxa_and_local_sites_kept_private":True,
        "original_island_study_species_names_or_coordinates_exported":False,
    }
    return output
