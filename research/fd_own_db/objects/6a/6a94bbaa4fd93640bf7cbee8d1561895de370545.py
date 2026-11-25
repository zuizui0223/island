"""EuPPollNet source original study/site/date × exact GIFT island plant bridge.

Positive European field interaction records are joined ONLY through source-
verified WGS84 GSHHG original polygon IDs and exact botanically accepted
binomials. Every original study, island, plant, coordinate, network and
original reference remains inside the protected GIFT workspace.
No source network nonedge becomes guild absence; site survey effort still
needs taxon-specific auditing for causal choice/effectiveness inference.
"""
from __future__ import annotations

import re
import numpy as np
import pandas as pd

from .independent_pollination_networks import original_binomial,visitor_guild
from .fullframe_globi_local_support import genuine_direct_motifs

REQUIRED=("Study_id","Network_id","Latitude","Longitude","Date",
          "Plant_accepted_name","Pollinator_accepted_name",
          "Pollinator_order","Pollinator_family")
REGIONS=("GLOBAL_ALL","northern_midlatitude","tropical",
         "southern_extratropical","northern_highlatitude",
         "geography_unresolved")
MOTIFS=("actinomorphic_deep","zygomorphic_deep","tube_bell_deep",
        "white_tube_bell","red_pink_tube_bell")
CATEGORIES=("bombus","other_bees","butterfly_family",
            "sphingidae_moth","other_lepidoptera",
            "flies_order","other_animals","unresolved")


def strict_original_accepted_plant(name:object)->str|None:
    """Do not reinterpret genus-only or infraspecific labels as source species.

    Botanical authorities after genuine accepted binomial are allowed, but
    subsp./var./f./cf./aff. qualifiers are NOT promoted to an observed
    accepted species-level visitor relationship without an original match.
    """
    if not isinstance(name,str) or not name.strip():
        return None
    tokens=name.strip().lower().replace(",", " ").split()
    if set(tokens) & {"subsp.","subsp","ssp.","ssp","var.","var",
                      "f.","forma","cf.","cf","aff.","aff","sp.","spp."}:
        return None
    return original_binomial(name)


def original_guild(order:str,family:str,animal:str)->str:
    raw=visitor_guild(order,family)
    animal=str(animal or "").strip()
    family=str(family or "").strip().lower()
    if raw=="bees_family_supported":
        return "bombus" if (animal.startswith("Bombus ") and
                           family=="apidae") else "other_bees"
    if raw=="butterfly_family_supported":return "butterfly_family"
    if raw=="sphingidae_moth":return "sphingidae_moth"
    if raw=="other_lepidoptera":return "other_lepidoptera"
    if raw=="flies_order_supported":return "flies_order"
    if raw=="unresolved_visitor_taxonomy":return "unresolved"
    return "other_animals"


def exact_source_network_site_union(source:pd.DataFrame,polygons):
    """Geographic original source site identity gate, no fuzzy island match."""
    import geopandas as gpd
    if not set(REQUIRED)<=set(source):
        raise ValueError("European original study-site positive interaction schema missing")
    if not {"island_id","geometry"}<=set(polygons):
        raise ValueError("Only GIFT verified original GSHHG polygon IDs allowed")
    if polygons.island_id.duplicated().any():
        raise ValueError("Ambiguous original geographic island source IDs")
    key=["Study_id","Network_id","Latitude","Longitude"]
    if source[key].isna().any().any():
        raise ValueError("Original public study/network/coordinate fields missing")
    network=source[key].drop_duplicates().copy()
    network["Latitude"]=pd.to_numeric(network.Latitude,errors="coerce")
    network["Longitude"]=pd.to_numeric(network.Longitude,errors="coerce")
    valid=network.Latitude.between(-90,90)&network.Longitude.between(-180,180)
    bad=int((~valid).sum())
    network=network.loc[valid].copy()
    unique=network[["Latitude","Longitude"]].drop_duplicates()
    pts=gpd.GeoDataFrame(
        unique,geometry=gpd.points_from_xy(unique.Longitude,unique.Latitude),
        crs="EPSG:4326")
    # The same GPS pair has one island or remains unresolved; no nearest guess.
    located=gpd.sjoin(pts,polygons[["island_id","geometry"]],
                      how="inner",predicate="within")
    bypoint=located.groupby(level=0).island_id.nunique()
    unique_points=set(bypoint.index[bypoint.eq(1)])
    matched=located.loc[located.index.isin(unique_points),
                        ["Latitude","Longitude","island_id"]].drop_duplicates()
    if matched[["Latitude","Longitude"]].duplicated().any():
        raise ValueError("One EuPPollNet source point matches multiple islands")
    network=network.merge(matched,on=["Latitude","Longitude"],
                          how="inner",validate="many_to_one")
    return network,dict(
        n_original_network_site_coordinate_combinations=int(
            len(source[key].drop_duplicates())),
        n_original_site_coordinate_combinations_missing_or_invalid=int(bad),
        n_original_distinct_source_coordinate_pairs=int(len(unique)),
        n_original_source_coordinate_pairs_inside_unique_GIFT_polygons=int(len(matched)),
        n_original_source_study_network_site_units_on_verified_GIFT_polygon=int(len(network)),
        n_source_original_study_ids_with_exact_GIFT_polygon=int(
            network.Study_id.nunique()),
        n_source_animal_site_polygon_ambiguous_geographic_points=int(
            bypoint.gt(1).sum()),
        no_unverified_nearest_island_substitution=True,
        source_site_positive_relation_not_pollination_effectiveness=True,
    )


def evaluate_original_european_island_visits(
    source:pd.DataFrame,polygons,flora:pd.DataFrame,traits:dict,
    regime:dict,*,original_parent_n:int=4505,
)->dict:
    if original_parent_n!=4505 or flora.island_id.nunique()!=4505:
        raise ValueError("Original full global GIFT botanical study population changed")
    need={"island_id","accepted_species","origin_status","endemic_status"}
    if not need<=set(flora) or flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original accepted GIFT flora source IDs not unique")
    if not set(REQUIRED)<=set(source):
        raise ValueError("Original accepted field-network source columns missing")
    # Keep WGS84 identity exactly aligned on source records and unique
    # network-site table; published source CSV coordinates are read as text.
    source=source.copy()
    source["Latitude"]=pd.to_numeric(source.Latitude,errors="coerce")
    source["Longitude"]=pd.to_numeric(source.Longitude,errors="coerce")
    valid=(source.Latitude.between(-90,90)&
           source.Longitude.between(-180,180))
    n_rejected=int((~valid).sum())
    source=source.loc[valid].copy()
    network,site_report=exact_source_network_site_union(source,polygons)
    site_report["n_original_source_rows_with_invalid_position_excluded"]=n_rejected
    if network.empty:
        return dict(status="SOURCE_FIELD_NETWORK_ISLAND_OPPORTUNITY_UNIDENTIFIABLE",
            n_original_GIFT_flora_population=4505,
            source_verified_site_opportunity=site_report,
            n_original_same_island_accepted_flower_visitor_units=0,
            positive_source_absence_not_biological_no_visitation=True)
    linked=source.merge(network,on=[
        "Study_id","Network_id","Latitude","Longitude"],
        how="inner",validate="many_to_one")
    normalize={}
    for value in linked.Plant_accepted_name.dropna().unique():
        name=strict_original_accepted_plant(value)
        normalize[value]=name
    linked["plant"]=linked.Plant_accepted_name.map(normalize)
    # No subspecific imputation, invalid name or synonym crosswalk.
    linked=linked.loc[linked.plant.notna()].copy()
    linked=linked.merge(
        flora[["island_id","accepted_species","origin_status",
               "endemic_status"]].rename(columns={"accepted_species":"plant"}),
        on=["island_id","plant"],how="inner",validate="many_to_one")
    linked["guild"]=[original_guild(o,f,a) for o,f,a in zip(
        linked.Pollinator_order,linked.Pollinator_family,
        linked.Pollinator_accepted_name)]
    linked["region"]=linked.island_id.map(
        lambda x:regime.get(x,"geography_unresolved"))
    linked["year"]=linked.Date.astype(str).str.slice(0,4)
    if len(linked) and not linked.year.str.fullmatch(r"\d{4}").all():
        raise ValueError("Original study date incomplete: do not manufacture year")
    # A row is a positive source edge/event. Multiple identical recording
    # interactions become ONE documented positive study×network×site×year
    # botanical plant×animal guild opportunity, not independent animal tests.
    key=["Study_id","Network_id","Latitude","Longitude","year",
         "island_id","plant","guild","origin_status","endemic_status","region"]
    events=linked[key].drop_duplicates()
    motifs=genuine_direct_motifs(traits)
    cells=[]
    for r in REGIONS:
        rr=events if r=="GLOBAL_ALL" else events.loc[events.region.eq(r)]
        for scope in ("all_recorded","original_native",
                      "original_native_endemic"):
            if scope=="all_recorded":ss=rr
            elif scope=="original_native":
                ss=rr.loc[rr.origin_status.eq("native")]
            else:ss=rr.loc[rr.origin_status.eq("native")&
                          rr.endemic_status.eq("endemic")]
            for motif in MOTIFS:
                annotated=motifs[motif]
                yes=set(sp for sp,val in annotated.items() if val)
                known=set(annotated)
                known_ss=ss.loc[ss.plant.isin(known)]
                with_trait=known_ss.loc[known_ss.plant.isin(yes)]
                for guild in CATEGORIES:
                    q=with_trait.loc[with_trait.guild.eq(guild)]
                    cells.append(dict(region=r,botanical_status=scope,
                        same_accepted_species_direct_flower_pair=motif,
                        original_european_source_guild=guild,
                        n_source_studies_with_positive_motif_visits=int(
                            q.Study_id.nunique()),
                        n_source_study_network_site_year_positive_units=int(len(q)),
                        n_source_independently_verified_GIFT_islands=int(
                            q.island_id.nunique()),
                        n_source_accepted_botanical_plant_species=int(
                            q.plant.nunique()),
                        positive_source_only_not_measured_effective_pollen=True))
    return dict(
        status="SOURCE_EUPPOLLNET_EUROPEAN_FIELD_SITE_EXACT_GIFT_FLORA_OPPORTUNITY",
        n_original_global_GIFT_botanical_floras=4505,
        n_original_source_european_rows=len(source),
        n_original_source_european_study_ids=int(source.Study_id.nunique()),
        n_original_source_study_network_date_positive_units=int(
            source[["Study_id","Network_id","Date"]].drop_duplicates().shape[0]),
        source_exact_site_polygon_eligibility=site_report,
        source_verified_site_opportunity=site_report,
        n_original_source_positive_recorded_rows_matching_same_island_GIFT_plant=int(len(linked)),
        n_original_source_distinct_islands_with_matched_flora_plant=int(
            linked.island_id.nunique()),
        n_original_source_study_ids_with_matched_flora_plant=int(
            linked.Study_id.nunique()),
        n_original_accepted_botanical_species_with_source_island_visits=int(
            linked.plant.nunique()),
        n_original_unique_study_network_site_year_plant_guild_positive_units=int(
            len(events)),
        true_source_direct_flower_and_animal_guild_documentation=cells,
        no_source_interaction_edge_implicit_negative=True,
        source_field_plant_phenotype_is_species_wide_NOT_measured_on_site=True,
        source_comparison_only_european_not_global_tropical=True,
        source_original_study_overlap_with_GloBI_WOL_DoPI_not_checked=True,
        no_true_effective_pollen_delivery_or_seed_set=True,
        source_real_guild_choices_not_causal_selection=True,
        private_island_study_plant_and_site_ids_never_in_aggregate=True,
    )
