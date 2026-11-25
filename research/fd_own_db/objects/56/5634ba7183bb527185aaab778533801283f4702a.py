"""Source-pinned EuPPollNet v1.3.0 × original island floral repertoire gate.

The source is *actual documented plant–visitor interactions*, not pollen
delivery or a historical loss experiment. Only original fixed island polygons
and exact botanical species names are joined. Aggregates suppress site,
network, provider, plant, and pollinator identifiers.
"""
from __future__ import annotations

from collections import defaultdict
import numpy as np
import pandas as pd

SOURCE_COLUMNS = (
    "Study_id", "Network_id", "Date", "Latitude", "Longitude",
    "Plant_accepted_name", "Plant_rank", "Pollinator_accepted_name",
    "Pollinator_order", "Pollinator_genus", "Pollinator_family"
)
SOURCE_GROUPS = ("Hymenoptera", "Diptera", "Lepidoptera", "Coleoptera")


def validate_original_visits(frame: pd.DataFrame) -> pd.DataFrame:
    if not set(SOURCE_COLUMNS).issubset(frame):
        raise ValueError("EuPPollNet v1.3.0 original interaction header changed")
    x = frame[list(SOURCE_COLUMNS)].copy()
    x["Latitude"] = pd.to_numeric(x.Latitude, errors="coerce")
    x["Longitude"] = pd.to_numeric(x.Longitude, errors="coerce")
    x["Date"] = pd.to_datetime(x.Date, errors="coerce")
    x = x.loc[
       x.Latitude.between(-90, 90) & x.Longitude.between(-180, 180) &
       x.Date.notna() &
       x.Study_id.notna() & x.Network_id.notna() &
       x.Plant_accepted_name.notna() &
       x.Pollinator_accepted_name.notna()
    ].copy()
    if x.empty:
        raise ValueError("No original dated/site georeferenced plant visits")
    x["Year"] = x.Date.dt.year
    x["Study_id"] = x.Study_id.astype(str).str.strip()
    x["Network_id"] = x.Network_id.astype(str).str.strip()
    x["Plant_accepted_name"] = x.Plant_accepted_name.astype(str).str.strip()
    x["Pollinator_order"] = x.Pollinator_order.fillna("unresolved").astype(str).str.strip()
    x["Pollinator_genus"] = x.Pollinator_genus.fillna("").astype(str).str.strip()
    x["Pollinator_family"] = x.Pollinator_family.fillna("").astype(str).str.strip()
    if (x.Study_id.eq("") | x.Network_id.eq("") |
        x.Plant_accepted_name.eq("")).any():
        raise ValueError("Blank original study, sampling network, or plant")
    return x


def _guild(row: pd.Series) -> str:
    order = row["Pollinator_order"]
    genus = row["Pollinator_genus"]
    family = row["Pollinator_family"]
    if genus == "Bombus":
        return "bombus"
    if order == "Hymenoptera" and family in {
        "Apidae", "Andrenidae", "Colletidae", "Megachilidae",
        "Halictidae", "Melittidae", "Stenotritidae"
    }:
        return "other_bees"
    if order == "Diptera":
        return "diptera"
    if order == "Lepidoptera":
        return "lepidoptera_pooled"
    if order == "Coleoptera":
        return "coleoptera"
    return "unresolved_or_other"


def audit_original_island_euppollnet(
    interaction: pd.DataFrame,
    site_to_island: pd.DataFrame,
    eligible_928_flora: pd.DataFrame,
    island_geo: pd.DataFrame,
    *, require_min_islands: int = 10,
    require_min_blocks: int = 5,
) -> dict:
    """All returned records are coarse anonymous sums, never original IDs."""
    if require_min_islands != 10 or require_min_blocks != 5:
        raise ValueError("Prespecified no-claim support gate changed")
    x = validate_original_visits(interaction)
    keys = ("Study_id", "Network_id", "Latitude", "Longitude")
    if not set((*keys, "island_id")).issubset(site_to_island):
        raise ValueError("Original polygon site match contract missing")
    sites = site_to_island[list((*keys, "island_id"))].copy()
    if sites[list(keys)].duplicated().any() or sites.isna().any().any():
        raise ValueError("One interaction site must match no more than one island")
    if not {"island_id", "accepted_species", "origin_status"}.issubset(eligible_928_flora):
        raise ValueError("Original Direct island flora data missing")
    if eligible_928_flora.duplicated(["island_id", "accepted_species"]).any():
        raise ValueError("Original flora is not one island × accepted plant species")
    if not {"island_id", "analysis_regime", "spatial_block"}.issubset(island_geo):
        raise ValueError("Original geographic block join missing")
    if island_geo.island_id.duplicated().any():
        raise ValueError("Duplicated original physical island")
    # Sites are assigned by actual polygon spatial join, not centroids/buffers.
    matched = x.merge(sites, on=list(keys), how="inner", validate="many_to_one")
    eligible = set(eligible_928_flora.island_id)
    matched = matched.loc[matched.island_id.isin(eligible)].copy()
    if len(matched) == 0:
        return {
            "status": "NO_SOURCE_VERIFIED_FD_ISLAND_NETWORK_INTERACTIONS",
            "original_euppollnet_rows_with_date_and_coordinates": int(len(x)),
            "original_euppollnet_geosites": int(x[list(keys)].drop_duplicates().shape[0]),
            "n_original_928_flora_islands_with_exact_network_site": 0,
            "n_same_flora_plant_source_visits": 0,
            "source_island_pollinator_shift_or_loss_identified": False,
            "original_taxon_network_island_or_coordinates_released": False,
        }
    flora = eligible_928_flora[["island_id", "accepted_species", "origin_status"]].copy()
    flora["origin_status"] = flora.origin_status.astype(str)
    matched = matched.merge(
        flora.rename(columns={"accepted_species":"Plant_accepted_name"}),
        on=["island_id", "Plant_accepted_name"],
        how="left", validate="many_to_one", indicator="floral_species_join"
    )
    matched = matched.merge(
        island_geo[["island_id", "analysis_regime", "spatial_block"]],
        on="island_id", how="left", validate="many_to_one"
    )
    if matched.analysis_regime.isna().any() or matched.spatial_block.isna().any():
        raise ValueError("Source geography missing for original island")
    matched["exact_plant_on_flora"] = matched.floral_species_join.eq("both")
    matched["strict_native_direct_plant"] = (
         matched.exact_plant_on_flora & matched.origin_status.eq("native")
    )
    matched["documented_guild"] = matched.apply(_guild, axis=1)
    # Both number of years and exact site reuse are distinct from number of
    # interactions. Network ID may represent an entire sampling site not visit.
    group_year = matched.groupby(
        ["island_id", "Study_id", "Network_id", "Latitude", "Longitude"]
    ).Year.nunique()
    repeated_site_years = group_year[group_year >= 2]
    networks = matched[["island_id","Study_id","Network_id","Year"]].drop_duplicates()
    island_network_summary = matched.groupby("island_id").agg(
        interaction_rows=("Study_id", "size"),
        studies=("Study_id", "nunique"),
        networks=("Network_id", "nunique"),
        years=("Year", "nunique"),
        same_flora=("exact_plant_on_flora", "sum"),
        same_native=("strict_native_direct_plant", "sum")
    )
    per_regime=[]
    for regime, v in matched.groupby("analysis_regime"):
        i=v.island_id.unique()
        blocks=v.spatial_block.nunique()
        admissible=(len(i)>=require_min_islands and
                    blocks>=require_min_blocks and v.Study_id.nunique()>=3)
        rec = dict(
            regime=str(regime),
            n_original_eligible_islands_with_visitor_network=int(len(i)),
            n_spatial_blocks=int(blocks),
            n_distinct_original_studies=int(v.Study_id.nunique()),
            n_distinct_study_network_year_units=int(v[
                 ["Study_id","Network_id","Year"]].drop_duplicates().shape[0]),
            n_observed_interaction_rows=int(len(v)),
            n_direct_species_joined_interaction_rows=int(v.exact_plant_on_flora.sum()),
            n_strict_native_direct_plant_interaction_rows=int(v.strict_native_direct_plant.sum()),
            n_distinct_original_plant_names_visiting_island=int(v.Plant_accepted_name.nunique()),
            n_exact_match_to_direct_species=int(v.loc[
                v.exact_plant_on_flora,"Plant_accepted_name"].nunique()),
            n_island_sites_repeated_across_years=int(
                len(set(repeated_site_years.index.get_level_values("island_id")) & set(i))),
            family_or_study_disjoint_network_comparison_eligible=bool(admissible),
            abundance_or_effective_pollen_delivery_measured=False
        )
        for guild in ("bombus","other_bees","diptera","lepidoptera_pooled","coleoptera"):
            rec["n_documented_"+guild+"_interaction_rows"]=int(
               (v.documented_guild == guild).sum())
        per_regime.append(rec)
    return {
        "status":"ORIGINAL_EUPPOLLNET_ISLAND_VISITOR_NETWORK_OPPORTUNITY_ONLY",
        "source_version":"v1.3.0",
        "source_revision_has_duplicates_corrected":True,
        "original_euppollnet_rows_with_date_and_coordinates":int(len(x)),
        "original_euppollnet_geosites":int(x[list(keys)].drop_duplicates().shape[0]),
        "n_original_928_flora_islands_with_exact_network_site":int(matched.island_id.nunique()),
        "n_original_island_study_year_network_units":int(networks.shape[0]),
        "n_direct_species_joined_interaction_rows":int(matched.exact_plant_on_flora.sum()),
        "n_strict_native_direct_species_interactions":int(matched.strict_native_direct_plant.sum()),
        "n_study_network_coordinate_sites_with_two_years":int(len(repeated_site_years)),
        "n_islands_with_two_years_site":int(len(set(
            repeated_site_years.index.get_level_values("island_id")))),
        "n_regions_meeting_10_islands_5blocks_3studies":int(sum(
            q["family_or_study_disjoint_network_comparison_eligible"] for q in per_regime)),
        "regions":per_regime,
        "source_island_pollinator_shift_or_loss_identified":False,
        "actual_visited_flower_not_pollen_delivery":True,
        "study_taxonomic_sampling_differences_not_pollinator_absence":True,
        "time_matched_flora_dynamics_measured":False,
        "original_taxon_network_island_or_coordinates_released":False,
    }
