"""Original GIFT full-botanical-flora GloBI exact-island source-visitor audit.

The 4,505-island sample is immutable. Only source SHA verified + species
accepted + exact independently GSHHG-hash-polygon positive relations count.
GloBI georeferenced records denote documentary metadata, NOT independently
verified observation, pollination service, guild absence or local evolution.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

GUILDS=("bombus","non_bombus_bees","lepidoptera",
        "flower_visiting_birds","diptera")
REGIONS=("GLOBAL_ALL","northern_midlatitude","tropical",
         "southern_extratropical","northern_highlatitude")
MOTIFS=("actinomorphic_deep","zygomorphic_deep","tube_bell_deep",
        "white_tube_bell","red_pink_tube_bell")
SCOPES=("all_recorded","original_native","original_native_endemic")
MODES=("all_cited","largest_archive_excluded")


def genuine_direct_motifs(traits:dict)->dict[str,dict[str,bool]]:
    """SAME accepted species must be scored on BOTH axes; non-recorded
    states remain UNKNOWN rather than a false botanical absence.
    """
    need={"floral_symmetry","floral_form","tube_depth_class",
          "flower_primary_color"}
    if not need<=set(traits):
        raise ValueError("Original Direct flower trait mapping absent")
    s,f,d,c=(traits[k] for k in ("floral_symmetry","floral_form",
             "tube_depth_class","flower_primary_color"))
    out={}
    for name in MOTIFS:
        if name in ("actinomorphic_deep","zygomorphic_deep"):
            set_names=set(s)&set(d)
            target=("actinomorphic" if name=="actinomorphic_deep"
                    else "zygomorphic")
            out[name]={sp:(s[sp]==target and d[sp]=="deep") for sp in set_names}
        elif name=="tube_bell_deep":
            out[name]={sp:(f[sp] in ("tubular","bell_campanulate")
                           and d[sp]=="deep") for sp in set(f)&set(d)}
        else:
            target="white" if name=="white_tube_bell" else "red_pink"
            out[name]={sp:(c[sp]==target and
                        f[sp] in ("tubular","bell_campanulate"))
                       for sp in set(c)&set(f)}
    return out


def strict_source_geography(georeferenced:list[tuple],
                           polygons,flora:pd.DataFrame):
    """Return private exact island×plant×reference positive links.

    Each input tuple has (plant, guild, reference, provider, relation,
    lat, lon). Original verified GSHHG polygons only. Ambiguous boundary
    matches are discarded. Caller must filter multi-provider bibliography.
    """
    import geopandas as gpd
    req={"island_id","geometry"}
    if not req<=set(polygons) or polygons.island_id.duplicated().any():
        raise ValueError("Verified unique GSHHG polygon source required")
    need={"island_id","accepted_species","origin_status","endemic_status"}
    if not need<=set(flora) or flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("Original canonical botanical source ambiguous")
    if not georeferenced:
        return pd.DataFrame(columns=[
            "island_id","accepted_species","origin_status","endemic_status",
            "guild","reference_id","provider_id","evidence_class",
            "latitude","longitude"]),dict(n_original_candidate_points=0,
            n_candidate_points_in_exact_island_polygons=0,
            n_candidate_points_with_multiple_island_match=0,
            n_original_accepted_island_plant_positive_source_rows=0)
    raw=pd.DataFrame.from_records(georeferenced,columns=[
        "accepted_species","guild","reference_id","provider_id",
        "evidence_class","latitude","longitude"])
    if (raw.isna().any().any() or not raw.guild.isin(GUILDS).all()
            or not raw.evidence_class.isin(("flower_visit","pollination")).all()):
        raise ValueError("Unqualified positive original GloBI source relation")
    raw["latitude"]=pd.to_numeric(raw.latitude,errors="coerce")
    raw["longitude"]=pd.to_numeric(raw.longitude,errors="coerce")
    if not (raw.latitude.between(-90,90)&raw.longitude.between(-180,180)).all():
        raise ValueError("Original source location out of range")
    raw=raw.drop_duplicates().reset_index(drop=True)
    points=gpd.GeoDataFrame(raw,
        geometry=gpd.points_from_xy(raw.longitude,raw.latitude),crs="EPSG:4326")
    hit=gpd.sjoin(points,polygons[["island_id","geometry"]],
                  how="inner",predicate="within")
    counts=hit.groupby(level=0).island_id.nunique()
    unique_points=set(counts.index[counts.eq(1)])
    ambiguous=int(counts.gt(1).sum())
    qualified=hit.loc[hit.index.isin(unique_points)].drop(
        columns=["geometry","index_right"],errors="ignore")
    qualified=qualified.merge(
        flora[["island_id","accepted_species","origin_status","endemic_status"]],
        on=["island_id","accepted_species"],how="inner",
        validate="many_to_one")
    qualified=qualified.drop_duplicates([
        "island_id","accepted_species","guild","reference_id",
        "provider_id","evidence_class","latitude","longitude"])
    return pd.DataFrame(qualified),dict(
        n_original_candidate_points=int(len(points)),
        n_candidate_points_in_exact_island_polygons=int(len(unique_points)),
        n_candidate_points_with_multiple_island_match=ambiguous,
        n_original_accepted_island_plant_positive_source_rows=int(len(qualified)),
        georeferenced_public_relation_not_verifiable_same_plot_visit=True,
        uncertain_geographical_sites_never_inferred_absent=True,
    )


def original_source_guild_morphology_availability(
    linked:pd.DataFrame,flora:pd.DataFrame,traits:dict,
    largest_archive_label:str,region:dict,
)->dict:
    """Report only anonymity-preserving positive source opportunity cells.

    All source/island/plant/provider/location IDs must be removed before any
    public JSON. Count independent *bibliographic labels* but do not assert
    those are independent research sites or original datasets.
    """
    if not isinstance(largest_archive_label,str) or not largest_archive_label:
        raise ValueError("Largest source selected globally before outcomes")
    need={"island_id","accepted_species","origin_status","endemic_status",
          "guild","reference_id","provider_id","evidence_class"}
    if not need<=set(linked) or not {"island_id","accepted_species"}<=set(flora):
        raise ValueError("GloBI exact-local plant evidence not certified")
    if (linked.guild.isna().any() or not linked.guild.isin(GUILDS).all()
            or linked[list(need)].isna().any().any()):
        raise ValueError("Original source-linked record missing taxonomy/status")
    if not set(zip(linked.island_id,linked.accepted_species))<=set(zip(
        flora.island_id,flora.accepted_species)):
        raise ValueError("Documentary records not on original true island flora")
    motifs=genuine_direct_motifs(traits)
    mapped=linked.copy()
    mapped["region"]=mapped.island_id.map(
        lambda island:region.get(island,"geography_unresolved"))
    cells=[]
    for mode in MODES:
        source=(mapped if mode=="all_cited" else
                mapped.loc[mapped.provider_id.ne(largest_archive_label)])
        for regime in REGIONS:
            indexed=(source if regime=="GLOBAL_ALL"
                     else source.loc[source.region.eq(regime)])
            for scope in SCOPES:
                if scope=="all_recorded":q=indexed
                elif scope=="original_native":
                    q=indexed.loc[indexed.origin_status.eq("native")]
                else:
                    q=indexed.loc[indexed.origin_status.eq("native") &
                                  indexed.endemic_status.eq("endemic")]
                for motif in MOTIFS:
                    known=motifs[motif]
                    known_names=set(known)
                    q_known=q.loc[q.accepted_species.isin(known_names)].copy()
                    if len(q_known):
                        q_known["source_direct_motif_true"]=q_known.accepted_species.map(known)
                        hits=q_known.loc[q_known.source_direct_motif_true.eq(True)]
                    else:
                        hits=q_known
                    for guild in GUILDS:
                        x=hits.loc[hits.guild.eq(guild)]
                        cells.append(dict(source_mode=mode,region=regime,
                            botanical_status=scope,real_species_direct_pair=motif,
                            original_source_visitor_guild=guild,
                            n_positive_cited_relation_units=int(len(x)),
                            n_documented_islands=int(x.island_id.nunique()),
                            n_documented_accepted_plant_species=int(x.accepted_species.nunique()),
                            n_documented_reference_labels=int(x.reference_id.nunique()),
                            n_documented_archive_first_labels=int(x.provider_id.nunique()),
                            n_explicit_pollination_relations=int(
                                x.evidence_class.eq("pollination").sum()),
                            n_documented_distinct_site_coordinates=int(
                                x[["latitude","longitude"]].drop_duplicates().shape[0])
                              if {"latitude","longitude"}<=set(x) else None,
                            no_recorded_positive_is_NOT_ecological_absence=True))
    return dict(status="FULL_FRAME_GLOBI_EXACT_SAME_ISLAND_SOURCE_MOTIF_AUDIT",
        n_original_parent_GIFT_botanical_floras=4505,
        n_source_positive_original_island_plant_pairs=int(
            linked[["island_id","accepted_species"]].drop_duplicates().shape[0]),
        n_source_original_islands_with_some_positive_relation=int(
            linked.island_id.nunique()),
        n_original_accepted_species_with_positive_reference=int(
            linked.accepted_species.nunique()),
        n_positive_cited_source_rows=int(len(linked)),
        n_original_source_bibliographic_keys=int(linked.reference_id.nunique()),
        n_original_source_provider_archive_first_labels=int(linked.provider_id.nunique()),
        cells=cells,
        original_source_provider_omission_not_independent_replication=True,
        cited_document_is_NOT_observed_effective_pollination=True,
        original_plant_traits_species_wide_not_site_measured=True,
        no_documented_guild_not_a_biological_nondetection=True)
