"""Independent plant–visitor network screening, without pseudo-absence.

Web of Life is a publication-curated source of network sites and sampled
plant/animal edges. Mangal archives heterogeneous interaction networks.
NEITHER is automatically an independent underlying study from GloBI, nor a
measure of effective pollen deposition/seed set. Source sample metadata and
original 4,505-island GIFT botanicals remain separate.

This module contains only pure parsing/aggregation; run original sources in
the island protected Actions workflow. No original private GIFT island IDs,
plant identities or coordinates are returned in public reports.
"""
from __future__ import annotations

import math
import re
from collections import defaultdict

import pandas as pd

MORPHOLOGY=("actinomorphic_deep","zygomorphic_deep","tube_bell_deep")
GUILDS=("bird_family_supported","butterfly_family_supported",
        "sphingidae_moth","other_lepidoptera","bees_family_supported",
        "flies_order_supported","unresolved_visitor_taxonomy")
BUTTERFLY_FAMILIES=frozenset({
    "nymphalidae","papilionidae","pieridae","lycaenidae",
    "hesperiidae","riodinidae",
})
BEE_FAMILIES=frozenset({
    "apidae","halictidae","andrenidae","colletidae",
    "megachilidae","melittidae","stenotritidae",
})
BIRD_FAMILIES=frozenset({
    "trochilidae","nectariniidae","meliphagidae",
    "zosteropidae","promeropidae","drepanididae",
    "dicaeidae","philepittidae",
})


def original_binomial(name:object)->str|None:
    """Exact unambiguous genus+epithet: no automatic taxon synonym join."""
    if name is None:return None
    q=str(name).strip().replace("_"," ")
    m=re.match(r"^([A-Z][a-zA-Z-]+)\s+([a-z][a-z-]+)(?:\b|$)",q)
    if not m:return None
    genus,epi=m.group(1),m.group(2)
    if epi in ("sp","spp","cf","aff"):return None
    return genus+" "+epi


def visitor_guild(order:object=None,family:object=None)->str:
    o=str(order or "").strip().lower()
    f=str(family or "").strip().lower()
    if f in BIRD_FAMILIES:return "bird_family_supported"
    if f in BUTTERFLY_FAMILIES:return "butterfly_family_supported"
    if f=="sphingidae":return "sphingidae_moth"
    if o=="lepidoptera":return "other_lepidoptera"
    if f in BEE_FAMILIES:return "bees_family_supported"
    if o=="diptera":return "flies_order_supported"
    return "unresolved_visitor_taxonomy"


def network_location_metadata(raw:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    """Web of Life metadata: accept published network IDs and explicit GPS."""
    if raw.empty:
        raise ValueError("Web of Life network catalogue empty")
    cols=set(raw.columns)
    ids=next((c for c in ("network","network_name","id") if c in cols),None)
    kind=next((c for c in ("interaction_type","network_type") if c in cols),None)
    if not ids or "latitude" not in cols or "longitude" not in cols:
        raise ValueError("Web of Life site coordinates and network ID absent")
    d=raw.copy()
    d["network_id"]=d[ids].astype(str).str.strip()
    d["lat"]=pd.to_numeric(d.latitude,errors="coerce")
    d["lon"]=pd.to_numeric(d.longitude,errors="coerce")
    poll=(d[kind].astype(str).str.contains("pollina",case=False,regex=False)
          if kind else pd.Series(False,index=d.index))
    # Historical M_PL convention independently documented by Web of Life.
    poll|=d.network_id.str.match(r"^M_PL_\d+$")
    selected=d.loc[poll].copy()
    valid=(selected.lat.between(-90,90)&selected.lon.between(-180,180)&
           selected.network_id.ne("")&
           selected.lat.notna()&selected.lon.notna())
    out=selected.loc[valid].copy()
    if out.network_id.duplicated().any():
        raise ValueError("same network published with conflicting metadata")
    return out,dict(
        n_all_web_of_life_network_records=int(len(raw)),
        n_declared_pollination_networks=int(len(selected)),
        n_pollination_networks_with_coordinates=int(len(out)),
        n_pollination_networks_without_supported_coordinates=int((~valid).sum()),
        n_unique_reference_metadata=int(
            out.reference.dropna().nunique()) if "reference" in out else None,
        site_coordinate_precision_provider_supplied_not_independent_GSHHG=True,
    )


def mangal_location_metadata(raw_networks:list[dict],datasets:list[dict]
    )->tuple[pd.DataFrame,dict]:
    byid={str(d.get("id")):d for d in datasets}
    points=[]
    candidates=0
    for network in raw_networks:
        descriptor=byid.get(str(network.get("dataset_id")),{})
        description=(" ".join(str(x or "") for x in [
            descriptor.get("name"),descriptor.get("description"),
            network.get("name"),network.get("description")
        ])).lower()
        # This is only a search heuristic; actual node/edge type MUST be
        # checked before calling the network a flower-visitation source.
        if not re.search(r"pollinat|flower.visitor|flower.visiting|"
                         r"nectar|anthophil|bees?",description):
            continue
        candidates+=1
        geom=network.get("geom") or {}
        if not isinstance(geom,dict) or geom.get("type")!="Point":continue
        coords=geom.get("coordinates") or []
        if len(coords)!=2:continue
        try:
            lon,lat=float(coords[0]),float(coords[1])
        except (TypeError,ValueError):continue
        if not (-180<=lon<=180 and -90<=lat<=90 and
                math.isfinite(lat) and math.isfinite(lon)):continue
        points.append(dict(network_id=str(network.get("id")),
                           dataset_id=str(network.get("dataset_id")),
                           lat=lat,lon=lon,
                           source_reference_id=str(descriptor.get("ref_id"))))
    q=pd.DataFrame.from_records(points,columns=[
        "network_id","dataset_id","lat","lon","source_reference_id"
    ])
    if not q.empty and q.network_id.duplicated().any():
        raise ValueError("duplicate Mangal network metadata ID")
    return q,dict(
        n_mangal_datasets_retrieved=len(datasets),
        n_mangal_networks_retrieved=len(raw_networks),
        n_mangal_keyword_candidate_networks=candidates,
        n_mangal_candidate_georeferenced_points=len(q),
        source_interaction_types_and_nodes_not_yet_verified=True,
    )


def _plants_role(v:object)->bool:
    q=str(v or "").lower().strip()
    return q in ("plant","plants","plantae","flower","flowers")


def identify_actual_plant_edges(edges:pd.DataFrame,
                                species:pd.DataFrame)->tuple[pd.DataFrame,dict]:
    """Only explicit plant-role edges; never infer a visit from matrix zeros.

    Web of Life edge CSV documents row/column species and row/column roles.
    Species CSV supplies taxonomy for each endpoint, when available. Never
    guess taxonomy if no family/order supplied by original curator.
    """
    fields={"row_species","column_species","row_role","column_role"}
    if missing:=fields-set(edges):
        raise ValueError(f"Web of Life edge role contract missing {missing}")
    if "species" not in species and "name" not in species:
        raise ValueError("Web of Life species taxonomy lacks species name")
    name_col="species" if "species" in species else "name"
    taxonomy={}
    for _,r in species.iterrows():
        key=str(r[name_col]).strip()
        if key in taxonomy:continue
        taxonomy[key]=(
            r.get("order",""),r.get("family","")
        )
    records=[]
    rejected=0
    for row in edges.itertuples(index=False):
        x,y=str(row.row_species).strip(),str(row.column_species).strip()
        left=_plants_role(row.row_role)
        right=_plants_role(row.column_role)
        if left==right:
            rejected+=1
            continue
        plant,animal=(x,y) if left else (y,x)
        binomial=original_binomial(plant)
        if binomial is None:
            rejected+=1
            continue
        order,family=taxonomy.get(animal,("",""))
        records.append(dict(
            plant=binomial,
            guild=visitor_guild(order,family),
            visitor_taxon_name=animal,
            source_original_plant_name=plant,
            evidence_type="network_source_recorded_edge",
        ))
    result=pd.DataFrame(records,columns=[
        "plant","guild","visitor_taxon_name",
        "source_original_plant_name","evidence_type",
    ])
    return result,dict(
        n_original_source_edge_records=int(len(edges)),
        n_explicit_plant_role_binomial_edge_records=int(len(result)),
        n_unresolved_plant_role_or_name=int(rejected),
        n_taxonomy_classified_guild_edge_records=int(
            result.guild.ne("unresolved_visitor_taxonomy").sum()
        ) if len(result) else 0,
        empty_edges_does_not_mean_no_pollination=True,
    )


def real_same_species_morphologies(plant:str,traits:dict)->dict[str,bool]|None:
    a=traits.get("floral_symmetry",{}).get(plant)
    b=traits.get("tube_depth_class",{}).get(plant)
    form=traits.get("floral_form",{}).get(plant)
    if b is None:
        return None
    return {
        "actinomorphic_deep":a=="actinomorphic" and b=="deep"
          if a is not None else None,
        "zygomorphic_deep":a=="zygomorphic" and b=="deep"
          if a is not None else None,
        "tube_bell_deep":form in ("tubular","bell_campanulate") and b=="deep"
          if form is not None else None,
    }


def source_bounded_site_plant_audit(
    eligible_network_edges:pd.DataFrame,
    original_flora:pd.DataFrame,
    traits:dict,*,global_original_islands:int=4505
)->dict:
    """One network/plant/guild per source edge, citations as network IDs.

    eligible_network_edges are PRIVATE already geospatially matched to
    published original GSHHG island polygons; do not return these IDs.
    """
    if global_original_islands!=4505:
        raise ValueError("full GIFT 4505 botanical frame reduced")
    req={"network_id","island_id","plant","guild"}
    if not req <= set(eligible_network_edges):
        raise ValueError("network polygon/name matched edge fields missing")
    if not {"island_id","accepted_species","origin_status","endemic_status"}<=set(original_flora):
        raise ValueError("original GIFT flora native status missing")
    if original_flora.duplicated(["island_id","accepted_species"]).any():
        raise ValueError("original GIFT accepted flora duplicated")
    names=original_flora.rename(columns={"accepted_species":"plant"})
    merged=eligible_network_edges.merge(
        names[["island_id","plant","origin_status","endemic_status"]],
        on=["island_id","plant"],how="inner",validate="many_to_one"
    ).drop_duplicates(["network_id","island_id","plant","guild"])
    if not merged.empty:
        for k in MORPHOLOGY:
            merged[k]=merged.plant.map(
                lambda s:real_same_species_morphologies(s,traits)[k]
                if real_same_species_morphologies(s,traits) is not None
                else None
            )
    cells=[]
    for scope in ("all_recorded","native_label","endemic_label","origin_unresolved"):
        if scope=="all_recorded":sub=merged
        elif scope=="native_label":sub=merged.loc[merged.origin_status.eq("native")]
        elif scope=="endemic_label":sub=merged.loc[merged.origin_status.eq("native") & merged.endemic_status.eq("endemic")]
        else:sub=merged.loc[merged.origin_status.eq("unresolved")]
        for shape in MORPHOLOGY:
            for guild in GUILDS:
                x=sub.loc[sub[shape].eq(True)&sub.guild.eq(guild)] if not sub.empty else sub
                cells.append(dict(
                    botanical_scope=scope,verified_real_joint_morphology=shape,
                    source_network_visitor_taxonomy=guild,
                    n_verified_site_networks_with_cited_plant_and_visitor=int(
                        x.network_id.nunique()) if len(x) else 0,
                    n_original_islands_with_matching_source_network=int(
                        x.island_id.nunique()) if len(x) else 0,
                    n_plant_guild_network_documented_units=int(len(x)),
                    effective_pollen_transfer_not_measured=True,
                    no_network_edge_not_a_biological_absence=True,
                ))
    return dict(
        status="SOURCE_NETWORK_EXACT_SITE_SAME_SPECIES_FLORAL_TRAIT_SUPPORT",
        n_original_global_GIFT_botanical_islands=4505,
        n_original_external_networks_with_plant_guild_overlap=int(
            merged.network_id.nunique()),
        n_original_external_network_islands_with_plant_guild_overlap=int(
            merged.island_id.nunique()),
        n_original_species_network_plant_guild_units=int(len(merged)),
        cells=cells,
        limitations=dict(
            original_full_4505_GIFT_botanical_universe_retained=True,
            candidate_independent_aggregator_not_verified_independent_original_study=True,
            true_species_level_exact_original_flora_plant_match_only=True,
            real_two_trait_profiles_never_imputed_from_marginals=True,
            source_visitor_guild_taxonomy_requires_family_or_order=True,
            only_original_spatially_verified_island_polygons_count_as_site=True,
            bird_butterfly_sphingidae_subgroups_taxonomic_not_effective_service=True,
            network_nonedge_never_recoded_as_true_absence=True,
            time_matched_mainland_visit_and_seed_set_unavailable=True,
        ),
    )
