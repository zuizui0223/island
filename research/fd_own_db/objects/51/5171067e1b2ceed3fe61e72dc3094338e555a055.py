"""Extend ORIGINAL GIFT island / Web of Life positive network evidence.

This geometry screen is a per-island EXACT-HASH intersection with a
SHA-pinned GSHHG GeoPackage, NOT a claim of global (8,265-island) geometry
equivalence: 9 historical original IDs were not reproducible. No position,
original island, accepted plant or network ID can be exported.
"""
from __future__ import annotations
import re
import pandas as pd

ORIGINAL_GIFT_N=4505
HISTORICAL_COMPLETE_N=928
HISTORICAL_MATCHED_N=927
SOURCE_GSHHG_SHA="8dbbe7e071e77e9e75f2d639239099ebca8d5c16d6a07df8169729d49f15cf41"
SOURCE_PREFIX="gshhg_2.3.7_h_"
HASH=re.compile(r"^[0-9a-f]{20}$")


def full_botanical_exact_geometry_gate(
    gpkg, flora, historical_928:set[str], receipt:dict,
):
    """Return source-private verified GeoDataFrame plus anonymous support counts.

    The intersection is based on the SAME original cryptographic island ID,
    and a second SHA-pinned geometry file. The fact that 9 IDs elsewhere
    failed historical rebuilding is preserved, not filled by centroids.
    """
    req={"island_id","geometry_sha256","geometry"}
    if not req<=set(gpkg):
        raise ValueError("GSHHG source polygons must contain real hash identities")
    if not {"island_id","accepted_species"}<=set(flora):
        raise ValueError("Original accepted GIFT species roster missing")
    if flora.island_id.nunique()!=ORIGINAL_GIFT_N or not len(historical_928)==HISTORICAL_COMPLETE_N:
        raise ValueError("Global GIFT or historical phenotype cohort changed")
    if gpkg.island_id.duplicated().any() or gpkg.island_id.isna().any():
        raise ValueError("Duplicate GSHHG original island identity")
    if not receipt.get("GSHHG_archive_SHA256")==SOURCE_GSHHG_SHA:
        raise ValueError("Original archive has been replaced")
    if not (receipt.get("GSHHG_version")=="2.3.7"
            and receipt.get("GSHHG_resolution")=="h"
            and receipt.get("n_original_focal_islands")==21
            and receipt.get("n_focal_islands_with_exact_historical_polygon")==21):
        raise ValueError("Historical GSHHG source provenance inconsistent")
    # This is a deliberate NEW per-ID audit. Never claim all 8,265 historical
    # polygons are exact if the archived full-universe receipt says otherwise.
    if receipt.get("global_full_universe_identity_verified") is not False:
        raise ValueError("Unexpected global GSHHG historical ID status")
    original=set(flora.island_id.astype(str))
    got=gpkg.island_id.astype(str)
    hashed=gpkg.geometry_sha256.astype(str)
    valid=got.str.startswith(SOURCE_PREFIX) & hashed.str.fullmatch(HASH.pattern)
    valid &= got.str.slice(len(SOURCE_PREFIX)).eq(hashed)
    # No wrong geometry hash is permitted among candidate botanical islands.
    if (~valid & got.isin(original)).any():
        raise ValueError("A plant-containing original polygon has incompatible source ID and hash")
    compatible=gpkg.loc[valid & got.isin(original)].copy()
    historical_covered=set(compatible.island_id)&set(historical_928)
    if len(historical_covered)!=HISTORICAL_MATCHED_N:
        raise ValueError("Original 927-of-928 exact polygon baseline not recovered")
    covered=set(compatible.island_id)
    if len(covered)>len(original):
        raise AssertionError("Geometry coverage exceeds original botanical universe")
    return compatible,dict(
        n_original_accepted_GIFT_botanical_island_floras=len(original),
        n_original_island_floras_with_exact_hash_rebuilt_polygon=len(covered),
        n_original_island_floras_without_verified_polygon=len(original-covered),
        n_historical_928_direct_floras_with_exact_polygon=len(historical_covered),
        n_historical_928_not_covered=HISTORICAL_COMPLETE_N-len(historical_covered),
        n_source_polygon_records=int(len(gpkg)),
        n_source_original_universe_ids_missing_from_rebuild=int(
            receipt.get("n_IDs_missing_from_geometry",0)),
        all_8265_original_GSHHG_hashes_rebuilt=False,
        only_source_exact_shared_hash_IDs_used=True,
        no_centroid_or_nearest_polygon_substitution=True,
        incomplete_polygons_not_a_biological_absence=True,
    )


def anonymous_network_extension_counts(
    located, historical_928:set, *, all_global_regions=None,
):
    """Only numerical candidate network-site support, not study/pollen efficacy.

    One network ID is required to have at most one polygon match; use
    geospatial 'within' and reject ambiguous boundaries before calling this.
    """
    if not {"network_id","island_id"}<=set(located):
        raise ValueError("No source network/site-to-island IDs")
    if (located.network_id.duplicated().any() or
            located[["network_id","island_id"]].isna().any().any()):
        raise ValueError("Source network point matches are ambiguous")
    old=located.loc[located.island_id.isin(historical_928)]
    extra=located.loc[~located.island_id.isin(historical_928)]
    if len(old)!=75 or old.island_id.nunique()!=16:
        raise ValueError("75 original site/16 GIFT island baseline not reproduced")
    return dict(
        n_source_point_networks_inside_verified_full_botanical_GIFT_polygons=int(len(located)),
        n_source_distinct_botanical_islands_with_network_sites=int(
            located.island_id.nunique()),
        n_original_927_polygon_network_sites_reproduced=int(len(old)),
        n_original_927_polygon_islands_with_network_sites=int(old.island_id.nunique()),
        n_additional_verified_botanical_island_network_sites=int(len(extra)),
        n_additional_verified_botanical_islands_with_network_sites=int(extra.island_id.nunique()),
        only_source_site_points_with_unique_exact_polygon_used=True,
        source_network_site_not_effective_pollen_delivery=True,
    )


def source_scoped_region_network_audit(
    edges, flora, traits:dict, regime:dict, *, source_audit,
):
    """Regional positive network × SAME-SPECIES observed morphology counts.

    Species-wide floral trait states need not be separately observed at the
    source network site. With no plant match, zero means NO DOCUMENTED
    source support, never an ecological nonvisit.
    """
    from .independent_pollination_networks import source_bounded_site_plant_audit
    if not {"island_id","accepted_species"}<=set(flora):
        raise ValueError("Original plant roster missing")
    if not {"network_id","island_id","plant","guild"}<=set(edges):
        raise ValueError("Original positive network edge schema missing")
    if not set(edges.island_id)<=set(flora.island_id):
        raise ValueError("An original network site has no verified floristic polygon")
    out=[]
    categories=("northern_midlatitude","tropical","southern_extratropical",
                "northern_highlatitude","geography_unresolved")
    for region in categories:
        ids={i for i in set(flora.island_id) if regime.get(i,"geography_unresolved")==region}
        scoped=flora[flora.island_id.isin(ids)]
        e=edges[edges.island_id.isin(ids)]
        if scoped.empty:continue
        audit=source_bounded_site_plant_audit(
            e,scoped,traits,global_original_islands=4505
        )
        out.append(dict(
            region=region,
            n_original_botanical_island_floras_in_verified_polygons=int(scoped.island_id.nunique()),
            n_original_published_networks_with_same_accepted_plant_and_guild=int(
                audit["n_original_external_networks_with_plant_guild_overlap"]),
            n_original_islands_with_source_positive_plant_guild_documentation=int(
                audit["n_original_external_network_islands_with_plant_guild_overlap"]),
            true_direct_paired_flower_guild_documentation_cells=audit["cells"],
        ))
    return dict(
        status="REGIONAL_VERIFIED_FULL_FRAME_WOL_SOURCE_POSITIVE_SITE_TRAITS",
        regional_original_source_opportunity=out,
        n_global_parent_botanical_islands=4505,
        site_specific_plant_phenotype_not_measured=True,
        absence_of_source_edges_not_ecological_nonvisit=True,
        independent_study_evidence_not_deduplicated_against_GloBI=True,
        effective_pollination_pollen_deposition_not_measured=True,
    )
