"""Rights-safe FCP photo evidence overlap with geologically typed island endemics.

FCP is an independently measured (iNaturalist) 4,475 accepted photographic
species ledger. Its four-state colour D/W labels are PHOTO-BASED, not genetic
morphs or botanical source-verified corolla colours, and the species-wide
ledger cannot localize photographs to *this particular island*. This program
audits SOURCE OPPORTUNITY ONLY, not oceanic/continental adaptive divergence.

The original GIFT native mainland source checklist is a potential donor
inventory, not actual propagules or the ancestral flowering population.
"""
from __future__ import annotations
from collections import defaultdict

import numpy as np
import pandas as pd


def _genus(species:str)->str:
    v=str(species).split()
    return v[0] if len(v)>=2 else ""


def audit_fcp_endemic_geology_opportunity(
    fcp:pd.DataFrame,status:pd.DataFrame,geology:pd.DataFrame,
    assignment:pd.DataFrame,source:pd.DataFrame
)->dict:
    fcp_fields={"species","n_source_photo","n_classified_photo",
                "D_photo_pairwise_unbiased","W_axis_sq_alignment",
                "n_local50km_colour_photo_pairs"}
    if not fcp_fields.issubset(fcp) or len(fcp)!=4475:
        raise ValueError("Original frozen FCP 4475 species photo ledger changed")
    if fcp.species.duplicated().any() or fcp.species.isna().any():
        raise ValueError("Cannot join ambiguous original photo taxon identities")
    photo=fcp.set_index("species",drop=False).copy()
    for col in ("n_source_photo","n_classified_photo",
                "n_local50km_colour_photo_pairs"):
        photo[col]=pd.to_numeric(photo[col],errors="raise")
    if int(photo.n_source_photo.sum())!=447500 or (
        int(photo.n_classified_photo.sum())!=227107):
        raise ValueError("FCP 447500/227107 original photo source counts changed")
    required={"island_id","accepted_species","origin_status",
              "endemic_status","floristic_status"}
    if not required.issubset(status) or status.duplicated(
        ["island_id","accepted_species"]).any():
        raise ValueError("Original source island/endemic status unavailable or duplicated")
    if ((status.floristic_status.eq("endemic")&
        (~status.origin_status.eq("native")|
         ~status.endemic_status.eq("endemic"))).any()):
        raise ValueError("Original accepted endemic status contradictory")
    if not {"island_id","geological_origin","analysis_regime","spatial_block"}.issubset(geology):
        raise ValueError("Original independently classified geology missing")
    if geology.island_id.duplicated().any():
        raise ValueError("Original geological island crosswalk duplicated")
    typed=geology.loc[geology.geological_origin.isin(("oceanic","continental"))]
    if len(typed)<4:
        raise ValueError("Insufficient independent oceanic/continental source labels")
    names=status.loc[status.island_id.isin(typed.island_id)].merge(
        typed[["island_id","geological_origin","analysis_regime","spatial_block"]],
        on="island_id",validate="many_to_one",how="inner")
    if not {"island_id","source_mode","source_rank","entity_ID"}.issubset(assignment):
        raise ValueError("Original closest GIFT mainland candidate reference missing")
    near=assignment.loc[assignment.source_mode.eq("geo_k5")].copy()
    near["source_rank"]=pd.to_numeric(near.source_rank,errors="raise")
    near=near.loc[near.source_rank.between(1,5)]
    if near.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Original candidate mainland source region duplicated")
    byisland=near.sort_values(["island_id","source_rank"]).groupby(
        "island_id").entity_ID.agg(tuple).to_dict()
    if not {"entity_ID","work_species"}.issubset(source):
        raise ValueError("Original GIFT native source species inventory missing")
    donor=source[["entity_ID","work_species"]].drop_duplicates().groupby(
        "entity_ID").work_species.agg(set).to_dict()
    fcp_taxa=set(photo.index)
    rows=[]
    for island,grp in names.groupby("island_id"):
        sources=byisland.get(str(island),())
        valid=(len(sources)==5 and len(set(sources))==5)
        own=set(grp.accepted_species)
        mainland=(set().union(*(donor.get(y,set()) for y in sources))-own
                  if valid else set())
        photo_source=(mainland&fcp_taxa)
        source_genera={_genus(x) for x in photo_source}
        for role in ("endemic","native_nonendemic"):
            eligible=grp.loc[grp.floristic_status.eq(role)]
            hits=eligible.loc[eligible.accepted_species.isin(fcp_taxa)].copy()
            for species in hits.accepted_species:
                rec=photo.loc[species]
                g=_genus(species)
                rows.append(dict(
                    island_id=str(island),region=str(grp.analysis_regime.iloc[0]),
                    origin=str(grp.geological_origin.iloc[0]),
                    block=str(grp.spatial_block.iloc[0]),
                    role=role,
                    taxon=species,genus=g,
                    exact_original_GIFT_source_available=valid,
                    independent_photo_source_same_genus_available=(valid and g in source_genera),
                    n_original_FCP_photo=int(rec["n_source_photo"]),
                    n_original_FCP_colour_classified=int(rec["n_classified_photo"]),
                    n_original_source_local50km_pairs=int(rec["n_local50km_colour_photo_pairs"]),
                    FCP_D_estimated=bool(pd.notna(rec["D_photo_pairwise_unbiased"])),
                    FCP_W_estimated=bool(pd.notna(rec["W_axis_sq_alignment"])),
                ))
    df=pd.DataFrame(rows)
    output=[]
    for region in sorted(names.analysis_regime.dropna().unique()):
        for origin in ("oceanic","continental"):
            all_islands=names.loc[
                names.analysis_regime.eq(region)&
                names.geological_origin.eq(origin)]
            for role in ("endemic","native_nonendemic"):
                elig=all_islands.loc[all_islands.floristic_status.eq(role)]
                x=df.loc[(df.region.eq(region))&
                         (df.origin.eq(origin))&
                         (df.role.eq(role))] if not df.empty else pd.DataFrame()
                src=x.loc[x.independent_photo_source_same_genus_available] if not x.empty else x
                enough=x.loc[x.n_original_FCP_colour_classified.ge(40)] if not x.empty else x
                output.append(dict(
                    region=str(region),geological_origin=origin,
                    native_endemic_status=role,
                    n_original_independently_geology_typed_islands=int(
                        all_islands.island_id.nunique()),
                    n_original_GIFT_status_island_species=int(len(elig)),
                    n_original_GIFT_status_distinct_taxa=int(elig.accepted_species.nunique()),
                    n_original_exact_FCP_measured_species_rows=int(len(x)),
                    n_original_exact_FCP_distinct_species=(
                        int(x.taxon.nunique()) if not x.empty else 0),
                    n_original_exact_FCP_distinct_genera=(
                        int(x.genus.nunique()) if not x.empty else 0),
                    n_independently_typed_islands_with_FCP_species=(
                        int(x.island_id.nunique()) if not x.empty else 0),
                    n_independent_geographic_blocks_with_FCP_species=(
                        int(x.block.nunique()) if not x.empty else 0),
                    n_FCP_photo_records_original=(
                        int(x.n_original_FCP_photo.sum()) if not x.empty else 0),
                    n_species_with_at_least_40_classifiable_FCP_photos=(
                        int(enough.taxon.nunique()) if not x.empty else 0),
                    n_species_with_FCP_D=(
                        int(x.loc[x.FCP_D_estimated].taxon.nunique()) if not x.empty else 0),
                    n_species_with_FCP_W=(
                        int(x.loc[x.FCP_W_estimated].taxon.nunique()) if not x.empty else 0),
                    n_FCP_species_with_local50km_photo_pair_opportunity=(
                        int(x.loc[x.n_original_source_local50km_pairs.gt(0)].taxon.nunique())
                        if not x.empty else 0),
                    n_species_with_GIFT_candidate_mainland_FCP_same_genus=(
                        int(src.taxon.nunique()) if not x.empty else 0),
                    n_distinct_genera_with_candidate_mainland_FCP_same_genus=(
                        int(src.genus.nunique()) if not x.empty else 0),
                ))
    return dict(
        status="SOURCE_VERIFIED_FCP_4475_GEOLOGY_ENDEMIC_PHOTOGRAPHIC_OPPORTUNITY",
        n_original_geologically_typed_islands=len(typed),
        n_FCP_species_original=4475,
        n_FCP_original_photo_measurements=447500,
        n_FCP_colour_classifiable_original=227107,
        original_regional_photo_opportunity=output,
        scope=dict(
            exact_species_name_join_only_no_inferred_synonyms=True,
            original_status_native_and_endemic_verified_independently=True,
            GIFT_native_region_donor_is_not_the_ancestor=True,
            FCP_species_wide_photo_colour_is_not_island_local_photo_trait=True,
            FCP_photo_colour_not_genetic_polymorphism_or_adaptive_selection=True,
            source_photo_record_numbers_reused_across_islands_not_unique_global_sum=True,
            original_taxa_or_island_IDs_never_exported=True,
        ),
    )
