"""Source-conservative FCP × original GIFT native/endemic CoL XR synonym opportunity.

Do NOT automatically correct botanical identity. A GBIF v2 matching
"automatic candidate" to a CoL XR name is not independent checklist proof
of synonymy, verified island local nativeness, or taxon establishment.
Treat as an ambiguity-flagged candidate only, never as an extra original
FCP photo observation or genetic evidence.

Original endemic labels and geologically typed island IDs are rights-safe
inside original CI. Only aggregate status counts leave CI.
"""
from __future__ import annotations
from collections import Counter
import numpy as np
import pandas as pd

COLXR_CHECKLIST="7ddf754f-d193-4cc9-b351-99906754a03b"


def original_candidate_fcp_name_map(
    fcp_species:set[str],colxr:pd.DataFrame
)->tuple[dict[str,str],dict[str,str],dict]:
    required={"submitted_name","candidate_accepted_name","resolution_status",
              "automatic_resolution_candidate","candidate_kingdom",
              "candidate_accepted_rank","candidate_genus","match_type",
              "matched_status","confidence","processing_flags","issues","checklist_key"}
    if not required.issubset(colxr) or len(colxr)!=115328:
        raise ValueError("Frozen original 115328 CoL XR schema or population changed")
    if colxr.submitted_name.duplicated().any() or colxr.submitted_name.isna().any():
        raise ValueError("CoL submitted name candidates duplicate")
    if colxr.checklist_key.ne(COLXR_CHECKLIST).any():
        raise ValueError("CoL XR checklist SHA/identity source changed")
    if len(fcp_species)!=4475:
        raise ValueError("Original measured FCP 4475 source taxonomy opportunity changed")
    original=set(colxr.submitted_name)
    if not original:
        raise ValueError("Original CoL species name frame empty")
    canon=lambda x:" ".join(str(x).strip().split())
    fcp_by_form={}
    for name in sorted(fcp_species):
        form=canon(name)
        if form in fcp_by_form:
            raise ValueError("FCP source names ambiguous under whitespace-only normalization")
        fcp_by_form[form]=name
    direct={}
    candidates={}
    selected=colxr.fillna("")
    for row in selected.itertuples(index=False):
        label=canon(row.submitted_name)
        target=canon(row.candidate_accepted_name)
        if label in fcp_by_form:
            direct[label]=fcp_by_form[label]
            continue
        if (target not in fcp_by_form or not target or target==label):
            continue
        # Provenance: the accepted usage pointer and synonym itself must be
        # independently exact, a species, and source explicitly a synonym.
        if (row.resolution_status!="exact_synonym_to_accepted_candidate" or
            row.automatic_resolution_candidate.lower()!="true" or
            row.match_type!="EXACT" or
            row.candidate_kingdom.lower()!="plantae" or
            row.candidate_accepted_rank.upper()!="SPECIES" or
            row.matched_status.upper()!="SYNONYM" or
            row.processing_flags.strip() or row.issues.strip()):
            continue
        try:
            if float(row.confidence)<95:
                continue
        except (ValueError,TypeError):
            continue
        candidates[label]=fcp_by_form[target]
    # Reused FCP same biological accepted taxon is a one-to-many original
    # source collision; no candidate can be promoted or counted as unique
    # botanical evidence when this cannot be resolved.
    alias_counts=Counter(candidates.values())
    one_to_one={src:dst for src,dst in candidates.items()
                if alias_counts[dst]==1 and dst not in direct.values()}
    collision=sum(alias_counts[dst]>1 or dst in direct.values()
                  for dst in candidates.values())
    return direct,one_to_one,dict(
        n_colxr_original_submitted_names=len(original),
        n_FCP_source_species=4475,
        n_whitespace_exact_source_candidates=len(direct),
        n_exact_synonym_to_original_photo_names_before_collision=len(candidates),
        n_one_to_one_original_photo_synonym_candidates=len(one_to_one),
        n_photo_name_synonym_collisions_quarantined=collision,
        CoL_match_is_candidate_not_taxonomic_review=True,
        no_synonym_candidate_promoted_to_accepted_name=True,
    )


def original_endemic_fcp_taxonomy_sensitivity(
    fcp_species:set[str],original_flora:pd.DataFrame,
    geology:pd.DataFrame,colxr:pd.DataFrame
)->dict:
    req={"island_id","accepted_species","origin_status","endemic_status","floristic_status"}
    if not req.issubset(original_flora) or original_flora.duplicated(
        ["island_id","accepted_species"]).any():
        raise ValueError("Original island GIFT native/endemic status incomplete")
    if ((original_flora.floristic_status.eq("endemic")&
        (~original_flora.origin_status.eq("native")|
         ~original_flora.endemic_status.eq("endemic"))).any()):
        raise ValueError("Original endemic status cannot be inferred from origin status")
    if not {"island_id","geological_origin","analysis_regime","spatial_block"}.issubset(geology):
        raise ValueError("Independent geologically typed island context missing")
    if geology.island_id.duplicated().any():
        raise ValueError("Independent geological type labels duplicated")
    typed=geology.loc[geology.geological_origin.isin(("oceanic","continental"))]
    if len(typed)<2 or set(typed.geological_origin)!={"oceanic","continental"}:
        raise ValueError("Both independent geological island types are required")
    direct,alias,audit=original_candidate_fcp_name_map(fcp_species,colxr)
    orig=original_flora.loc[
        original_flora.island_id.isin(typed.island_id)&
        original_flora.floristic_status.isin(("endemic","native_nonendemic"))
    ].merge(typed[["island_id","geological_origin","analysis_regime","spatial_block"]],
        how="inner",on="island_id",validate="many_to_one")
    orig["literal_FCP_record"]=orig.accepted_species.isin(direct)
    orig["one_to_one_CoL_candidate_FCP_record"]=orig.accepted_species.isin(alias)
    # The exact and synonym categories are DISJOINT.
    if (orig.literal_FCP_record & orig.one_to_one_CoL_candidate_FCP_record).any():
        raise AssertionError("Direct species named the same as a synonym candidate")
    groups=[]
    for (regime,origin,status),sub in orig.groupby(
        ["analysis_regime","geological_origin","floristic_status"],sort=True):
        a=sub.loc[sub.literal_FCP_record]
        b=sub.loc[sub.one_to_one_CoL_candidate_FCP_record]
        groups.append(dict(
            region=regime,geological_origin=origin,
            verified_plant_origin_status=status,
            n_original_islands_with_status=len(set(sub.island_id)),
            n_original_distinct_accepted_species=int(sub.accepted_species.nunique()),
            n_original_exact_FCP_measured_species=int(a.accepted_species.nunique()),
            n_original_exact_FCP_species_island_records=int(len(a)),
            n_one_to_one_CoL_synonym_candidate_FCP_species=int(b.accepted_species.nunique()),
            n_one_to_one_CoL_synonym_candidate_FCP_island_records=int(len(b)),
            n_islands_with_one_to_one_FCP_candidate=int(b.island_id.nunique()),
            n_geography_10degree_blocks_with_candidate=int(b.spatial_block.nunique()),
            n_distinct_FCP_labels_linked_by_candidate=int(len(set(alias.get(name)
                for name in b.accepted_species if name in alias))),
            candidate_requires_independent_botanical_review=True,
        ))
    return dict(
        status="ORIGINAL_GEOLOGICAL_NATIVE_ENDEMIC_FCP_COLXR_CANDIDATE_NAME_AUDIT",
        n_original_geologically_typed_islands=int(len(typed)),
        original_taxonomy_candidate_audit=audit,
        geological_native_endemic_candidate_overlap=groups,
        no_source_taxon_island_photo_or_candidate_accepted_name_exported=True,
        no_fcp_synonym_promoted_to_accepted_native_endemic_or_local_photo=True,
    )
