"""Strict, source-pinned CoL XR family map and 928-island family opportunity audit.

Source: original 2.0.0-alpha2-taxonomy, signed schema/candidate receipt
from zuizui0223/island, Col XR checklist 7ddf754f-d193-4cc9-b351-99906754a03b.

Crucial: alpha2 automatic-resolution records are CANDIDATES not promoted
taxonomic classifications. Only EXACT accepted-species candidate with the
submitted original family matching candidate family, exact accepted plant
name matching, and genus matching is accepted for a restricted sensitivity.
No COL XR name-match candidate is promoted to accepted taxonomy or phylogeny.
Never use genus prefix to GUESS a family.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import hashlib
import numpy as np
import pandas as pd

from .rare_profile_decomposition import genus

COLXR_EXPECTED_CHECKLIST="7ddf754f-d193-4cc9-b351-99906754a03b"
REQUIRED=("submitted_name","submitted_family","candidate_family",
          "candidate_genus","candidate_kingdom",
          "candidate_accepted_name","resolution_status","match_type",
          "automatic_resolution_candidate","checklist_key")


def strict_candidate_family_map(crosswalk:pd.DataFrame, names:set[str],*,
                                expected_rows:int=115328):
    if not set(REQUIRED).issubset(crosswalk):
        raise ValueError("Col XR original source taxonomy schema unexpected")
    if len(crosswalk)!=expected_rows or expected_rows<5:
        raise ValueError("Original pinned alpha2 source species universe drifted")
    raw=crosswalk[list(REQUIRED)].fillna("").astype(str)
    if raw.submitted_name.duplicated().any():
        raise ValueError("Taxonomic source has duplicate submitted names")
    if raw.checklist_key.ne(COLXR_EXPECTED_CHECKLIST).any():
        raise ValueError("Original source CoL XR checklist key changed")
    relevant=raw.loc[raw.submitted_name.isin(names)].copy()
    if relevant.empty:
        return {},{"n_original_species":len(names),"n_source_colxr_name_matched":0,
                  "n_strict_family_candidates":0,"n_genus_family_conflicts":0}
    # Original candidates—not independently reviewed accepted taxonomy.
    family=relevant.submitted_family.str.strip()
    accepted_family=relevant.candidate_family.str.strip()
    col_genus=relevant.candidate_genus.str.strip()
    reported_genus=relevant.submitted_name.map(genus)
    ok=(relevant.match_type.eq("EXACT")&
        relevant.automatic_resolution_candidate.str.lower().eq("true")&
        relevant.resolution_status.eq("exact_accepted_candidate")&
        relevant.candidate_kingdom.str.lower().eq("plantae")&
        relevant.candidate_accepted_name.eq(relevant.submitted_name)&
        family.ne("")&accepted_family.ne("")&
        family.str.casefold().eq(accepted_family.str.casefold())&
        col_genus.eq(reported_genus)&reported_genus.ne("")&
        ~family.str.lower().isin(("unknown","unresolved","incertae sedis")))
    selected=relevant.loc[ok].copy()
    # Fail closed if the same genus has conflicting original family
    # assignments even when all individual species pass exact CoL matching.
    selected["source_genus"]=selected.submitted_name.map(genus)
    conflict=(selected.groupby("source_genus").candidate_family.nunique()>1)
    bad=set(conflict[conflict].index)
    selected=selected.loc[~selected.source_genus.isin(bad)]
    mapping=selected.set_index("submitted_name").candidate_family.to_dict()
    audits={
        "n_original_species":len(names),
        "n_source_colxr_name_matched":len(relevant),
        "n_strict_family_candidates_before_genus_conflict_quarantine":int(ok.sum()),
        "n_strict_family_candidates":len(mapping),
        "n_original_direct_species_not_exact_family_candidate":len(names)-len(mapping),
        "n_candidate_families":int(selected.candidate_family.nunique()),
        "n_candidate_genera_with_family":int(selected.source_genus.nunique()),
        "n_genus_family_conflicts":len(bad),
        "n_candidate_species_discarded_for_genus_family_conflicts":int(
            ok.sum())-len(mapping),
        "candidate_taxonomy_not_verified_reviewed_final":True,
        "cannot_infer_ancestry_from_family_rank":True,
    }
    return mapping,audits


def group_by_family(names:set[str], mapping:dict[str,str]):
    fam=defaultdict(list)
    for name in sorted(names):
        family=mapping.get(name)
        if family:
            fam[family].append(name)
    return {k:tuple(v) for k,v in fam.items()}


def family_source_support(triple:pd.DataFrame,flora:pd.DataFrame,geo:pd.DataFrame,
                          assignments:pd.DataFrame,mainland:pd.DataFrame,
                          family_map:dict[str,str])->dict:
    needed=[
      (flora,{"island_id","accepted_species","origin_status"}),
      (geo,{"island_id","analysis_regime","spatial_block"}),
      (assignments,{"island_id","source_mode","source_rank","entity_ID"}),
      (mainland,{"entity_ID","work_species"})]
    for frame,columns in needed:
        if not columns.issubset(frame):
            raise ValueError("Missing source family opportunity input columns")
    valid=flora.loc[flora.accepted_species.isin(triple.index)].copy()
    counts=valid.groupby("island_id").accepted_species.nunique()
    eligible=set(counts[counts>=30].index)
    if len(eligible)!=928:
        raise ValueError("Original full 928 island floral scope changed")
    valid=valid.loc[valid.island_id.isin(eligible)]
    geo=geo[["island_id","analysis_regime","spatial_block"]].drop_duplicates()
    if geo.island_id.duplicated().any():
        raise ValueError("Duplicate island geographical block")
    g=geo.set_index("island_id")
    nearby=assignments.loc[assignments.island_id.isin(eligible)&
       assignments.source_mode.eq("geo_k5")&
       assignments.source_rank.between(1,5),
       ["island_id","source_rank","entity_ID"]].copy()
    if nearby.duplicated(["island_id","source_rank"]).any():
        raise ValueError("Candidate GIFT mainland source ranks duplicate")
    source=nearby.sort_values(["island_id","source_rank"]).groupby(
      "island_id").entity_ID.agg(tuple).to_dict()
    if any(len(source.get(i,()))!=5 for i in eligible):
        raise ValueError("Original geo_k5 five nearest mainland sources missing")
    pool=(mainland.loc[mainland.work_species.isin(triple.index),
       ["entity_ID","work_species"]].drop_duplicates().groupby(
       "entity_ID").work_species.agg(set).to_dict())
    by={i:z for i,z in valid.groupby("island_id")}
    rows=[]
    for island in sorted(eligible):
        islandset=set(by[island].accepted_species)
        donor=(set().union(*(pool.get(e,set()) for e in source[island]))&
               set(triple.index))-islandset
        b=group_by_family(donor,family_map)
        for scope,own in (
           ("all_recorded",islandset),
           ("source_native",set(by[island].loc[
                by[island].origin_status.eq("native"),"accepted_species"]))):
            if len(own)<30:continue
            a=group_by_family(own,family_map)
            both=set(a)&set(b)
            # Family can be shared even when only the SAME genus is
            # represented in donor and island. For a stricter n-family
            # alternative-genus counter-check, require at least one original
            # source plant pair with different GENUS inside the family.
            different_genus={
               family for family in both
               if len({genus(n) for n in a[family]}|
                      {genus(n) for n in b[family]})>=2
               and any(genus(n)!=genus(m) for n in a[family] for m in b[family])
            }
            rows.append(dict(
               island_id=island,scope=scope,
               regime=str(g.at[island,"analysis_regime"]),
               spatial_block=str(g.at[island,"spatial_block"]),
               source_signature=source[island][0],
               n_local_family_mapped_species=int(sum(map(len,a.values()))),
               n_source_family_mapped_species=int(sum(map(len,b.values()))),
               n_local_families=len(a),n_donor_families=len(b),
               n_shared_families=len(both),
               n_shared_families_with_different_genus_options=len(different_genus),
               supports_equal_five_family_identity=len(both)>=5,
               supports_equal_five_family_different_genus=len(different_genus)>=5,
               supports_equal_eight_family_identity=len(both)>=8,
               supports_equal_15_distinct_family_richness=len(a)>=15 and len(b)>=15,
               supports_equal_ten_distinct_family_richness=len(a)>=10 and len(b)>=10,
               eligible_original_geographic_candidate_source=len(donor)>=30
            ))
    private=pd.DataFrame(rows)
    result=[]
    for (scope,regime),q in private.groupby(["scope","regime"]):
        if len(q)<10:continue
        rec={
           "scope":scope,"regime":regime,
           "n_original_eligible_flora_islands":len(q),
           "n_source_disjoint_FD_islands":int(q.eligible_original_geographic_candidate_source.sum()),
           "n_family_5_shared_identity":int(q.supports_equal_five_family_identity.sum()),
           "n_family_5_different_genus_identity":int(q.supports_equal_five_family_different_genus.sum()),
           "n_family_8_shared_identity":int(q.supports_equal_eight_family_identity.sum()),
           "n_family_15_equal_richness":int(q.supports_equal_15_distinct_family_richness.sum()),
           "n_family_10_equal_richness":int(q.supports_equal_ten_distinct_family_richness.sum()),
           "median_family_mapped_island_Direct_species":float(q.n_local_family_mapped_species.median()),
           "median_family_mapped_source_Direct_species":float(q.n_source_family_mapped_species.median()),
           "median_shared_candidate_family_count":float(q.n_shared_families.median()),
           "n_distinct_physical_blocks_with_family5":int(q.loc[
              q.supports_equal_five_family_identity,"spatial_block"].nunique()),
           "n_distinct_first_GIFT_source_groups_with_family5":int(q.loc[
              q.supports_equal_five_family_identity,"source_signature"].nunique()),
           "n_distinct_physical_blocks_with_different_genus_family5":int(q.loc[
              q.supports_equal_five_family_different_genus,"spatial_block"].nunique()),
           "n_first_GIFT_sources_with_different_genus_family5":int(q.loc[
              q.supports_equal_five_family_different_genus,"source_signature"].nunique()),
        }
        result.append(rec)
    return {
       "status":"COLXR_STRICT_FAMILY_CANDIDATE_928_ISLAND_POSITIVITY_ONLY",
       "n_original_eligible_islands":928,
       "n_direct_flora_scope_records":len(private),
       "n_taxonomy_candidate_mapped_species_in_triple":len(family_map),
       "region_family_support":result,
       "family_accepted_only_as_unreviewed_original_source_taxonomy_candidate":True,
       "no_family_or_phylogenetic_effect_estimated":True,
       "original_plant_family_map_and_taxon_IDs_not_published":True,
    }
