"""Exact GSHHG island × accepted-plant × original GloBI citation versus real two-trait flowers.

This is an evidence-availability audit, NOT a pollinator or fitness comparison.
Missing animal records are never classified as guild absence.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

GUILDS=("bombus","non_bombus_bees","lepidoptera","flower_visiting_birds","diptera")
PAIRS={
    "actinomorphic_deep":("floral_symmetry",frozenset({"actinomorphic"})),
    "zygomorphic_deep":("floral_symmetry",frozenset({"zygomorphic"})),
    "tubular_bell_deep":("floral_form",frozenset({"tubular","bell_campanulate"})),
}
SCOPES=("all_recorded","strict_native","unresolved_origin")
MODES=("original_citations","largest_archive_omitted")


def real_same_species_island_visitor_audit(
    flora:pd.DataFrame,linked:pd.DataFrame,
    direct:dict[str,dict[str,str]],globally_dominant_archive:str,
    *,original_global_flora_n:int=4505,verified_polygons_n:int=927,
)->dict:
    """Return ONLY non-identifying aggregate counts, no private IDs."""
    if flora.island_id.nunique()!=928 or original_global_flora_n!=4505 or verified_polygons_n!=927:
        raise ValueError("original 4505 population / 928 selected / 927 polygon cohort drift")
    if flora.duplicated(["island_id","accepted_species"]).any() or flora.origin_status.isna().any():
        raise ValueError("duplicate or unknown origin-status source flora")
    need={"island_id","accepted_species","origin_status","guild","reference_id",
          "provider_id","evidence_class","latitude","longitude"}
    if need-set(linked):
        raise ValueError("missing source-verified exact polygon GloBI relation keys")
    if not isinstance(globally_dominant_archive,str) or not globally_dominant_archive:
        raise ValueError("dominant source label must be selected globally first")
    if not linked.empty:
        if linked[["island_id","accepted_species","guild","reference_id","provider_id"]].isna().any().any():
            raise ValueError("source references and plant identities missing")
        if (not linked.guild.isin(GUILDS).all() or
            not linked.evidence_class.isin(("flower_visit","pollination")).all()):
            raise ValueError("source interaction predicates or guilds invalid")
        if not set(zip(linked.island_id,linked.accepted_species)) <= set(zip(flora.island_id,flora.accepted_species)):
            raise ValueError("local animal relation is NOT an exact original island flora species")
    if not {"floral_form","floral_symmetry","tube_depth_class"} <= set(direct):
        raise ValueError("species-level Direct both axes not available")
    result=[]
    for label,(axis,positive) in PAIRS.items():
        q=flora[["island_id","accepted_species","origin_status"]].copy()
        q["a"]=q.accepted_species.map(direct[axis])
        q["b"]=q.accepted_species.map(direct["tube_depth_class"])
        q=q.loc[q.a.notna()&q.b.notna()].copy()
        q["pair_state"]=np.where(q.a.isin(positive)&q.b.eq("deep"),
                                  "joint_true","both_axes_known_not_joint")
        visit=linked.merge(q,on=["island_id","accepted_species","origin_status"],
                           how="inner",validate="many_to_one")
        if len(visit)>len(linked):
            raise AssertionError("GloBI rows multiplied by botanical trait merge")
        for scope in SCOPES:
            if scope=="strict_native":
                sub=q.loc[q.origin_status.eq("native")]
                v=visit.loc[visit.origin_status.eq("native")]
            elif scope=="unresolved_origin":
                sub=q.loc[q.origin_status.eq("unresolved")]
                v=visit.loc[visit.origin_status.eq("unresolved")]
            else:
                sub=q;v=visit
            for state in ("joint_true","both_axes_known_not_joint"):
                eligible=sub.loc[sub.pair_state.eq(state)]
                x=v.loc[v.pair_state.eq(state)]
                for mode in MODES:
                    z=x if mode==MODES[0] else x.loc[x.provider_id.ne(globally_dominant_archive)]
                    for guild in GUILDS:
                        g=z.loc[z.guild.eq(guild)]
                        by=g.groupby(["island_id","accepted_species"])
                        pair=g[["island_id","accepted_species"]].drop_duplicates()
                        result.append({
                            "real_same_accepted_species_morphological_pair":label,
                            "original_plant_origin_scope":scope,
                            "actual_direct_pair_state":state,
                            "GloBI_source_variant":mode,
                            "documented_original_visitor_guild":guild,
                            "n_original_island_plant_pairs_with_both_direct_traits":int(len(eligible)),
                            "n_eligible_original_island_floras":int(eligible.island_id.nunique()),
                            "n_cited_island_plant_pairs":int(len(pair)),
                            "n_islands_with_source_cited_exact_plant_and_guild":int(g.island_id.nunique()),
                            "n_plant_species_with_source_cited_locality":int(g.accepted_species.nunique()),
                            "n_source_bibliographic_keys":int(g.reference_id.nunique()),
                            "n_source_archive_labels":int(g.provider_id.nunique()),
                            "n_island_plant_pairs_with_3plus_distinct_references":int(
                                by.reference_id.nunique().ge(3).sum()) if len(g) else 0,
                            "n_island_plant_pairs_with_2plus_archive_labels":int(
                                by.provider_id.nunique().ge(2).sum()) if len(g) else 0,
                            "n_island_plant_pairs_with_pollination_predicate":int(
                                len(g.loc[g.evidence_class.eq("pollination"),
                                    ["island_id","accepted_species"]].drop_duplicates())),
                            "source_locality_is_not_effective_plant_pollination":True,
                            "unrecorded_guild_is_not_absent":True,
                        })
    if len(result)!=len(PAIRS)*len(SCOPES)*2*len(MODES)*len(GUILDS):
        raise AssertionError("source-guild joint morphology support lattice incomplete")
    return {
        "status":"ORIGINAL_GLOBI_EXACT_SITE_SAME_SPECIES_TRUE_FLOWER_PAIR_COVERAGE",
        "n_original_GIFT_botanical_island_population":4505,
        "n_original_selected_colour_form_symmetry_complete_islands":928,
        "n_original_geometrically_verifiable_islands":927,
        "n_actual_true_species_direct_phenotype_combinations":len(PAIRS),
        "n_source_guild_trait_plant_origin_support_cells":len(result),
        "support":result,
        "limitations":{
            "both_morphological_traits_required_on_same_accepted_plant_species":True,
            "GloBI_source_site_within_exact_original_GSHHG_polygon_required":True,
            "actual_original_flora_plant_name_at_same_island_required":True,
            "animal_nondetection_is_never_inferred":True,
            "citation_predicate_not_actual_pollination_efficiency":True,
            "lepidoptera_butterfly_and_moth_separation_not_supported":True,
            "no_mainland_island_matched_effective_visits":True,
            "source_archive_omission_determined_before_local_trait_counts":True,
            "selected_928_not_universal_4505_population":True,
            "source_island_species_citation_and_exact_site_ids_private":True,
        },
    }
