"""Source-availability (NOT ecological-effect) gate for original 4,505 GIFT floras.

Takes original SHA-locked ANONYMOUS GloBI positive citation support cells.
It never loads, infers, exports, or imputes accepted plants/islands/sites.
Gates are exploratory operational thresholds chosen after source inventory:
>=10 documented islands, >=10 distinct plants, >=10 bibliography labels,
>=2 original archive/provider labels. Same thresholds MUST hold after removing
the globally largest provider before an availability signal is called robust.

Passing this gate is still NOT sufficient for independent field studies,
genuine pollinator nondetection or measured pollen delivery/seed set.
"""
from __future__ import annotations
from itertools import product

REGIONS=("GLOBAL_ALL","northern_midlatitude","tropical",
         "southern_extratropical","northern_highlatitude")
SCOPES=("all_recorded","original_native","original_native_endemic")
MOTIFS=("actinomorphic_deep","zygomorphic_deep","tube_bell_deep",
        "white_tube_bell","red_pink_tube_bell")
GUILDS=("bombus","non_bombus_bees","lepidoptera",
        "flower_visiting_birds","diptera")
MODES=("all_cited","largest_archive_excluded")
UNITS=("n_positive_cited_relation_units","n_documented_islands",
       "n_documented_accepted_plant_species","n_documented_reference_labels",
       "n_documented_archive_first_labels","n_explicit_pollination_relations")
MIN_ISLANDS=10
MIN_PLANTS=10
MIN_REFERENCE_LABELS=10
MIN_ARCHIVE_LABELS=2


def _eligible(q:dict)->bool:
    return (q["n_documented_islands"]>=MIN_ISLANDS
        and q["n_documented_accepted_plant_species"]>=MIN_PLANTS
        and q["n_documented_reference_labels"]>=MIN_REFERENCE_LABELS
        and q["n_documented_archive_first_labels"]>=MIN_ARCHIVE_LABELS)


def audit_original_guild_documentation(source:dict,*,wol:dict|None=None)->dict:
    """Return anonymous source support, never infer animal absence/effects."""
    if source.get("status")!="ORIGINAL_FULL_4505_GIFT_TRUE_PAIR_GLOBI_CITED_LOCAL_GUILD_OPPORTUNITY":
        raise ValueError("Not original full-parent SHA-locked GloBI source audit")
    if (source.get("n_original_parent_GIFT_botanical_island_floras")!=4505
            or source.get("n_original_GloBI_rows_checked")!=24577183):
        raise ValueError("Original GIFT population or GloBI source was altered")
    if source.get("old_928_locality_species_source_replay",{}).get("replay_passed") is not True:
        raise ValueError("Original triple-Direct 928/927/552 source replay absent")
    a=source.get("original_site_plant_real_morphological_guild_support")
    if not isinstance(a,dict) or a.get("status")!="FULL_FRAME_GLOBI_EXACT_SAME_ISLAND_SOURCE_MOTIF_AUDIT":
        raise ValueError("Original source-matched phenotype guild cells missing")
    if a.get("n_original_parent_GIFT_botanical_floras")!=4505:
        raise ValueError("Observed source botanical parent drift")
    cells=a.get("cells")
    if not isinstance(cells,list) or len(cells)!=len(REGIONS)*len(SCOPES)*len(MOTIFS)*len(GUILDS)*len(MODES):
        raise ValueError("Original full anonymous source support tensor incomplete")
    indexed={}
    for q in cells:
        if not isinstance(q,dict):
            raise ValueError("Original anonymized source cell malformed")
        key=(q.get("source_mode"),q.get("region"),q.get("botanical_status"),
             q.get("real_species_direct_pair"),q.get("original_source_visitor_guild"))
        if key in indexed or any(k is None for k in key):
            raise ValueError("Source support cells duplicate or absent classification")
        for field in UNITS:
            v=q.get(field)
            if not isinstance(v,int) or v<0:
                raise ValueError("Source support counts invalid: "+field)
        if (q["n_documented_reference_labels"]>q["n_positive_cited_relation_units"]
            or q["n_documented_archive_first_labels"]>q["n_positive_cited_relation_units"]
            or q["n_explicit_pollination_relations"]>q["n_positive_cited_relation_units"]):
            raise ValueError("Original reference/provider/pollination subset exceeds citations")
        indexed[key]=q
    expected=set(product(MODES,REGIONS,SCOPES,MOTIFS,GUILDS))
    if set(indexed)!=expected:
        raise ValueError("Original positive-source cross-table changed")
    output=[]
    for region,scope in product(REGIONS,SCOPES):
        preliminary=0
        provider_stable=0
        trait_guild=[]
        for motif,guild in product(MOTIFS,GUILDS):
            before=indexed[("all_cited",region,scope,motif,guild)]
            leaveout=indexed[("largest_archive_excluded",region,scope,motif,guild)]
            for f in UNITS:
                if leaveout[f]>before[f]:
                    raise ValueError("Largest-archive omission gained original documentary support")
            pre=_eligible(before)
            post=_eligible(leaveout)
            if post and not pre:
                raise ValueError("Omission creates false source support")
            preliminary+=int(pre)
            provider_stable+=int(post)
            if pre or post or (
                guild in ("flower_visiting_birds","lepidoptera") and
                motif in ("actinomorphic_deep","zygomorphic_deep")):
                trait_guild.append(dict(motif=motif,guild=guild,
                    n_cited_islands=int(before["n_documented_islands"]),
                    n_cited_plants=int(before["n_documented_accepted_plant_species"]),
                    n_provider_labels=int(before["n_documented_archive_first_labels"]),
                    n_cited_bibliography_labels=int(before["n_documented_reference_labels"]),
                    after_largest_archive_omitted_islands=int(leaveout["n_documented_islands"]),
                    after_largest_archive_omitted_plants=int(leaveout["n_documented_accepted_plant_species"]),
                    after_largest_archive_omitted_provider_labels=int(
                        leaveout["n_documented_archive_first_labels"]),
                    preliminary_documentary_support=pre,
                    stable_documentary_support_after_provider_removal=post))
        output.append(dict(region=region,botanical_status=scope,
            n_preselected_animal_guild_x_true_morphology_endpoints=25,
            n_preliminary_documentation_gates_passing=preliminary,
            n_provider_removal_documentation_gates_passing=provider_stable,
            pre_vs_post_candidate_supported_profile_guilds=trait_guild))
    more={}
    if wol is not None:
        if (wol.get("status")!="ORIGINAL_GIFT_4505_EXACT_HASH_FULLFRAME_WOL_LOCAL_POLLINATOR_SUPPORT"
            or wol.get("parent_botanical_floras_original")!=4505):
            raise ValueError("WOL original source-protected botanical parent drift")
        w=wol["full_original_botanical_flora_island_visitor_evidence"]
        more=dict(
            n_WOL_source_geographically_supported_botanical_islands=int(
                w["n_original_external_network_islands_with_plant_guild_overlap"]),
            n_WOL_original_positive_network_botanical_units=int(
                w["n_original_species_network_plant_guild_units"]),
            WOL_and_GloBI_underlying_study_overlap_unknown=True,
            no_sum_of_WOL_and_GloBI_sites_or_species_or_citations=True,
        )
    return dict(
        status="ORIGINAL_4505_GIFT_SOURCE_POLLINATION_DOCUMENTATION_CAPACITY",
        n_original_global_botanical_island_floras=4505,
        n_original_exact_GSHHG_polygon_flora=int(
            source["geometry_original_biological_source"][
                "n_original_island_floras_with_exact_hash_rebuilt_polygon"]),
        n_original_GloBI_positive_documented_islands=int(
            a["n_source_original_islands_with_some_positive_relation"]),
        n_original_GloBI_positive_plant_x_island_pairs=int(
            a["n_source_positive_original_island_plant_pairs"]),
        n_original_GloBI_source_bibliographic_label_strings=int(
            a["n_original_source_bibliographic_keys"]),
        minimum_source_islands=MIN_ISLANDS,
        minimum_distinct_botanical_plant_species=MIN_PLANTS,
        minimum_citation_labels=MIN_REFERENCE_LABELS,
        minimum_original_archive_provider_labels=MIN_ARCHIVE_LABELS,
        source_qualification_panels=output,
        other_provider_source_opportunity=more,
        conclusions=dict(
            threshold_is_exploratory_posthoc_feasibility_not_preregistered_inference=True,
            passing_documentary_gate_not_evidence_of_independent_original_publications=True,
            missing_source_guild_document_not_animal_absence=True,
            original_species_wide_flower_character_not_local_phenotype=True,
            source_label_pollination_is_not_measured_pollen_transfer=True,
            documented_animal_guild_not_effective_pollinator_selection=True,
            no_pvalue_or_effect_coefficient_calculated=True,
            original_macro_populations_not_narrowed_to_928_or_232_islands=True,
            no_source_island_genus_coordinates_or_publication_ids_exported=True,
        ),
    )
