"""GBIF-name-usage taxonomic rescue for Mangal original visitor nodes.

Only ORIGINAL Mangal node.taxonomy.gbif keys on actual positive source
network interactions of SAME original GIFT flower-morphology plant pairs may
be looked up. Use immutable input source nodes and record GBIF response
snapshot; no fuzzy lookup on arbitrary visitor name, no inference from
non-observation. GBIF taxonomy does NOT independently confirm a visit or
effective pollination and source network metadata may duplicate GloBI.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from .independent_pollination_networks import (
    visitor_guild,real_same_species_morphologies,
)


def gbif_taxonomy_visitor_guild(usage:dict,key:int)->str:
    if not isinstance(usage,dict):
        return "unresolved_visitor_taxonomy"
    try: returned=int(usage.get("key",0))
    except (TypeError,ValueError):return "unresolved_visitor_taxonomy"
    if returned!=int(key):
        return "unresolved_visitor_taxonomy"
    if str(usage.get("kingdom") or "").lower()!="animalia":
        return "unresolved_visitor_taxonomy"
    return visitor_guild(usage.get("order"),usage.get("family"))


def target_original_morphology_partner_keys(
    linked:pd.DataFrame,source_direct:dict[str,dict[str,str]]
)->pd.DataFrame:
    required={"network_id","island_id","plant",
              "original_partner_gbif_key","original_partner_scientific_name",
              "interaction_id"}
    if needed:=required-set(linked):
        raise ValueError(f"missing original source node+real flower link {sorted(needed)}")
    records=[]
    for q in linked.itertuples(index=False):
        phen=real_same_species_morphologies(q.plant,source_direct)
        if phen is None or not any(x is True for x in phen.values()):
            continue
        key=q.original_partner_gbif_key
        if pd.isna(key): gbifkey=None
        else:
            try: gbifkey=int(key)
            except (TypeError,ValueError):gbifkey=None
        records.append(dict(
            network_id=q.network_id,
            island_id=q.island_id,
            plant=q.plant,
            partner_gbif_key=gbifkey,
            source_interaction_id=q.interaction_id,
            **{key+"_true":v for key,v in phen.items()},
        ))
    return pd.DataFrame.from_records(records,columns=[
        "network_id","island_id","plant","partner_gbif_key",
        "source_interaction_id","actinomorphic_deep_true",
        "zygomorphic_deep_true","tube_bell_deep_true"
    ])


def summarize_gbif_resolved_source(
    rows:pd.DataFrame,gbif_responses:dict[int,dict],
    *,total_original_islands:int=4505,
)->dict:
    if total_original_islands!=4505:
        raise ValueError("original full botanical island population changed")
    required={"network_id","island_id","plant","partner_gbif_key",
              "origin_status","endemic_status",
              "actinomorphic_deep_true","zygomorphic_deep_true",
              "tube_bell_deep_true"}
    if required-set(rows):
        raise ValueError("original true pair partner source not supplied")
    out=[]
    all_sources=[
        ("all_recorded",rows),
        ("strict_native",rows.loc[rows.origin_status.eq("native")]),
        ("strict_native_endemic",rows.loc[
            rows.origin_status.eq("native")&rows.endemic_status.eq("endemic")
        ]),
        ("origin_unresolved",rows.loc[rows.origin_status.eq("unresolved")]),
    ]
    for scope,eligible in all_sources:
        for pair in ("actinomorphic_deep","zygomorphic_deep",
                     "tube_bell_deep"):
            matches=eligible.loc[eligible[pair+"_true"].eq(True)].copy()
            matches["guild"]=matches.partner_gbif_key.map(
                lambda q:gbif_taxonomy_visitor_guild(
                    gbif_responses.get(int(q),{}),int(q)
                ) if pd.notna(q) else "unresolved_visitor_taxonomy"
            )
            for guild in (
                "bird_family_supported","butterfly_family_supported",
                "sphingidae_moth","other_lepidoptera",
                "bees_family_supported","flies_order_supported",
                "unresolved_visitor_taxonomy",
            ):
                q=matches.loc[matches.guild.eq(guild)].copy()
                out.append(dict(
                    botanical_scope=scope,
                    same_accepted_plant_true_flower_pair=pair,
                    GBIF_verified_partner_guild=guild,
                    n_source_positive_network_island_plant_partner_units=int(len(q)),
                    n_source_networks_with_positive_guild=int(q.network_id.nunique()),
                    n_original_islands_with_source_positive_guild=int(q.island_id.nunique()),
                    n_source_network_plant_pairs_with_guild=int(
                        len(q[["network_id","plant"]].drop_duplicates())),
                    no_unobserved_taxon_or_absence_imputed=True,
                ))
    return dict(
        status="ORIGINAL_MANGAL_SAME_PLANT_FLORAL_PAIR_GBIF_TAXONOMY_EVIDENCE",
        n_original_GIFT_botanical_islands=4505,
        n_original_direct_flower_pair_partner_source_interactions=int(len(rows)),
        n_distinct_original_GBIF_visitor_taxonomy_keys_in_pairs=int(
            rows.partner_gbif_key.dropna().nunique()
        ) if len(rows) else 0,
        n_GBIF_keys_retrieved_from_public_backbone=len(gbif_responses),
        n_source_guild_support_cells=len(out),
        support=out,
        limitations=dict(
            GBIF_taxonomy_retrieved_from_original_Mangal_source_node_identifier=True,
            original_Mangal_node_gbif_taxonomy_is_not_independent_visitor_observation=True,
            morphology_requires_both_original_Direct_species_axes=True,
            no_morphological_visitor_unrecorded_state_called_negative=True,
            these_positive_links_not_effective_pollen_or_seed_set=True,
            current_GBIF_backbone_may_be_newer_than_original_network_taxonomy=True,
            all_original_4505_botanical_islands_remain_full_population=True,
        ),
    )
