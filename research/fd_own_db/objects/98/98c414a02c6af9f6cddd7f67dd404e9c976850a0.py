"""Strict Mangal pollination network node/interactions to GIFT source support.

Mangal keyword-selected sites alone do NOT establish plant-pollinator edges.
Only original published network node IDs, original-name accepted plant and
positive reported interactions of compatible types are counted. Other
endpoint's guild defaults to UNRESOLVED absent order/family provenance.
These are literature-reported interactions, NOT independently measured
effective pollen delivery, and not guaranteed original-study independence
from either GloBI or Web of Life.
"""
from __future__ import annotations

import pandas as pd
from .independent_pollination_networks import (
    original_binomial,visitor_guild
)

ALLOWED_RELATION_TYPES=frozenset({
    "mutualism","unspecified","consumption","symbiosis",
})


def strict_mangal_pollination_candidate_metadata(network:dict,dataset:dict)->bool:
    """Require source text about floral/pollination study, not just 'bees'."""
    description=" ".join(str(x or "") for x in (
        network.get("name"),network.get("description"),
        dataset.get("name"),dataset.get("description"),
    )).lower()
    return any(x in description for x in (
        "pollinat","flower visit","flower-visit","flower_vis",
        "anthophil","nectar",
    ))


def extract_mangal_positive_plant_links(
    nodes:list[dict],interactions:list[dict],
    original_focal_island_plants:set[str],
    *,network_id:str,
)->tuple[pd.DataFrame,dict]:
    """Network-specific source links with zero-source/pseudoabsence safeguards.

    Node taxonomy may supply only rank/name (not full lineage), so a partner
    with unverified family/order is retained strictly as unknown visitor guild.
    A known original GIFT plant name at the exact island is necessary.
    """
    if not isinstance(network_id,str) or not network_id:
        raise ValueError("original Mangal network ID missing")
    if not isinstance(nodes,list) or not isinstance(interactions,list):
        raise ValueError("original Mangal nodes/interactions not arrays")
    byid={}
    duplicate_node=0
    for n in nodes:
        key=str(n.get("id",""))
        if not key or key=="None":
            raise ValueError("original Mangal node without ID")
        if key in byid:
            duplicate_node+=1
            continue
        if str(n.get("network_id"))!=network_id:
            raise ValueError("Mangal node belongs to another network")
        byid[key]=n
    if duplicate_node:
        raise ValueError("original Mangal node ID duplicates")
    records=[]
    rejected_type=0
    unknown_nodes=0
    possible_edges=0
    for t in interactions:
        if str(t.get("network_id"))!=network_id:
            raise ValueError("Mangal interaction belongs to another original network")
        k=str(t.get("type") or "").lower().strip()
        if k not in ALLOWED_RELATION_TYPES:
            rejected_type+=1
            continue
        # Mangal sometimes represents observed 0 interaction values; this
        # IS NOT a documented visit and must never enter a positive count.
        if t.get("value") is not None:
            try:
                value=float(t.get("value"))
            except (TypeError,ValueError):
                unknown_nodes+=1
                continue
            if value<=0:
                continue
        a=byid.get(str(t.get("node_from")))
        b=byid.get(str(t.get("node_to")))
        if a is None or b is None:
            unknown_nodes+=1
            continue
        names=[
            original_binomial(a.get("original_name")),
            original_binomial(b.get("original_name")),
        ]
        matched=[i for i,x in enumerate(names)
                 if x is not None and x in original_focal_island_plants]
        if len(matched)!=1:
            continue
        p=matched[0]
        other=b if p==0 else a
        partner_taxonomy=other.get("taxonomy") or {}
        if not isinstance(partner_taxonomy,dict):
            partner_taxonomy={}
        rank=str(partner_taxonomy.get("rank","")).lower()
        family=(partner_taxonomy.get("name") if rank=="family" else "")
        order=(partner_taxonomy.get("name") if rank=="order" else "")
        # If source metadata explicitly provides lineage fields, use them;
        # otherwise DO NOT guess from a Latin binomial.
        family=other.get("family",family)
        order=other.get("order",order)
        guild=visitor_guild(order,family)
        gbif_key=partner_taxonomy.get("gbif")
        if gbif_key is not None and str(gbif_key).strip():
            try:
                gbif_key=int(gbif_key)
            except (TypeError,ValueError):
                gbif_key=None
            if gbif_key is not None and gbif_key<=0:
                gbif_key=None
        else:
            gbif_key=None
        records.append(dict(
            plant=names[p],
            guild=guild,
            network_id=network_id,
            interaction_id=str(t.get("id")),
            interaction_type=k,
            original_partner_rank=rank,
            original_method=str(t.get("method") or ""),
            original_partner_gbif_key=gbif_key,
            original_partner_scientific_name=str(
                partner_taxonomy.get("name") or
                other.get("original_name") or ""
            ),
        ))
        possible_edges+=1
    table=pd.DataFrame(records,columns=[
        "plant","guild","network_id","interaction_id",
        "interaction_type","original_partner_rank","original_method",
        "original_partner_gbif_key","original_partner_scientific_name",
    ])
    if len(table) and table.interaction_id.duplicated().any():
        # Distinct positive source interactions are uniquely identified.
        raise ValueError("Mangal source interaction ID duplicated")
    return table,dict(
        n_original_mangal_source_nodes=len(nodes),
        n_original_mangal_source_interactions=len(interactions),
        n_original_source_positive_original_plant_links=int(possible_edges),
        n_original_interactions_with_incompatible_types=int(rejected_type),
        n_original_interactions_with_unresolved_endpoints_or_strength=int(unknown_nodes),
        n_partner_taxonomy_classified_to_family_or_order=int(
            table.guild.ne("unresolved_visitor_taxonomy").sum()
        ),
        n_original_positive_animal_links_with_GBIF_key=int(
            table.original_partner_gbif_key.notna().sum()
        ),
        n_distinct_original_partner_GBIF_keys=int(
            table.original_partner_gbif_key.dropna().nunique()
        ),
        negative_or_zero_mangal_strength_not_a_positive_visit=True,
        source_interaction_type_not_effective_pollination=True,
    )
