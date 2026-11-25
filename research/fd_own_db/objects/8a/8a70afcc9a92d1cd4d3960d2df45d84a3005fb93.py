"""Mangal actual pollination network nodes and positive edges, no invented guild."""
import pytest
from fd.independent_mangal_pollination_links import (
    strict_mangal_pollination_candidate_metadata,
    extract_mangal_positive_plant_links,
)


def test_keyword_bees_alone_not_measured_flower_visitation():
    assert strict_mangal_pollination_candidate_metadata(
        {"name":"island_net"}, {"description":"Flower visitation by insects"}
    )
    assert not strict_mangal_pollination_candidate_metadata(
        {"name":"bees"}, {"description":"Bee larval parasitism"}
    )


def source():
    nodes=[
        dict(id=1,network_id=7,original_name="Impatiens balsamina",
             taxonomy=dict(name="Impatiens balsamina",rank="species")),
        dict(id=2,network_id=7,original_name="Hyles gallii",
             taxonomy=dict(name="Sphingidae",rank="family")),
        dict(id=3,network_id=7,original_name="visitor unidentified",
             taxonomy=dict(name="Animalia",rank="kingdom")),
        dict(id=4,network_id=7,original_name="Camellia japonica",
             taxonomy=dict(name="Camellia japonica",rank="species")),
    ]
    edges=[
        dict(id=10,network_id=7,node_from=1,node_to=2,type="mutualism",
             method="observation",value=1),
        dict(id=11,network_id=7,node_from=1,node_to=3,type="mutualism",
             method="observation",value=3),
        dict(id=12,network_id=7,node_from=1,node_to=2,type="predation",
             method="observation",value=10),
        dict(id=13,network_id=7,node_from=1,node_to=2,type="mutualism",
             method="observation",value=0),
        dict(id=14,network_id=7,node_from=4,node_to=3,type="mutualism",
             method="biblio",value=1),
    ]
    return nodes,edges


def test_source_nodes_and_interactions_must_be_positive_and_plant_exact():
    nodes,edges=source()
    raw,counts=extract_mangal_positive_plant_links(
        nodes,edges,{"Impatiens balsamina"},network_id="7"
    )
    assert len(raw)==2
    assert set(raw.guild)=={"sphingidae_moth","unresolved_visitor_taxonomy"}
    assert counts["n_original_interactions_with_incompatible_types"]==1
    assert counts["n_original_source_positive_original_plant_links"]==2
    assert counts["n_partner_taxonomy_classified_to_family_or_order"]==1
    assert counts["negative_or_zero_mangal_strength_not_a_positive_visit"]


def test_unmatched_island_plant_is_not_source_matched():
    nodes,edges=source()
    a,meta=extract_mangal_positive_plant_links(
        nodes,edges,{"Trifolium repens"},network_id="7")
    assert a.empty and meta["n_original_source_positive_original_plant_links"]==0


def test_cross_network_node_or_interaction_is_fail_closed():
    nodes,edges=source()
    nodes[0]["network_id"]=9
    with pytest.raises(ValueError,match="another network"):
        extract_mangal_positive_plant_links(
            nodes,edges,{"Impatiens balsamina"},network_id="7"
        )
    nodes,edges=source()
    edges[0]["network_id"]=9
    with pytest.raises(ValueError,match="another original network"):
        extract_mangal_positive_plant_links(
            nodes,edges,{"Impatiens balsamina"},network_id="7"
        )


def test_unknown_source_guild_must_remain_unclassified():
    nodes,edges=source()
    nodes[1]["taxonomy"]=dict(name="Hyles gallii",rank="species")
    a,_=extract_mangal_positive_plant_links(
        nodes,edges,{"Impatiens balsamina"},network_id="7"
    )
    assert all(x=="unresolved_visitor_taxonomy" for x in a.guild)
