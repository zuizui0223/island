"""Strict positive citation availability gates only; no species/site source IDs."""
import copy
from itertools import product
import pytest
from fd.visitor_documentation_support import (
    audit_original_guild_documentation,MODES,REGIONS,SCOPES,MOTIFS,GUILDS,
)


def original_source_fixture():
    cells=[]
    for mode,reg,scope,motif,guild in product(MODES,REGIONS,SCOPES,MOTIFS,GUILDS):
        row=dict(source_mode=mode,region=reg,botanical_status=scope,
            real_species_direct_pair=motif,original_source_visitor_guild=guild,
            n_positive_cited_relation_units=0,
            n_documented_islands=0,n_documented_accepted_plant_species=0,
            n_documented_reference_labels=0,
            n_documented_archive_first_labels=0,
            n_explicit_pollination_relations=0)
        if (reg=="GLOBAL_ALL" and scope=="all_recorded"
            and motif=="actinomorphic_deep" and guild=="lepidoptera"):
            v=12 if mode=="all_cited" else 2
            row.update(n_positive_cited_relation_units=v,
                n_documented_islands=v,
                n_documented_accepted_plant_species=v,
                n_documented_reference_labels=v,
                n_documented_archive_first_labels=3 if mode=="all_cited" else 1)
        cells.append(row)
    return dict(
        status="ORIGINAL_FULL_4505_GIFT_TRUE_PAIR_GLOBI_CITED_LOCAL_GUILD_OPPORTUNITY",
        n_original_parent_GIFT_botanical_island_floras=4505,
        n_original_GloBI_rows_checked=24577183,
        old_928_locality_species_source_replay=dict(replay_passed=True),
        geometry_original_biological_source=dict(
            n_original_island_floras_with_exact_hash_rebuilt_polygon=4499),
        original_site_plant_real_morphological_guild_support=dict(
            status="FULL_FRAME_GLOBI_EXACT_SAME_ISLAND_SOURCE_MOTIF_AUDIT",
            n_original_parent_GIFT_botanical_floras=4505,
            n_source_original_islands_with_some_positive_relation=232,
            n_source_positive_original_island_plant_pairs=1799,
            n_original_source_bibliographic_keys=8646,
            cells=cells))


def test_exploratory_source_support_fails_after_global_provider_removal():
    q=audit_original_guild_documentation(original_source_fixture())
    assert q["n_original_global_botanical_island_floras"]==4505
    assert q["n_original_exact_GSHHG_polygon_flora"]==4499
    global_all=next(x for x in q["source_qualification_panels"]
                    if x["region"]=="GLOBAL_ALL"
                    and x["botanical_status"]=="all_recorded")
    assert global_all["n_preliminary_documentation_gates_passing"]==1
    assert global_all["n_provider_removal_documentation_gates_passing"]==0
    assert q["conclusions"]["no_pvalue_or_effect_coefficient_calculated"]
    assert q["conclusions"]["original_macro_populations_not_narrowed_to_928_or_232_islands"]


def test_incomplete_duplicated_and_nonmonotone_citation_provenance_rejected():
    source=original_source_fixture()
    bad=copy.deepcopy(source);bad["original_site_plant_real_morphological_guild_support"]["cells"].pop()
    with pytest.raises(ValueError,match="incomplete"):
        audit_original_guild_documentation(bad)
    bad=copy.deepcopy(source)
    c=bad["original_site_plant_real_morphological_guild_support"]["cells"]
    c[0]=copy.deepcopy(c[1])
    with pytest.raises(ValueError,match="duplicate"):
        audit_original_guild_documentation(bad)
    bad=copy.deepcopy(source)
    c=bad["original_site_plant_real_morphological_guild_support"]["cells"]
    target=next(x for x in c if x["source_mode"]=="largest_archive_excluded" and
        x["region"]=="GLOBAL_ALL" and x["botanical_status"]=="all_recorded"
        and x["real_species_direct_pair"]=="actinomorphic_deep"
        and x["original_source_visitor_guild"]=="lepidoptera")
    target["n_documented_islands"]=100
    with pytest.raises(ValueError,match="gained original"):
        audit_original_guild_documentation(bad)
    bad=copy.deepcopy(source)
    bad["old_928_locality_species_source_replay"]["replay_passed"]=False
    with pytest.raises(ValueError,match="replay absent"):
        audit_original_guild_documentation(bad)


def test_cross_provider_unknown_original_study_identity_never_becomes_independent_evidence():
    sample=original_source_fixture()
    wol=dict(
        status="ORIGINAL_GIFT_4505_EXACT_HASH_FULLFRAME_WOL_LOCAL_POLLINATOR_SUPPORT",
        parent_botanical_floras_original=4505,
        full_original_botanical_flora_island_visitor_evidence=dict(
            n_original_external_network_islands_with_plant_guild_overlap=18,
            n_original_species_network_plant_guild_units=2529))
    q=audit_original_guild_documentation(sample,wol=wol)
    other=q["other_provider_source_opportunity"]
    assert other["WOL_and_GloBI_underlying_study_overlap_unknown"]
    assert other["n_WOL_source_geographically_supported_botanical_islands"]==18
    assert other["no_sum_of_WOL_and_GloBI_sites_or_species_or_citations"]
    assert not any("total_independent_studies" in x for x in q)
    with pytest.raises(ValueError,match="parent drift"):
        audit_original_guild_documentation(sample,wol=dict(wol,parent_botanical_floras_original=927))
