"""NHM constituent text-source provenance and strict same-species agreement."""
import pandas as pd
import pytest

from fd_nhm_original_provider_audit import (
    consensus_provider_profiles,validate_external_labels,
    read_provider_file,analyze,
)
from fd_nhm_external_flora_preflight import EXTERNAL_COLUMNS


def test_corpus_consensus_requires_actual_distinct_corpora_and_never_votes_conflict():
    data={
        "BHL_text":{"primary":{"A first":"white","B second":"yellow_orange","C third":"red_pink"},
                    "expanded":{},"form":{},"symmetry":{}},
        "FloraChina_text":{"primary":{"A first":"white","B second":"red_pink"},
                    "expanded":{},"form":{},"symmetry":{}},
        "Ecoflora_text":{"primary":{"A first":"white"},
                    "expanded":{},"form":{},"symmetry":{}},
    }
    consensus,two,conflicts,report=consensus_provider_profiles(data)
    assert consensus=={"A first":"white","C third":"red_pink"}
    assert two=={"A first":"white"}
    assert conflicts==1
    assert len(report["A first"])==3


def test_original_island_name_excluded_from_external_source_donors_even_if_consistent():
    data={
        "BHL_text":{"primary":{"A first":"white","B second":"yellow_orange",
                               "C third":"red_pink"},
                    "expanded":{"A first":"white"},"form":{},"symmetry":{}},
        "FloraChina_text":{"primary":{"A first":"white","B second":"yellow_orange"},
                    "expanded":{"A first":"white"},"form":{},"symmetry":{}},
    }
    original={"primary":{"A first":"white","C third":"white"},
              "expanded":{"A first":"white"},"form":{},"symmetry":{}}
    r=validate_external_labels(
        data,original,{"A first"},{"A first":1,"B second":1,"C third":1})
    q=r["primary"]
    assert q["n_independent_NHM_provider_unequivocal_names"]==3
    assert q["n_consistent_two_or_more_source_corpus_names"]==2
    assert q["n_in_GIFT_mainland_not_in_original_4505_islands"]==2
    assert q["n_two_or_more_corpora_in_GIFT_outside_original_4505"]==1
    assert q["n_exact_original_Direct_species_names_comparable"]==2
    assert q["n_agree_original_direct_category"]==1


def test_small_source_file_allows_optional_missing_columns_without_faking_source(tmp_path):
    p=tmp_path/"original.csv"
    p.write_text("taxon,flower colour\nA first,red\nB second,white\n")
    r=read_provider_file(p,"BHL_text")
    assert r["source"].unique().tolist()==["BHL_text"]
    assert set(EXTERNAL_COLUMNS)<=set(r.columns)
    assert r["flower shape"].eq("").all()
    assert len(r)==2


def test_public_construction_fails_if_private_source_islands_change():
    raw={"BHL_text":pd.DataFrame([{
        **{x:"" for x in EXTERNAL_COLUMNS},
        "taxon":"A first","flower colour":"white","source":"BHL_text"
    }])}
    original={"flower_primary_color":{"A first":"white"},
              "floral_form":{},"floral_symmetry":{}}
    flora=pd.DataFrame([{
        "island_id":"IslandA","accepted_species":"A first"
    }])
    gift=pd.DataFrame([{"entity_ID":"X","work_species":"A first"}])
    with pytest.raises(ValueError):
        analyze(raw,original,flora,gift)
