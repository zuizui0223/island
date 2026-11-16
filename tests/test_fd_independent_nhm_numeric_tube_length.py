"""Synthetic strictly numeric corolla tube length validation and provenance."""
import numpy as np
import pandas as pd
import pytest

from run_fd_independent_nhm_numeric_tube_length import (
    strict_single_length,parse_nhm_numeric_file,
    combine_external_provider_lengths
)


def test_strict_cm_numeric_range_and_OCR_ambiguity():
    assert np.isclose(strict_single_length("1.0","3.0"),2.)
    assert np.isclose(strict_single_length("1.5",""),1.5)
    assert strict_single_length("3.0","1.0") is None
    assert strict_single_length("1.0-3.0","") is None
    assert strict_single_length("NA","") is None
    assert strict_single_length("nan","nan") is None
    assert strict_single_length("0","") is None
    assert strict_single_length("101","") is None
    assert np.isclose(strict_single_length("1e-1",""),.1)


def test_parser_exact_binomial_not_aggregate_species(tmp_path):
    path=tmp_path/"fake_provider.csv"
    pd.DataFrame([
        {"taxon":"Alpha beta","corolla tube min. length [cm]":"1",
         "corolla tube max. length [cm]":"3"},
        {"taxon":"Alpha beta","corolla tube min. length [cm]":"2",
         "corolla tube max. length [cm]":"4"},
        {"taxon":"Alpha unknown","corolla tube min. length [cm]":"5",
         "corolla tube max. length [cm]":"3"},
        {"taxon":"Genus spp.","corolla tube min. length [cm]":"1",
         "corolla tube max. length [cm]":"1"},
        {"taxon":"Beta gamma","corolla tube min. length [cm]":"",
         "corolla tube max. length [cm]":"2"},
    ]).to_csv(path,index=False)
    d,a=parse_nhm_numeric_file(path,"synthetic_corpus")
    assert d=={"Alpha beta":2.5,"Beta gamma":2.}
    assert a["n_raw_rows_with_numeric_tube_length"]==3
    assert a["n_unique_exact_names_with_numeric_length"]==2
    assert a["n_invalid_numeric_range_rows"]==1
    assert a["n_names_with_multiple_numeric_entries"]==1


def test_cross_provider_agreement_is_strict_and_source_independent():
    pooled,consistent,summary=combine_external_provider_lengths({
        "A":{"Alpha beta":1.0,"Beta gamma":2.0},
        "B":{"Alpha beta":1.2,"Beta gamma":8.0},
        "C":{"Alpha beta":1.1},
    })
    assert np.isclose(pooled["Alpha beta"],1.1)
    assert np.isclose(consistent["Alpha beta"],1.1)
    assert "Beta gamma" not in consistent
    assert summary["n_original_names_recorded_by_two_consistent_corpora"]==1
    assert summary["n_multi_corpus_gt_1p5_max_to_min_conflicts"]==1


def test_invalid_external_numeric_cm_fails_closed():
    with pytest.raises(ValueError):
        combine_external_provider_lengths({"A":{"Alpha beta":float("nan")}})
