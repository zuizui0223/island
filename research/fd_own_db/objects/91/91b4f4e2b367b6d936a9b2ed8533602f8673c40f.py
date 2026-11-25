import pandas as pd
import pytest
from fd.trait_direct_evidence_priority import prospective_sc_donor_coverage


def fixture():
    names = [
        {"accepted_species":f"Genus{i//10} test{i}","genus":f"Genus{i//10}"}
        for i in range(100)
    ]
    donors = [{"accepted_species":"Genus0 donor0","genus":"Genus0"}]
    candidates = [
        {"accepted_species":f"Genus{j} prospect0","genus":f"Genus{j}"}
        for j in range(1,8)
    ]
    return (pd.DataFrame(names),pd.DataFrame(donors),
            pd.DataFrame(candidates),
            {t["accepted_species"] for t in names})


def test_upper_bound_only_no_sc_imputation():
    test, donor, candidate, heldout=fixture()
    z=prospective_sc_donor_coverage(
        test,donor,candidate,all_island_species=heldout)
    assert len(z)==1
    assert z.iloc[0].current_same_genus_coverage==pytest.approx(.1)
    assert z.iloc[0].hypothetical_coverage_after_5_new_genera==pytest.approx(.6)
    assert z.iloc[0].max_candidate_supported_coverage==pytest.approx(.8)
    assert not {"genus","accepted_species","state","source_id"} & set(z.columns)


def test_global_species_name_holdout_and_no_reused_sc_donor():
    test, donor, candidate, heldout=fixture()
    with pytest.raises(ValueError,match="holdout"):
        prospective_sc_donor_coverage(
            test,donor,candidate,all_island_species=set())
    candidate=pd.concat([candidate,donor],ignore_index=True)
    with pytest.raises(ValueError,match="already has"):
        prospective_sc_donor_coverage(
            test,donor,candidate,all_island_species=heldout)


def test_privacy_and_duplicate_guards():
    test, donor, candidate, heldout=fixture()
    with pytest.raises(ValueError,match="Privacy"):
        prospective_sc_donor_coverage(
            test,donor,candidate,all_island_species=heldout,
            minimum_targets=2)
    duplicate=pd.concat([candidate,candidate.iloc[[0]]])
    with pytest.raises(ValueError,match="duplicate"):
        prospective_sc_donor_coverage(
            test,donor,duplicate,all_island_species=heldout)
