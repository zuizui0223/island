import pandas as pd
import pytest

from fd.trait_pair_opportunity import matched_trait_gains


def fixture():
    train = pd.DataFrame([{
        "accepted_species": f"Genus{i%3} donor{i}",
        "genus": f"Genus{i%3}",
        "self_incompatibility": i%2,
        "floral_symmetry": 1 if i%3==0 else 0,
    } for i in range(60)])
    target = pd.DataFrame([{
        "accepted_species": f"Genus{i%5} test{i}",
        "genus": f"Genus{i%5}",
        "self_incompatibility": i%2,
        "floral_symmetry": 1 if i%5==0 else 0,
    } for i in range(55)])
    return train, target, set(target.accepted_species)


def test_same_train_test_names_and_two_traits():
    donor, test, excluded = fixture()
    z = matched_trait_gains(donor, test, globally_excluded_species=excluded)
    assert set(z.trait)=={"self_incompatibility","floral_symmetry"}
    assert (z.n_train_species==60).all()
    assert (z.n_test_species==55).all()
    assert (z.donor_coverage<1).all()
    assert not {"accepted_species","genus","island_id","entity_ID"} & set(z.columns)


def test_donor_overlap_fail_closed():
    donor,test,excluded=fixture()
    excluded.add(donor.iloc[0].accepted_species)
    with pytest.raises(ValueError,match="leakage"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded)


def test_missing_and_ambiguous_states_fail():
    donor,test,excluded=fixture()
    donor.loc[0,"self_incompatibility"]=None
    with pytest.raises(ValueError,match="measured binary"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded)
    donor,test,excluded=fixture()
    test["floral_symmetry"]=test["floral_symmetry"].astype(object)
    test.loc[0,"floral_symmetry"]="multiple"
    with pytest.raises(ValueError,match="measured binary"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded)


def test_genus_consistency_and_duplicate_rejected():
    donor,test,excluded=fixture()
    donor.loc[0,"genus"]="Incorrect"
    with pytest.raises(ValueError,match="disagree"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded)
    donor,test,excluded=fixture()
    donor.loc[1,"accepted_species"]=donor.loc[0,"accepted_species"]
    with pytest.raises(ValueError,match="duplicate"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded)


def test_coverage_and_smoothing_guards():
    donor,test,excluded=fixture()
    with pytest.raises(ValueError,match="threshold"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded,
                            min_train_species=2)
    with pytest.raises(ValueError,match="positive"):
        matched_trait_gains(donor,test,globally_excluded_species=excluded,
                            tau=0)
