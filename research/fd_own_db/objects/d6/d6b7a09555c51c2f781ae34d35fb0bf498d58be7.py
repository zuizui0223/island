"""Pure synthetic COL XR source taxonomy strict family validation.

Never invent family from genus or treat a candidate match as a reviewed
accepted taxonomy classification or a phylogeny.
"""
import pandas as pd
import pytest
from fd.colxr_family_opportunity import (
    strict_candidate_family_map, group_by_family,COLXR_EXPECTED_CHECKLIST
)


def _crosswalk():
    rows=[]
    for s,g,f,candidate_f,status in (
        ("Alpha one","Alpha","Aaceae","Aaceae","exact_accepted_candidate"),
        ("Alpha two","Alpha","Aaceae","Aaceae","exact_accepted_candidate"),
        ("Beta one","Beta","Baceae","Baceae","exact_accepted_candidate"),
        ("Gamma one","Gamma","Caceae","Caceae","review_required"),
        ("Delta one","Delta","Daceae","Qaceae","exact_accepted_candidate"),
        ("Epsilon one","Epsilon","Eaceae","Eaceae","exact_synonym_to_accepted_candidate"),
    ):
        rows.append({
            "submitted_name":s,"submitted_family":f,"candidate_family":candidate_f,
            "candidate_genus":g,"candidate_kingdom":"Plantae",
            "candidate_accepted_name":s,"resolution_status":status,
            "match_type":"EXACT","automatic_resolution_candidate":"true",
            "checklist_key":COLXR_EXPECTED_CHECKLIST,
        })
    return pd.DataFrame(rows)


def test_strict_family_not_genus_guessed_and_only_col_exact_accepted_candidates():
    source=_crosswalk()
    result,audit=strict_candidate_family_map(source,set(source.submitted_name),expected_rows=6)
    assert result=={"Alpha one":"Aaceae","Alpha two":"Aaceae","Beta one":"Baceae"}
    assert audit["n_strict_family_candidates"]==3
    assert audit["n_candidate_families"]==2
    assert audit["candidate_taxonomy_not_verified_reviewed_final"] is True
    assert group_by_family({"Alpha one","Beta one","Gamma one"},result)=={
        "Aaceae":("Alpha one",),"Baceae":("Beta one",)
    }


def test_conflicting_family_per_genus_quarantines_all_its_species():
    source=_crosswalk()
    source.loc[source.submitted_name.eq("Alpha two"),"candidate_family"]="Baceae"
    source.loc[source.submitted_name.eq("Alpha two"),"submitted_family"]="Baceae"
    result,audit=strict_candidate_family_map(source,set(source.submitted_name),expected_rows=6)
    assert "Alpha one" not in result and "Alpha two" not in result
    assert audit["n_genus_family_conflicts"]==1
    assert result=={"Beta one":"Baceae"}


def test_duplicate_and_changed_taxonomy_source_fail_closed():
    source=_crosswalk()
    duplicated=pd.concat([source,source.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="pinned"):
        strict_candidate_family_map(duplicated,set(source.submitted_name),expected_rows=6)
    source.loc[0,"checklist_key"]="different"
    with pytest.raises(ValueError,match="checklist"):
        strict_candidate_family_map(source,set(source.submitted_name),expected_rows=6)
