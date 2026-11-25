"""Independent exact name versus unpromoted CoL XR synonym-candidate checks."""
import pandas as pd
import pytest
from fd.fcp_endemic_colxr_synonym_gate import (
    original_candidate_fcp_name_map,
    original_endemic_fcp_taxonomy_sensitivity,
    COLXR_CHECKLIST
)


def _source():
    fcp={"Photo exact","Photo accepted","Photo collision"}
    fcp|={f"Extra{i} specimen" for i in range(4472)}
    assert len(fcp)==4475
    names=[f"Source{i} taxon" for i in range(115322)]+[
        "Photo exact","Original synonym","Collision synonym1",
        "Collision synonym2","Inexact candidate","Unmatched record"]
    assert len(names)==115328
    col=pd.DataFrame({
        "submitted_name":names,
        "candidate_accepted_name":[""]*len(names),
        "resolution_status":["unmatched"]*len(names),
        "automatic_resolution_candidate":["false"]*len(names),
        "candidate_kingdom":["Plantae"]*len(names),
        "candidate_accepted_rank":["SPECIES"]*len(names),
        "candidate_genus":[""]*len(names),
        "match_type":["NONE"]*len(names),
        "matched_status":[""]*len(names),
        "confidence":["0"]*len(names),
        "processing_flags":[""]*len(names),
        "issues":[""]*len(names),
        "checklist_key":[COLXR_CHECKLIST]*len(names),
    })
    for name,target in (
        ("Original synonym","Photo accepted"),
        ("Collision synonym1","Photo collision"),
        ("Collision synonym2","Photo collision"),
        ("Inexact candidate","Photo accepted"),
    ):
        idx=col.index[col.submitted_name.eq(name)][0]
        col.loc[idx,["candidate_accepted_name","resolution_status",
            "automatic_resolution_candidate","match_type",
            "matched_status","confidence"]]=[
            target,"exact_synonym_to_accepted_candidate",
            "true","EXACT","SYNONYM","99"]
    col.loc[col.submitted_name.eq("Inexact candidate"),"match_type"]="VARIANT"
    return fcp,col


def test_exact_vs_unpromoted_one_to_one_colxr_candidates_and_collision():
    fcp,col=_source()
    exact,candidate,audit=original_candidate_fcp_name_map(fcp,col)
    assert exact=={"Photo exact":"Photo exact"}
    assert candidate=={"Original synonym":"Photo accepted"}
    assert audit["n_photo_name_synonym_collisions_quarantined"]==2
    assert audit["no_synonym_candidate_promoted_to_accepted_name"]


def test_endemism_counts_by_origin_without_name_export():
    fcp,col=_source()
    flora=[]
    for origin,island in (("oceanic","isl1"),("continental","isl2")):
        for name in ("Original synonym","Photo exact","Collision synonym1"):
            flora.append(dict(island_id=island,accepted_species=name,
                              origin_status="native",endemic_status="endemic",
                              floristic_status="endemic"))
    geo=pd.DataFrame([
        dict(island_id="isl1",geological_origin="oceanic",
             analysis_regime="tropical",spatial_block="block1"),
        dict(island_id="isl2",geological_origin="continental",
             analysis_regime="tropical",spatial_block="block2"),
        dict(island_id="isl3",geological_origin="unresolved",
             analysis_regime="tropical",spatial_block="block3"),
        dict(island_id="isl4",geological_origin="unresolved",
             analysis_regime="tropical",spatial_block="block4"),
    ])
    out=original_endemic_fcp_taxonomy_sensitivity(fcp,pd.DataFrame(flora),geo,col)
    for group in out["geological_native_endemic_candidate_overlap"]:
        assert group["n_original_exact_FCP_measured_species"]==1
        assert group["n_one_to_one_CoL_synonym_candidate_FCP_species"]==1
        assert "Original synonym" not in str(group)
        assert group["candidate_requires_independent_botanical_review"]


def test_acceptance_never_promoted_without_strict_match():
    fcp,col=_source()
    col.loc[col.submitted_name.eq("Original synonym"),"issues"]="ambiguous"
    _,cand,_=original_candidate_fcp_name_map(fcp,col)
    assert not cand
    col.loc[col.submitted_name.eq("Original synonym"),"issues"]=""
    col.loc[col.submitted_name.eq("Original synonym"),"confidence"]="94"
    _,cand,_=original_candidate_fcp_name_map(fcp,col)
    assert not cand
