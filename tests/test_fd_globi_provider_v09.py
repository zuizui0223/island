import pandas as pd
import pytest

from island_v2.fd_globi_provider_v09 import (
    RAW, canonicalize_chunk, ProviderReferenceEvidence
)

def raw(interaction="http://purl.obolibrary.org/obo/RO_0002622",
        inverted=False, reference="10.5555/test", provider="archive:original"):
    r={c:"" for c in RAW}
    for role in ("source","target"):
        is_visitor=(role=="source") != inverted
        r[role+"TaxonGenusName"]="Bombus" if is_visitor else "Campanula"
        r[role+"TaxonKingdomName"]="Animalia" if is_visitor else "Plantae"
        r[role+"TaxonSpeciesName"]="Bombus terrestris" if is_visitor else "Campanula punctata"
        r[role+"TaxonFamilyName"]="Apidae" if is_visitor else "Campanulaceae"
        r[role+"TaxonOrderName"]="Hymenoptera" if is_visitor else "Asterales"
        r[role+"TaxonClassName"]="Insecta" if is_visitor else "Magnoliopsida"
    r["interactionTypeId"]=interaction
    r["sourceArchiveURI"]=provider
    r["referenceDoi"]=reference
    return r


def test_direction_and_reference_source_semantics():
    a=raw()
    b=raw(interaction="http://purl.obolibrary.org/obo/RO_0002623",inverted=True)
    c=raw(interaction="http://purl.obolibrary.org/obo/RO_0002455",reference="10.5555/poll")
    frame,receipt=canonicalize_chunk(pd.DataFrame([a,b,c]),{"Campanula"})
    assert len(frame)==2
    assert set(frame.visitor_channel)=={"bombus"}
    assert set(frame.evidence_class)=={"flower_visit","pollination"}
    assert receipt["raw_input_rows"]==3


def test_provider_missing_cannot_be_synthesized_from_citation():
    a=raw();a["sourceArchiveURI"]=""
    b=raw(reference="");b["referenceCitation"]="";b["referenceUrl"]=""
    d,rec=canonicalize_chunk(pd.DataFrame([a,b]))
    assert d.empty
    assert rec["missing_source_provider"]==1
    assert rec["missing_bibliographic_reference"]==1


def test_source_dedup_and_citation_collision_excluded():
    rows=[]
    for j in range(3):
        for k in range(3 if j==2 else 2):
            rows.append(raw(provider=f"provider{j}",reference=f"10.1111/{j}-{k}"))
    rows.append(raw(provider="provider2",reference="10.1111/0-0"))
    d,_=canonicalize_chunk(pd.DataFrame(rows))
    a=ProviderReferenceEvidence()
    a.add(d)
    original, after, rec=a.summarize()
    assert rec["n_distinct_source_providers"]==3
    assert rec["n_crossprovider_repeated_references_excluded"]==1
    assert rec["n_genera_with_3plus_refs_before"]==1
    assert rec["n_genera_with_3plus_refs_after"]==1
    assert after.to_numpy().sum()<original.to_numpy().sum()
    assert "provider_id" not in rec


def test_nonplant_and_unqualified_visits_are_excluded():
    x=raw();x["targetTaxonKingdomName"]="Fungi"
    q,receipt=canonicalize_chunk(pd.DataFrame([x]))
    assert q.empty and receipt["nonplant_or_unknown_kingdom"]==1
    x=raw(interaction="http://purl.obolibrary.org/obo/RO_9999999")
    q,_=canonicalize_chunk(pd.DataFrame([x]))
    assert q.empty


def test_missing_original_column_fails_closed():
    with pytest.raises(ValueError,match="Missing original"):
        canonicalize_chunk(pd.DataFrame({"referenceDoi":["A"]}))
