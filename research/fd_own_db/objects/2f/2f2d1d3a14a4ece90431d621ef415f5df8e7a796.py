"""Original endemicity is a separate trait-opportunity axis, not a native synonym."""
import pandas as pd
import pytest

from fd.geological_endemic_trait_opportunity import (
    original_endemic_trait_opportunity, TRAITS
)

def fixture():
    islands=["A1","A2","B1","B2"]
    names=["Endema alba","Nativea rubra","Introducia foo"]
    rows=[]
    for i in islands:
        for n in names:
            isend=n.startswith("Endema")
            native=n!="Introducia foo"
            rows.append(dict(
                island_id=i,accepted_species=n,
                origin_status="native" if native else "introduced",
                endemic_status="endemic" if isend else
                    ("nonendemic" if native else "unresolved"),
                floristic_status=("endemic" if isend else
                    "native_nonendemic" if native else "introduced")))
    flora=pd.DataFrame(rows)
    geology=pd.DataFrame([
        {"island_id":i,"geological_origin":
            "oceanic" if i.startswith("A") else "continental",
         "analysis_regime":"tropical","spatial_block":i}
        for i in islands])
    near=pd.DataFrame([dict(
        island_id=i,source_mode="geo_k5",source_rank=k,
        entity_ID=f"S{k}") for i in islands for k in range(1,6)])
    source=pd.DataFrame([dict(entity_ID=f"S{k}",work_species=n) for k in range(1,6)
        for n in ("Endema minor","Nativea minor")])
    ledger=pd.DataFrame([dict(
        accepted_species=n,trait_name=t,normalized_value=value,
        evidence_scope="species_direct",evidence_quality="high"
    ) for t in ("flower_primary_color","self_incompatibility") for n,value in [
        ("Endema alba","white"),("Nativea rubra","red"),
        ("Endema minor","white"),("Nativea minor","red")
    ]])
    return ledger,flora,geology,near,source

def test_no_triple_complete_requirement_and_real_endemic_congener():
    data=fixture()
    r=original_endemic_trait_opportunity(*data)
    assert r["n_original_geological_typed_islands"]==4
    assert r["measurement_only"] and r[
        "no_actual_colonization_adaptive_selection_or_pollinator_loss_inferred"]
    assert len(r["trait_names"])==10
    for orig in ("oceanic","continental"):
        row=next(x for x in r["strict_original_status_counts_by_region_geology_trait"]
            if x["geological_origin"]==orig and x["trait"]=="flower_primary_color"
            and x["original_status_tier"]=="strict_endemic")
        assert row["n_original_high_medium_Direct_island_taxon_trait_records"]==2
        assert row["n_Direct_taxon_rows_with_comparable_same_genus_source"]==2
        assert row["n_geologically_typed_islands_with_Direct_trait"]==2
        # Different traits are NOT silently inferred from recorded colour.
        missing=next(x for x in r["strict_original_status_counts_by_region_geology_trait"]
            if x["geological_origin"]==orig and x["trait"]=="floral_symmetry"
            and x["original_status_tier"]=="strict_endemic")
        assert missing["n_original_high_medium_Direct_island_taxon_trait_records"]==0

def test_source_endemism_conflict_fails_closed():
    data=list(fixture())
    data[1].loc[0,"endemic_status"]="nonendemic"
    with pytest.raises(ValueError, match="Native/endemic"):
        original_endemic_trait_opportunity(*data)
