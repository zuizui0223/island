import pandas as pd
import pytest
from island_v2.fd_globi_exact_species_syndrome_cv_v09 import (
    species_family_documentation_syndrome_transfer as audit,
)

def evidence(n=90):
    ledger=[]; taxa=[]; units=set(); byref={}; flora=set()
    for i in range(n):
        genus=f"Plant{i}"
        plant=f"{genus} species{i}"
        flora.add(plant); taxa.append({"genus":genus,"family":f"Family{i//3}"})
        flower=i%2==0
        for trait,state in (
            ("flower_primary_color","white" if i%3 else "green_brown_inconspicuous"),
            ("floral_form","open_radial" if flower else "tubular"),
            ("floral_symmetry","actinomorphic" if flower else "zygomorphic"),
        ):
            ledger.append({"accepted_species":plant,"trait_name":trait,
                "normalized_value":state,"evidence_scope":"species_direct",
                "evidence_quality":"high"})
        for k in range(3):
            ref=f"secondary/{i}/{k}"
            guild=("flower_visiting_birds" if k<(2 if flower else 1)
                   else "non_bombus_bees")
            units.add((plant,guild,ref,"secondary","flower_visit"))
            byref[ref]={"secondary"}
        for k in range(6):
            ref=f"dominant/{i}/{k}"
            units.add((plant,"lepidoptera" if k%2 else "non_bombus_bees",
                       ref,"dominant","pollination"))
            byref[ref]={"dominant"}
    return units,byref,pd.DataFrame(ledger),pd.DataFrame(taxa),flora

def test_whole_family_source_matched_real_compositional_improvement():
    args=evidence()
    r=audit(*args,bootstrap=299)
    assert r["status"].startswith("COMPLETED_")
    assert r["n_strict_common_species_with_3plus_documentation"]==90
    assert r["n_strict_common_botanical_families"]==30
    assert r["same_species_same_family_before_after_sample"]
    assert r["models_after_source_label_omission"]["plus_open_flower"][
        "species_equal_gain_over_source_effort"]>0
    assert not r["source_location_visitor_effectiveness_tested"]
    assert not r["evolutionary_pollinator_shift_confirmed"]
    assert not r["original_plant_family_archive_or_citation_identifiers_released"]

def test_documentation_loss_not_falsely_zero_biological_effect():
    units,refs,ledger,taxonomy,flora=evidence()
    units={u for u in units if u[3]=="dominant"}
    refs={ref:p for ref,p in refs.items() if p=={"dominant"}}
    r=audit(units,refs,ledger,taxonomy,flora,bootstrap=299)
    assert r["status"]=="NOT_IDENTIFIABLE_AFTER_ORIGINAL_SOURCE_OMISSION"
    assert r["n_strict_common_species_with_3plus_documentation"]==0

def test_conflicting_direct_evidence_rejected():
    units,refs,ledger,taxonomy,flora=evidence()
    ledger=pd.concat([ledger,ledger.iloc[[0]]],ignore_index=True)
    with pytest.raises(ValueError,match="conflict"):
        audit(units,refs,ledger,taxonomy,flora,bootstrap=299)

def test_conflicting_botanical_family_discarded_without_majority_assumption():
    units,refs,ledger,taxonomy,flora=evidence()
    taxonomy=pd.concat([taxonomy,pd.DataFrame([{
        "genus":"Plant0","family":"AnotherFamily",
    }])],ignore_index=True)
    r=audit(units,refs,ledger,taxonomy,flora,bootstrap=299)
    assert r["n_strict_common_species_with_3plus_documentation"]==89

def test_deduplicate_original_citation_shared_across_archives():
    units,refs,ledger,taxonomy,flora=evidence()
    old=next(u for u in units if u[0]=="Plant0 species0"
             and u[2]=="secondary/0/0")
    units.add((old[0],old[1],old[2],"dominant",old[4]))
    refs["secondary/0/0"]={"dominant","secondary"}
    r=audit(units,refs,ledger,taxonomy,flora,bootstrap=299)
    assert r["n_crossarchive_bibliographic_key_collisions_removed"]==1
    assert r["n_strict_common_species_with_3plus_documentation"]==89
