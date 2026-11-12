import pandas as pd
import pytest

from island_v2.fd_globi_provider_v09 import RAW, RO
from island_v2.fd_globi_species_shift_v09 import (
    ExactSpeciesGloBICoverage, select_native_source_names,
    GEO_FIELDS,
)


def original_fixture():
    species=[f"Genus{k%2} species{k}" for k in range(35)]
    ledger,flora,states,sources,gift=[],[],[],[],[]
    for name in species:
        for trait,state in [
            ("flower_primary_color","white"),
            ("floral_form","open_radial"),
            ("floral_symmetry","actinomorphic")]:
            ledger.append(dict(accepted_species=name,trait_name=trait,
                               normalized_value=state,
                               evidence_scope="species_direct",
                               evidence_quality="high"))
    for island in ("I0","I1"):
        for name in species[:30]:
            flora.append(dict(island_id=island,accepted_species=name,
                              origin_status="native"))
        states.extend([
            dict(island_id=island,channel_id="bombus",
                 channel_state="disrupted",source_state="available"),
            dict(island_id=island,channel_id="flower_visiting_birds",
                 channel_state="retained",source_state="available"),
        ])
        for rank in range(1,6):
            sources.append(dict(island_id=island,source_mode="geo_k5",
                                source_rank=rank,entity_ID=rank))
    for rank in range(1,6):
        for name in species[:12]:
            gift.append(dict(entity_ID=rank,work_species=name))
    return [pd.DataFrame(t) for t in (ledger,flora,states,sources,gift)]


def row(plant="Genus0 species0",ref="10.1/a",prov="src0",
        visitor_in_source=True,
        relation="http://purl.obolibrary.org/obo/RO_0002622",
        geo=True):
    r={k:"" for k in (*RAW,*GEO_FIELDS)}
    v="source" if visitor_in_source else "target"
    p="target" if visitor_in_source else "source"
    r["interactionTypeId"]=relation
    for prefix in (v,p):
        is_visitor=prefix==v
        r[prefix+"TaxonSpeciesName"]="Bombus terrestris" if is_visitor else plant
        r[prefix+"TaxonGenusName"]="Bombus" if is_visitor else "Genus0"
        r[prefix+"TaxonFamilyName"]="Apidae" if is_visitor else "Rosaceae"
        r[prefix+"TaxonOrderName"]="Hymenoptera" if is_visitor else "Rosales"
        r[prefix+"TaxonClassName"]="Insecta" if is_visitor else "Magnoliopsida"
        r[prefix+"TaxonKingdomName"]="Animalia" if is_visitor else "Plantae"
    r["referenceDoi"]=ref
    r["sourceArchiveURI"]=prov
    if geo:
        r["decimalLatitude"]="-4.0"
        r["decimalLongitude"]="102.1"
        r["eventDateUnixEpoch"]="1629000000"
    return r


def test_native_same_species_selection_is_not_an_evolutionary_shift():
    names,receipt=select_native_source_names(*original_fixture())
    assert len(names)==12
    assert receipt["n_native_same_name_GIFT_source_island_plant_rows"]==24
    assert receipt["n_candidate_islands_qualified_Bombus_nondetection_bird_recorded"]==2
    assert not receipt["focal_species_names_exported"]


def test_missing_one_of_five_GIFT_source_regions_fails_closed():
    frames=original_fixture()
    frames[3]=frames[3].iloc[1:]
    with pytest.raises(ValueError,match="Five candidate mainland"):
        select_native_source_names(*frames)


def test_visitor_source_and_target_roles_both_work_on_exact_species():
    q=ExactSpeciesGloBICoverage({"Genus0 species0"})
    r1=row()
    r2=row(ref="10.1/b",prov="src1",visitor_in_source=False,
           relation="http://purl.obolibrary.org/obo/RO_0002623")
    r3=row(ref="10.1/c",prov="src2")
    q.add_chunk(pd.DataFrame([r1,r2,r3,r1]))
    receipt=q.summarize()
    assert receipt["n_species_with_3plus_exact_distinct_bibliography"]==1
    assert receipt["n_species_with_2plus_provider_archive_labels"]==1
    assert receipt["n_species_with_any_georeferenced_relation"]==1
    assert receipt["n_confirmed_evolutionary_pollinator_shifts"]==0
    assert not {"plant_genus","plant_species","island_id","provider_id","reference_id"} & set(receipt)


def test_no_exact_species_name_does_not_mean_no_actual_pollinator():
    q=ExactSpeciesGloBICoverage({"Genus0 species0"})
    q.add_chunk(pd.DataFrame([row(plant="Genus0 species99")]))
    receipt=q.summarize()
    assert receipt["n_exact_species_predicate_rows"]==0
    assert receipt["n_species_with_at_least_one_exact_cited_flower_interaction"]==0
    assert receipt["n_confirmed_evolutionary_pollinator_shifts"]==0


def test_cross_archive_same_bibliography_is_not_independent():
    q=ExactSpeciesGloBICoverage({"Genus0 species0"})
    q.add_chunk(pd.DataFrame([
        row(ref="10.1/shared",prov="A"),
        row(ref="10.1/shared",prov="B"),
        row(ref="10.1/other",prov="C"),
    ]))
    receipt=q.summarize()
    assert receipt["n_crossprovider_shared_reference_keys_excluded"]==1
    assert receipt["n_species_with_at_least_one_exact_cited_flower_interaction"]==1
    assert receipt["n_species_with_3plus_exact_distinct_bibliography"]==0


def test_qualified_source_nondetection_cannot_be_imputed():
    frames=original_fixture()
    states=frames[2]
    states.loc[(states.island_id=="I0") &
               (states.channel_id=="bombus"),"source_state"]="unresolved"
    with pytest.raises(ValueError,match="Unqualified"):
        select_native_source_names(*frames)


def test_unavailable_GloBI_geography_fields_fail_closed():
    q=ExactSpeciesGloBICoverage({"Genus0 species0"})
    with pytest.raises(ValueError,match="geography schema"):
        q.add_chunk(pd.DataFrame([{k:"" for k in RAW}]))
