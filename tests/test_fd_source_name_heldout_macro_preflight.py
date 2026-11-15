"""Synthetic taxon/name holdout and proper-source scoring gates; NO original botanical IDs."""
import numpy as np
import pandas as pd
import pytest

from fd_source_name_heldout_macro_preflight import (
    genus,donor_probability,score,exact_join_and_source,predict_one_cohort,
    short_distance_block_disjoint_endpoint_union
)


def _data(n_islands=40):
    flowers={}
    all_flora={}
    native_flora={}
    geog=[]
    original=[]
    assigns=[]
    # Test names are distinct from all eligible additional Direct donor plants.
    for i in range(n_islands):
        island=f"I{i}"
        names={f"Gen{k%3} island{i}_{k}" for k in range(5)}
        all_flora[island]=names
        native_flora[island]=names
        for k in range(5):
            flowers[f"Gen{k%3} island{i}_{k}"]=(
                "white" if k%2 else "red","open" if i%2 else "tube")
        geog.append(dict(island_id=island,spatial_block=f"B{i//2}"))
        assigns += [
            dict(island_id=island,source_mode="geo_k5",source_rank=str(j),
                 entity_ID=f"S{(i+j)%5}")
            for j in range(1,6)]
    # The synthetic mainland candidate source contains all heldout target
    # names, plus 30 genuine different original plant names.
    for e in range(5):
        for s in flowers:
            original.append(dict(entity_ID=f"S{e}",work_species=s))
        for j in range(30):
            n=f"Gen{j%3} donor{j}"
            flowers[n]=("red" if j%2 else "white",
                        "tube" if j%3 else "open")
            original.append(dict(entity_ID=f"S{e}",work_species=n))
    mainland=pd.DataFrame(original).drop_duplicates()
    assigned=pd.DataFrame(assigns)
    return flowers,all_flora,native_flora,mainland,assigned,pd.DataFrame(geog)


def test_full_heldout_names_cannot_be_donors_at_any_source_region():
    labels,all_f,native_f,main,asgn,geo=_data()
    donor,near=exact_join_and_source(main,asgn,labels)
    q=predict_one_cohort("test",sorted(all_f),all_f,native_f,
                         labels,donor,near,geo,test_scope="all_recorded")
    assert q["n_original_focal_islands"]==40
    assert q["n_original_focal_full_flora_names"]==200
    assert q["n_unique_source_direct_pair_donor_species_after_exact_name_holdout"]==30
    assert q["n_unique_source_direct_pair_donor_species_before_holdout"]==230
    assert q["n_source_entities_with_some_direct_pair_donor"]==5
    assert q["n_evaluated_original_test_names"]==200
    assert q["ratio_test_carriers_with_same_genus_local_five_source_donor"]==1
    assert q["status"]=="POSTHOC_OUT_OF_TARGET_NAME_SOURCE_PREDICTION"
    assert q["n_scored_physical_blocks"]==20
    assert q["source_local_minus_global_genus_improvement_island_log_loss"]==pytest.approx(0)
    assert q["source_local_minus_global_genus_improvement_island_brier"]==pytest.approx(0)


def test_no_post_holdout_donors_abstains_instead_of_fabricating_profiles():
    labels,all_f,native_f,main,asgn,geo=_data()
    main=main.loc[~main.work_species.str.contains("donor")]
    donor,near=exact_join_and_source(main,asgn,labels)
    q=predict_one_cohort("no-donors",sorted(all_f),all_f,native_f,
                         labels,donor,near,geo,test_scope="all_recorded")
    assert q["status"]=="NOT_IDENTIFIABLE_NO_DIRECT_PAIR_GIFT_DONORS_AFTER_HOLDOUT"
    assert q["n_unique_source_direct_pair_donor_species_after_exact_name_holdout"]==0


def test_exact_native_scope_does_not_promote_unresolved_plants():
    labels,all_f,native_f,main,asgn,geo=_data()
    donor,near=exact_join_and_source(main,asgn,labels)
    for island in native_f:
        native_f[island]=set(sorted(native_f[island])[:3])
    q=predict_one_cohort("native",sorted(all_f),all_f,native_f,
                         labels,donor,near,geo,test_scope="original_native")
    assert q["n_original_focal_full_flora_names"]==200
    assert q["n_test_direct_colour_form_species_carriers"]==120
    assert q["n_evaluated_original_test_names"]==120


def test_proper_positive_joint_probability_and_identity_gates():
    prior=np.array([.25,.75])
    q=donor_probability(np.array([10,0]),prior)
    assert q[0]>.25 and q[1]>.0
    l,b=score(q,0)
    assert l>0 and b>=0
    with pytest.raises(ValueError): donor_probability(np.array([1,2]),np.array([1,0]))
    with pytest.raises(ValueError): score(np.array([1,0]),0)
    with pytest.raises(ValueError): score(q,2)
    assert genus("A b")=="A"
    assert genus("SingleSpecies")=="UNRESOLVED:SingleSpecies"


def test_duplicate_source_and_rank_leakage_fail_closed():
    labels,all_f,native_f,main,asgn,geo=_data()
    with pytest.raises(ValueError):
        exact_join_and_source(pd.concat([main,main.iloc[:1]]),asgn,labels)
    with pytest.raises(ValueError):
        exact_join_and_source(main,pd.concat([asgn,asgn.iloc[:1]]),labels)


def test_unrelated_work_species_does_not_inflate_profile_donors():
    labels,all_f,native_f,main,asgn,geo=_data()
    donor,near=exact_join_and_source(
        pd.concat([main,pd.DataFrame([{
            "entity_ID":"S0","work_species":"Mystery unknown"}])],
            ignore_index=True),
        asgn,labels)
    q=predict_one_cohort("unknown-outcomes",sorted(all_f),
                         all_f,native_f,labels,donor,near,geo,
                         test_scope="all_recorded")
    assert q["n_unique_source_direct_pair_donor_species_before_holdout"]==230


def test_outcome_blind_original_distance_source_mask_selects_exact_pair_endpoints():
    full={
        'A':{'One a','One b'},
        'B':{'Two a','Two b'},
        'C':{'Three a','Three b'},
        'D':{'Four a','Four b'},
    }
    geo=pd.DataFrame([
        dict(island_id=i,analysis_regime='northern_midlatitude',
             spatial_block=f'block{i}',island_latitude=0,
             island_longitude=j*10)
        for j,i in enumerate(['A','B','C','D'])])
    endpoints,n_pairs,n_blocks=short_distance_block_disjoint_endpoint_union(
        list(full),full,geo)
    assert (endpoints,n_pairs,n_blocks)==(['A','B','C','D'],3,4)
    full['A'].add('Shared source plant')
    full['B'].add('Shared source plant')
    endpoints,n_pairs,n_blocks=short_distance_block_disjoint_endpoint_union(
        list(full),full,geo)
    assert n_pairs==2
    assert endpoints==['B','C','D']
    with pytest.raises(ValueError):
        short_distance_block_disjoint_endpoint_union(
            list(full),full,geo.drop(columns=['island_latitude']))
