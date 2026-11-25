import pandas as pd
import pytest
from fd.source_disjoint_trait_cv import holdout_genus_transfer


def fixtures():
    regions=[];sources=[];flora=[];gift=[];trait=[]
    for i in range(80):
        group="A" if i<40 else "B"
        regions.append({"island_id":f"I{i}","analysis_regime":"north"})
        for rank in range(1,6):
            sources.append({"island_id":f"I{i}","source_mode":"geo_k5",
                            "source_rank":rank,"entity_ID":f"{group}{rank}"})
        for j in range(3):
            genus="FlowerA" if j==0 else "FlowerB"
            name=f"{genus} test{i}_{j}"
            flora.append({"island_id":f"I{i}","accepted_species":name,
                          "origin_status":"native" if j<2 else "unresolved"})
            trait.append({"accepted_species":name,"trait_name":"floral_symmetry",
                          "normalized_value":"actinomorphic" if j==0 else "zygomorphic"})
    for i in range(100):
        genus="FlowerA" if i%2==0 else "FlowerB"
        name=f"{genus} reference{i}"
        gift.append({"entity_ID":"External","work_species":name})
        trait.append({"accepted_species":name,"trait_name":"floral_symmetry",
                      "normalized_value":"actinomorphic" if i%2==0 else "zygomorphic"})
    # A separate source may also report the evaluation species: strict
    # global name holdout must exclude them even from that OTHER region.
    for f in flora:
        gift.append({"entity_ID":"External","work_species":f["accepted_species"]})
    return [pd.DataFrame(x) for x in (regions,sources,flora,gift,trait)]


def run(d,policy="all_g0"):
    return holdout_genus_transfer(*d,ledger_label="synthetic",
                                  scope_label="Direct",holdout_policy=policy)


def test_strict_all_island_names_excluded_from_donors():
    z=run(fixtures())
    assert len(z)==4  # two donor components x all/native targets
    assert (z.n_islands==40).all()
    assert set(z.holdout_policy)=={"all_g0"}
    assert (z.brier_gain_tau_1>0).all()
    assert (z.n_train_source_disjoint_trait_species==100).all()
    assert not {"island_id","accepted_species","entity_ID","genus"} & set(z.columns)
    # Weaker component-only name holdout admits test species from the
    # OTHER component as training donors, and must be clearly distinguished.
    weak=run(fixtures(),policy="component")
    assert (weak.n_train_source_disjoint_trait_species==220).all()


def test_five_source_regions_required():
    f=fixtures();f[1]=f[1].iloc[1:]
    with pytest.raises(ValueError,match="Incomplete"):
        run(f)


def test_duplicate_island_taxon_fails():
    f=fixtures();f[2]=pd.concat([f[2],f[2].iloc[[0]]])
    with pytest.raises(ValueError,match="Duplicate"):
        run(f)


def test_ambiguous_trait_is_missing_not_zero():
    f=fixtures()
    f[4].loc[f[4].accepted_species.str.startswith("FlowerA test"),
             "normalized_value"]="unknown"
    z=run(f)
    assert z.empty or (z.n_test_species<=80).all()


def test_privacy_and_holdout_policy_guards():
    with pytest.raises(ValueError,match="Privacy"):
        holdout_genus_transfer(*fixtures(),ledger_label="x",
                              scope_label="All",min_test_species=2)
    with pytest.raises(ValueError,match="Unknown"):
        run(fixtures(),policy="invented")


def test_repeated_dataframe_index_is_not_repeated_island_identity():
    f=fixtures()
    x=f[0]
    f[0]=pd.concat([x.iloc[:40],x.iloc[40:].reset_index(drop=True)])
    assert len(run(f)) == 4
