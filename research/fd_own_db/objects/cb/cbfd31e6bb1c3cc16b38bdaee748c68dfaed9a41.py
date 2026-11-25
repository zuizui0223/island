import numpy as np
import pandas as pd
import pytest
from fd.globi_channel_composition_cv import (
    CHANNELS, TRAITS, direct_complete_traits, prepare_genus_frame,
    simplex_project, heldout_channel_scores, _folds,
)


def fixture(n=65):
    raw, globi = [], []
    for i in range(n):
        genus = f"Genus{i:03d}"
        state = i % 2
        for j in range(3):
            species = f"{genus} specimen{j}"
            for t, v in [
                ("flower_primary_color", "red_pink" if state else "green_brown_inconspicuous"),
                ("floral_form", "open_radial" if state else "tubular"),
                ("floral_symmetry", "actinomorphic" if state else "zygomorphic"),
            ]:
                raw.append(dict(accepted_species=species, trait_name=t,
                                normalized_value=v, evidence_scope="species_direct",
                                evidence_quality="high"))
        refs = [3+10*state, 7, 3+5*state, 2, 5]
        globi.append(dict(genus=genus, n_independent_references=22+15*state,
                          n_matched_plant_species=2,
                          n_reference_channel_units=sum(refs),
                          **dict(zip(("refs__"+c for c in CHANNELS), refs))))
    return pd.DataFrame(raw), pd.DataFrame(globi)


def test_direct_only_never_promotes_genus_consensus_or_mixed_state():
    raw, g = fixture()
    assert len(direct_complete_traits(raw)) == 195
    x = raw.iloc[0].copy()
    x["accepted_species"] = "Genus999 specimen0"
    x["evidence_scope"] = "genus_consensus"
    raw = pd.concat([raw, pd.DataFrame([x])], ignore_index=True)
    assert len(direct_complete_traits(raw)) == 195
    raw.loc[(raw.accepted_species == "Genus000 specimen0") &
            (raw.trait_name == "flower_primary_color"), "normalized_value"] = "red_pink|white"
    assert len(direct_complete_traits(raw)) == 194


def test_conflicting_direct_evidence_fails_closed():
    raw, _ = fixture()
    with pytest.raises(ValueError, match="Duplicate"):
        direct_complete_traits(pd.concat([raw, raw.iloc[[0]]], ignore_index=True))


def test_reference_units_and_genus_uniqueness():
    raw, g = fixture()
    g.loc[0, "refs__bombus"] += 1
    with pytest.raises(ValueError, match="do not sum"):
        prepare_genus_frame(direct_complete_traits(raw), g)
    raw, g = fixture()
    with pytest.raises(ValueError, match="repeated GloBI"):
        prepare_genus_frame(direct_complete_traits(raw),
                            pd.concat([g,g.iloc[[0]]], ignore_index=True))


def test_synthetic_morphology_predicts_citation_composition_outofgenus():
    raw, g = fixture()
    X, models, Y, info = prepare_genus_frame(direct_complete_traits(raw), g)
    assert info["n_source_direct_species"] == 195
    assert info["n_genera"] == 65
    assert np.allclose(Y.sum(axis=1), 1)
    out = heldout_channel_scores(X, models, Y,
                                 seeds=(20261111,), ridge_penalties=(100.,))
    assert len(out) == 24
    assert out.n_genera.nunique() == 1
    assert not {"genus", "accepted_species", "island_id", "entity_ID"} & set(out.columns)
    gain = out[(out.channel == "all_five_channels") &
               (out.model == "simple_morphology")].MSE_gain_vs_reference_only.iloc[0]
    assert gain > 0


def test_simplex_and_crossfold_safety():
    p = simplex_project(np.array([[-.5, 1.4, .1, 0, 0],
                                  [.2, .2, .2, .2, .2]]))
    assert (p>=0).all() and np.allclose(p.sum(axis=1), 1)
    pairs = list(_folds(65, seed=123))
    assert len(pairs) == 5
    assert sorted(np.concatenate([test for _,test in pairs]).tolist()) == list(range(65))
    for train, test in pairs:
        assert not set(train) & set(test)
