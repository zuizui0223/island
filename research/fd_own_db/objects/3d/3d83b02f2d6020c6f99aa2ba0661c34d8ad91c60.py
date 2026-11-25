import numpy as np
import pandas as pd
import pytest

from fd.globi_family_holdout_cv import (
    family_disjoint_scores, family_folds, within_family_morphology_null,
    coarse_morphology_ablation,
)

MODELS = {
    "evidence_only": ["log_refs", "log_plant", "log_n_morph"],
    "simple_morphology": ["log_refs", "log_plant", "log_n_morph",
                          "prop_actinomorphic", "prop_open_radial", "prop_inconspicuous"],
    "all_marginals": ["log_refs", "log_plant", "log_n_morph",
                      "prop_actinomorphic", "prop_open_radial", "prop_inconspicuous"],
    "marginals_plus_joints": ["log_refs", "log_plant", "log_n_morph",
                              "prop_actinomorphic", "prop_open_radial", "prop_inconspicuous"],
}


def synthetic(n=100):
    rng = np.random.default_rng(42)
    families = pd.Series([f"Fam{i//5}" for i in range(n)],
                         index=[f"Genus{i}" for i in range(n)])
    t = rng.uniform(0, 1, n)
    x = pd.DataFrame({
        "log_refs": rng.normal(3, 1, n),
        "log_plant": rng.normal(3, .6, n),
        "log_n_morph": rng.normal(2, .4, n),
        "prop_actinomorphic": t,
        "prop_open_radial": t*.8,
        "prop_inconspicuous": 1-t,
    }, index=families.index)
    y = np.column_stack([.15+.1*t, .1+.02*t, .25+.1*(1-t), .2+.1*t, .3-.02*t])
    return x, y/y.sum(axis=1)[:, None], families


def test_family_splits_are_disjoint_and_exhaustive():
    x, y, families = synthetic()
    folds = family_folds(families.to_numpy(), seed=20261111)
    assert len(folds) == 5
    assert sorted(np.concatenate([test for _,test in folds]).tolist()) == list(range(len(families)))
    assert all(not (set(families.iloc[train]) & set(families.iloc[test]))
               for train, test in folds)


def test_blocked_morphology_gain_and_private_aggregate():
    x, y, families = synthetic()
    summary = family_disjoint_scores(
        x, MODELS, y, families, seeds=(20261111,), alphas=(100.,),
        bootstrap_draws=199)
    assert len(summary) == 4
    row = summary.set_index("model").loc["simple_morphology"]
    assert row.n_genera == 100 and row.n_family_blocks == 20
    assert row.MSE_gain_over_reference_only > 0
    assert row.family_equal_bootstrap_lower <= row.family_equal_gain
    assert row.family_equal_bootstrap_upper >= row.family_equal_gain
    assert (row.conditional_family_bootstrap_lower <= row.MSE_gain_over_reference_only
            <= row.conditional_family_bootstrap_upper)
    assert not {"genus", "family", "island_id", "accepted_species"} & set(summary)


def test_unknown_family_and_missing_trait_fail_closed():
    x, y, families = synthetic()
    families.iloc[0] = "unknown"
    with pytest.raises(ValueError, match="Unresolved"):
        family_disjoint_scores(x, MODELS, y, families, seeds=(1,), alphas=(100.,),
                               bootstrap_draws=199)
    x, y, families = synthetic()
    x.loc[x.index[0], "prop_actinomorphic"] = np.nan
    with pytest.raises(ValueError, match="Unknown direct morphology"):
        family_disjoint_scores(x, MODELS, y, families, seeds=(1,), alphas=(100.,),
                               bootstrap_draws=199)


def test_within_family_whole_triplet_placebo_not_causal():
    x, y, families = synthetic()
    result = within_family_morphology_null(
        x, y, families, base_features=MODELS["evidence_only"],
        morph_features=["prop_actinomorphic", "prop_open_radial",
                        "prop_inconspicuous"], permutations=99)
    assert result["n_genera"] == 100
    assert result["n_families"] == 20
    assert result["n_genera_in_shufflable_families"] == 100
    assert not result["actual_island_visitors_or_pollination_effectiveness_tested"]



def test_morphology_ablation_disaggregates_colour_and_structure():
    x, y, families = synthetic()
    z = coarse_morphology_ablation(
        x, y, families, base_features=MODELS["evidence_only"],
        symmetry="prop_actinomorphic", form="prop_open_radial",
        colour="prop_inconspicuous", seeds=(20261111,), alpha=100.)
    assert len(z) == 8
    scores = z.set_index("model").MSE_gain_over_evidence_only
    assert scores["all_three"] > 0
    assert scores["symmetry_form"] > 0
    assert scores["evidence_only"] == pytest.approx(0)
    assert not {"family", "genus", "accepted_species"} & set(z.columns)


def test_competing_family_labels_must_be_excluded_before_blocked_cv():
    x, y, families = synthetic()
    families.iloc[0] = "AMBIG:Apocynaceae&Asclepiadaceae"
    with pytest.raises(ValueError, match="conflicting"):
        family_disjoint_scores(x, MODELS, y, families,
                               seeds=(20261111,), alphas=(100.),
                               bootstrap_draws=199)
