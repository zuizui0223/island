import numpy as np
import pandas as pd
import pytest

from fd.globi_family_lofo import leave_one_family_out_scores


def fixture(n=80):
    rng = np.random.default_rng(20)
    families = pd.Series([f"Fam{i//4}" for i in range(n)],
                         index=[f"Genus{i}" for i in range(n)])
    form = rng.uniform(0, 1, n)
    features = pd.DataFrame({
        "reference_effort": rng.uniform(1, 4, n),
        "open_form": form, "colour": rng.uniform(0, 1, n),
    }, index=families.index)
    composition = np.column_stack([
        .15+.25*form, .25-.1*form, .2+.1*form,
        .2-.1*form, .2-.15*form,
    ])
    composition /= composition.sum(axis=1)[:, None]
    models = {
        "evidence_only": ["reference_effort"],
        "open_form_only": ["reference_effort", "open_form"],
        "with_colour": ["reference_effort", "open_form", "colour"],
    }
    return features, families, composition, models


def test_every_family_excluded_and_no_source_id_exported():
    f, families, y, models = fixture()
    x = leave_one_family_out_scores(f, families, y, models, bootstrap_draws=299)
    assert len(x) == 3
    assert set(x.n_families) == {20}
    a = x.set_index("model")
    assert a.loc["open_form_only", "genus_equal_MSE_gain"] > 0
    assert a.loc["open_form_only", "genus_equal_JSD"] < a.loc["evidence_only", "genus_equal_JSD"]
    assert not {"genus", "family", "island_id", "source_provider"} & set(x)


def test_unresolved_family_or_misaligned_target_rejected():
    features, fam, y, models = fixture()
    fam.iloc[0] = "AMBIG:FamilyA&FamilyB"
    with pytest.raises(ValueError, match="ambiguous"):
        leave_one_family_out_scores(features, fam, y, models, bootstrap_draws=299)
    features, fam, y, models = fixture()
    with pytest.raises(ValueError, match="misaligned"):
        leave_one_family_out_scores(features, fam.iloc[::-1], y, models, bootstrap_draws=299)


def test_nonbinary_source_proportions_and_unrecorded_morphology_rejected():
    features, fam, y, models = fixture()
    features.loc[features.index[0], "open_form"] = np.nan
    with pytest.raises(ValueError, match="Missing Direct"):
        leave_one_family_out_scores(features, fam, y, models, bootstrap_draws=299)
    features, fam, y, models = fixture()
    y[0, 0] = .95
    with pytest.raises(ValueError, match="Invalid GloBI"):
        leave_one_family_out_scores(features, fam, y, models, bootstrap_draws=299)
