import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from fd.nectarivore_flora_test import (
    NectarTestResult,
    add_controls,
    bird_design,
    clade_counts,
    decision,
    group_index,
    permute_within_blocks,
    plant_block_expectation,
    run_test,
    within_fit,
)
from fd.nectarivore_sim import fit_e1_e2, simulate
from fd.range_map_nectarivore import OTHER, PARTIAL, SPECIALIST, diet_counts


def test_block_expectation_is_exact_block_mean():
    plants = pd.DataFrame(
        {
            "island": ["i1", "i1", "i2", "i2", "i2"],
            "species": ["a1", "b1", "a2", "b1", "b2"],
            "genus": ["A", "B", "A", "B", "B"],
            "grade": ["h"] * 5,
            "trait": [1, 0, 0, 0, 1],
        }
    )
    out = plant_block_expectation(plants)
    # genus A mean over distinct species = .5, genus B (b1=0, b2=1) = .5
    assert out.loc["i1", "share"] == 0.5
    assert out.loc["i1", "expected"] == pytest.approx(0.5)
    assert out.loc["i2", "share"] == pytest.approx(1 / 3)
    assert out.loc["i2", "residual"] == pytest.approx(1 / 3 - 0.5)


def test_block_expectation_rejects_inconsistent_species():
    bad = pd.DataFrame(
        {"island": ["i1", "i2"], "species": ["a", "a"], "genus": ["A", "A"], "grade": ["h", "h"], "trait": [1, 0]}
    )
    with pytest.raises(ValueError, match="conflicting"):
        plant_block_expectation(bad)


def _linear_islands(rng, n_groups, per=6, b_spec=0.3, b_part=-0.1):
    rows = []
    for g in range(n_groups):
        fe = rng.normal()
        for k in range(per):
            s, p, o = rng.poisson(4), rng.poisson(4), rng.poisson(20)
            rows.append(
                dict(group=g, specialist=s, partial=p, other=o, area_km2=np.exp(rng.normal(2, 1)),
                     dist_km=rng.uniform(1, 900), n=int(rng.integers(5, 40)), lat=rng.uniform(-20, 20), fe=fe)
            )
    d = add_controls(pd.DataFrame(rows))
    counts = np.c_[d["other"], d["partial"], d["specialist"]].astype(float)
    y = d["fe"] + b_spec * np.log1p(d["specialist"]) + b_part * np.log1p(d["partial"]) + 0.01 * rng.normal(size=len(d))
    return d, counts, y.to_numpy()


def test_within_fit_recovers_coefficients_despite_group_effects():
    rng = np.random.default_rng(0)
    d, counts, y = _linear_islands(rng, 40)
    beta = within_fit(bird_design(counts, d), y, group_index(d["group"]))
    assert beta[0] == pytest.approx(0.3, abs=0.01)
    assert beta[1] == pytest.approx(-0.1, abs=0.01)


def test_permutation_preserves_block_composition():
    rng = np.random.default_rng(1)
    labels = np.array([0, 1, 2, 2, 0, 1, 1])
    blocks = np.array([0, 0, 0, 1, 1, 1, 2])
    out = permute_within_blocks(labels, blocks, rng)
    for b in np.unique(blocks):
        assert sorted(out[blocks == b]) == sorted(labels[blocks == b])


def _random_species_world(rng, n_groups, per=5, n_species=60):
    inc = sparse.csr_matrix(rng.random((n_groups * per, n_species)) < 0.15)
    diet = rng.choice([OTHER, PARTIAL, SPECIALIST], n_species, p=[0.7, 0.2, 0.1])
    genus = np.arange(n_species) // 3
    family = np.arange(n_species) // 12
    isl = pd.DataFrame(
        dict(group=np.repeat(np.arange(n_groups), per), area_km2=np.exp(rng.normal(2, 1, n_groups * per)),
             dist_km=rng.uniform(1, 900, n_groups * per), n=rng.integers(5, 40, n_groups * per),
             lat=rng.uniform(-20, 20, n_groups * per))
    )
    counts = diet_counts(inc, diet)
    isl["share"] = np.clip(0.2 + 0.05 * np.log1p(counts[:, SPECIALIST]) + rng.normal(0, 0.02, len(isl)), 0, 1)
    isl["residual"] = rng.normal(0, 0.02, len(isl))
    return isl, inc, diet, {"genus": genus, "family": family}


def test_run_test_gate_and_outputs():
    rng = np.random.default_rng(2)
    isl, inc, diet, blocks = _random_species_world(rng, 8)
    r = run_test(isl, inc, diet, blocks, n_boot=49, n_perm=49)
    assert r.status == "insufficient_support"
    assert decision(r) == "insufficient_support"
    isl, inc, diet, blocks = _random_species_world(rng, 45)
    r = run_test(isl, inc, diet, blocks, n_boot=99, n_perm=99)
    assert r.status == "ok" and r.n_islands == 225
    assert r.estimates["share"]["E"] > 0
    assert r.bootstrap["share"]["q025"] > 0
    assert set(r.permutation["share"]) == {"genus", "family"}
    assert 0 < r.permutation["share"]["family"]["p_one_sided"] <= 1
    assert decision(r) == "lineage_cooccurrence"


def test_run_test_rejects_small_floras():
    rng = np.random.default_rng(3)
    isl, inc, diet, blocks = _random_species_world(rng, 3)
    isl.loc[0, "n"] = 4
    with pytest.raises(ValueError, match="fewer than"):
        run_test(isl, inc, diet, blocks, n_boot=9, n_perm=9)


def _result(q1, q2, p):
    return NectarTestResult(
        "ok", 300, 40, 35, {},
        {"share": {"q025": q1[0], "q975": q1[1]}, "residual": {"q025": q2[0], "q975": q2[1]}},
        {"residual": {"genus": {"p_one_sided": p}}},
    )


@pytest.mark.parametrize(
    "q1,q2,p,expected",
    [
        ((0.1, 0.3), (0.05, 0.2), 0.01, "beyond_plant_and_bird_lineage"),
        ((0.1, 0.3), (0.05, 0.2), 0.20, "beyond_plant_lineage_only"),
        ((0.1, 0.3), (-0.05, 0.1), 0.01, "lineage_cooccurrence"),
        ((-0.1, 0.3), (-0.05, 0.1), 0.50, "no_generalisation"),
        ((-0.3, -0.1), (-0.05, 0.1), 0.50, "negative_association"),
    ],
)
def test_decision_table(q1, q2, p, expected):
    assert decision(_result(q1, q2, p)) == expected


def test_clade_counts_split_specialists():
    inc = sparse.csr_matrix(np.array([[1, 1, 1, 1], [0, 1, 0, 1]], bool))
    diet = np.array([SPECIALIST, SPECIALIST, PARTIAL, OTHER])
    clade = np.array([0, 1, -1, -1])
    counts, extra = clade_counts(inc, diet, clade, focal=0)
    assert counts.tolist() == [[1, 1, 1], [1, 0, 0]]
    assert extra.tolist() == [1, 1]


@pytest.mark.parametrize("scenario,e1_pos,e2_pos", [("cooccurrence", True, False), ("matching", True, True)])
def test_simulation_separates_cooccurrence_from_matching(scenario, e1_pos, e2_pos):
    d = simulate(scenario, 40, seed=3)
    f = fit_e1_e2(d, n_boot=99)
    assert (f["share"][1] > 0) == e1_pos
    assert (f["residual"][1] > 0) == e2_pos
