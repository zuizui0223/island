import numpy as np
import pandas as pd
import pytest

from fd.floral_syndrome_coupling import (
    _entropy, _excess_one_sample, audit_island_flora,
    marginal_rao_equal_traits, template_overlap, TEMPLATES,
)


def test_equal_marginal_diversity_can_hide_different_joint_pairing():
    aligned = np.array([[0, 0, 0]] * 4 + [[1, 1, 1]] * 4)
    independent = np.array([
        [c, f, s] for c in (0, 1) for f in (0, 1) for s in (0, 1)
    ])
    assert marginal_rao_equal_traits(aligned) == pytest.approx(.5)
    assert marginal_rao_equal_traits(independent) == pytest.approx(.5)
    def tc(x):
        c, f, s = x.T
        return (_entropy(c) + _entropy(f) + _entropy(s)
                - _entropy((c * 2 + f) * 2 + s))
    assert tc(aligned) > 1.3
    assert tc(independent) == pytest.approx(0, abs=1e-12)


def test_colour_architecture_excess_and_genus_null_identifiability():
    arr = np.array([[0, 0, 0]] * 15 + [[1, 1, 1]] * 15)
    unique_genus = np.array([f"G{i}" for i in range(30)])
    r = _excess_one_sample(arr, unique_genus,
                           np.random.default_rng(42), permutations=99)
    assert r["excess_ca"] > 0
    assert r["excess_fs"] > 0
    assert r["fraction_genus_color_informative"] == 0
    assert r["excess_ca_genus"] == pytest.approx(0, abs=1e-12)


def _test_islands(n=24, per=30):
    records = []
    for isl in range(n):
        for i in range(per):
            g = "GenusA" if i < per // 2 else "GenusB"
            state = 0 if i < per // 2 else 1
            records.append(dict(
                island_id=f"I{isl}", accepted_species=f"{g} species{i}",
                genus=g, analysis_regime="northern_midlatitude",
                spatial_block=f"BLOCK{isl//3}",
                origin_status="unresolved",
                flower_primary_color="white" if state == 0 else "blue",
                floral_form="open_radial" if state == 0 else "tubular",
                floral_symmetry="actinomorphic" if state == 0 else "zygomorphic",
            ))
    return pd.DataFrame(records)


def test_aggregate_never_exports_species_or_island_identifiers():
    x = audit_island_flora(_test_islands(), rarefactions=1,
                           permutations=29, block_draws=199)
    assert len(x) == 1
    assert x.iloc[0].n_islands == 24
    assert x.iloc[0].n_spatial_blocks == 8
    assert x.iloc[0].median_color_architecture_excess_nats > 0
    assert not {"accepted_species", "genus", "island_id", "spatial_block"} & set(x.columns)


def test_mixed_colour_states_and_conflicting_species_evidence_fail_closed():
    f = _test_islands()
    f.loc[0, "flower_primary_color"] = "red|white"
    with pytest.raises(ValueError, match="Mixed-state"):
        audit_island_flora(f, rarefactions=1, permutations=29, block_draws=199)
    f = _test_islands()
    f.loc[30, "floral_form"] = "new_form"
    with pytest.raises(ValueError, match="Conflicting"):
        audit_island_flora(f, rarefactions=1, permutations=29, block_draws=199)


def test_status_and_island_support_gates():
    f = _test_islands(n=22, per=30)
    f["spatial_block"] = "SINGLE_BLOCK"
    assert audit_island_flora(f, rarefactions=1,
                             permutations=29, block_draws=199).empty
    with pytest.raises(ValueError, match="Privacy"):
        audit_island_flora(f, min_islands=2)


def test_candidate_syndrome_overlap_is_not_realized_pollinator():
    rows = []
    for i in range(55):
        for tpl in TEMPLATES:
            rows.append(dict(
                accepted_species=f"Spec{i}", syndrome=tpl,
                soft_membership=.9, n_informative_traits=3))
    d = pd.DataFrame(rows)
    s = template_overlap(d, high_cutoff=.75)
    assert s["n_comparable_species"] == 55
    assert s["n_high_two_or_more"] == 55
    assert s["n_high_all_three"] == 55
    assert s["pollinator_identity_established"] is False
    d.loc[d.syndrome.eq("large_bee_like"), "n_informative_traits"] = 2
    with pytest.raises(ValueError, match="Too few"):
        template_overlap(d)
