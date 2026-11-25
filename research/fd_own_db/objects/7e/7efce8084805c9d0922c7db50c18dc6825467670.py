"""Synthetic, independent finite-sample checks for native island floral FD."""

from itertools import combinations, permutations

import numpy as np
import pytest

from fd.regional_native_flora_fd import exact_rarefied_flora_fd


def _direct_q(sample):
    n = len(sample)
    marginal = []
    for j in range(3):
        _, counts = np.unique(sample[:, j], return_counts=True)
        marginal.append(1.0 - np.square(counts/n).sum())
    _, joint_counts = np.unique(sample, axis=0, return_counts=True)
    return float(np.mean(marginal)), float(1.0 - np.square(joint_counts/n).sum())


def test_finite_raoq_matches_exhaustive_without_replacement():
    a = np.array([
        ["A", "open_radial", "actinomorphic"],
        ["A", "tubular", "zygomorphic"],
        ["B", "open_radial", "actinomorphic"],
        ["B", "tubular", "zygomorphic"],
        ["B", "tubular", "actinomorphic"],
    ], dtype=str)
    observed = np.array([
        _direct_q(a[list(idx)]) for idx in combinations(range(len(a)), 3)
    ])
    exact = exact_rarefied_flora_fd(a, n=3)
    assert np.allclose(
        observed.mean(axis=0),
        [exact["marginal_rao"], exact["joint_rao"]],
        atol=1e-14,
    )


def test_pairing_shuffled_null_is_exact_colour_permutation_average():
    a = np.array([
        ["A", "open_radial", "actinomorphic"],
        ["A", "open_radial", "actinomorphic"],
        ["B", "open_radial", "zygomorphic"],
        ["B", "tubular", "zygomorphic"],
        ["C", "tubular", "actinomorphic"],
    ], dtype=str)
    original = exact_rarefied_flora_fd(a, n=3)
    shuffled_q = []
    for idx in permutations(range(len(a))):
        b = a.copy()
        b[:, 0] = a[np.asarray(idx), 0]
        shuffled_q.append(exact_rarefied_flora_fd(b, n=3)["joint_rao"])
    assert np.isclose(
        original["joint_rao"] - float(np.mean(shuffled_q)),
        original["colour_architecture_pairing_departure"],
        atol=1e-14,
    )


def test_marginal_diversity_same_but_joint_and_pairing_differ():
    x = [["A", "open_radial", "actinomorphic"]] * 4 + [
        ["B", "tubular", "zygomorphic"]] * 4
    y = (
        [["A", "open_radial", "actinomorphic"]] * 2 +
        [["A", "tubular", "zygomorphic"]] * 2 +
        [["B", "open_radial", "actinomorphic"]] * 2 +
        [["B", "tubular", "zygomorphic"]] * 2
    )
    one = exact_rarefied_flora_fd(np.array(x), n=4)
    two = exact_rarefied_flora_fd(np.array(y), n=4)
    assert np.isclose(one["marginal_rao"], two["marginal_rao"])
    assert one["joint_rao"] < two["joint_rao"]
    assert one["colour_architecture_pairing_departure"] < two[
        "colour_architecture_pairing_departure"]


def test_fail_closed_incomplete_ambiguous_and_undersupported():
    a = np.array([["A", "open_radial", "actinomorphic"]] * 5)
    with pytest.raises(ValueError):
        exact_rarefied_flora_fd(a, n=30)
    bad = a.copy()
    bad[0, 1] = "open|tubular"
    with pytest.raises(ValueError):
        exact_rarefied_flora_fd(bad, n=3)
    with pytest.raises(ValueError):
        exact_rarefied_flora_fd(a, n=1)
