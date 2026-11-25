"""Synthetic-only validations of mainland source opportunity and 17D KDE."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import chi2
from scipy.special import gammaln

from fd.core import ProfileSpace
from fd.source_genus_opportunity import (
    ranking_metrics, source_candidate_scores, _paired_uncertainty,
    exclude_heldout_species,
)
from fd.hypervolume17 import (
    DIMENSIONS, kde_hdr_17, pcoa_embedding_17,
)


def test_rank_metrics_with_tied_opportunity_scores():
    y = np.array([1, 0, 1, 0])
    scores = np.array([3.0, 3.0, 2.0, 1.0])
    result = ranking_metrics(y, scores)
    assert result["average_precision"] == pytest.approx(
        0.5 * 0.5 + 0.5 * (2 / 3))
    assert result["roc_auc"] == pytest.approx(0.625)
    same = ranking_metrics(y, np.ones(4))
    assert same["average_precision"] == pytest.approx(0.5)
    assert same["roc_auc"] == pytest.approx(0.5)
    with pytest.raises(ValueError):
        ranking_metrics(np.ones(4), scores)


def test_mainland_genus_source_scores_only_use_native_reference():
    pool = {
        "A": {"GenA": 6, "GenB": 2, "GenC": 1},
        "B": {"GenA": 3, "GenD": 1, "GenE": 2},
        "C": {"GenB": 9, "GenF": 5},
    }
    catalog_sizes = {"A": 9, "B": 6, "C": 14}
    global_occ = {"GenA": 2, "GenB": 2, "GenC": 1,
                  "GenD": 1, "GenE": 1, "GenF": 1}
    ids, score = source_candidate_scores(
        [("A", 10), ("B", 80), ("C", 500)], pool, catalog_sizes,
        global_occ, 3)
    assert ids == sorted(global_occ)
    assert len(score) == 4
    assert score["global"][ids.index("GenA")] == pytest.approx(2 / 3)
    assert score["local_count"][ids.index("GenA")] == 2
    assert score["local_density"][ids.index("GenA")] > 0
    with pytest.raises(ValueError):
        source_candidate_scores([("A", 10)], pool, catalog_sizes, global_occ, 3)


def test_paired_opportunity_uncertainty_is_deterministic():
    d = np.array([0.1, 0.4, -0.2, 0.3, -0.05, 0.2])
    a = _paired_uncertainty(d, seed=99, draws=999)
    b = _paired_uncertainty(d, seed=99, draws=999)
    assert a == b
    assert a[0] < np.mean(d) < a[1]
    assert 0 < a[2] <= 1


def test_pcoa17_rejects_ill_defined_reference():
    profiles = {
        i: (frozenset({f"c{i}"}), frozenset({"bell"}))
        for i in range(25)
    }
    space = ProfileSpace(profiles)
    coordinates, audit = pcoa_embedding_17(space)
    assert coordinates.shape == (25, DIMENSIONS)
    assert 0.6 < audit["positive_inertia_17axes"] < 1.0
    assert audit["n_positive_axes"] >= 17
    tiny = ProfileSpace({0: profiles[0], 1: profiles[1]})
    with pytest.raises(ValueError):
        pcoa_embedding_17(tiny)


def test_projected_kde17_integrator_against_single_gaussian():
    # 17 identical centers produce exactly one isotropic 17D Gaussian.
    # The analytical volume of its 90% HDR is a sphere of radius
    # h * sqrt(chi2_17(.90)).
    points = np.zeros((17, 17))
    h = 0.55
    log_volume, rel_mc = kde_hdr_17(
        points, bandwidth=h, n_draws=22000, seed=2026)
    radius = h * math.sqrt(chi2.ppf(.9, df=17))
    log_volume_true = (
        DIMENSIONS/2 * math.log(math.pi)
        - gammaln(DIMENSIONS/2+1)
        + DIMENSIONS * math.log(radius))
    assert log_volume == pytest.approx(log_volume_true, abs=0.30)
    assert rel_mc < 0.07
    with pytest.raises(ValueError):
        kde_hdr_17(points, bandwidth=h, n_draws=100, seed=2026)


def test_strict_heldout_taxon_names_are_removed_from_all_mainland_sources():
    mainland = pd.DataFrame({
        "entity_ID": ["E1", "E2", "E2", "E3", "E3"],
        "work_species": ["X test", "X test", "Y unrelated",
                         "Z test", "W unrelated"],
    })
    removed = exclude_heldout_species(mainland, {"X test", "Z test"})
    assert len(removed) == 2
    assert set(removed.work_species) == {"Y unrelated", "W unrelated"}
    assert not {"X test", "Z test"} & set(removed.work_species)
    with pytest.raises(ValueError):
        exclude_heldout_species(mainland, set(mainland.work_species))
