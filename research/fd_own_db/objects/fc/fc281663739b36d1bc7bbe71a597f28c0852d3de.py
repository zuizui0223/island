"""Synthetic-only N2 five-geography-fold prediction and cluster uncertainty."""
import numpy as np
import pandas as pd
import pytest

from fd.n2_crossfit_macro import (
    METRICS, block_bootstrap_mean, summarize, taxon_balanced_source_summary,
    summarize_trait_axes,
)


def test_block_bootstrap_clusters_not_individual_islands_and_is_deterministic():
    # 10 dependent islands in one block + one island in each other block.
    groups = ["A"] * 10 + ["B", "C", "D", "E", "F"]
    values = [2.] * 10 + [-1.] * 5
    data = pd.DataFrame({"spatial_block": groups, "delta": values})
    mu, low, high = block_bootstrap_mean(
        data, response="delta", draws=999, seed=42)
    assert mu == pytest.approx(1.)
    assert low < mu < high
    assert (mu, low, high) == block_bootstrap_mean(
        data, response="delta", draws=999, seed=42)
    with pytest.raises(ValueError):
        block_bootstrap_mean(
            data.assign(spatial_block=None), response="delta", draws=999)
    with pytest.raises(ValueError):
        block_bootstrap_mean(
            data.loc[data.spatial_block.isin(["A", "B"])],
            response="delta", draws=999)


def test_grouped_summary_keeps_the_five_fold_source_contrast_sign():
    rows = []
    for scope in ("all", "direct"):
        for block in range(6):
            rows.append({
                "scope": scope, "analysis_regime": "tropical",
                "spatial_block": f"B{block}", "island_id": f"secret{block}",
                "global_coverage": .8, "local_coverage": .3,
                "prior_logloss": 6., "global_logloss": 5.,
                "local_logloss": 5.1,
                "genus_gain": 1., "source_gain": -.1,
                "source_distance_gain": .02,
                "genus_brier_gain": .3, "source_brier_gain": .004,
                "genus_kernel_gain": .5, "source_kernel_gain": .01,
            })
    output = summarize(rows, draws=299, min_cell=5)
    assert len(output) == 4  # one global and one tropical row per scope
    assert all(r["mean_source_gain"] == pytest.approx(-.1) for r in output)
    assert all(r["source_gain_block_q975"] < 0 for r in output)
    assert all(r["genus_gain_block_q025"] > 0 for r in output)
    assert all(r["mean_source_distance_gain"] == pytest.approx(.02)
               for r in output)
    assert all(r["mean_source_brier_gain"] == pytest.approx(.004) for r in output)
    assert all(r["mean_source_kernel_gain"] == pytest.approx(.01) for r in output)
    assert all(r["source_kernel_gain_block_q025"] > 0 for r in output)
    assert not any("island_id" in row or "spatial_block" in row for row in output)


def test_scant_spatial_replication_suppresses_predictive_claim():
    rows = [{
        "scope": "all", "analysis_regime": "southern_extratropical",
        "spatial_block": "ONE", "global_coverage": .1,
        "local_coverage": .05, "prior_logloss": 5.0,
        "global_logloss": 4.0, "local_logloss": 4.5,
        "genus_gain": 1.0, "source_gain": -.5,
        "source_distance_gain": .01,
    } for _ in range(8)]
    out = summarize(rows, draws=99)
    assert len(out) == 2
    assert all(r["status"] == "insufficient_independent_geographic_support"
               for r in out)
    assert all("mean_genus_gain" not in r for r in out)
    assert METRICS == (
        "genus_gain", "source_gain", "source_distance_gain",
        "genus_brier_gain", "source_brier_gain",
        "genus_kernel_gain", "source_kernel_gain",
    )


def test_taxon_balanced_source_gain_does_not_multiply_repeated_island_names():
    # One broadly recorded taxon (five folds) versus 19 unique taxa.
    # A record-weighted mean is 5/24, but correct one-name mean = 1/20.
    rows = []
    for g in range(20):
        folds = range(5) if g == 0 else [0]
        for fold in folds:
            rows.append({
                "accepted_species": f"Gen{g} taxon",
                "genus": f"Gen{g}", "scope": "all", "fold": fold,
                "source_gain_kernel": float(g == 0),
                "source_gain_brier": .001 * (g + 1),
                "source_gain_logloss": 0.,
                "source_gain_trait_distance_loss": 0.,
                "has_local_genus_donor": float(g == 0),
            })
    out = taxon_balanced_source_summary(
        [pd.DataFrame(rows)], draws=399)
    assert len(out) == 1
    assert out[0]["n_unique_names"] == 20
    assert out[0]["n_unique_genera"] == 20
    assert out[0]["n_name_fold_units"] == 24
    assert out[0]["mean_source_gain_kernel"] == pytest.approx(1/20)
    assert out[0]["mean_local_genus_donor_support"] == pytest.approx(1/20)
    assert out[0]["mean_source_gain_brier"] == pytest.approx(.0105)
    assert all("accepted_species" not in x and "genus" not in x
               for x in out)


def test_taxon_balance_rejects_duplicated_name_fold_rows():
    rows = pd.DataFrame([{
        "accepted_species": "Gen one", "genus": "Gen", "scope": "all",
        "fold": 0, "source_gain_kernel": .1, "source_gain_brier": .1,
        "source_gain_logloss": .1, "source_gain_trait_distance_loss": .1,
        "has_local_genus_donor": 1.,
    }] * 2)
    with pytest.raises(ValueError, match="deduplicated"):
        taxon_balanced_source_summary([rows], draws=299)



def test_trait_axis_macro_summary_checks_additive_identity_and_block_uncertainty():
    rows = []
    for scope in ("all", "direct"):
        for j in range(6):
            rows.append(dict(
                scope=scope, spatial_block=f"block{j}",
                analysis_regime="tropical",
                genus_colour=.3, genus_form=.1,
                genus_additive=.2, genus_joint_product=.12,
                source_colour=.01, source_form=-.002,
                source_additive=.004, source_joint_product=.006,
            ))
    report = summarize_trait_axes(rows, draws=199)
    assert len(report) == 4
    assert all(r["mean_source_additive"] == pytest.approx(.004) for r in report)
    assert all(r["mean_source_joint_product"] == pytest.approx(.006) for r in report)
    assert all(r["source_joint_product_block_q025"] > 0 for r in report)
    assert all("spatial_block" not in r for r in report)
    rows[0]["source_additive"] = .9
    with pytest.raises(AssertionError, match="Additive"):
        summarize_trait_axes(rows, draws=199)
