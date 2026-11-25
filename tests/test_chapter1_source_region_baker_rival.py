"""Purely synthetic tests for the post-hoc source-region/Baker rival audit."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "audit_chapter1_source_region_baker_rival.py"
)
spec = importlib.util.spec_from_file_location("source_region_rival", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.fixture
def synthetic_data():
    rng = np.random.default_rng(20261117)
    scores = []
    covariates = []
    assignments = []
    for region, higher in (
        ("northern_midlatitude", True),
        ("tropical", False),
    ):
        for k in range(84):
            island = f"{region}_{k}"
            x = rng.uniform(0.3, 5.6)
            covariates.append({
                "island_id": island,
                "analysis_regime": region,
                "spatial_block": f"{region}_block_{k % 14}",
                "log_distance_to_continent_km": x,
                "log_island_area_km2": rng.normal(),
                **{f"climate_pc{j}": rng.normal() for j in range(1, 5)},
            })
            for mode in module.MODES:
                assignments.append({
                    "island_id": island,
                    "source_mode": mode,
                    "entity_ID": f"donor_{k % 8}",
                })
                scv = 0.1 * x + rng.normal(scale=0.12)
                for outcome in module.TRAITS:
                    genus = rng.normal(scale=0.05)
                    within = (
                        (0.14 * x if higher else -0.02 * x)
                        if outcome == "generalized_form"
                        else 0.01 * x
                    ) + rng.normal(scale=0.06)
                    if outcome == "self_compatibility":
                        within = scv
                    source = 0.2 + rng.normal(scale=0.02)
                    scores.append({
                        "island_id": island,
                        "context": region,
                        "evidence_scope": "direct_only",
                        "source_mode": mode,
                        "outcome": outcome,
                        "identity_error": 0.0,
                        "raw_h1_mean": source + genus + within,
                        "source_species_expectation": source,
                        "genus_species_structure_enrichment": genus,
                        "total_species_sorting_enrichment": genus + within,
                        "within_genus_species_sorting_enrichment": within,
                    })
    return (
        pd.DataFrame(scores), pd.DataFrame(covariates),
        pd.DataFrame(assignments),
    )


def test_direct_contrast_keeps_source_and_sc_adjustment(synthetic_data):
    result = module.run(*synthetic_data)
    focal = result.loc[
        result.source_mode.eq("geo_k10")
        & result.outcome.eq("generalized_form")
        & result.sc_adjusted
        & result.window.eq("common_05_95")
        & result.component.eq("within_genus_species_sorting_enrichment")
    ]
    assert len(focal) == 1
    row = focal.iloc[0]
    assert row.north_status == row.tropical_status == "fit"
    assert row.contrast_north_minus_tropical > 0.07
    assert row.contrast_p < 0.05
    assert result.supports_novel_rule.eq(False).all()


def test_identity_and_geo_uniqueness_guards(synthetic_data):
    scores, covariates, assignments = synthetic_data
    corrupted = scores.copy()
    corrupted.loc[corrupted.index[0], "identity_error"] = 1e-4
    with pytest.raises(ValueError, match="identity"):
        module.prepare(corrupted, covariates, assignments)
    with pytest.raises(ValueError, match="duplicate corrected geography"):
        module.prepare(
            scores, pd.concat([covariates, covariates.head(1)]),
            assignments,
        )


def test_region_support_gate(synthetic_data):
    scores, covariates, assignments = synthetic_data
    covariates = covariates.copy()
    covariates["spatial_block"] = "one_block"
    result = module.run(scores, covariates, assignments)
    assert result.contrast_p.isna().all()
    assert result.north_status.eq("insufficient_support").all()
