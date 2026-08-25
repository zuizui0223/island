import math

import numpy as np
import pandas as pd

from island_v2.chapter1_trait_resolution_bias import (
    build_trait_resolution_coverage,
    run_coverage_adjusted_omnibus,
)


def test_coverage_keeps_zero_trait_resolution_as_evidence_missingness():
    flora = pd.DataFrame(
        {
            "island_id": ["i1", "i1", "i2", "i2"],
            "accepted_species": ["a", "b", "c", "d"],
            "origin_status": ["native"] * 4,
            "floristic_status": ["native_nonendemic"] * 4,
        }
    )
    audit = pd.DataFrame(
        {
            "accepted_species": ["a"],
            "trait_name": ["self_incompatibility"],
            "resolved_for_primary": [True],
        }
    )
    coverage = build_trait_resolution_coverage(flora, audit)
    si = coverage.loc[
        coverage["stratum"].eq("all_native")
        & coverage["trait_name"].eq("self_incompatibility")
    ].set_index("island_id")
    assert si.loc["i1", "n_direct_trait_species"] == 1
    assert si.loc["i2", "n_direct_trait_species"] == 0
    assert si.loc["i2", "direct_trait_fraction"] == 0


def _synthetic():
    rng = np.random.default_rng(59)
    comp = []
    coverage = []
    cov = []
    idx = 0
    responses = [
        ("floral_form", "bell", 1.15),
        ("floral_form", "salver", -1.0),
        ("flower_primary_color", "red", 0.9),
    ]
    for context, multiplier in [("northern_midlatitude", 1.0), ("tropical", -0.55)]:
        for j in range(65):
            island_id = f"c{idx}"
            distance = -1.5 + 3.0 * j / 64
            cov.append(
                {
                    "island_id": island_id,
                    "analysis_regime": context,
                    "spatial_block": f"b{idx // 3}",
                    "log_distance_to_continent_km": distance,
                    "log_island_area_km2": 2.0 + rng.normal(0, 0.2),
                    "climate_pc1": rng.normal(),
                    "climate_pc2": rng.normal(),
                    "climate_pc3": rng.normal(),
                    "climate_pc4": rng.normal(),
                }
            )
            for stratum in ["all_native", "native_nonendemic"]:
                coverage_by_trait = {}
                for trait in ["floral_form", "flower_primary_color"]:
                    # Coverage is geographically biased but not deterministically
                    # equal to distance. Independent island-level variation keeps
                    # the ecological distance effect statistically identifiable.
                    coverage_logit = (
                        -0.2
                        + 0.45 * distance
                        + (0.2 if trait == "floral_form" else 0)
                        + rng.normal(0, 0.55)
                    )
                    total = 100
                    covered = max(1, min(99, int(total / (1 + math.exp(-coverage_logit)))))
                    observed_coverage_logit = math.log(
                        (covered + 0.5) / (total - covered + 0.5)
                    )
                    coverage_by_trait[trait] = observed_coverage_logit
                    coverage.append(
                        {
                            "island_id": island_id,
                            "stratum": stratum,
                            "trait_name": trait,
                            "n_observed_stratum_species": total,
                            "n_direct_trait_species": covered,
                            "direct_trait_fraction": covered / total,
                            "coverage_logit_smoothed": observed_coverage_logit,
                        }
                    )
                for trait, state, slope in responses:
                    eta = (
                        -0.3
                        + multiplier * slope * distance
                        + 0.15 * coverage_by_trait[trait]
                    )
                    p = 1 / (1 + math.exp(-eta))
                    trials = 60
                    comp.append(
                        {
                            "island_id": island_id,
                            "stratum": stratum,
                            "trait_name": trait,
                            "trait_state": state,
                            "successes": int(rng.binomial(trials, p)),
                            "trials": trials,
                        }
                    )
            idx += 1
    return pd.DataFrame(comp), pd.DataFrame(coverage), pd.DataFrame(cov)


def test_coverage_adjusted_omnibus_recovers_headline_structure():
    comp, coverage, cov = _synthetic()
    config = {
        "baseline_covariates": [
            "log_island_area_km2",
            "climate_pc1",
            "climate_pc2",
            "climate_pc3",
            "climate_pc4",
        ],
        "distance_gradient": "log_distance_to_continent_km",
        "context_column": "analysis_regime",
        "cluster_column": "spatial_block",
        "contexts": ["northern_midlatitude", "tropical"],
        "headline_contexts": ["northern_midlatitude", "tropical"],
        "headline_strata": ["all_native", "native_nonendemic"],
        "confirmatory_islands_per_response": 50,
    }
    within, between, summary = run_coverage_adjusted_omnibus(
        comp, coverage, cov, config
    )
    assert len(within) == 4
    assert len(between) == 2
    assert summary["headline_replication"].all()
