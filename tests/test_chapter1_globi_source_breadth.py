from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_globi_source_breadth import (
    PlantBreadthAccumulator,
    PlantTaxonomyMatcher,
    analyse_source_breadth,
)


def _config() -> dict:
    return {
        "contract": "chapter1_globi_source_breadth_v1",
        "sampling_effort_gate": {
            "primary_min_independent_references": 3,
            "sensitivity_min_independent_references": [2, 5],
        },
        "source_matching": {
            "source_modes": ["geo_k5"],
            "primary_matching": "prevalence_richness",
            "sensitivity_matching": "prevalence_only",
            "minimum_represented_genera": 5,
            "strata": ["all_native", "native_nonendemic"],
        },
        "island_enrichment": {
            "outcomes": ["entry_enrichment", "species_enrichment", "loading_increment"],
            "primary_outcome": "entry_enrichment",
            "weights": {
                "entry_enrichment": "n_represented_genera",
                "species_enrichment": "n_represented_species",
                "loading_increment": "n_represented_genera",
            },
        },
        "model": {
            "response_metrics": ["sampled_channel_count", "effective_channel_number"],
            "predictors": [
                "log_distance_to_continent_km", "log_island_area_km2",
                "climate_pc1", "climate_pc2", "climate_pc3", "climate_pc4",
            ],
            "cluster_column": "spatial_block",
            "contexts": ["northern_midlatitude"],
            "confirmatory_min_islands": 50,
            "pilot_min_islands": 30,
        },
        "primary_prediction": {
            "metric": "effective_channel_number",
            "outcome": "entry_enrichment",
            "source_matching": "prevalence_richness",
            "effort_threshold": 3,
        },
        "claim_ceiling": "test D3 boundary",
    }


def test_breadth_accumulator_collapses_reference_channel_units() -> None:
    taxonomy = pd.DataFrame(
        {
            "accepted_species": ["Alpha alba", "Alpha rubra", "Beta minor"],
            "genus": ["Alpha", "Alpha", "Beta"],
        }
    )
    matcher = PlantTaxonomyMatcher(taxonomy)
    accumulator = PlantBreadthAccumulator(
        matcher,
        ["bombus", "non_bombus_bees", "lepidoptera", "flower_visiting_birds", "diptera"],
    )
    evidence = pd.DataFrame(
        [
            {"plant_taxon": "Alpha alba", "channel_id": "bombus", "reference_key": "r1", "evidence_strength": "pollination"},
            {"plant_taxon": "Alpha alba", "channel_id": "bombus", "reference_key": "r1", "evidence_strength": "pollination"},
            {"plant_taxon": "Alpha rubra Smith", "channel_id": "lepidoptera", "reference_key": "r2", "evidence_strength": "flower_visit"},
            {"plant_taxon": "Alpha rubra", "channel_id": "lepidoptera", "reference_key": "r3", "evidence_strength": "flower_visit"},
            {"plant_taxon": "Unknown plant", "channel_id": "bombus", "reference_key": "r4", "evidence_strength": "flower_visit"},
        ]
    )
    accumulator.add(evidence)
    result = accumulator.to_frame().set_index("genus")
    alpha = result.loc["Alpha"]
    assert alpha["n_independent_references"] == 3
    assert alpha["refs__bombus"] == 1
    assert alpha["refs__lepidoptera"] == 2
    assert alpha["sampled_channel_count"] == 2
    assert 1.0 < alpha["effective_channel_number"] <= 2.0
    assert accumulator.n_unmatched_evidence_rows == 1
    assert accumulator.n_binomial_rows >= 1


def test_source_breadth_analysis_detects_positive_synthetic_enrichment(tmp_path: Path) -> None:
    genera = ["Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta"]
    breadth = pd.DataFrame(
        {
            "genus": genera,
            "n_independent_references": [6] * 6,
            "sampled_channel_count": [1, 1, 2, 2, 3, 5],
            "effective_channel_number": [1.0, 1.2, 1.8, 2.2, 3.0, 4.8],
        }
    )
    breadth_path = tmp_path / "breadth.csv"
    breadth.to_csv(breadth_path, index=False)

    gift_rows = []
    for genus in genera:
        gift_rows.append({"entity_ID": 1, "work_species": f"{genus} alba"})
    gift = pd.DataFrame(gift_rows)
    gift_path = tmp_path / "gift.csv"
    gift.to_csv(gift_path, index=False)

    status_rows = []
    assignment_rows = []
    cov_rows = []
    for i in range(60):
        island = f"island_{i:02d}"
        x = i / 10.0
        represented = genera[:5] if i < 30 else genera[1:]
        for genus in represented:
            status_rows.append(
                {
                    "island_id": island,
                    "accepted_species": f"{genus} alba",
                    "origin_status": "native",
                    "floristic_status": "native_nonendemic",
                }
            )
        assignment_rows.append({"island_id": island, "entity_ID": 1, "source_mode": "geo_k5"})
        cov_rows.append(
            {
                "island_id": island,
                "analysis_regime": "northern_midlatitude",
                "spatial_block": f"block_{i // 5}",
                "log_distance_to_continent_km": x,
                "log_island_area_km2": np.sin(x) + x / 20.0,
                "climate_pc1": np.cos(x),
                "climate_pc2": np.sin(x / 2.0),
                "climate_pc3": np.cos(x / 2.0),
                "climate_pc4": np.sin(x / 3.0),
            }
        )
    status = pd.DataFrame(status_rows)
    assignments = pd.DataFrame(assignment_rows)
    covariates = pd.DataFrame(cov_rows)
    status_path = tmp_path / "status.csv"
    assignments_path = tmp_path / "assignments.csv"
    covariates_path = tmp_path / "covariates.csv"
    status.to_csv(status_path, index=False)
    assignments.to_csv(assignments_path, index=False)
    covariates.to_csv(covariates_path, index=False)

    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(_config()), encoding="utf-8")
    output = tmp_path / "out"
    manifest = analyse_source_breadth(
        genus_breadth_csv=breadth_path,
        gift_flora_csv=gift_path,
        assignments_csv=assignments_path,
        status_flora_csv=status_path,
        covariates_csv=covariates_path,
        config_path=config_path,
        output_dir=output,
        predictor_sha256="",
    )
    slopes = pd.read_csv(output / "globi_source_breadth_context_slopes.csv")
    primary = slopes.loc[
        slopes["metric"].eq("effective_channel_number")
        & slopes["outcome"].eq("entry_enrichment")
        & slopes["source_matching"].eq("prevalence_richness")
        & slopes["min_independent_references"].eq(3)
        & slopes["stratum"].eq("all_native")
    ]
    assert len(primary) == 1
    assert primary.iloc[0]["distance_slope"] > 0
    assert primary.iloc[0]["support_class"] == "confirmatory"
    assert manifest["N1_rescued"] is False
    assert manifest["N2_opened"] is False
