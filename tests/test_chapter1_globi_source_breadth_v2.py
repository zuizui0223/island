from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from island_v2.chapter1_globi_source_breadth_v2 import (
    analyse_source_breadth_v2,
    compute_effort_matched_enrichment,
    reference_effort_bins,
)


def _config() -> dict:
    return {
        "contract": "chapter1_globi_source_breadth_v2",
        "sampling_effort_gate": {
            "primary_min_independent_references": 3,
            "sensitivity_min_independent_references": [2, 5],
        },
        "source_matching": {
            "source_modes": ["geo_k5"],
            "primary_matching": "prevalence_richness_effort",
            "sensitivity_matching": "prevalence_richness",
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
            "source_matching": "prevalence_richness_effort",
            "effort_threshold": 3,
        },
        "claim_ceiling": "test v2 boundary",
    }


def test_reference_effort_bins_and_matching() -> None:
    assert reference_effort_bins(np.array([1, 2, 3, 4, 7, 8])).tolist() == [0, 1, 1, 2, 2, 3]
    prevalence = np.array([1, 1, 1, 1, 1, 1])
    richness = np.ones(6)
    counts = np.array([1, 1, 1, 1, 1, 0], dtype=float)
    positions = np.array([1, 2, 3, 4, 5, 9], dtype=float)
    effort = np.array([1, 1, 1, 1, 1, 3], dtype=np.int16)
    matched = compute_effort_matched_enrichment(
        prevalence, richness, counts, positions, effort,
        matching="prevalence_richness_effort", minimum_represented_genera=5,
    )
    unmatched = compute_effort_matched_enrichment(
        prevalence, richness, counts, positions, effort,
        matching="prevalence_richness", minimum_represented_genera=5,
    )
    assert matched is not None and unmatched is not None
    assert abs(matched["entry_enrichment"]) < 1e-12
    assert unmatched["entry_enrichment"] < 0


def test_v2_analysis_runs_with_frozen_predictor(tmp_path: Path) -> None:
    genera = ["Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta"]
    breadth = pd.DataFrame(
        {
            "genus": genera,
            "n_independent_references": [6, 6, 6, 6, 6, 6],
            "sampled_channel_count": [1, 1, 2, 2, 3, 5],
            "effective_channel_number": [1.0, 1.2, 1.8, 2.2, 3.0, 4.8],
        }
    )
    breadth_path = tmp_path / "breadth.csv"
    breadth.to_csv(breadth_path, index=False)
    pd.DataFrame(
        [{"entity_ID": 1, "work_species": f"{genus} alba"} for genus in genera]
    ).to_csv(tmp_path / "gift.csv", index=False)

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
                "log_island_area_km2": np.sin(x) + x / 20,
                "climate_pc1": np.cos(x),
                "climate_pc2": np.sin(x / 2),
                "climate_pc3": np.cos(x / 2),
                "climate_pc4": np.sin(x / 3),
            }
        )
    pd.DataFrame(status_rows).to_csv(tmp_path / "status.csv", index=False)
    pd.DataFrame(assignment_rows).to_csv(tmp_path / "assignments.csv", index=False)
    pd.DataFrame(cov_rows).to_csv(tmp_path / "covariates.csv", index=False)
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(_config()), encoding="utf-8")

    manifest = analyse_source_breadth_v2(
        genus_breadth_csv=breadth_path,
        gift_flora_csv=tmp_path / "gift.csv",
        assignments_csv=tmp_path / "assignments.csv",
        status_flora_csv=tmp_path / "status.csv",
        covariates_csv=tmp_path / "covariates.csv",
        config_path=config_path,
        output_dir=tmp_path / "out",
    )
    slopes = pd.read_csv(tmp_path / "out/globi_source_breadth_v2_context_slopes.csv")
    primary = slopes.loc[
        slopes["metric"].eq("effective_channel_number")
        & slopes["outcome"].eq("entry_enrichment")
        & slopes["source_matching"].eq("prevalence_richness_effort")
        & slopes["min_independent_references"].eq(3)
        & slopes["stratum"].eq("all_native")
    ]
    assert len(primary) == 1
    assert primary.iloc[0]["distance_slope"] > 0
    assert manifest["N1_rescued"] is False
    assert manifest["N2_opened"] is False
