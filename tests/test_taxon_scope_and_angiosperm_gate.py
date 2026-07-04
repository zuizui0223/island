import importlib

import pandas as pd

observation_diagnostics = importlib.import_module("island_v2.observation_diagnostics")
taxon_scope = importlib.import_module("island_v2.taxon_scope")


CONFIG = {
    "version": "test",
    "min_species_for_proportions": 2,
    "min_records": 4,
    "flags": {
        "single_dataset_max": 1,
        "low_specimen_support_fraction": 0.2,
        "high_cultivation_fraction": 0.1,
        "coarse_coordinate_uncertainty_m": 5000,
        "stale_last_record_year": 1980,
    },
    "species_threshold_sensitivity": [1, 2],
}


def _scope_decisions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": "Flower one",
                "taxonomic_group": "angiosperm",
                "decision_basis": "reviewed backbone taxonomy",
                "source_citation": "GBIF backbone",
                "source_url": "https://example.org/flower-one",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
            {
                "accepted_species": "Flower two",
                "taxonomic_group": "angiosperm",
                "decision_basis": "reviewed backbone taxonomy",
                "source_citation": "GBIF backbone",
                "source_url": "https://example.org/flower-two",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
            {
                "accepted_species": "Fern one",
                "taxonomic_group": "non_seed_vascular",
                "decision_basis": "reviewed backbone taxonomy",
                "source_citation": "GBIF backbone",
                "source_url": "https://example.org/fern-one",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
        ]
    )


def test_taxon_scope_excludes_ferns_and_unreviewed_species_from_angiosperm_universe():
    island_species = pd.DataFrame(
        {
            "island_id": ["flowers", "flowers", "ferns", "ferns", "unreviewed"],
            "species": ["Flower one", "Flower two", "Fern one", "Fern one", "Unknown one"],
            "n_records": [2, 3, 20, 1, 15],
        }
    )
    scoped, angiosperms = taxon_scope.build_taxon_scope(island_species, _scope_decisions())
    coverage = taxon_scope.island_angiosperm_coverage(scoped, island_species).set_index("island_id")

    assert set(angiosperms["species"]) == {"Flower one", "Flower two"}
    assert coverage.loc["flowers", "n_accepted_angiosperm_species"] == 2
    assert coverage.loc["flowers", "n_accepted_angiosperm_records"] == 5
    assert coverage.loc["ferns", "n_accepted_angiosperm_species"] == 0
    assert coverage.loc["unreviewed", "n_accepted_angiosperm_species"] == 0


def test_primary_observation_gate_uses_accepted_angiosperm_not_raw_vascular_counts():
    effort = pd.DataFrame(
        {
            "island_id": ["flowers", "ferns", "unreviewed"],
            "n_records": [5, 21, 15],
            "n_species": [2, 1, 1],
            "n_datasets": [2, 2, 2],
            "n_preserved_specimen": [5, 21, 15],
            "n_establishment_cultivated": [0, 0, 0],
            "median_coordinate_uncertainty_m": [50, 50, 50],
            "year_max": [2020, 2020, 2020],
        }
    )
    island_species = pd.DataFrame(
        {
            "island_id": ["flowers", "flowers", "ferns", "unreviewed"],
            "species": ["Flower one", "Flower two", "Fern one", "Unknown one"],
            "n_records": [2, 3, 21, 15],
        }
    )
    scoped, _ = taxon_scope.build_taxon_scope(island_species, _scope_decisions())
    coverage = taxon_scope.island_angiosperm_coverage(scoped, island_species)
    diagnostics = observation_diagnostics.compute_island_diagnostics(
        effort, CONFIG, coverage
    ).set_index("island_id")

    assert diagnostics.loc["flowers", "taxonomic_scope_for_inclusion"] == "accepted_angiosperms"
    assert bool(diagnostics.loc["flowers", "analysis_included"]) is True
    assert bool(diagnostics.loc["ferns", "analysis_included"]) is False
    assert bool(diagnostics.loc["unreviewed", "analysis_included"]) is False
    assert "sparse_flora" in diagnostics.loc["ferns", "data_process_flags"]


def test_raw_vascular_screening_remains_explicit_when_scope_table_is_missing():
    effort = pd.DataFrame(
        {
            "island_id": ["raw"],
            "n_records": [10],
            "n_species": [5],
            "n_datasets": [2],
            "n_preserved_specimen": [10],
            "n_establishment_cultivated": [0],
            "median_coordinate_uncertainty_m": [50],
            "year_max": [2020],
        }
    )
    diagnostics = observation_diagnostics.compute_island_diagnostics(effort, CONFIG)

    assert diagnostics.loc[0, "taxonomic_scope_for_inclusion"] == "raw_tracheophyta_screening_only"
    assert bool(diagnostics.loc[0, "analysis_included"]) is True
