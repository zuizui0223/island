import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon

from island_v2.gbif_collect import (
    EFFORT_COLUMNS,
    assign_occurrences_to_islands,
    summarize_observation_effort,
)
from island_v2.observation_diagnostics import (
    compute_island_diagnostics,
    coverage_sensitivity,
    load_config,
)

CONFIG = {
    "version": "test",
    "min_species_for_proportions": 5,
    "min_records": 10,
    "flags": {
        "single_dataset_max": 1,
        "low_specimen_support_fraction": 0.2,
        "high_cultivation_fraction": 0.1,
        "coarse_coordinate_uncertainty_m": 5000.0,
        "stale_last_record_year": 1980,
    },
    "species_threshold_sensitivity": [3, 5, 10],
}


def test_shipped_config_loads_and_matches_module_contract():
    config = load_config(__import__("pathlib").Path("config/observation_diagnostics.yml"))
    assert "flags" in config
    assert config["min_species_for_proportions"] >= 1


def test_inclusion_gate_uses_realised_coverage_not_area():
    effort = pd.DataFrame(
        {
            "island_id": ["rich", "sparse"],
            "n_records": [40, 4],
            "n_species": [12, 2],
            "n_datasets": [4, 1],
            "n_preserved_specimen": [30, 0],
            "n_establishment_cultivated": [0, 0],
            "median_coordinate_uncertainty_m": [50, 100],
            "year_max": [2020, 2019],
        }
    )
    diagnostics = compute_island_diagnostics(effort, CONFIG).set_index("island_id")

    # A data-rich island is included; a sparse one is not -- regardless of area.
    assert bool(diagnostics.loc["rich", "analysis_included"]) is True
    assert bool(diagnostics.loc["sparse", "analysis_included"]) is False
    # The sparse, single-source, observation-only, no-voucher island is flagged.
    sparse_flags = set(diagnostics.loc["sparse", "data_process_flags"].split("|"))
    assert {"sparse_flora", "single_dataset", "low_specimen_support"} <= sparse_flags


def test_specimen_and_cultivation_fractions_and_flags():
    effort = pd.DataFrame(
        {
            "island_id": ["cultivated"],
            "n_records": [100],
            "n_species": [20],
            "n_datasets": [5],
            "n_preserved_specimen": [90],
            "n_establishment_cultivated": [40],
            "median_coordinate_uncertainty_m": [10000],
            "year_max": [1975],
        }
    )
    row = compute_island_diagnostics(effort, CONFIG).iloc[0]

    assert row["specimen_fraction"] == 0.9
    assert row["cultivated_fraction"] == 0.4
    flags = set(row["data_process_flags"].split("|"))
    assert {"high_cultivation", "coarse_coordinates", "stale_records"} <= flags
    assert "low_specimen_support" not in flags  # 90% vouchered
    assert bool(row["analysis_included"]) is True  # coverage gate is separate from flags


def test_coverage_sensitivity_counts_decline_with_threshold():
    effort = pd.DataFrame(
        {
            "island_id": ["a", "b", "c"],
            "n_records": [50, 50, 50],
            "n_species": [3, 6, 12],
            "n_datasets": [2, 2, 2],
            "n_preserved_specimen": [10, 10, 10],
            "n_establishment_cultivated": [0, 0, 0],
            "median_coordinate_uncertainty_m": [10, 10, 10],
            "year_max": [2020, 2020, 2020],
        }
    )
    table = coverage_sensitivity(effort, CONFIG).set_index("min_species")
    assert table.loc[3, "n_islands_included"] == 3
    assert table.loc[5, "n_islands_included"] == 2
    assert table.loc[10, "n_islands_included"] == 1


def test_consumes_real_collector_effort_schema_end_to_end():
    # Build a real effort table via the collector, then diagnose it: guards the
    # collect -> diagnostics hand-off against schema drift.
    islands = gpd.GeoDataFrame(
        {"island_id": ["isl_a"]},
        geometry=[Polygon([(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)])],
        crs=4326,
    )
    occ = pd.DataFrame(
        {
            "gbif_id": [str(i) for i in range(12)],
            "dataset_key": ["d1"] * 12,
            "family": ["Fabaceae"] * 12,
            "genus": ["Aaa"] * 12,
            "species": [f"Aaa sp{i}" for i in range(12)],
            "scientific_name": [f"Aaa sp{i} L." for i in range(12)],
            "decimal_longitude": [0.5] * 12,
            "decimal_latitude": [0.5] * 12,
            "coordinate_uncertainty_m": ["100"] * 12,
            "year": ["2020"] * 12,
            "basis_of_record": ["PRESERVED_SPECIMEN"] * 12,
            "establishment_means": [""] * 12,
        }
    )
    assigned = assign_occurrences_to_islands(occ, islands)
    effort = summarize_observation_effort(assigned)

    assert {
        "island_id", "n_records", "n_species", "n_datasets", "n_preserved_specimen",
        "n_establishment_cultivated", "median_coordinate_uncertainty_m", "year_max",
    } <= set(EFFORT_COLUMNS)

    diagnostics = compute_island_diagnostics(effort, CONFIG)
    assert len(diagnostics) == 1
    assert bool(diagnostics.iloc[0]["analysis_included"]) is True  # 12 species, 12 records
