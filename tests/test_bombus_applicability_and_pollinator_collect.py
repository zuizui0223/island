import importlib
from pathlib import Path

import pandas as pd

applicability = importlib.import_module("island_v2.bombus_applicability")
pollinator_collect = importlib.import_module("island_v2.gbif_pollinator_collect")


def _regions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source_region_id": "north",
                "biogeographic_region": "temperate_north",
                "floristic_source_pool": "north_pool",
                "plausible_dispersal_connection_or_gateway": "nearest continental temperate source",
                "native_bombus_status": "native_present",
                "applicability": "applicable",
                "applicability_rationale": "Native Bombus is documented in the source pool.",
                "source_region_evidence_id": "ev_north",
                "review_status": "accepted",
                "regional_analysis_only": "false",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
            {
                "source_region_id": "south",
                "biogeographic_region": "tropical_oceanic",
                "floristic_source_pool": "south_pool",
                "plausible_dispersal_connection_or_gateway": "nearest tropical source",
                "native_bombus_status": "native_absent",
                "applicability": "structurally_not_applicable",
                "applicability_rationale": "No native Bombus channel in the predeclared source pool.",
                "source_region_evidence_id": "ev_south",
                "review_status": "accepted",
                "regional_analysis_only": "false",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
        ]
    )


def _evidence() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "source_region_evidence_id": "ev_north",
                "source_region_id": "north",
                "evidence_type": "native_bombus_checklist",
                "source_citation": "North checklist",
                "source_url": "https://example.org/north",
                "evidence_excerpt": "Native Bombus present.",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
            {
                "source_region_evidence_id": "ev_south",
                "source_region_id": "south",
                "evidence_type": "native_bombus_checklist",
                "source_citation": "South checklist",
                "source_url": "https://example.org/south",
                "evidence_excerpt": "Native Bombus absent.",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
        ]
    )


def _assignments() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "island_id": "island_north",
                "source_region_id": "north",
                "assignment_rationale": "Temperate source-pool connection.",
                "assignment_evidence_id": "assignment_north",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
            {
                "island_id": "island_south",
                "source_region_id": "south",
                "assignment_rationale": "Tropical source-pool connection.",
                "assignment_evidence_id": "assignment_south",
                "review_status": "accepted",
                "reviewer": "reviewer",
                "decision_date": "2026-07-04",
            },
        ]
    )


def test_shipped_applicability_config_loads():
    config = applicability.load_config(Path("config/bombus_applicability.yml"))
    assert set(config["categories"]) == {
        "applicable",
        "structurally_not_applicable",
        "unresolved",
    }


def test_applicability_registry_keeps_unassigned_islands_explicitly_unresolved():
    islands = pd.DataFrame({"island_id": ["island_north", "island_south", "island_unassigned"]})
    table = applicability.build_applicability_table(islands, _assignments(), _regions(), _evidence()).set_index(
        "island_id"
    )

    assert table.loc["island_north", "applicability"] == "applicable"
    assert table.loc["island_south", "applicability"] == "structurally_not_applicable"
    assert table.loc["island_unassigned", "applicability"] == "unresolved"
    assert bool(table.loc["island_unassigned", "source_region_assignment_missing"]) is True
    summary = applicability.applicability_summary(table.reset_index())
    assert summary["n_core_confirmatory_islands"] == 1
    assert summary["n_out_of_domain_islands"] == 1
    assert summary["n_unresolved_islands"] == 1


def test_applicability_registry_rejects_prohibited_outcome_columns():
    regions = _regions().assign(floral_traits="not_allowed")
    try:
        applicability.build_applicability_table(
            pd.DataFrame({"island_id": ["island_north"]}),
            _assignments().iloc[[0]],
            regions,
            _evidence(),
        )
    except Exception as exc:
        assert "prohibited" in str(exc)
    else:
        raise AssertionError("Expected prohibited outcome column to be rejected")


def test_exact_pollinator_records_keep_only_exact_island_rows_but_do_not_silently_drop_taxonomy_mismatch():
    assigned = pd.DataFrame(
        {
            "island_id": ["a", "a", pd.NA, "b"],
            "gbif_id": ["1", "2", "3", "4"],
            "block_id": ["block_1", "block_1", "block_1", "block_2"],
            "family": ["Apidae", "Halictidae", "Apidae", "Apidae"],
            "genus": ["Bombus", "Lasioglossum", "Bombus", "Apis"],
            "dataset_key": ["d1", "d1", "d1", "d2"],
            "year": [2020, 2020, 2020, 2021],
        }
    )
    records = pollinator_collect.exact_island_pollinator_records(assigned, "Apidae")

    assert set(records["gbif_id"]) == {"1", "2", "4"}
    assert int(records["expected_family_match"].sum()) == 2
    coverage = pollinator_collect.summarize_pollinator_records(records, "Apidae").set_index("island_id")
    assert coverage.loc["a", "n_exact_island_records"] == 2
    assert coverage.loc["a", "n_expected_family_records"] == 1
    assert coverage.loc["a", "n_bombus_records"] == 1
    assert coverage.loc["b", "n_bombus_records"] == 0
