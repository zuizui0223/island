import importlib
from pathlib import Path

import pandas as pd

attrition_audit = importlib.import_module("island_v2.attrition_audit")
bombus_diagnostics = importlib.import_module("island_v2.bombus_diagnostics")


BOMBUS_CONFIG = {
    "version": "test",
    "primary_target_group": {"name": "Apidae"},
    "primary_thresholds": {
        "min_background_records": 5,
        "min_background_spatial_units": 3,
        "min_background_temporal_units": 2,
        "min_distinct_datasets": 2,
        "required_bombus_records_for_detection": 1,
    },
    "spatial_grid_degrees": 0.1,
}

ATTRITION_CONFIG = {
    "version": "test",
    "primary_evidence_track": "direct_conservative",
    "trait_evidence_tracks": {
        "direct_conservative": {
            "accepted_scopes": ["species_direct", "species_indirect"],
            "accepted_reliability": ["A_primary_or_monograph", "B_curated_database_or_institution"],
        },
        "direct_broad_web": {
            "accepted_scopes": ["species_direct", "species_indirect"],
            "accepted_reliability": [
                "A_primary_or_monograph",
                "B_curated_database_or_institution",
                "C_curated_specialist_web",
            ],
        },
    },
    "domain_rules": {
        "pollen_vector": {
            "trait_names": ["pollen_vector_mode"],
            "min_covered_species": 2,
            "min_coverage_fraction": 0.5,
        },
        "floral_signal": {
            "trait_names": ["flower_primary_color"],
            "min_covered_species": 2,
            "min_coverage_fraction": 0.5,
        },
        "floral_architecture": {
            "trait_names": ["floral_form", "floral_symmetry"],
            "min_covered_species": 2,
            "min_coverage_fraction": 0.5,
        },
        "reproductive_assurance": {
            "trait_names": ["self_incompatibility", "autonomous_selfing_capacity"],
            "min_covered_species": 2,
            "min_coverage_fraction": 0.5,
        },
    },
    "attrition_rules": {
        "confirmatory_applicability": "applicable",
        "evaluable_bombus_evidence": ["detected", "adequate_non_detection"],
    },
    "continuation_thresholds": {
        "m1_pilot_min_islands": 1,
        "m1_confirmatory_target_islands": 2,
        "m2_m3_exploratory_min_islands": 1,
        "m2_m3_confirmatory_target_islands": 2,
    },
}


def _record(
    island_id: str,
    genus: str,
    *,
    family: str = "Apidae",
    year: int = 2020,
    dataset: str = "d1",
    lon: float = 0.01,
    lat: float = 0.01,
    gbif_id: str = "",
) -> dict[str, object]:
    return {
        "island_id": island_id,
        "family": family,
        "genus": genus,
        "dataset_key": dataset,
        "year": year,
        "decimal_latitude": lat,
        "decimal_longitude": lon,
        "basis_of_record": "PRESERVED_SPECIMEN",
        "gbif_id": gbif_id,
    }


def test_shipped_bombus_and_attrition_configs_load():
    assert (
        bombus_diagnostics.load_config(Path("config/bombus_observation_diagnostics.yml"))["primary_target_group"]["name"]
        == "Apidae"
    )
    assert (
        attrition_audit.load_config(Path("config/attrition_audit.yml"))["primary_evidence_track"]
        == "direct_conservative"
    )


def test_adequate_non_detection_requires_background_spatial_temporal_and_dataset_coverage():
    records = [
        _record("adequate", "Apis", year=2020, dataset="d1", lon=0.01, lat=0.01, gbif_id="1"),
        _record("adequate", "Andrena", year=2020, dataset="d1", lon=0.21, lat=0.01, gbif_id="2"),
        _record("adequate", "Osmia", year=2021, dataset="d2", lon=0.41, lat=0.01, gbif_id="3"),
        _record("adequate", "Eucera", year=2021, dataset="d2", lon=0.41, lat=0.21, gbif_id="4"),
        _record("adequate", "Anthophora", year=2021, dataset="d2", lon=0.21, lat=0.21, gbif_id="5"),
        _record("sparse", "Apis", year=2020, dataset="d1", lon=1.01, lat=1.01, gbif_id="6"),
    ]
    result = bombus_diagnostics.compute_bombus_diagnostics(
        pd.DataFrame(records), BOMBUS_CONFIG
    ).set_index("island_id")

    assert result.loc["adequate", "bombus_occurrence_evidence"] == "adequate_non_detection"
    assert bool(result.loc["adequate", "bombus_channel_evaluable"]) is True
    assert result.loc["sparse", "bombus_occurrence_evidence"] == "insufficient_effort"
    sparse_flags = set(result.loc["sparse", "observation_diagnostic_flags"].split("|"))
    assert {"below_background_records", "low_spatial_dispersion", "low_temporal_dispersion"} <= sparse_flags


def test_detected_bombus_is_not_reclassified_as_non_detection_and_duplicate_gbif_ids_do_not_inflate_background():
    records = [
        _record("detected", "Bombus", gbif_id="same"),
        _record("detected", "Bombus", gbif_id="same"),
        _record("detected", "Apis", gbif_id="other"),
    ]
    result = bombus_diagnostics.compute_bombus_diagnostics(
        pd.DataFrame(records), BOMBUS_CONFIG
    ).iloc[0]

    assert result["bombus_occurrence_evidence"] == "detected"
    assert result["bombus_record_count"] == 1
    assert result["target_group_record_count"] == 2
    assert bool(result["bombus_channel_evaluable"]) is True


def test_unresolved_target_taxonomy_is_not_silently_called_absent():
    records = [_record("unknown", "", gbif_id="1")]
    result = bombus_diagnostics.compute_bombus_diagnostics(
        pd.DataFrame(records), BOMBUS_CONFIG
    ).iloc[0]

    assert result["bombus_occurrence_evidence"] == "unresolved"
    assert "unresolved_target_taxonomy" in result["observation_diagnostic_flags"]


def _accepted_evidence(species: list[str]) -> pd.DataFrame:
    rows = []
    for taxon in species:
        for trait_name in [
            "pollen_vector_mode",
            "flower_primary_color",
            "floral_form",
            "self_incompatibility",
        ]:
            rows.append(
                {
                    "accepted_species": taxon,
                    "trait_name": trait_name,
                    "evidence_scope": "species_direct",
                    "source_reliability": "B_curated_database_or_institution",
                    "review_status": "accepted",
                }
            )
    return pd.DataFrame(rows)


def test_trait_coverage_and_attrition_keep_applicability_non_detection_and_trait_coverage_separate():
    island_species = pd.DataFrame(
        {
            "island_id": ["a", "a", "b", "b", "c", "c"],
            "species": ["A one", "A two", "B one", "B two", "C one", "C two"],
        }
    )
    evidence = _accepted_evidence(["A one", "A two", "B one", "B two", "C one", "C two"])
    coverage = attrition_audit.compute_trait_coverage(island_species, evidence, ATTRITION_CONFIG)

    flora = pd.DataFrame({"island_id": ["a", "b", "c"], "analysis_included": [True, True, True]})
    applicability = pd.DataFrame(
        {
            "island_id": ["a", "b", "c"],
            "applicability": ["applicable", "applicable", "structurally_not_applicable"],
        }
    )
    bombus = pd.DataFrame(
        {
            "island_id": ["a", "b", "c"],
            "bombus_occurrence_evidence": ["detected", "insufficient_effort", "detected"],
        }
    )

    attrition = attrition_audit.build_attrition_table(
        flora, applicability, bombus, coverage, ATTRITION_CONFIG
    ).set_index("island_id")
    assert bool(attrition.loc["a", "m1_eligible"]) is True
    assert bool(attrition.loc["a", "m2_m3_eligible"]) is True
    assert bool(attrition.loc["b", "m1_eligible"]) is False  # observation effort, not floral coverage
    assert bool(attrition.loc["c", "m1_eligible"]) is False  # out of Bombus domain, not Bombus deficit
    assert attrition.loc["c", "applicability"] == "structurally_not_applicable"

    stages = attrition_audit.attrition_stage_summary(
        attrition.reset_index(), ATTRITION_CONFIG
    ).set_index("stage")
    assert stages.loc["all_islands", "n_islands"] == 3
    assert stages.loc["m1_eligible", "n_islands"] == 1
    summary = attrition_audit.continuation_summary(attrition.reset_index(), ATTRITION_CONFIG)
    assert summary["m1_status"] == "exploratory_or_pilot_only"
    assert summary["m2_m3_status"] == "exploratory_or_pilot_only"
