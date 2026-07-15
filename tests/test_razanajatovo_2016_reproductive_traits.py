from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import yaml

from island_v2.razanajatovo_2016_reproductive_traits import (
    classify_autofertility,
    classify_self_compatibility,
    prepare_razanajatovo_2016,
)


def test_index_mappings_do_not_turn_low_compatibility_into_si() -> None:
    assert classify_self_compatibility([0.9, 1.0], lower=0.2, upper=0.8) == "SC"
    assert classify_self_compatibility([0.0, 0.2], lower=0.2, upper=0.8) == ""
    assert (
        classify_self_compatibility([0.1, 0.6], lower=0.2, upper=0.8)
        == "mixed_or_variable"
    )
    assert classify_autofertility([0.0, 0.2], lower=0.2, upper=0.8) == "absent"
    assert classify_autofertility([0.8, 1.0], lower=0.2, upper=0.8) == "autonomous"
    assert (
        classify_autofertility([0.1, 0.9], lower=0.2, upper=0.8)
        == "mixed_or_variable"
    )


def test_prepare_writes_separate_compatibility_and_autofertility_traits(tmp_path: Path) -> None:
    workbook = tmp_path / "source.xlsx"
    workbook.write_bytes(b"test workbook placeholder")
    sheet = tmp_path / "sheet.csv"
    pd.DataFrame(
        [
            {
                "Data_source": "Ref_A",
                "Species": "Alpha_a",
                "Family": "FamA",
                "Self-compatibility_index_FS": "0.9",
                "Self-compatibility_index_SFL": "1.0",
                "Autofertility_index_FS": "0.1",
                "Autofertility_index_SFL": "",
            },
            {
                "Data_source": "Ref_B",
                "Species": "Alpha_b",
                "Family": "FamA",
                "Self-compatibility_index_FS": "0.1",
                "Self-compatibility_index_SFL": "",
                "Autofertility_index_FS": "0.9",
                "Autofertility_index_SFL": "",
            },
            {
                "Data_source": "Ref_C",
                "Species": "Alpha_c",
                "Family": "FamA",
                "Self-compatibility_index_FS": "0.5",
                "Self-compatibility_index_SFL": "",
                "Autofertility_index_FS": "0.5",
                "Autofertility_index_SFL": "",
            },
        ]
    ).to_csv(sheet, index=False)
    master = tmp_path / "master.csv"
    pd.DataFrame(
        {
            "accepted_species": ["Alpha a", "Alpha b", "Alpha c"],
            "genus": ["Alpha", "Alpha", "Alpha"],
            "family": ["FamA", "FamA", "FamA"],
            "n_islands": [1, 1, 1],
            "n_records": [1, 1, 1],
        }
    ).to_csv(master, index=False)
    config = {
        "source": {
            "name": "razanajatovo_test",
            "citation": "Test citation",
            "article_url": "https://example.test/article",
            "license": "CC-BY-4.0",
            "expected_workbook_sha256": hashlib.sha256(workbook.read_bytes()).hexdigest(),
            "expected_sheet_csv_sha256": hashlib.sha256(sheet.read_bytes()).hexdigest(),
            "expected_rows": 3,
        },
        "columns": {
            "scientific_name": "Species",
            "family": "Family",
            "source_reference": "Data_source",
            "self_compatibility_indices": [
                "Self-compatibility_index_FS",
                "Self-compatibility_index_SFL",
            ],
            "autofertility_indices": [
                "Autofertility_index_FS",
                "Autofertility_index_SFL",
            ],
        },
        "thresholds": {"lower": 0.2, "upper": 0.8},
        "family_aliases": {},
        "gbif_reconciliation": {
            "endpoint": "https://api.gbif.org/v1/species/match",
            "minimum_confidence": 95,
            "allowed_match_type": "EXACT",
            "allowed_rank": "SPECIES",
            "allowed_kingdom": "Plantae",
            "allowed_status": ["SYNONYM"],
        },
    }
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    scope_config = tmp_path / "scope.yml"
    scope_config.write_text(
        yaml.safe_dump(
            {
                "family_column": "family",
                "eligible_group": "angiosperm",
                "unresolved_group": "family_unresolved",
                "non_angiosperm_families": {},
            }
        ),
        encoding="utf-8",
    )

    report = prepare_razanajatovo_2016(
        workbook_path=workbook,
        sheet_csv=sheet,
        master_csv=master,
        output_dir=tmp_path / "output",
        config_path=config_path,
        scope_config_path=scope_config,
        resolve_gbif=False,
    )
    candidates = pd.read_csv(report["outputs"]["trait_candidates"], dtype=str).fillna("")
    observed = set(
        zip(
            candidates["accepted_species"],
            candidates["trait_name"],
            candidates["standardized_value"],
            strict=True,
        )
    )
    assert observed == {
        ("Alpha a", "self_incompatibility", "SC"),
        ("Alpha a", "autonomous_selfing_capacity", "absent"),
        ("Alpha b", "autonomous_selfing_capacity", "autonomous"),
        ("Alpha c", "self_incompatibility", "mixed_or_variable"),
        ("Alpha c", "autonomous_selfing_capacity", "mixed_or_variable"),
    }
    assert "SI" not in set(candidates["standardized_value"])
    assert report["n_direct_candidate_cells"] == 5
