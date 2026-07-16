from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

from island_v2.meyer_2026_reproductive_traits import (
    CANDIDATE_COLUMNS,
    build_genus_inferences,
    prepare_meyer_2026,
    strict_gbif_resolution,
)


def _config() -> dict:
    return {
        "source": {
            "name": "meyer_test",
            "citation": "Meyer test source",
            "dataset_url": "https://doi.org/10.5061/dryad.cc2fqz6hr",
            "license": "CC0-1.0",
            "expected_output_sha256": "",
            "expected_output_rows": 0,
        },
        "value_map": {"sc": "SC", "si": "SI"},
        "excluded_dio_only_names": ["Bad dioica"],
        "gbif_reconciliation": {
            "endpoint": "https://api.gbif.org/v1/species/match",
            "minimum_confidence": 95,
            "allowed_match_type": "EXACT",
            "allowed_rank": "SPECIES",
            "allowed_kingdom": "Plantae",
            "allowed_status": ["SYNONYM"],
        },
        "genus_inference": {
            "enabled": True,
            "target_trait": "self_incompatibility",
            "minimum_direct_species": 5,
            "minimum_modal_purity": 1.0,
            "confidence": "low",
            "analysis_role": "sensitivity_only_genus_inference",
        },
    }


def test_strict_gbif_resolution_accepts_only_exact_family_consistent_synonym() -> None:
    master = {
        "feijoa sellowiana": {
            "accepted_species": "Feijoa sellowiana",
            "family": "Myrtaceae",
        }
    }
    reconciliation = _config()["gbif_reconciliation"]
    payload = {
        "matchType": "EXACT",
        "rank": "SPECIES",
        "kingdom": "Plantae",
        "status": "SYNONYM",
        "confidence": 98,
        "acceptedUsageKey": 3181977,
        "species": "Feijoa sellowiana",
        "family": "Myrtaceae",
    }
    accepted, reason = strict_gbif_resolution(
        source_name="Acca sellowiana",
        source_family="Myrtaceae",
        payload=payload,
        master_by_name=master,
        reconciliation=reconciliation,
    )
    assert (accepted, reason) == ("Feijoa sellowiana", "accepted_exact_gbif_synonym")

    fuzzy = {**payload, "matchType": "FUZZY"}
    assert strict_gbif_resolution(
        source_name="Acca sellowiana",
        source_family="Myrtaceae",
        payload=fuzzy,
        master_by_name=master,
        reconciliation=reconciliation,
    )[1] == "reject_non_exact_match"

    wrong_family = {**payload, "family": "Rosaceae"}
    assert strict_gbif_resolution(
        source_name="Acca sellowiana",
        source_family="Myrtaceae",
        payload=wrong_family,
        master_by_name=master,
        reconciliation=reconciliation,
    )[1] == "reject_family_conflict"


def test_genus_inference_requires_unanimous_five_species_and_is_not_direct() -> None:
    master = pd.DataFrame(
        [
            *[(f"Alpha {name}", "Alpha", "FamA") for name in "abcdefg"],
            ("Beta a", "Beta", "FamB"),
            ("Beta b", "Beta", "FamB"),
            ("Beta target", "Beta", "FamB"),
        ],
        columns=["accepted_species", "genus", "family"],
    )
    direct = pd.DataFrame(
        [
            *[
                {
                    "accepted_species": f"Alpha {name}",
                    "genus": "Alpha",
                    "family": "FamA",
                    "standardized_value": "SC",
                }
                for name in "abcdef"
            ],
            {
                "accepted_species": "Beta a",
                "genus": "Beta",
                "family": "FamB",
                "standardized_value": "SC",
            },
            {
                "accepted_species": "Beta b",
                "genus": "Beta",
                "family": "FamB",
                "standardized_value": "SI",
            },
        ]
    )
    inference, report = build_genus_inferences(master, direct, _config())
    assert inference["accepted_species"].tolist() == ["Alpha g"]
    assert inference["candidate_kind"].eq("hierarchical_inference").all()
    assert inference["evidence_scope"].eq("genus_inference").all()
    assert inference["confidence"].eq("low").all()
    assert report["n_qualifying_genera"] == 1
    assert report["leave_one_out_n"] == 6
    assert report["leave_one_out_accuracy"] == 1.0


def test_prepare_excludes_dio_and_writes_direct_plus_genus_layers(tmp_path: Path) -> None:
    source = tmp_path / "Output_Data_MS.csv"
    pd.DataFrame(
        [
            *[
                {
                    "Unnamed: 0": index,
                    "genus_species": f"Alpha {name}",
                    "Family": "FamA",
                    "mating_system": "sc",
                }
                for index, name in enumerate("abcde", start=1)
            ],
            {
                "Unnamed: 0": 6,
                "genus_species": "Bad dioica",
                "Family": "FamA",
                "mating_system": "si",
            },
            {
                "Unnamed: 0": 7,
                "genus_species": "Unused name",
                "Family": "FamA",
                "mating_system": "na",
            },
        ]
    ).to_csv(source, index=False)

    master = tmp_path / "master.csv"
    pd.DataFrame(
        [
            *[
                {
                    "accepted_species": f"Alpha {name}",
                    "genus": "Alpha",
                    "family": "FamA",
                    "n_islands": 1,
                    "n_records": 1,
                }
                for name in "abcdef"
            ],
            {
                "accepted_species": "Bad dioica",
                "genus": "Bad",
                "family": "FamA",
                "n_islands": 1,
                "n_records": 1,
            },
        ]
    ).to_csv(master, index=False)

    config = _config()
    config["source"]["expected_output_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    config["source"]["expected_output_rows"] = 7
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

    report = prepare_meyer_2026(
        source_csv=source,
        master_csv=master,
        output_dir=tmp_path / "output",
        config_path=config_path,
        scope_config_path=scope_config,
        resolve_gbif=False,
    )
    direct = pd.read_csv(report["outputs"]["trait_candidates"], dtype=str).fillna("")
    inferred = pd.read_csv(
        report["outputs"]["genus_inference_candidates"], dtype=str
    ).fillna("")
    assert list(direct.columns) == CANDIDATE_COLUMNS
    assert len(direct) == 5
    assert set(direct["standardized_value"]) == {"SC"}
    assert direct["candidate_kind"].eq("source_backed").all()
    assert inferred["accepted_species"].tolist() == ["Alpha f"]
    assert report["n_dio_only_rows_excluded"] == 1
    assert report["n_source_rows"] == 7
    assert json.loads(Path(report["manifest"]).read_text(encoding="utf-8"))["source_license"] == (
        "CC0-1.0"
    )
