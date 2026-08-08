from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import yaml

from island_v2.ferrer_2024_si_traits import (
    CANDIDATE_COLUMNS,
    _safe_source_rows,
    prepare_ferrer_2024,
)


def _write_workbook(path: Path) -> None:
    references = pd.DataFrame(
        [
            {
                "SI_orden": "1",
                "Family": "FamA",
                "Genera": "Alpha",
                "species": "one",
                "work_species": "Alpha one",
                "SI analyses": "SC",
                "Review1": "Review A, 2001",
            },
            {
                "SI_orden": "2",
                "Family": "FamB",
                "Genera": "Beta",
                "species": "two",
                "work_species": "Beta two",
                "SI analyses": "SI",
                "Reference1": "Study B, 2002",
            },
            {
                "SI_orden": "3",
                "Family": "FamC",
                "Genera": "Gamma",
                "species": "three",
                "work_species": "Gamma three",
                "SI analyses": "SC-SI",
                "Reference1": "Study C, 2003",
            },
            {
                "SI_orden": "4",
                "Family": "FamD",
                "Genera": "Delta",
                "species": "four",
                "work_species": "Delta four",
                "SI analyses": "SI",
                "Reference1": "Review-only D, 2004",
            },
            {
                "SI_orden": "5",
                "Family": "FamE",
                "Genera": "Epsilon",
                "species": "five",
                "work_species": "Epsilon five",
                "SI analyses": "SS",
                "Reference1": "Study E, 2005",
            },
            {
                "SI_orden": "6",
                "Family": "WrongFamily",
                "Genera": "Zeta",
                "species": "six",
                "work_species": "Zeta six",
                "SI analyses": "SC",
                "Reference1": "Study F, 2006",
            },
            {
                "SI_orden": "7",
                "Family": "FamA",
                "Genera": "Alpha",
                "species": "one",
                "work_species": "Alpha one",
                "SI analyses": "SI",
                "Reference1": "Study G, 2007",
            },
        ]
    )
    phenotype = pd.DataFrame(
        [
            {
                "SI_orden": "2",
                "work_species": "Beta two",
                "SI status": "SI",
                "SI reaction": "Stylar SI",
                "Notes BS": "Pollen tubes arrested in the style",
                "Ref BS": "Study B, 2002",
                "Ref Gen": "",
            },
            {
                "SI_orden": "3",
                "work_species": "Gamma three",
                "SI status": "SC-SI",
                "SI reaction": "Late acting SI",
                "Notes BS": "Variable fruit set",
                "Ref BS": "Study C, 2003",
                "Ref Gen": "",
            },
            {
                "SI_orden": "7",
                "work_species": "Alpha one",
                "SI status": "SI",
                "SI reaction": "Stigmatic SI",
                "Notes BS": "Pollen arrested",
                "Ref BS": "Study G, 2007",
                "Ref Gen": "",
            },
        ]
    )
    with pd.ExcelWriter(path, engine="openpyxl") as writer:
        references.to_excel(writer, sheet_name="References", index=False)
        phenotype.to_excel(writer, sheet_name="SI Phenotype data", index=False)


def _config(path: Path, workbook: Path) -> None:
    path.write_text(
        yaml.safe_dump(
            {
                "source": {
                    "name": "ferrer_2024_si_database",
                    "citation": "Ferrer test database",
                    "dataset_url": "https://doi.org/10.1093/aob/mcae056",
                    "license": "test",
                    "expected_workbook_sha256": hashlib.sha256(
                        workbook.read_bytes()
                    ).hexdigest(),
                    "expected_source_rows": 7,
                },
                "gbif_reconciliation": {
                    "endpoint": "https://api.gbif.org/v1/species/match",
                    "minimum_confidence": 95,
                    "allowed_match_type": "EXACT",
                    "allowed_rank": "SPECIES",
                    "allowed_kingdom": "Plantae",
                    "allowed_status": ["SYNONYM"],
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )


def _scope_config(path: Path) -> None:
    path.write_text(
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


def test_safe_subset_does_not_treat_genus_si_or_self_sterility_as_direct(
    tmp_path: Path,
) -> None:
    workbook = tmp_path / "source.xlsx"
    _write_workbook(workbook)
    safe, exclusions = _safe_source_rows(workbook)

    values = set(zip(safe["source_record_id"], safe["standardized_value"], strict=True))
    assert ("1", "SC") in values
    assert ("2", "SI") in values
    assert ("3", "mixed_or_variable") in values
    assert ("4", "SI") not in values
    assert ("5", "SC") not in values
    reasons = set(exclusions["reason"])
    assert "dataset_level_SI_without_confirming_phenotype" in reasons
    assert "self_sterility_not_self_incompatibility" in reasons


def test_prepare_quarantines_conflicts_family_mismatch_and_emits_no_genus_low(
    tmp_path: Path,
) -> None:
    workbook = tmp_path / "source.xlsx"
    _write_workbook(workbook)
    config = tmp_path / "config.yml"
    _config(config, workbook)
    scope = tmp_path / "scope.yml"
    _scope_config(scope)
    master = tmp_path / "master.csv"
    pd.DataFrame(
        [
            ("Alpha one", "Alpha", "FamA"),
            ("Beta two", "Beta", "FamB"),
            ("Gamma three", "Gamma", "FamC"),
            ("Delta four", "Delta", "FamD"),
            ("Epsilon five", "Epsilon", "FamE"),
            ("Zeta six", "Zeta", "FamZ"),
            ("Gamma target", "Gamma", "FamC"),
        ],
        columns=["accepted_species", "genus", "family"],
    ).to_csv(master, index=False)
    meyer = tmp_path / "meyer.csv.gz"
    pd.DataFrame(
        [{"accepted_species": "Beta two", "standardized_value": "SI"}]
    ).to_csv(meyer, index=False)

    report = prepare_ferrer_2024(
        source_xlsx=workbook,
        master_csv=master,
        output_dir=tmp_path / "output",
        config_path=config,
        scope_config_path=scope,
        meyer_candidates_csv=meyer,
        resolve_gbif=False,
    )
    candidates = pd.read_csv(report["outputs"]["trait_candidates"], dtype=str).fillna("")
    audit = pd.read_csv(report["outputs"]["name_match_audit"], dtype=str).fillna("")

    assert list(candidates.columns) == CANDIDATE_COLUMNS
    assert set(candidates["accepted_species"]) == {"Beta two", "Gamma three"}
    assert dict(zip(candidates["accepted_species"], candidates["standardized_value"])) == {
        "Beta two": "SI",
        "Gamma three": "mixed_or_variable",
    }
    assert candidates["inference_rule"].eq("").all()
    assert candidates["evidence_scope"].eq("species_direct").all()
    assert candidates.loc[
        candidates["accepted_species"].eq("Beta two"), "source_lineage"
    ].item() == "doi:10.5061/dryad.cc2fqz6hr"
    assert report["genus_inference_emitted"] is False
    assert report["conflicting_master_species"] == ["Alpha one"]
    assert audit.loc[audit["source_scientific_name"].eq("Zeta six"), "reason"].item() == (
        "reject_family_conflict"
    )
    assert not (tmp_path / "output" / "genus_inference_candidates.csv.gz").exists()
