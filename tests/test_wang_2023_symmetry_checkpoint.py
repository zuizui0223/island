from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.wang_2023_symmetry_checkpoint import build_checkpoint

REPO_ROOT = Path(__file__).resolve().parents[1]


def _source(species: str, family: str, value: str = "RadSymmetry") -> dict[str, str]:
    return {
        "source_row": species,
        "order": "Exampleales",
        "family": family,
        "genus": species.split()[0],
        "accepted_species": species,
        "raw_symmetry": value,
        "note1": "",
        "note2": "",
    }


def test_family_validation_selects_only_exact_high_precision_species_rows() -> None:
    source = pd.DataFrame(
        [
            _source("Alpha one", "Passaceae"),
            _source("Alpha two", "Passaceae"),
            _source("Alpha three", "Passaceae"),
            _source("Beta one", "Failaceae"),
            _source("Beta two", "Failaceae", "BiSymmetry"),
            _source("Beta three", "Failaceae"),
            _source("Gamma one", "Asteraceae"),
            _source("Gamma two", "Asteraceae"),
            _source("Delta one", "Passaceae"),
            _source("Outside one", "Passaceae"),
            _source("Wrong one", "Wrongaceae"),
        ]
    )
    master = pd.DataFrame(
        {
            "accepted_species": [
                "Alpha one",
                "Alpha two",
                "Alpha three",
                "Beta one",
                "Beta two",
                "Beta three",
                "Gamma one",
                "Gamma two",
                "Delta one",
                "Wrong one",
            ],
            "family": [
                "Passaceae",
                "Passaceae",
                "Passaceae",
                "Failaceae",
                "Failaceae",
                "Failaceae",
                "Asteraceae",
                "Asteraceae",
                "Passaceae",
                "Passaceae",
            ],
        }
    )
    direct = pd.DataFrame(
        {
            "accepted_species": [
                "Alpha one",
                "Alpha two",
                "Beta one",
                "Beta two",
                "Gamma one",
                "Gamma two",
                "Delta one",
            ],
            "trait_name": ["floral_symmetry"] * 7,
            "normalized_value": [
                "actinomorphic",
                "actinomorphic",
                "zygomorphic",
                "zygomorphic",
                "actinomorphic",
                "actinomorphic",
                "asymmetric",
            ],
            "source_lineages": [
                "origin:one",
                "origin:two",
                "origin:three",
                "origin:four",
                "origin:five",
                "origin:six",
                "origin:seven",
            ],
        }
    )

    evidence, audit, selection, families = build_checkpoint(
        source,
        direct,
        master,
        min_family_audit=2,
        min_family_precision=0.95,
    )

    assert evidence["accepted_species"].tolist() == [
        "Alpha one",
        "Alpha two",
        "Alpha three",
        "Delta one",
    ]
    assert evidence["quality"].eq("medium").all()
    assert evidence["trait_name"].eq("floral_symmetry").all()
    assert evidence["axis"].eq("floral_structural_complexity").all()
    assert evidence.set_index("accepted_species").loc["Alpha one", "source_lineage"] == "origin:one"
    assert (
        "upstream-row-lineage-unmapped"
        in evidence.set_index("accepted_species").loc["Alpha three", "source_lineage"]
    )
    assert len(audit) == 3
    assert audit["accepted_correct"].eq("true").sum() == 2
    assert audit.set_index("accepted_species").loc["Delta one", "accepted_correct"] == "false"
    assert families.set_index("family").loc["Passaceae", "production_approved"]
    assert not families.set_index("family").loc["Asteraceae", "production_approved"]

    reasons = selection.set_index("accepted_species")["selection_reason"].to_dict()
    assert reasons["Beta one"] == "family_precision_below_threshold"
    assert reasons["Beta three"] == "family_precision_below_threshold"
    assert reasons["Gamma one"] == "controversial_family_definition"
    assert reasons["Delta one"] == "selected_family_validated_species_row"
    assert reasons["Outside one"] == "not_in_strict_island_universe"
    assert reasons["Wrong one"] == "family_conflict"


def test_missing_source_columns_fail_closed() -> None:
    source = pd.DataFrame([_source("Alpha one", "Passaceae")]).drop(columns=["raw_symmetry"])
    direct = pd.DataFrame(
        columns=[
            "accepted_species",
            "trait_name",
            "normalized_value",
            "source_lineages",
        ]
    )
    master = pd.DataFrame({"accepted_species": ["Alpha one"], "family": ["Passaceae"]})
    try:
        build_checkpoint(source, direct, master, min_family_audit=1)
    except ValueError as error:
        assert "raw_symmetry" in str(error)
    else:
        raise AssertionError("missing source columns must fail closed")


def test_committed_wang_checkpoint_is_exact_and_reproducible() -> None:
    checkpoint_dir = (
        REPO_ROOT
        / "data/v2/staging/traits/open_web_pilot"
        / "wang_2023_symmetry_checkpoint_20260811"
    )
    combined_dir = (
        REPO_ROOT / "data/v2/staging/traits/open_web_pilot" / "bsdb_csiro_wang_checkpoint_20260811"
    )

    manifest = json.loads(
        (checkpoint_dir / "wang_2023_symmetry_checkpoint_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = pd.read_csv(
        checkpoint_dir / "wang_2023_symmetry_evidence_20260811.csv.gz",
        low_memory=False,
    )
    audit = pd.read_csv(
        checkpoint_dir / "wang_2023_symmetry_independent_audit_20260811.csv.gz",
        low_memory=False,
    )
    families = pd.read_csv(
        checkpoint_dir / "wang_2023_symmetry_family_gate_20260811.csv",
        low_memory=False,
    )
    combined_manifest = json.loads(
        (combined_dir / "bsdb_csiro_wang_symmetry_source_package_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )

    assert len(evidence) == evidence["accepted_species"].nunique() == 28_083
    assert evidence["trait_name"].eq("floral_symmetry").all()
    assert evidence["quality"].eq("medium").all()
    assert evidence["source_lineage"].nunique() > 1
    assert len(audit) == 1_636
    assert audit["accepted_correct"].astype(str).str.lower().eq("true").sum() == 1_601
    assert families["production_approved"].sum() == 38
    assert manifest["guardrails"]["axis_only_join"] is False
    assert manifest["guardrails"]["family_inference"] is False
    assert manifest["gate"]["precision"] >= 0.95
    assert combined_manifest["candidate_evidence_rows"] == 30_200
    assert combined_manifest["selected_evidence_rows"] == 30_165
    assert combined_manifest["source_package_gate"]["package_gate_passed"] is True
