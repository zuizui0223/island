from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

from analysis.acquire_pngtrees_source_package_20260813 import (
    SOURCE_LINEAGE,
    build,
    extract_traits,
)
from island_v2.open_web_common import reviewed_source_package_evidence

ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "pngtrees_source_package_20260813"
)


def test_structured_fields_remain_trait_specific_and_multistate() -> None:
    excerpt = (
        "Inflorescence terminal, flowers single or flowers on a branched axis; "
        "flowers with many planes of symmetry or with one plane of symmetry, "
        "diameter small (up to10 mm diam.); inner perianth white, pink, or mauve; "
        "5, free; stamens 10."
    )

    assert extract_traits(excerpt) == [
        (
            "floral_symmetry",
            "actinomorphic|zygomorphic",
            "explicit symmetry field",
        ),
        ("flower_size_class", "small", "diameter small (up to 10 mm)"),
        (
            "flower_primary_color",
            "blue_purple|red_pink|white",
            "inner perianth white, pink, or mauve",
        ),
        (
            "inflorescence_display",
            "raceme_spike_panicle|solitary",
            "explicit flowers arrangement field",
        ),
    ]

    syconium = (
        "Inflorescence axillary, flowers single (syconium/fig solitary), "
        "flowers with many planes of symmetry (syconium/fig), with one plane "
        "of symmetry (male and female flowers), 2 mm long, diameter small "
        "(up to10 mm diam.); inner perianth red; 3, free; stamens 1."
    )
    assert extract_traits(syconium) == [
        ("flower_size_class", "small", "diameter small (up to 10 mm)"),
        ("flower_primary_color", "red_pink", "inner perianth red"),
    ]


def test_committed_pngtrees_package_passes_the_shared_review_gate() -> None:
    acquisition = json.loads(
        (ROOT / "pngtrees_acquisition_summary_20260813.json").read_text(
            encoding="utf-8"
        )
    )
    evidence = pd.read_csv(
        ROOT / "pngtrees_source_package_evidence_20260813.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "pngtrees_source_package_audit_200_20260813.csv", dtype=str
    ).fillna("")
    pages = pd.read_csv(
        ROOT / "pngtrees_page_manifest_20260813.csv.gz", dtype=str
    ).fillna("")

    selected, scopes, gate = reviewed_source_package_evidence(evidence, audit)

    assert len(pages) == 312
    assert acquisition["credential_free"] is True
    assert acquisition["http_requests"] == 313
    assert acquisition["search_api_queries"] == 0
    assert acquisition["search_cost_usd"] == 0.0
    assert acquisition["candidate_rows"] == 1_143
    assert pages["fetch_status"].value_counts().to_dict() == {
        "accepted_page": 301,
        "family_rejected": 11,
    }
    assert len(evidence) == len(selected) == 1_143
    assert evidence["source_record_id"].is_unique
    assert evidence[["accepted_species", "trait_name"]].duplicated().sum() == 0
    assert evidence["accepted_species"].nunique() == 301
    assert evidence["source_lineage"].eq(SOURCE_LINEAGE).all()
    assert evidence["quality"].eq("medium").all()
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert set(evidence["trait_name"]) == {
        "floral_symmetry",
        "flower_primary_color",
        "flower_size_class",
        "inflorescence_display",
    }
    assert not set(evidence["trait_name"]).intersection(
        {"pollen_vector_mode", "reward_type"}
    )

    assert len(audit) == 200
    assert audit["accepted_species"].nunique() == 200
    assert audit["trait_name"].value_counts().to_dict() == {
        "floral_symmetry": 50,
        "flower_primary_color": 50,
        "flower_size_class": 50,
        "inflorescence_display": 50,
    }
    assert gate["package_gate_passed"] is True
    assert gate["precision"] == 1.0
    assert gate["cultivar_contamination_rate"] == 0.0
    assert scopes["production_approved"].all()


def test_family_drift_is_audited_not_silently_reconciled() -> None:
    pages = pd.read_csv(
        ROOT / "pngtrees_page_manifest_20260813.csv.gz", dtype=str
    ).fillna("")
    rejected = pages.loc[pages["fetch_status"].eq("family_rejected")]

    assert len(rejected) == 11
    assert rejected["expected_family"].ne(rejected["displayed_family"]).all()
    assert set(rejected["accepted_species"]).issuperset(
        {"Calophyllum inophyllum", "Pangium edule", "Platea latifolia"}
    )


def test_pngtrees_package_rebuild_is_byte_reproducible(tmp_path: Path) -> None:
    manifest = json.loads(
        (ROOT / "pngtrees_source_package_manifest_20260813.json").read_text(
            encoding="utf-8"
        )
    )
    for filename in manifest["input_hashes"]:
        shutil.copy2(ROOT / filename, tmp_path / filename)

    rebuilt = build(tmp_path)

    assert rebuilt["input_hashes"] == manifest["input_hashes"]
    assert rebuilt["output_hashes"] == manifest["output_hashes"]
    assert rebuilt["review"] == manifest["review"]
    assert rebuilt["selected"] == manifest["selected"]
