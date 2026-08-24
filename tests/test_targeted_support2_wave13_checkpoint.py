import json
from pathlib import Path

import pandas as pd

from island_v2.open_web_finalize import (
    validate_individually_reviewed_evidence,
)

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_support2_wave13_checkpoint_20260814"
)
SOURCE_GROUP = "targeted_support2_wave13_checkpoint_20260814"


def test_wave13_rows_are_species_direct_and_trait_specific() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "targeted_support2_wave13_evidence_20260814.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "targeted_support2_wave13_manual_audit_20260814.csv", dtype=str
    ).fillna("")

    assert len(evidence) == len(audit) == 2
    assert set(evidence["accepted_species"]) == {
        "Anthemis rigida",
        "Decaspermum fruticosum",
    }
    assert set(evidence["trait_name"]) == {
        "flower_size_class",
        "floral_symmetry",
    }
    assert evidence["evidence_scope"].eq("species_direct").all()
    assert evidence["evidence_quality"].eq("medium").all()
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_url"].str.startswith("https://").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].str.casefold().eq("true").all()
    assert audit["value_correct"].str.casefold().eq("true").all()


def test_wave13_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = pd.read_csv(
        CHECKPOINT / "combined_curated_evidence_20260814.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        CHECKPOINT / "combined_curated_manual_audit_20260814.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )
    assert len(evidence) == len(audit) == len(accepted) == 2_751
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 2


def test_wave13_manifest_records_fail_closed_rejections() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave13.json").read_text(
            encoding="utf-8"
        )
    )
    assert manifest["manual_open_web_queries"] == 48
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["accepted_species_trait_rows"] == 2
    assert len(manifest["withheld_or_rejected"]) == 3
    assert manifest["guardrails"]["search_snippet_as_evidence"] is False
    assert manifest["guardrails"]["cross_trait_substitution"] is False
