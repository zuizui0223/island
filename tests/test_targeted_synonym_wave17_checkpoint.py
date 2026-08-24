import gzip
import json
from pathlib import Path

import pandas as pd

from island_v2.bsdb_two_backbone_synonym_checkpoint import (
    _mapping_reason,
    build_mapping_audit,
    read_response_records,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

CHECKPOINT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "targeted_synonym_wave17_checkpoint_20260820"
)
SOURCE_GROUP = "targeted_synonym_wave17_checkpoint_20260820"


def _read(name: str) -> pd.DataFrame:
    return pd.read_csv(CHECKPOINT / name, dtype=str).fillna("")


def test_frozen_responses_reproduce_exact_two_backbone_gate() -> None:
    response_path = CHECKPOINT / "bsdb_name_resolution_responses_20260820.jsonl.gz"
    records = read_response_records(response_path)
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")
    audit = build_mapping_audit(records, master)
    committed = _read("bsdb_two_backbone_mapping_audit_20260820.csv.gz")

    assert len(records) == len(audit) == len(committed) == 1_383
    assert audit["mapping_reason"].value_counts().to_dict() == committed[
        "mapping_reason"
    ].value_counts().to_dict()
    accepted = audit.loc[
        audit["mapping_reason"].eq("accepted_strict_two_backbone")
    ]
    assert len(accepted) == 21
    assert accepted["accepted_species"].eq(accepted["wfo_accepted_species"]).all()
    assert accepted["accepted_species"].eq(accepted["gbif_accepted_species"]).all()
    assert accepted["gbif_match_type"].eq("EXACT").all()
    assert accepted["gbif_rank"].eq("SPECIES").all()
    assert accepted["gbif_kingdom"].eq("Plantae").all()
    with gzip.open(response_path, "rt", encoding="utf-8") as handle:
        assert sum(1 for _ in handle) == 1_383


def test_mapping_gate_rejects_one_backbone_fuzzy_and_family_conflicts() -> None:
    base = {
        "wfo_status": 200,
        "gbif_status": 200,
        "source_family": "Testaceae",
        "wfo": {
            "match": {
                "placement": "Code/Plantae/Testaceae/Testus/accepted$Oldus/name",
                "wfo_id": "wfo-test",
            },
            "parsedName": {"rank": "species"},
            "params": {"fuzzyNameParts": 0},
        },
        "gbif": {
            "matchType": "EXACT",
            "rank": "SPECIES",
            "kingdom": "Plantae",
            "species": "Testus accepted",
            "family": "Testaceae",
        },
    }
    families = {"Testus accepted": "Testaceae"}
    assert _mapping_reason(base, families) == (
        "accepted_strict_two_backbone",
        "Testus accepted",
    )

    one_backbone = json.loads(json.dumps(base))
    one_backbone["wfo"]["match"] = None
    assert _mapping_reason(one_backbone, families)[0] == "wfo_no_exact_match"
    fuzzy = json.loads(json.dumps(base))
    fuzzy["wfo"]["params"]["fuzzyNameParts"] = 1
    assert _mapping_reason(fuzzy, families)[0] == "wfo_fuzzy_match_rejected"
    family_conflict = json.loads(json.dumps(base))
    family_conflict["gbif"]["family"] = "Wrongaceae"
    assert _mapping_reason(family_conflict, families)[0] == "family_conflict"


def test_wave17_has_reviewed_direct_rows_and_retains_direct_conflict() -> None:
    evidence = _read("targeted_synonym_wave17_evidence_20260820.csv")
    audit = _read("targeted_synonym_wave17_manual_audit_20260820.csv")

    assert len(evidence) == len(audit) == 12
    assert evidence[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 10
    assert evidence["source_group"].eq(SOURCE_GROUP).all()
    assert evidence["inference_rule"].eq("").all()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["source_lineage"].ne("").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert audit["decision"].eq("accept").all()
    assert audit["species_identity_correct"].eq("true").all()
    assert audit["value_correct"].eq("true").all()
    assert audit["cultivar_contamination"].eq("false").all()

    doona = evidence.loc[
        evidence["accepted_species"].eq("Doona cordifolia")
        & evidence["trait_name"].eq("self_incompatibility")
    ]
    assert len(doona) == 3
    assert set(doona["normalized_value"]) == {"SC", "SI"}
    assert doona["source_lineage"].nunique() == 3
    monopsis = evidence.loc[evidence["accepted_species"].eq("Monopsis stellarioides")]
    assert set(monopsis["trait_name"]) == {"floral_form", "flower_primary_color"}
    assert set(monopsis["evidence_quality"]) == {"medium"}


def test_wave17_combined_checkpoint_passes_common_review_gate() -> None:
    evidence = _read("combined_curated_evidence_20260820.csv")
    audit = _read("combined_curated_manual_audit_20260820.csv")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")
    accepted = validate_individually_reviewed_evidence(
        evidence, audit, master, Path("config/trait_ontology.yml")
    )

    assert len(evidence) == len(audit) == len(accepted) == 2_897
    assert len(accepted.loc[accepted["source_group"].eq(SOURCE_GROUP)]) == 12


def test_wave17_manifest_records_real_query_cost_and_fail_closed_rules() -> None:
    manifest = json.loads(
        (CHECKPOINT / "source_acquisition_manifest_wave17.json").read_text(
            encoding="utf-8"
        )
    )

    assert manifest["accepted_evidence_rows"] == 12
    assert manifest["accepted_species_trait"] == 10
    assert manifest["bsdb_acquisition"]["unresolved_unique_names_queried"] == 1_383
    assert manifest["bsdb_acquisition"]["total_web_api_requests"] == 2_766
    assert manifest["bsdb_acquisition"]["strict_two_backbone_agreements"] == 21
    assert manifest["bsdb_acquisition"]["new_direct_rows"] == 10
    assert manifest["bsdb_acquisition"]["direct_conflict_species_trait"] == 1
    assert manifest["regional_flora_resume"]["new_species_pages_fetched"] == 333
    assert manifest["search_cost_usd"] == 0.0
    assert manifest["guardrails"]["fuzzy_or_one_backbone_match"] is False
    assert manifest["guardrails"]["direct_conflict_forced_resolved"] is False
    assert manifest["guardrails"]["genus_axis_only_join"] is False
