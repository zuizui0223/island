from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.nparks_near_rule_morphology_checkpoint import (
    _canonical_file_size,
    normalize_field,
    page_record_kind,
)
from island_v2.open_web_finalize import validate_individually_reviewed_evidence

ROOT = Path(
    "data/v2/staging/traits/open_web_pilot/"
    "nparks_near_rule_checkpoint_20260811"
)


def _sha256(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.lower() in {".csv", ".json"}:
        payload = payload.replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


def test_species_gate_rejects_forms_but_keeps_botanical_authorship() -> None:
    assert page_record_kind("Begonia rex-cultorum", "Begonia") == (
        "cultivar_hybrid_or_infraspecific"
    )
    assert page_record_kind("Pinanga coronata (Blume) Blume", "Pinanga") == (
        "species_binomial"
    )
    assert page_record_kind("Nepenthes ampullaria (Brunei Red)", "Nepenthes") == (
        "cultivar_hybrid_or_infraspecific"
    )
    assert page_record_kind("Heliconia indica 'Spectabilis'", "Heliconia") == (
        "cultivar_hybrid_or_infraspecific"
    )
    assert page_record_kind("Heliconia bihai × caribaea", "Heliconia") == (
        "cultivar_hybrid_or_infraspecific"
    )


def test_field_mapping_preserves_states_and_fails_closed() -> None:
    mapping = {
        "red": "red_pink",
        "yellow / golden": "yellow_orange",
    }
    assert normalize_field("Red, Yellow / Golden", mapping) == (
        "red_pink|yellow_orange",
        [],
    )
    assert normalize_field("Red, chartreuse", mapping) == ("", ["chartreuse"])


def test_committed_checkpoint_is_complete_reviewed_and_hash_pinned() -> None:
    manifest = json.loads(
        (ROOT / "nparks_near_rule_checkpoint_manifest_20260811.json").read_text(
            encoding="utf-8"
        )
    )
    inventory = pd.read_csv(
        ROOT / "nparks_near_rule_inventory_20260811.csv.gz", dtype=str
    ).fillna("")
    evidence = pd.read_csv(
        ROOT / "nparks_near_rule_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "nparks_near_rule_manual_audit_20260811.csv", dtype=str
    ).fillna("")

    assert manifest["inventory_rows"] == len(inventory) == 501
    assert manifest["evidence_rows"] == len(evidence) == len(audit) == 82
    assert manifest["species"] == evidence["accepted_species"].nunique() == 54
    assert manifest["species_trait"] == 82
    assert manifest["trait_counts"] == {
        "floral_form": 8,
        "floral_symmetry": 20,
        "flower_primary_color": 37,
        "inflorescence_display": 17,
    }
    assert manifest["inventory_status_counts"]["rejected_non_species_record"] == 89
    assert evidence["candidate_id"].is_unique
    assert not evidence.duplicated(["accepted_species", "trait_name"]).any()
    assert evidence["source_excerpt"].ne("").all()
    assert evidence["content_sha256"].str.fullmatch(r"[0-9a-f]{64}").all()
    assert evidence["name_match_method"].eq("accepted_name_exact").all()
    assert evidence["evidence_quality"].eq("high").all()
    assert evidence["domain"].eq("nparks.gov.sg").all()
    assert not set(evidence["trait_name"]).intersection(
        {"reward_type", "pollen_vector_mode"}
    )
    assert not evidence["page_title"].str.contains(
        r"Brunei Red|'|Dwarf|Variegated|Hybrid", case=False, regex=True
    ).any()
    assert audit["decision"].eq("accept").all()
    assert audit["cultivar_contamination"].eq("false").all()

    for filename, metadata in manifest["files"].items():
        path = ROOT / filename
        assert _canonical_file_size(path) == metadata["size_bytes"]
        assert _sha256(path) == metadata["sha256"]


def test_combined_curated_checkpoint_passes_the_formal_gate() -> None:
    evidence = pd.read_csv(
        ROOT / "combined_curated_evidence_20260811.csv", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        ROOT / "combined_curated_manual_audit_20260811.csv", dtype=str
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv", dtype=str
    ).fillna("")

    accepted = validate_individually_reviewed_evidence(
        evidence,
        audit,
        master,
        Path("config/trait_ontology.yml"),
    )

    assert len(evidence) == len(audit) == len(accepted) == 170
    assert evidence["candidate_id"].is_unique
    assert audit["candidate_id"].is_unique
