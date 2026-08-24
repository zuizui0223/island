import hashlib
import json
from pathlib import Path

import pandas as pd

from island_v2.florml_colour_recovery_checkpoint import high_precision_candidate
from island_v2.open_web_common import reviewed_source_package_evidence


def _sha256(path: Path) -> str:
    payload = path.read_bytes()
    if path.suffix.casefold() in {".csv", ".json"}:
        payload = payload.replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


def test_florml_colour_gate_rejects_coloured_coverings() -> None:
    base = {
        "accepted_species": "Testus exactus",
        "trait_name": "flower_primary_color",
        "normalized_value": "yellow_orange",
        "source_provider": "flora_malesiana",
        "evidence_scope": "species_direct",
        "name_match_method": "wfo_june_2026_exact_accepted_xml_species_family_consistent",
    }
    accepted, _ = high_precision_candidate(
        {**base, "source_excerpt": "Flowers yellow; petals glabrous."}
    )
    covered, reason = high_precision_candidate(
        {**base, "source_excerpt": "Flowers with yellow puberulous hairs."}
    )
    assert accepted
    assert not covered
    assert reason == "coloured_hair_scale_or_pubescence_context"


def test_committed_florml_package_passes_shared_review_gate() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "combined_florml_source_package_20260811"
    )
    evidence_path = root / "combined_source_package_with_florml_20260811_evidence.csv.gz"
    audit_path = root / "combined_source_package_with_florml_20260811_audit.csv.gz"
    manifest_path = root / "combined_source_package_with_florml_20260811_manifest.json"
    evidence = pd.read_csv(evidence_path, dtype=str).fillna("")
    audit = pd.read_csv(audit_path, dtype=str).fillna("")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert len(evidence) == len(selected) == 1_524
    assert len(audit) == 800
    assert summary["package_gate_passed"]
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert scopes["production_approved"].all()
    assert _sha256(evidence_path) == manifest["output_hashes"][evidence_path.name]
    assert _sha256(audit_path) == manifest["output_hashes"][audit_path.name]
