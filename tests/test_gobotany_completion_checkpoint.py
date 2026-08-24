import json
from pathlib import Path

import pandas as pd

from island_v2.gobotany_completion_checkpoint import (
    _sample_audit,
    finalize_audit,
    parse_page,
    resolved_direct_pairs,
    species_links,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS


def page_html(*, species: str = "Solanum example", family: str = "Solanaceae") -> bytes:
    return f"""
    <html><head><title>{species} (test): Go Botany</title></head><body>
      <h1>{species} — test plant</h1>
      <dl><dt>Cleistogamous flowers</dt><dd>
        <ul><li>the plant has some cleistogamous flower</li>
        <li>there are no cleistogamous flowers on the plan</li></ul>
      </dd></dl>
      <dl><dt>Flower petal color</dt><dd><ul><li>blue to purple</li><li>white</li></ul></dd></dl>
      <dl><dt>Flower symmetry</dt><dd>
        there is only one way to evenly divide the flower (the flower is bilaterally symmetrical)
      </dd></dl>
      <div class="section"><h2>Family</h2><p><a href="/family/test/">{family}</a></p></div>
    </body></html>
    """.encode()


def test_species_links_accepts_only_binomial_species_pages() -> None:
    html = """
    <a href="/species/solanum/example/">Solanum example</a>
    <a href="/species/solanum/example/">not a binomial title</a>
    <a href="https://example.test/species/solanum/wrong/">Solanum wrong</a>
    """
    assert species_links(html) == {
        "Solanum example": "https://gobotany.nativeplanttrust.org/species/solanum/example/"
    }


def test_resume_excludes_only_resolved_direct_pairs() -> None:
    ledger = pd.DataFrame(
        [
            {
                "accepted_species": "Solanum resolved",
                "trait_name": "cleistogamy",
                "resolution_status": "resolved",
            },
            {
                "accepted_species": "Solanum conflict",
                "trait_name": "cleistogamy",
                "resolution_status": "unresolved",
            },
        ]
    )
    assert resolved_direct_pairs(ledger) == {("Solanum resolved", "cleistogamy")}


def test_parse_page_keeps_trait_specific_multistate_values() -> None:
    audit, rows = parse_page(
        expected_species="Solanum example",
        expected_family="Solanaceae",
        source_url="https://gobotany.nativeplanttrust.org/species/solanum/example/",
        payload=page_html(),
    )
    assert audit.status == "accepted_identity"
    by_trait = {row["trait_name"]: row for row in rows}
    assert by_trait["cleistogamy"]["normalized_value"] == "absent|facultative"
    assert by_trait["flower_primary_color"]["normalized_value"] == "blue_purple|white"
    assert by_trait["floral_symmetry"]["normalized_value"] == "zygomorphic"
    assert all(row["evidence_scope"] == "species_direct" for row in rows)
    assert all(row["quality"] == "medium" for row in rows)
    assert all(row["source_excerpt"] for row in rows)
    assert all(row["source_lineage"].startswith("url:https://") for row in rows)


def test_parse_page_fails_closed_on_name_and_family() -> None:
    name_audit, name_rows = parse_page(
        expected_species="Solanum target",
        expected_family="Solanaceae",
        source_url="https://gobotany.nativeplanttrust.org/species/solanum/target/",
        payload=page_html(species="Solanum other"),
    )
    assert name_audit.status == "name_mismatch"
    assert name_rows == []

    family_audit, family_rows = parse_page(
        expected_species="Solanum example",
        expected_family="Solanaceae",
        source_url="https://gobotany.nativeplanttrust.org/species/solanum/example/",
        payload=page_html(family="Convolvulaceae"),
    )
    assert family_audit.status == "family_mismatch"
    assert family_rows == []


def test_deterministic_audit_and_review_provenance(tmp_path: Path) -> None:
    rows = []
    for trait, axis, count in (
        ("cleistogamy", "reproductive_assurance", 120),
        ("flower_primary_color", "flower_colour", 60),
        ("floral_symmetry", "floral_structural_complexity", 60),
    ):
        for index in range(count):
            row = {column: "x" for column in EVIDENCE_COLUMNS}
            row.update(
                {
                    "accepted_species": f"Testus {trait.replace('_', '')}{index}",
                    "axis": axis,
                    "trait_name": trait,
                    "normalized_value": "absent" if trait == "cleistogamy" else "white",
                    "quality": "medium",
                    "source_url": f"https://example.test/{trait}/{index}",
                    "source_record_id": f"{trait}-{index}",
                    "source_excerpt": f"{trait}: evidence {index}",
                    "page_title": f"title {index}",
                    "content_fingerprint": f"fingerprint-{trait}-{index}",
                }
            )
            rows.append(row)
    evidence = pd.DataFrame(rows)
    first = _sample_audit(evidence)
    second = _sample_audit(evidence.sample(frac=1, random_state=4))
    assert len(first) == 200
    assert first["source_url"].nunique() == 200
    assert set(first["candidate_id"]) == set(second["candidate_id"])
    assert first["trait_name"].value_counts().to_dict() == {
        "cleistogamy": 100,
        "flower_primary_color": 50,
        "floral_symmetry": 50,
    }

    candidates = tmp_path / "audit_candidates.csv"
    reviewed = tmp_path / "audit_reviewed.csv"
    manifest = tmp_path / "manifest.json"
    first.to_csv(candidates, index=False)
    manifest.write_text('{"files": {}}', encoding="utf-8")
    summary = finalize_audit(
        candidates_csv=candidates,
        output_csv=reviewed,
        reviewer="Codex evidence audit",
        reviewed_at_utc="2026-08-11T03:00:00Z",
        manifest_json=manifest,
    )
    assert summary["reviewed"] == 200
    assert summary["precision"] == 1.0
    result = pd.read_csv(reviewed, dtype=str).fillna("")
    assert result["accepted_correct"].eq("True").all()
    assert result["reviewer"].eq("Codex evidence audit").all()
    assert result["reviewed_at_utc"].eq("2026-08-11T03:00:00Z").all()
    manifest_payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert manifest_payload["completed_review"]["reviewed"] == 200
    assert manifest_payload["files"][reviewed.name]["size_bytes"] == reviewed.stat().st_size
