from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd

from island_v2.iospe_orchid_colour_checkpoint import (
    deterministic_audit_template,
    extract_page_candidates,
    observed_orchid_colour_states,
)
from island_v2.open_web_common import reviewed_source_package_evidence


def test_orchid_colour_states_keep_multistate_flower_markings() -> None:
    quote = (
        "The flowers have green tepals with purple veins, white petals, "
        "and a yellow lip spotted red."
    )
    assert observed_orchid_colour_states(quote) == {
        "blue_purple",
        "green_brown_inconspicuous",
        "red_pink",
        "white",
        "yellow_orange",
    }


def test_orchid_colour_states_reject_nearer_non_target_organs() -> None:
    quote = (
        "The plant has purple pseudobulbs and green leaves and carries "
        "white flowers with a yellow column and brown anther."
    )
    assert observed_orchid_colour_states(quote) == {"white"}


def test_orchid_colour_states_handle_source_typos_and_long_predicates() -> None:
    quote = (
        "The floral bracts are brown and carrying withe, pirple, carmine "
        "with terracotta striped flowers."
    )
    assert observed_orchid_colour_states(quote) == {
        "blue_purple",
        "red_pink",
        "white",
        "yellow_orange",
    }
    assert observed_orchid_colour_states(
        "The flowers become apricot with age."
    ) == {"yellow_orange"}


def test_extract_page_requires_exact_heading_and_description_text() -> None:
    html = b"""
      <html><head><title>Alpha beta</title></head><body>
      <p>Alpha beta Author 1900</p>
      <p>Common Name</p><p>The Red Alpha</p>
      <p>Found in forests and carrying white petals and a purple lip.</p>
      <p>Synonyms</p><p>Alpha gamma</p>
      </body></html>
    """
    digest = hashlib.sha256(html).hexdigest()
    rows, audit = extract_page_candidates(
        species="Alpha beta",
        source_url="https://example.org/alpha.htm?tracking=1",
        page_sha256=digest,
        html=html,
    )
    assert audit["identity_exact"] is True
    assert audit["candidate_nodes"] == 1
    assert rows[0]["normalized_value"] == "blue_purple|white"
    assert rows[0]["quality"] == "medium"
    assert rows[0]["retrieval_date"] == "2026-08-12"
    assert rows[0]["source_url"] == "https://example.org/alpha.htm"
    assert rows[0]["source_excerpt"] == (
        "Found in forests and carrying white petals and a purple lip."
    )

    cf_html = html.replace(b"Alpha beta Author", b"Alpha cf. beta Author")
    rejected, rejected_audit = extract_page_candidates(
        species="Alpha beta",
        source_url="https://example.org/cf.htm",
        page_sha256=hashlib.sha256(cf_html).hexdigest(),
        html=cf_html,
    )
    assert rejected == []
    assert rejected_audit["identity_exact"] is False


def test_deterministic_audit_covers_genera_before_fill() -> None:
    rows = []
    for genus in ("Alpha", "Beta", "Gamma"):
        for index in range(3):
            species = f"{genus} species{index}"
            rows.append(
                {
                    "accepted_species": species,
                    "trait_name": "flower_primary_color",
                    "normalized_value": "white",
                    "source_record_id": f"record:{species}",
                    "source_url": f"https://example.org/{genus}/{index}",
                    "page_title": species,
                    "source_excerpt": "The flowers are white.",
                    "content_fingerprint": f"fingerprint:{species}",
                    "source_lineage": f"lineage:{species}",
                    "source_provider": "IOSPE",
                    "source_citation": species,
                }
            )
    evidence = pd.DataFrame(rows)
    first = deterministic_audit_template(evidence, audit_size=5)
    second = deterministic_audit_template(evidence, audit_size=5)
    assert first.equals(second)
    assert len(first) == 5
    assert set(first["accepted_species"].str.split().str[0]) == {
        "Alpha",
        "Beta",
        "Gamma",
    }


def test_committed_checkpoint_passes_independent_production_gate() -> None:
    root = Path(
        "data/v2/staging/traits/open_web_pilot/"
        "iospe_orchid_colour_checkpoint_20260812"
    )
    evidence = pd.read_csv(
        root / "iospe_orchid_colour_evidence_20260812.csv.gz", dtype=str
    ).fillna("")
    audit = pd.read_csv(
        root / "iospe_orchid_colour_audit_200_20260812.csv", dtype=str
    ).fillna("")
    selected, trait_audit, summary = reviewed_source_package_evidence(
        evidence, audit
    )

    assert len(evidence) == 668
    assert len(audit) == 200
    assert summary["accepted_correct"] == 193
    assert summary["precision"] == 0.965
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"] is True
    assert trait_audit.to_dict("records") == [
        {
            "trait_name": "flower_primary_color",
            "reviewed": 200,
            "accepted_correct": 193,
            "precision": 0.965,
            "cultivar_contamination_rate": 0.0,
            "production_approved": True,
        }
    ]
    assert len(selected) == 661
    assert selected[["accepted_species", "trait_name"]].drop_duplicates().shape[0] == 631
