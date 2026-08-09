import hashlib
from pathlib import Path

import pandas as pd

from analysis.acquire_plantnet_nsw_traits_20260809 import (
    _source_lineage,
    build_evidence,
    deterministic_audit_queue,
    parse_profile,
)
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS
from island_v2.open_web_common import reviewed_source_package_evidence
from island_v2.trait_measurements import load_rules

COMMITTED_PACKAGE = (
    Path("data/v2/staging/traits/direct_llm_pilot")
    / "20260809_plantnet_nsw_source_acquisition"
)


def _rules():
    return load_rules(Path("config/measurement_classification.yml"))


def test_profile_parser_limits_text_to_official_description() -> None:
    html = """
    <html><body><table><tr>
      <td class="head2"><i>Abelmoschus manihot</i> (L.) Medik.</td>
      <td align="right">Family <b><a>Malvaceae</a></b></td>
    </tr></table>
    <p><b>Description:</b> Herb. Flowers solitary or in racemes.
    <p>Petals 5-7 cm long, white or yellow with a purple spot.
    <p><b>Distribution and occurrence:</b> Recorded from NSW.
    <div>Text by A. Botanist<br>Taxon concept: current</div>
    </body></html>
    """

    parsed = parse_profile(html)

    assert parsed["profile_scientific_name"] == "Abelmoschus manihot"
    assert parsed["profile_family"] == "Malvaceae"
    assert "Flowers solitary or in racemes" in parsed["description"]
    assert "Petals 5-7 cm long" in parsed["description"]
    assert "Distribution" not in parsed["description"]
    assert parsed["source_attribution"] == "Text by A. Botanist"


def test_exact_identity_and_family_gate_produces_six_trait_ledger_rows() -> None:
    inventory = pd.DataFrame(
        [
            {
                "plantnet_species": "Testa floribunda",
                "profile_scientific_name": "Testa floribunda",
                "profile_family": "Testaceae",
                "description": (
                    "Flowers solitary, zygomorphic and tubular, 12 mm long, red; "
                    "corolla tube 8 mm long. Inflorescences also in racemes."
                ),
                "source_attribution": "Text by A. Botanist",
                "profile_url": "https://example.test/profile",
                "content_sha256": "abc",
                "fetch_status": "cached",
            }
        ]
    )
    master = pd.DataFrame(
        [{"accepted_species": "Testa floribunda", "family": "Testaceae"}]
    )

    evidence, treatments = build_evidence(inventory, master, set(), _rules())

    assert list(evidence.columns) == EVIDENCE_COLUMNS
    assert set(evidence["trait_name"]) == {
        "flower_primary_color",
        "floral_form",
        "floral_symmetry",
        "flower_size_class",
        "tube_depth_class",
        "inflorescence_display",
    }
    assert evidence["quality"].eq("high").all()
    assert evidence["source_excerpt"].ne("").all()
    assert treatments.loc[0, "identity_gate_passed"]


def test_lineage_uses_citation_and_normalized_excerpt_fingerprint() -> None:
    first = _source_lineage(
        "PlantNET NSW; Text by A. Botanist",
        "Flowers white.",
        "Testa alba",
    )
    second = _source_lineage(
        "  plantnet NSW TEXT BY A botanist ",
        " flowers WHITE ",
        "Othera alba",
    )

    assert first == second
    assert first[1] == "citation_plus_normalized_excerpt_fingerprint"


def test_audit_queue_is_trait_stratified_and_deterministic() -> None:
    rows = []
    traits = {
        "flower_primary_color": "flower_colour",
        "floral_form": "floral_structural_complexity",
        "floral_symmetry": "floral_structural_complexity",
        "flower_size_class": "floral_structural_complexity",
        "tube_depth_class": "floral_structural_complexity",
        "inflorescence_display": "floral_structural_complexity",
    }
    for trait, axis in traits.items():
        for index in range(40):
            row = {column: "x" for column in EVIDENCE_COLUMNS}
            row.update(
                {
                    "accepted_species": f"Genus species{index}",
                    "axis": axis,
                    "trait_name": trait,
                    "source_record_id": f"{trait}-{index}",
                }
            )
            rows.append(row)
    evidence = pd.DataFrame(rows)

    first = deterministic_audit_queue(evidence)
    second = deterministic_audit_queue(evidence.sample(frac=1, random_state=9))

    assert len(first) == 200
    assert first["candidate_id"].tolist() == second["candidate_id"].tolist()
    assert first.groupby("trait_name").size().min() >= 20


def test_committed_source_package_approves_only_well_audited_traits() -> None:
    evidence = pd.read_csv(
        COMMITTED_PACKAGE / "plantnet_nsw_evidence.csv.gz",
        dtype=str,
    ).fillna("")
    audit = pd.read_csv(
        COMMITTED_PACKAGE / "plantnet_nsw_manual_audit_200.csv",
        dtype=str,
    ).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 200
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["candidate_evidence_rows"] == 758
    assert summary["selected_evidence_rows"] == 745
    assert summary["approved_traits"] == [
        "flower_size_class",
        "inflorescence_display",
        "tube_depth_class",
    ]
    assert len(selected) == 745
    assert not scopes.loc[
        scopes["trait_name"].isin(
            {"floral_form", "floral_symmetry", "flower_primary_color"}
        ),
        "production_approved",
    ].any()


def test_committed_source_package_manifest_is_exact() -> None:
    entries = {}
    for line in (COMMITTED_PACKAGE / "file_manifest.sha256").read_text(
        encoding="ascii"
    ).splitlines():
        digest, filename = line.split("  ", maxsplit=1)
        entries[filename] = digest

    files = {
        path.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in COMMITTED_PACKAGE.iterdir()
        if path.is_file() and path.name != "file_manifest.sha256"
    }
    assert entries == files
