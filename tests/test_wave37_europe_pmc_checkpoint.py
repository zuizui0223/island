from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from island_v2.wave37_europe_pmc_checkpoint import MANUAL_RECORDS, build_checkpoint

ROOT = Path(__file__).resolve().parents[1]
DECISIONS = (
    ROOT
    / "data/v2/staging/traits/wave37_europe_pmc_reproduction"
    / "wave37_candidate_review_decisions.csv"
)
BACKBONE = (
    ROOT
    / "data/v2/staging/traits/wave37_europe_pmc_reproduction"
    / "albizia_saman_two_backbone_snapshot.json"
)

PROMOTED = {
    "EPMC-2a3a05e020ca0d28185e": (
        "Albizia saman",
        "mating_system",
        "predominantly_outcrossing",
        "PMC7042685",
        "10.1002/ece3.6005",
    ),
    "EPMC-b804b22bb08fbf7aa75a": (
        "Galanthus nivalis",
        "self_incompatibility",
        "SC",
        "PMC12903926",
        "10.3389/fpls.2026.1760796",
    ),
    "EPMC-8c81fd8f578725ab749b": (
        "Guettarda scabra",
        "autonomous_selfing_capacity",
        "autonomous",
        "PMC6137094",
        "10.1038/s41598-018-32143-5",
    ),
    "EPMC-d6c18ab50bc787a3cda6": (
        "Hedyotis acutangula",
        "self_incompatibility",
        "SC",
        "PMC9572539",
        "10.3390/plants11192640",
    ),
    "EPMC-e941ce1fb462a58cec9b": (
        "Metrosideros excelsa",
        "mating_system",
        "mixed_mating",
        "PMC5395442",
        "10.1002/ece3.2867",
    ),
    "EPMC-61d3fdc7fe33f29a7da5": (
        "Ophiorrhiza japonica",
        "self_incompatibility",
        "SC",
        "PMC9572539",
        "10.3390/plants11192640",
    ),
    "EPMC-0d799145c5b8ed8befd4": (
        "Symphonia globulifera",
        "mating_system",
        "predominantly_outcrossing",
        "PMC11069024",
        "10.1111/eva.13691",
    ),
    "EPMC-4e93d77bb66b7f4325bd": (
        "Vitex rotundifolia",
        "mating_system",
        "mixed_mating",
        "PMC10897527",
        "10.1002/ece3.10927",
    ),
}


def _candidate_fixture(path: Path) -> None:
    decisions = pd.read_csv(DECISIONS, dtype=str).fillna("")
    rows = []
    for index, candidate_id in enumerate(decisions["candidate_id"]):
        species, trait, value, pmcid, doi = PROMOTED.get(
            candidate_id,
            ("Rejected example", "self_incompatibility", "SI", f"PMC{index}", ""),
        )
        rows.append(
            {
                "candidate_id": candidate_id,
                "accepted_species": species,
                "trait_name": trait,
                "normalized_value": value,
                "pmcid": pmcid,
                "doi": doi,
                "article_title": f"Reviewed article {index}",
                "source_url": f"https://europepmc.org/articles/{pmcid}",
                "exact_supporting_quote": f"Exact supporting quote {index}.",
                "name_match_method": "accepted_name_exact",
                "source_lineage": f"doi:{doi}" if doi else f"pmcid:{pmcid}",
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)


def _queue_fixture(path: Path) -> None:
    pd.DataFrame(
        {
            "genus": [f"Genus{index}" for index in range(500)],
            "trait_name": "self_incompatibility",
            "current_support_actual": "1",
            "state_count_actual": "1",
        }
    ).to_csv(path, index=False)


def _xml_fixtures(path: Path) -> None:
    path.mkdir()
    by_pmcid: dict[str, list[str]] = {}
    for record in MANUAL_RECORDS:
        by_pmcid.setdefault(record["pmcid"], []).append(record["source_excerpt"])
    for pmcid, quotes in by_pmcid.items():
        payload = "<article><body><p>" + " </p><p>".join(quotes) + "</p></body></article>"
        with gzip.open(path / f"{pmcid}.xml.gz", "wt", encoding="utf-8") as handle:
            handle.write(payload)


def test_wave37_checkpoint_is_complete_and_fail_closed(tmp_path: Path) -> None:
    candidates = tmp_path / "candidates.csv"
    queue = tmp_path / "queue.csv"
    raw_xml = tmp_path / "raw_xml"
    master = tmp_path / "master.csv"
    output = tmp_path / "output"
    _candidate_fixture(candidates)
    _queue_fixture(queue)
    _xml_fixtures(raw_xml)
    master_species = sorted(
        {
            "Samanea saman",
            *(row[0] for row in PROMOTED.values() if row[0] != "Albizia saman"),
            *(record["accepted_species"] for record in MANUAL_RECORDS),
        }
    )
    pd.DataFrame({"accepted_species": master_species}).to_csv(master, index=False)

    summary = build_checkpoint(
        candidate_csv=candidates,
        raw_xml_dir=raw_xml,
        decisions_csv=DECISIONS,
        backbone_json=BACKBONE,
        frozen_queue_csv=queue,
        master_csv=master,
        source_run_id="test-run",
        source_artifact="test-artifact",
        output_dir=output,
        expected_species=len(master_species),
    )

    evidence = pd.read_csv(output / "wave37_reviewed_direct_evidence.csv.gz")
    audit = pd.read_csv(output / "wave37_candidate_review_audit.csv.gz")
    assert summary["machine_candidates"] == 44
    assert summary["machine_candidates_promoted"] == 8
    assert summary["manual_corrections_promoted"] == 3
    assert len(evidence) == summary["reviewed_direct_evidence_rows"] == 11
    assert len(audit) == 47
    assert "Samanea saman" in set(evidence["accepted_species"])
    assert "Albizia saman" not in set(evidence["accepted_species"])
    assert "Saussurea polylepis" in set(evidence["accepted_species"])
    assert "Saussurea was" not in set(evidence["accepted_species"])
    vitex = evidence.loc[evidence["accepted_species"].eq("Vitex doniana")]
    assert set(vitex["trait_name"]) == {
        "mating_system",
        "autonomous_selfing_capacity",
    }
    assert evidence["source_excerpt"].str.len().gt(0).all()
    assert evidence["source_run_id"].eq("test-run").all()


def test_wave37_frozen_review_covers_every_candidate_once() -> None:
    decisions = pd.read_csv(DECISIONS, dtype=str).fillna("")
    assert len(decisions) == decisions["candidate_id"].nunique() == 44
    assert int(decisions["promote_to_ledger"].str.casefold().eq("true").sum()) == 8
    assert decisions["review_decision"].value_counts().to_dict() == {
        "exclude": 20,
        "accept_existing_not_promoted": 7,
        "accept_direct": 6,
        "exclude_false_closed_conflict": 3,
        "accept_counterevidence": 2,
        "superseded_by_manual_correction": 2,
        "accept_rule_support_only_not_integrated": 1,
        "exclude_infraspecific_scope": 1,
        "superseded_by_primary_source": 1,
        "exclude_name_backbone_disagreement": 1,
    }
