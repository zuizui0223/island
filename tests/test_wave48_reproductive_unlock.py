from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

import island_v2.wave48_reproductive_checkpoint as checkpoint_module
import island_v2.wave48_reproductive_overlay as overlay_module
from island_v2.all_evidence_trait_audit import build_rule_audit
from island_v2.wave48_incremental_all_evidence import build_incremental_audit

AXES = (
    "floral_structural_complexity",
    "flower_colour",
    "reproductive_assurance",
)


def _write(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, index=False, compression="gzip" if path.suffix == ".gz" else None)


def _coverage(species: tuple[str, ...]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": accepted_species,
                "axis": axis,
                "trait_composition": "",
                "trait_names": "",
                "source_groups": "",
                "source_lineages": "",
                "quality": "",
            }
            for accepted_species in species
            for axis in AXES
        ]
    )


def _evidence(
    species: str,
    value: str,
    lineage: str,
    *,
    scope: str = "species_direct",
    record_id: str | None = None,
) -> dict[str, str]:
    external = scope == "external_congener_species_direct"
    return {
        "accepted_species": species,
        "axis": "reproductive_assurance",
        "trait_name": "self_incompatibility",
        "normalized_value": value,
        "quality": "high",
        "source_group": "primary",
        "source_provider": "Original source",
        "source_url": "https://example.test/source",
        "source_record_id": record_id or f"record:{species}:{value}",
        "source_citation": "Original source 2026",
        "source_excerpt": f"{species} is {value}.",
        "evidence_scope": scope,
        "name_match_method": "strict_wfo_gbif_two_backbone",
        "source_lineage": lineage,
        "lineage_method": "original_article_bibliographic_identity",
        "source_run_id": "web:retrieved-20260831",
        "source_artifact": "test-source",
        "source_file": "source.html",
        "acceptance_contract": (
            "external_congener_species_direct_strict_two_backbone_v1"
            if external
            else "species_direct_original_source_v1"
        ),
    }


def _resolved_cells(rows: list[dict[str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": row["accepted_species"],
                "axis": row["axis"],
                "trait_name": row["trait_name"],
                "classification": "single_independent_lineage",
                "resolution_status": "resolved",
                "quality": row["quality"],
                "state_set": json.dumps([row["normalized_value"]]),
                "normalized_value": row["normalized_value"],
                "source_groups": row["source_group"],
                "source_lineages": row["source_lineage"],
                "genus": row["accepted_species"].split()[0],
            }
            for row in rows
        ]
    )


def _lineages(cells: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "accepted_species": row.accepted_species,
                "genus": row.genus,
                "axis": row.axis,
                "trait_name": row.trait_name,
                "state_set": row.state_set,
                "source_lineage": row.source_lineages,
                "source_group": row.source_groups,
                "ontology_valid": True,
                "lineage_internal_conflict": False,
            }
            for row in cells.itertuples(index=False)
        ]
    )


def test_wave48_checkpoint_verifies_source_quote_and_two_backbone_identity(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(checkpoint_module, "DIRECT_ROWS", 1)
    monkeypatch.setattr(checkpoint_module, "EXTERNAL_ROWS", 1)
    monkeypatch.setattr(checkpoint_module, "REVIEW_ROWS", 2)
    monkeypatch.setattr(checkpoint_module, "IDENTITY_ROWS", 2)
    monkeypatch.setattr(checkpoint_module, "REJECTED_ROWS", 1)
    species = ("Pelargonium alpha", "Pelargonium beta")
    target = tmp_path / "target.csv.gz"
    _write(_coverage(species), target)

    direct = _evidence("Pelargonium alpha", "SC", "doi:test")
    external = _evidence(
        "Pelargonium gamma",
        "SC",
        "doi:test",
        scope="external_congener_species_direct",
    )
    quote = "Pelargonium alpha and Pelargonium gamma are SC."
    fingerprint = hashlib.sha256(quote.encode()).hexdigest()
    direct["source_excerpt"] = quote
    direct["content_fingerprint"] = fingerprint
    external["source_excerpt"] = quote
    external["content_fingerprint"] = fingerprint

    packet = tmp_path / "packet"
    packet.mkdir()
    pd.DataFrame([direct]).to_csv(
        packet / "wave48_reviewed_direct_evidence.csv", index=False
    )
    pd.DataFrame([external]).to_csv(
        packet / "wave48_external_congener_evidence.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "record_id": row["source_record_id"],
                "exact_supporting_quote": quote,
                "accepted_correct": "true",
                "source_lineage": row["source_lineage"],
                "content_fingerprint": fingerprint,
            }
            for row in (direct, external)
        ]
    ).to_csv(packet / "wave48_source_review_audit.csv", index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": name,
                "target_universe_status": status,
                "family": "Geraniaceae",
                "name_match_method": "strict_wfo_gbif_two_backbone",
                "wfo_match_id": f"wfo-{index:010d}",
                "wfo_status": "accepted",
                "wfo_rank": "species",
                "gbif_usage_key": str(index),
                "gbif_status": "ACCEPTED",
                "gbif_rank": "SPECIES",
                "gbif_match_type": "EXACT",
                "gbif_confidence": "99",
                "gbif_family": "Geraniaceae",
            }
            for index, (name, status) in enumerate(
                (
                    ("Pelargonium alpha", "fixed_target"),
                    ("Pelargonium gamma", "external_congener_only"),
                ),
                start=1,
            )
        ]
    ).to_csv(packet / "wave48_identity_audit.csv", index=False)
    pd.DataFrame([{"candidate": "unsafe"}]).to_csv(
        packet / "wave48_rejected_candidates.csv", index=False
    )
    source_dir = tmp_path / "sources"
    source_dir.mkdir()
    (source_dir / "source.html").write_text(f"<p>{quote}</p>", encoding="utf-8")
    (packet / "source_manifest.json").write_text(
        json.dumps(
            {
                "fixed_target_species": 2,
                "formal_wave33_baseline": {"run_id": 32932103226},
                "immediate_formal_baseline": {"run_id": 33380486845},
                "evidence_counts": {
                    "direct_rows": 1,
                    "external_rows": 1,
                    "review_rows": 2,
                    "identity_rows": 2,
                    "rejected_rows": 1,
                },
                "selection": {"source_pages_retrieved": 1},
                "sources": [
                    {
                        "source_id": "test",
                        "retrieved_filename": "source.html",
                        "source_lineage": "doi:test",
                        "content_anchor_sha256": fingerprint,
                        "reviewed_rows": 2,
                    }
                ],
                "inference_policy": {
                    "join_key": "genus x axis x trait_name",
                    "minimum_species": 3,
                    "family_inference": False,
                    "global_fallback": False,
                    "reproductive_traits_interchangeable": False,
                    "single_source_lineage_can_unlock_rule": False,
                },
            }
        ),
        encoding="utf-8",
    )

    summary = checkpoint_module.validate_packet(
        packet_dir=packet,
        target_coverage_csv=target,
        retrieved_source_dir=source_dir,
        output_dir=tmp_path / "output",
        output_json=tmp_path / "output" / "summary.json",
        expected_species=2,
    )
    assert summary["evidence"]["new_direct_species_axis"] == 1
    assert summary["checks"]["retrieved_sources_verified"] is True
    assert summary["checks"]["content_fingerprints_verified"] is True


def test_wave48_incremental_rebuild_requires_trait_and_independent_lineages(
    tmp_path: Path,
) -> None:
    species = (
        "Pelargonium alpha",
        "Pelargonium beta",
        "Pelargonium gamma",
        "Callicarpa alpha",
    )
    master = pd.DataFrame(
        [
            {
                "accepted_species": name,
                "genus": name.split()[0],
                "family": "Testaceae",
                "n_islands": "1",
                "n_records": "1",
            }
            for name in species
        ]
    )
    master_path = tmp_path / "master.csv"
    master.to_csv(master_path, index=False)
    baseline_path = tmp_path / "baseline.csv.gz"
    _write(_coverage(species), baseline_path)

    prior_direct_rows = [_evidence("Pelargonium alpha", "SC", "doi:p-alpha")]
    prior_external_rows = [
        _evidence(
            "Pelargonium delta",
            "SC",
            "doi:p-delta",
            scope="external_congener_species_direct",
        ),
        _evidence(
            "Callicarpa beta",
            "SI",
            "doi:c-beta",
            scope="external_congener_species_direct",
        ),
        _evidence(
            "Callicarpa gamma",
            "SI",
            "doi:c-gamma",
            scope="external_congener_species_direct",
        ),
    ]
    prior_direct = _resolved_cells(prior_direct_rows)
    prior_external = _resolved_cells(prior_external_rows)
    old_low = pd.DataFrame(
        columns=["accepted_species", "genus", "axis", "trait_name", "state_set"]
    )
    previous_rules = build_rule_audit(
        pd.concat([prior_direct, prior_external], ignore_index=True),
        _lineages(pd.concat([prior_direct, prior_external], ignore_index=True)),
        old_low,
    )
    direct_path = tmp_path / "direct.csv.gz"
    external_path = tmp_path / "external.csv.gz"
    rules_path = tmp_path / "rules.csv.gz"
    low_path = tmp_path / "low.csv.gz"
    conflicts_path = tmp_path / "conflicts.csv.gz"
    _write(prior_direct, direct_path)
    _write(prior_external, external_path)
    _write(previous_rules, rules_path)
    _write(old_low, low_path)
    _write(
        pd.DataFrame(
            columns=[
                "accepted_species",
                "axis",
                "trait_name",
                "classification",
                "resolution_status",
            ]
        ),
        conflicts_path,
    )

    new_direct = pd.DataFrame(
        [
            _evidence("Pelargonium beta", "SC", "doi:p-beta"),
            _evidence("Callicarpa alpha", "SC", "doi:c-alpha"),
        ]
    )
    new_external = pd.DataFrame(
        [
            _evidence(
                "Pelargonium epsilon",
                "SC",
                "doi:p-epsilon",
                scope="external_congener_species_direct",
            )
        ]
    )
    new_direct_path = tmp_path / "new-direct.csv.gz"
    new_external_path = tmp_path / "new-external.csv.gz"
    _write(new_direct, new_direct_path)
    _write(new_external, new_external_path)

    summary = build_incremental_audit(
        master_csv=master_path,
        ontology_yaml=Path(__file__).resolve().parents[1] / "config" / "trait_ontology.yml",
        baseline_coverage_csv=baseline_path,
        previous_rule_audit_csv=rules_path,
        previous_resolved_direct_csv=direct_path,
        previous_external_resolved_csv=external_path,
        previous_external_conflicts_csv=conflicts_path,
        previous_rebuilt_low_csv=low_path,
        new_direct_evidence_csv=new_direct_path,
        new_external_evidence_csv=new_external_path,
        output_dir=tmp_path / "audit",
        expected_species=4,
        expected_direct_rows=2,
        expected_external_rows=1,
        expected_blocked_rules=frozenset(),
    )
    assert summary["new_eligible_rules"] == ["Pelargonium x self_incompatibility"]
    assert summary["new_rule"] == {
        "n_direct_species": 4,
        "dominance": 1.0,
        "species_loo_accuracy": 1.0,
        "lineage_loo_accuracy": 1.0,
        "inferred_value": "SC",
    }
    assert summary["counterexample_blocked_rules"] == [
        "Callicarpa x self_incompatibility"
    ]
    assert summary["new_validated_low_species_trait"] == 1


def test_wave48_overlay_relabels_outputs_and_pins_runs(tmp_path: Path, monkeypatch) -> None:
    old_names = (
        "wave45_species_axis_coverage.csv.gz",
        "wave45_resolved_direct_species_trait.csv.gz",
        "wave45_new_validated_low_species_trait.csv.gz",
        "wave45_new_trait_specific_genus_rule_audit.csv.gz",
        "wave45_change_audit.csv.gz",
        "wave45_external_congener_resolved_species_trait.csv.gz",
        "wave45_external_congener_source_conflicts.csv.gz",
    )

    def fake_build_wave45_overlay(**kwargs):
        output = kwargs["output_dir"]
        output.mkdir(parents=True, exist_ok=True)
        for name in old_names:
            _write(pd.DataFrame([{"value": "ok"}]), output / name)
        return {
            "contract": "wave45",
            "delta": {
                "resolved_wave45_direct_species_trait": 9,
                "resolved_wave45_direct_species_axis": 9,
            },
            "checks": {},
            "artifact_sha256": {},
        }

    monkeypatch.setattr(overlay_module, "build_wave45_overlay", fake_build_wave45_overlay)
    checkpoint = tmp_path / "checkpoint.json"
    checkpoint.write_text(
        json.dumps(
            {
                "checks": {
                    "retrieved_sources_verified": True,
                    "content_fingerprints_verified": True,
                }
            }
        ),
        encoding="utf-8",
    )
    output = tmp_path / "output"
    summary = overlay_module.build_wave48_overlay(
        baseline_csv=tmp_path / "unused-baseline.csv.gz",
        previous_rule_audit_csv=tmp_path / "unused-rules.csv.gz",
        all_evidence_dir=tmp_path / "unused-audit",
        direct_evidence_csv=tmp_path / "unused-direct.csv.gz",
        external_evidence_csv=tmp_path / "unused-external.csv.gz",
        checkpoint_summary_json=checkpoint,
        output_dir=output,
        expected_species=3,
    )
    assert summary["formal_wave33_run_id"] == 32932103226
    assert summary["baseline_formal_run_id"] == 33380486845
    assert summary["delta"]["resolved_wave48_direct_species_axis"] == 9
    assert summary["checks"]["all_source_receipts_verified"] is True
    assert (output / "wave48_species_axis_coverage.csv.gz").is_file()
    assert not list(output.glob("wave45_*"))
