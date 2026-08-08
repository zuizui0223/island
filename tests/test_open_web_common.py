from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.open_web_checkpoint import finalize_shard, prepare_shard
from island_v2.open_web_common import (
    NON_AXIS_TRAITS,
    STRICT_TRAIT_AXIS,
    compare_incremental_validated_low,
    coverage_change_counts,
    production_candidate_ids,
    promoted_to_common_evidence,
    reviewed_audit_metrics,
    reviewed_source_package_evidence,
)
from island_v2.open_web_finalize import (
    combine_public_web_ledgers,
    validate_individually_reviewed_evidence,
)


def _review_rows(domain: str, trait: str, count: int) -> list[dict[str, object]]:
    return [
        {
            "candidate_id": f"{domain}-{trait}-{index}",
            "domain": domain,
            "trait_name": trait,
            "decision": "accept",
            "species_identity_correct": True,
            "value_correct": True,
            "provenance_complete": True,
            "cultivar_contamination": False,
            "reviewer": "reviewer",
            "reviewed_at_utc": "2026-07-29T00:00:00Z",
        }
        for index in range(count)
    ]


def test_strict_axes_exclude_reward_and_pollen_vector_everywhere() -> None:
    assert set(NON_AXIS_TRAITS) == {"reward_type", "pollen_vector_mode"}
    assert "reward_type" not in STRICT_TRAIT_AXIS
    assert "pollen_vector_mode" not in STRICT_TRAIT_AXIS
    config = yaml.safe_load(
        Path("config/open_web_trait_acquisition.yml").read_text(encoding="utf-8")
    )
    strict = {trait for traits in config["axis_traits"].values() for trait in traits}
    assert "reward_type" not in strict
    assert "pollen_vector_mode" not in strict
    assert set(config["non_axis_traits"]) == set(NON_AXIS_TRAITS)


def test_pr132_has_no_independent_validated_low_implementation() -> None:
    pilot = Path("src/island_v2/open_web_pilot.py").read_text(encoding="utf-8")
    finalize = Path("src/island_v2/open_web_finalize.py").read_text(encoding="utf-8")
    common = Path("src/island_v2/open_web_common.py").read_text(encoding="utf-8")
    assert "def rebuild_all_direct_validated_low" not in pilot
    assert "def rebuild_validated_low" not in finalize
    assert "from island_v2.all_evidence_trait_audit import" in common
    assert "build_rule_audit(" in common
    assert "apply_genus_rules(" in common


def test_review_promotion_can_run_on_pr_with_exact_pinned_artifacts() -> None:
    workflow = Path(".github/workflows/open-web-review-promote.yml").read_text(encoding="utf-8")
    assert "pull_request:" in workflow
    assert "BASELINE_RUN_ID:" in workflow
    assert "'31245813999'" in workflow
    assert "'30434418380'" in workflow
    assert "open-web-multidomain-pilot-30434418380" in workflow
    assert "inputs.prior_public_web_run_id || '31257281120'" in workflow
    assert "reviewed-open-web-evidence-31257281120" in workflow
    assert "broad_web_medium_evidence.csv.gz" in workflow
    assert "prior_public_web_run_id:" in workflow
    assert "prior_public_web_artifact_name:" in workflow
    assert "prior_public_web_file_name:" in workflow
    assert '"/tmp/open-web-prior/$PRIOR_PUBLIC_WEB_FILE_NAME"' in workflow
    assert "--prior-public-web-csv" in workflow
    assert "source_package_evidence_csv_path:" in workflow
    assert "source_package_audit_csv_path:" in workflow
    assert "20260808_scale_source_acquisition/scale_source_incremental_evidence.csv.gz" in workflow
    assert "20260808_scale_source_acquisition/scale_source_manual_audit.csv" in workflow
    assert "scale-source-acquisition-20260808-v2" in workflow
    assert "--source-package-evidence-csv" in workflow
    assert ".source_package.novelty_rate >= 0.50" in workflow
    assert "coverage_change_species_axis.net_change" in workflow


def test_precision_is_accepted_correct_over_every_reviewed_row() -> None:
    rows = _review_rows("example.org", "floral_form", 10)
    rows.append(
        {
            **rows[0],
            "candidate_id": "bad-accept",
            "value_correct": False,
        }
    )
    _, summary = reviewed_audit_metrics(pd.DataFrame(rows))
    assert summary["accepted_correct"] == 10
    assert summary["reviewed_pages_or_candidates"] == 11
    assert summary["precision"] == pytest.approx(10 / 11)


def test_domain_trait_gate_requires_ten_reviews_and_thresholds() -> None:
    rows = _review_rows("good.example", "floral_form", 10)
    rows.extend(_review_rows("thin.example", "floral_form", 9))
    scopes, _ = reviewed_audit_metrics(pd.DataFrame(rows))
    approved = scopes.set_index(["domain", "trait_name"])["production_approved"]
    assert bool(approved.loc[("good.example", "floral_form")])
    assert not bool(approved.loc[("thin.example", "floral_form")])


def test_domain_trait_gate_requires_95_percent_precision() -> None:
    passing = _review_rows("pass.example", "flower_size_class", 20)
    passing[0]["value_correct"] = False
    failing = _review_rows("fail.example", "flower_size_class", 20)
    failing[0]["value_correct"] = False
    failing[1]["value_correct"] = False
    scopes, _ = reviewed_audit_metrics(pd.DataFrame([*passing, *failing]))
    approved = scopes.set_index(["domain", "trait_name"])["production_approved"]
    assert bool(approved.loc[("pass.example", "flower_size_class")])
    assert not bool(approved.loc[("fail.example", "flower_size_class")])


def test_production_approval_scales_scope_but_excludes_reviewed_failures() -> None:
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "good-reviewed",
                "domain": "flora.example",
                "trait_name": "flower_size_class",
            },
            {
                "candidate_id": "bad-reviewed",
                "domain": "flora.example",
                "trait_name": "flower_size_class",
            },
            {
                "candidate_id": "unreviewed",
                "domain": "flora.example",
                "trait_name": "flower_size_class",
            },
            {
                "candidate_id": "blocked-scope",
                "domain": "flora.example",
                "trait_name": "floral_form",
            },
        ]
    )
    audit = pd.DataFrame(
        [
            {
                **_review_rows("flora.example", "flower_size_class", 1)[0],
                "candidate_id": "good-reviewed",
            },
            {
                **_review_rows("flora.example", "flower_size_class", 1)[0],
                "candidate_id": "bad-reviewed",
                "decision": "reject",
                "value_correct": False,
            },
        ]
    )
    scopes = pd.DataFrame(
        [
            {
                "domain": "flora.example",
                "trait_name": "flower_size_class",
                "production_approved": True,
            },
            {"domain": "flora.example", "trait_name": "floral_form", "production_approved": False},
        ]
    )
    assert production_candidate_ids(candidates, audit, scopes) == {
        "good-reviewed",
        "unreviewed",
    }


def test_source_package_gate_scales_only_passing_traits_and_excludes_failures() -> None:
    audit_rows: list[dict[str, object]] = []
    evidence_rows: list[dict[str, object]] = []
    traits = [
        ("floral_form", 20, {0}),
        ("flower_size_class", 20, {0, 1}),
        ("flower_primary_color", 160, set()),
    ]
    for trait, count, failures in traits:
        axis = {
            "floral_form": "floral_structural_complexity",
            "flower_size_class": "floral_structural_complexity",
            "flower_primary_color": "flower_colour",
        }[trait]
        value = {
            "floral_form": "tubular",
            "flower_size_class": "small",
            "flower_primary_color": "white",
        }[trait]
        for index in range(count):
            candidate_id = f"{trait}-{index}"
            audit_rows.append(
                {
                    "candidate_id": candidate_id,
                    "trait_name": trait,
                    "accepted_correct": index not in failures,
                    "cultivar_status": "wild_or_flora_treatment",
                    "reviewer": "reviewer",
                    "reviewed_at_utc": "2026-08-08T00:00:00Z",
                    "audit_reason": "reviewed",
                }
            )
            evidence_rows.append(
                {
                    "accepted_species": f"Plantus {trait}{index}",
                    "axis": axis,
                    "trait_name": trait,
                    "normalized_value": value,
                    "quality": "high",
                    "source_group": "source_package",
                    "source_provider": "official_flora",
                    "source_url": f"https://flora.example/{candidate_id}",
                    "source_record_id": candidate_id,
                    "source_citation": "Official flora treatment",
                    "source_excerpt": f"Flowers are {value}.",
                    "evidence_scope": "species_direct",
                    "name_match_method": "accepted_name_exact",
                    "source_lineage": f"treatment:{candidate_id}",
                    "lineage_method": "provider_treatment_id",
                    "source_run_id": "package-run",
                    "source_artifact": "package-artifact",
                    "source_file": "package.csv.gz",
                    "acceptance_contract": "source_package_audit_v1",
                }
            )

    selected, scopes, summary = reviewed_source_package_evidence(
        pd.DataFrame(evidence_rows),
        pd.DataFrame(audit_rows),
    )

    approval = scopes.set_index("trait_name")["production_approved"]
    assert bool(approval.loc["floral_form"])
    assert bool(approval.loc["flower_primary_color"])
    assert not bool(approval.loc["flower_size_class"])
    assert summary["package_gate_passed"]
    assert summary["precision"] == pytest.approx(197 / 200)
    assert len(selected) == 179
    assert "floral_form-0" not in set(selected["source_record_id"])
    assert "flower_size_class-2" not in set(selected["source_record_id"])
    assert set(selected["trait_name"]) == {"floral_form", "flower_primary_color"}


def test_committed_wfo_source_package_audit_has_only_expected_approved_traits() -> None:
    root = Path(
        "data/v2/staging/traits/direct_llm_pilot/"
        "20260807_wfo_round2_context_gated_integrated"
    )
    evidence = pd.read_csv(root / "accepted_wfo_common_direct_evidence.csv.gz", dtype=str)
    audit = pd.read_csv(root / "stratified_trait_audit_200_20260808.csv", dtype=str)
    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 195
    assert summary["precision"] == pytest.approx(0.975)
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"]
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "floral_form",
        "floral_symmetry",
        "flower_primary_color",
        "inflorescence_display",
        "tube_depth_class",
    }
    assert set(selected["trait_name"]) == {
        "floral_form",
        "floral_symmetry",
        "flower_primary_color",
        "inflorescence_display",
        "tube_depth_class",
    }
    assert "WFOB-64bc03515616f3dc05b7" not in set(selected["source_record_id"])


def test_coverage_change_reports_gross_loss_and_net_separately() -> None:
    before = {
        ("A a", "flower_colour"),
        ("B b", "flower_colour"),
        ("C c", "reproductive_assurance"),
    }
    after = {
        ("B b", "flower_colour"),
        ("D d", "flower_colour"),
        ("E e", "floral_structural_complexity"),
    }
    change = coverage_change_counts(before, after)
    assert change["gross_gain"] == 2
    assert change["loss"] == 2
    assert change["net_change"] == 0
    assert change["by_axis"]["flower_colour"] == {
        "gross_gain": 1,
        "loss": 1,
        "net_change": 0,
    }
    assert change["by_axis"]["reproductive_assurance"]["net_change"] == -1
    assert change["by_axis"]["floral_structural_complexity"]["net_change"] == 1


def test_incremental_low_change_separates_add_invalidate_upgrade_and_value_change() -> None:
    before = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "state_set": '["tubular"]',
            },
            {
                "accepted_species": "Beta two",
                "trait_name": "flower_primary_color",
                "state_set": '["white"]',
            },
            {
                "accepted_species": "Gamma three",
                "trait_name": "floral_form",
                "state_set": '["rotate"]',
            },
        ]
    )
    after = pd.DataFrame(
        [
            {
                "accepted_species": "Delta four",
                "trait_name": "floral_form",
                "state_set": '["tubular"]',
            },
            {
                "accepted_species": "Gamma three",
                "trait_name": "floral_form",
                "state_set": '["campanulate"]',
            },
        ]
    )
    incremental_direct = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
            }
        ]
    )
    detail, summary = compare_incremental_validated_low(
        before,
        after,
        incremental_direct,
    )
    assert summary == {
        "added": 1,
        "invalidated": 1,
        "upgraded_to_direct": 1,
        "inferred_value_changed": 1,
        "retained": 0,
    }
    assert set(detail["pilot_change_status"]) == {
        "added",
        "invalidated",
        "upgraded_to_incremental_direct",
        "inferred_value_changed",
    }


def test_prior_formal_web_provenance_survives_common_schema_conversion() -> None:
    converted = promoted_to_common_evidence(
        pd.DataFrame(
            [
                {
                    "accepted_species": "Alpha one",
                    "trait_name": "floral_form",
                    "normalized_value": "tubular",
                    "evidence_quality": "medium",
                    "source_record_id": "prior-record",
                    "source_lineage": "citation:prior",
                    "source_run_id": "30141880859",
                    "source_artifact": "broad-web-medium-full-30141880859",
                }
            ]
        )
    )
    assert converted.loc[0, "source_record_id"] == "prior-record"
    assert converted.loc[0, "lineage_method"] == "artifact_source_lineage"
    assert converted.loc[0, "source_run_id"] == "30141880859"


def test_reviewed_web_ledger_appends_to_prior_instead_of_replacing_it() -> None:
    prior = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "source_lineage": "page:alpha",
                "source_url": "https://prior.example/alpha",
            },
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "source_lineage": "page:alpha",
                "source_url": "https://prior.example/alpha",
            },
        ]
    )
    promoted = pd.DataFrame(
        [
            {
                "accepted_species": "Beta two",
                "trait_name": "flower_primary_color",
                "normalized_value": "white",
                "source_lineage": "page:beta",
                "source_url": "https://new.example/beta",
                "source_run_id": "new-run",
                "source_artifact": "new-artifact",
            }
        ]
    )
    combined = combine_public_web_ledgers(
        prior,
        promoted,
        prior_run_id="prior-run",
        prior_artifact="prior-artifact",
    )
    assert len(combined) == 2
    by_species = combined.set_index("accepted_species")
    assert by_species.loc["Alpha one", "source_run_id"] == "prior-run"
    assert by_species.loc["Alpha one", "source_artifact"] == "prior-artifact"
    assert by_species.loc["Beta two", "source_run_id"] == "new-run"


def _individual_candidate() -> dict[str, object]:
    excerpt = "Alpha one has actinomorphic flowers."
    lineage = "doi:10.1234/example"
    return {
        "candidate_id": "66a9d783fb17b8510f7af25d",
        "accepted_species": "Alpha one",
        "trait_name": "floral_symmetry",
        "normalized_value": "actinomorphic",
        "evidence_quality": "high",
        "evidence_scope": "species_direct",
        "source_url": "https://journal.example/article.pdf",
        "source_citation": "Original article, p. 3",
        "source_excerpt": excerpt,
        "source_lineage": lineage,
        "name_match_method": "accepted_name_exact",
        "source_tier": "A",
        "content_sha256": "a" * 64,
        "content_sha256_basis": "downloaded_pdf_bytes",
        "retrieved_at_utc": "2026-08-08T00:00:00Z",
        "inference_rule": "",
    }


def _individual_audit() -> dict[str, object]:
    excerpt = "Alpha one has actinomorphic flowers."
    return {
        "candidate_id": "66a9d783fb17b8510f7af25d",
        "decision": "accept",
        "species_identity_correct": True,
        "value_correct": True,
        "provenance_complete": True,
        "cultivar_contamination": False,
        "reviewer": "reviewer",
        "reviewed_at_utc": "2026-08-08T00:05:00Z",
        "supporting_excerpt": excerpt,
        "normalized_excerpt_sha256": (
            "0ffb879a5307a170c0ffd3df61975cfb545c12141cd845d232b2f3d3b1f2e576"
        ),
        "content_fingerprint": ("93672467f7819aea0f5121b3d362199cc714eceb922308f075c7fa66d3a74c24"),
    }


def test_individual_manual_review_accepts_source_backed_direct_evidence(
    tmp_path: Path,
) -> None:
    ontology = tmp_path / "ontology.yml"
    ontology.write_text(
        "traits:\n  floral_symmetry:\n    allowed_values:\n      - actinomorphic\n",
        encoding="utf-8",
    )
    result = validate_individually_reviewed_evidence(
        pd.DataFrame([_individual_candidate()]),
        pd.DataFrame([_individual_audit()]),
        pd.DataFrame([{"accepted_species": "Alpha one"}]),
        ontology,
    )
    assert result["candidate_id"].tolist() == ["66a9d783fb17b8510f7af25d"]


def test_committed_high_leverage_batch_is_fully_reviewed_and_reproducible() -> None:
    root = Path("data/v2/staging/traits/open_web_pilot")
    candidates = pd.read_csv(
        root / "high_leverage_species_direct_evidence_20260808.csv",
        dtype=str,
    ).fillna("")
    audit = pd.read_csv(
        root / "high_leverage_species_direct_manual_audit_20260808.csv",
        dtype=str,
    ).fillna("")
    master = pd.read_csv(
        "data/v2/staging/gbif/collected/island_taxa.csv",
        dtype=str,
    ).fillna("")
    manifest = json.loads(
        (root / "high_leverage_source_manifest_20260808.json").read_text(
            encoding="utf-8"
        )
    )

    accepted = validate_individually_reviewed_evidence(
        candidates,
        audit,
        master,
        Path("config/trait_ontology.yml"),
    )
    assert len(accepted) == len(candidates) == len(audit) == 5
    assert set(accepted["candidate_id"]) == {
        record["candidate_id"] for record in manifest["records"]
    }
    assert "034e50e8e8d6fccc5109e223" in set(accepted["candidate_id"])


@pytest.mark.parametrize(
    ("candidate_patch", "audit_patch", "match"),
    [
        ({"source_tier": "B", "evidence_quality": "high"}, {}, "high_without_tier_a"),
        ({"name_match_method": "fuzzy"}, {}, "invalid_name_match"),
        ({"normalized_value": "unknown"}, {}, "invalid_ontology_value"),
        ({}, {"cultivar_contamination": True}, "accepted and correct"),
    ],
)
def test_individual_manual_review_fails_closed(
    tmp_path: Path,
    candidate_patch: dict[str, object],
    audit_patch: dict[str, object],
    match: str,
) -> None:
    ontology = tmp_path / "ontology.yml"
    ontology.write_text(
        "traits:\n  floral_symmetry:\n    allowed_values:\n      - actinomorphic\n",
        encoding="utf-8",
    )
    candidate = {**_individual_candidate(), **candidate_patch}
    audit = {**_individual_audit(), **audit_patch}
    with pytest.raises(ValueError, match=match):
        validate_individually_reviewed_evidence(
            pd.DataFrame([candidate]),
            pd.DataFrame([audit]),
            pd.DataFrame([{"accepted_species": "Alpha one"}]),
            ontology,
        )


def test_registry_has_at_least_five_generic_inventory_domains() -> None:
    registry = yaml.safe_load(
        Path("config/open_web_domain_registry_v2.yml").read_text(encoding="utf-8")
    )
    approved = [
        row
        for row in registry["domains"]
        if row["review_status"] == "approved" and row["inventory_urls"]
    ]
    assert len(approved) >= 5
    assert {row["domain"] for row in approved} >= {
        "mikawanoyasou.org",
        "gobotany.nativeplanttrust.org",
        "nzpcn.org.nz",
    }


def test_checkpoint_rejects_overlap_and_skips_completed_tasks(
    tmp_path: Path,
) -> None:
    tasks = pd.DataFrame(
        [
            {
                "task_id": f"task-{index}",
                "accepted_species": f"Plantus {index}",
                "trait_name": "floral_form",
                "query": f"query {index}",
            }
            for index in range(6)
        ]
    )
    tasks_csv = tmp_path / "tasks.csv"
    tasks.to_csv(tasks_csv, index=False)
    first = prepare_shard(
        tasks_csv=tasks_csv,
        output_dir=tmp_path / "first",
        offset=0,
        limit=3,
        source_run_id="100",
    )
    discovery = tmp_path / "discovery.json"
    discovery.write_text(
        json.dumps({"attempted_task_ids": first["pending_task_ids"]}),
        encoding="utf-8",
    )
    manifest_path = tmp_path / "first" / "checkpoint_manifest.json"
    finalize_shard(
        checkpoint_plan_json=tmp_path / "first" / "checkpoint_plan.json",
        discovery_report_json=discovery,
        output_json=manifest_path,
    )

    with pytest.raises(ValueError, match="overlaps prior run"):
        prepare_shard(
            tasks_csv=tasks_csv,
            output_dir=tmp_path / "overlap",
            offset=2,
            limit=3,
            source_run_id="101",
            prior_root=tmp_path / "first",
        )
    with pytest.raises(ValueError, match="already complete"):
        prepare_shard(
            tasks_csv=tasks_csv,
            output_dir=tmp_path / "resume",
            offset=0,
            limit=3,
            source_run_id="102",
            prior_root=tmp_path / "first",
            resume_run_id="100",
        )
