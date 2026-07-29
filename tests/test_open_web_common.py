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
    coverage_change_counts,
    reviewed_audit_metrics,
)
from island_v2.open_web_finalize import combine_public_web_ledgers


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
    workflow = Path(".github/workflows/open-web-review-promote.yml").read_text(
        encoding="utf-8"
    )
    assert "pull_request:" in workflow
    assert "BASELINE_RUN_ID:" in workflow
    assert "'30433986432'" in workflow
    assert "'30434418380'" in workflow
    assert "open-web-multidomain-pilot-30434418380" in workflow
    assert 'PRIOR_PUBLIC_WEB_RUN_ID: "30141880859"' in workflow
    assert "broad-web-medium-full-30141880859" in workflow
    assert "--prior-public-web-csv" in workflow
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
