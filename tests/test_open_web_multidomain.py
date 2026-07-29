from __future__ import annotations

from pathlib import Path

import pandas as pd

from island_v2.open_web_common import STRICT_TRAIT_AXIS
from island_v2.open_web_multidomain import _hash_sample


def test_multidomain_audit_sample_is_hash_based_trait_stratified_and_unique() -> None:
    rows = []
    for trait in ("flower_primary_color", "floral_form"):
        for index in range(20):
            rows.append(
                {
                    "candidate_id": f"{trait}-{index}",
                    "trait_name": trait,
                    "domain": ("one.example" if index % 2 else "two.example"),
                    "source_url": f"https://example.test/{trait}/{index}",
                }
            )
    frame = pd.DataFrame(rows)
    first = _hash_sample(frame, sample_pages=20, seed="fixed")
    second = _hash_sample(frame, sample_pages=20, seed="fixed")
    assert first["candidate_id"].tolist() == second["candidate_id"].tolist()
    assert first["source_url"].nunique() == 20
    assert first.groupby(["domain", "trait_name"]).size().nunique() == 1
    assert first.groupby("trait_name").size().to_dict() == {
        "floral_form": 10,
        "flower_primary_color": 10,
    }


def test_pollen_vector_does_not_fill_structure_axis() -> None:
    assert "pollen_vector_mode" not in STRICT_TRAIT_AXIS
    assert "reward_type" not in STRICT_TRAIT_AXIS


def test_multidomain_workflow_invokes_runner_subcommand_and_fails_closed() -> None:
    workflow = Path(
        ".github/workflows/open-web-trait-multidomain-pilot.yml"
    ).read_text(encoding="utf-8")
    assert "python -m island_v2.open_web_multidomain_runner run \\" in workflow
    assert 'if report["domain_artifacts"] < 5:' in workflow
    assert 'if report["audit_pages_sampled"] < 200:' in workflow
    assert 'if report["audit_domains_sampled"] < 5:' in workflow
