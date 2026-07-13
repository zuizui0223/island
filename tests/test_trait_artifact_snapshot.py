from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.trait_artifact_snapshot import build_snapshot


def test_snapshot_separates_machine_candidate_layers(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    wave = campaign / "waves" / "wave_00001_reproductive_wikimedia"
    wave.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "accepted_species": "Plant alpha",
                "trait_name": "pollen_vector_mode",
                "candidate_value": "biotic",
                "source_url": "https://example.org/a",
                "source_citation": "A",
                "source_excerpt": "pollinated by bees",
                "source_type": "wikimedia",
                "evidence_scope": "species_direct",
                "matched_term": "pollinated by bees",
                "confidence": "machine_extracted",
                "campaign_task": "reproductive_wikimedia",
                "campaign_phase": "reproductive_and_pollen_vector",
                "target_for_task": True,
            },
            {
                "accepted_species": "Plant beta",
                "trait_name": "floral_form",
                "candidate_value": "tubular",
                "source_url": "https://example.org/b",
                "source_citation": "B",
                "source_excerpt": "tubular flowers",
                "source_type": "wikimedia",
                "evidence_scope": "species_direct",
                "matched_term": "tubular",
                "confidence": "machine_extracted",
                "campaign_task": "floral_access_wikimedia",
                "campaign_phase": "floral_access",
                "target_for_task": True,
            },
        ]
    ).to_csv(wave / "machine_candidates.csv.gz", index=False, compression="gzip")
    pd.DataFrame([{"accepted_species": "Plant alpha"}, {"accepted_species": "Plant beta"}]).to_csv(
        campaign / "campaign_ledger.csv.gz", index=False, compression="gzip"
    )
    (campaign / "campaign_status.json").write_text(
        json.dumps({"n_global_species": 2, "active_task": "reproductive_wikimedia", "tasks": {}}),
        encoding="utf-8",
    )

    output = tmp_path / "snapshot"
    summary = build_snapshot(campaign, output)

    assert summary["n_target_candidate_rows"] == 2
    reproductive = pd.read_csv(output / "reproductive_and_pollen_vector_candidates.csv.gz")
    floral = pd.read_csv(output / "floral_access_candidates.csv.gz")
    assert reproductive["trait_name"].tolist() == ["pollen_vector_mode"]
    assert floral["trait_name"].tolist() == ["floral_form"]
    assert not reproductive["analysis_eligible"].astype(bool).any()
    assert set(reproductive["review_status"]) == {"machine_candidate_unreviewed"}
