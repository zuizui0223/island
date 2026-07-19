from pathlib import Path

import pandas as pd

from island_v2.full_coverage_trait_campaign import aggregate


def test_aggregate_preserves_unresolved_and_reports_provider_state(tmp_path: Path):
    root = tmp_path / "campaign" / "shard-0"
    root.mkdir(parents=True)
    pd.DataFrame(
        [
            {"species": "Alpha beta", "status": "completed"},
            {"species": "Gamma delta", "status": "completed"},
        ]
    ).to_csv(root / "checkpoint.csv", index=False)
    pd.DataFrame(
        [
            {"species": "Alpha beta", "provider": "flora", "status": "hit"},
            {"species": "Gamma delta", "provider": "flora", "status": "retry"},
        ]
    ).to_csv(root / "provider_checkpoint.csv", index=False)
    cumulative = root / "cumulative"
    cumulative.mkdir()
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "matched_source_name": "Alpha beta",
                "trait_name": "flower_primary_color",
                "reported_value": "red",
                "raw_candidate_value": "red",
                "source_url": "https://example.test/a",
                "provider": "flora",
                "source_title": "Flora",
                "excerpt": "Flowers red.",
                "retrieval_date": "2026-07-19",
                "stable_source_id": "record:a",
            }
        ]
    ).to_csv(cumulative / "evidence.csv", index=False)

    summary = aggregate(tmp_path / "campaign", tmp_path / "out")
    assert summary["total_species_processed"] == 2
    assert summary["species_successfully_queried"] == 1
    assert summary["direct_evidence_records_added"] == 1
    assert summary["direct_coverage"]["flower_colour"]["species"] == 1
    assert summary["species_with_all_three_axes_direct"] == 0
    assert summary["unresolved_species_count"] == 2
    assert summary["provider_errors_and_blocked_sources"] == {"flora": 1}
    assert summary["next_continuation_state"]["remaining_provider_tasks"] == 1
    assert not summary["global_fallback_used"]
