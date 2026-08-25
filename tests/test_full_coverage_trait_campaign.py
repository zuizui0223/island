from pathlib import Path

import pandas as pd

from island_v2.full_coverage_trait_campaign import aggregate


def test_aggregate_applies_high_medium_low_and_rejects_family(tmp_path: Path):
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
            {"species": "Alpha beta", "provider": "europe_pmc", "status": "hit"},
            {"species": "Gamma delta", "provider": "web_descriptions", "status": "retry"},
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
                "source_url": "https://literature.test/a",
                "provider": "europe_pmc",
                "source_title": "Journal article",
                "excerpt": "Flowers red.",
                "retrieval_date": "2026-07-19",
                "stable_source_id": "doi:10.1/a",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Alpha beta",
                "matched_source_name": "Alpha beta",
                "trait_name": "floral_form",
                "reported_value": "tubular",
                "source_url": "https://horticulture.test/a",
                "provider": "horticultural_web",
                "source_title": "Plant description",
                "excerpt": "Tubular flowers.",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Alpha beta",
                "matched_source_name": "Alpha",
                "trait_name": "mating_system",
                "reported_value": "mixed",
                "provider": "genus_consensus",
                "evidence_scope": "genus_inference",
            },
            {
                "accepted_species": "Gamma delta",
                "matched_source_name": "Gamma family",
                "trait_name": "flower_primary_color",
                "reported_value": "white",
                "provider": "family_consensus",
                "evidence_scope": "family_inference",
            },
        ]
    ).to_csv(cumulative / "evidence.csv", index=False)

    out = tmp_path / "out"
    summary = aggregate(tmp_path / "campaign", out)
    evidence = pd.read_csv(out / "all_master_evidence.csv.gz")
    coverage = pd.read_csv(out / "three_axis_coverage.csv.gz")

    assert summary["total_species_processed"] == 2
    assert summary["species_successfully_queried"] == 1
    assert summary["evidence_records_by_quality"] == {"high": 1, "low": 1, "medium": 1}
    assert summary["direct_evidence_records_added"] == 2
    assert summary["family_inference_used"] is False
    assert set(evidence["evidence_quality"]) == {"high", "medium", "low"}
    assert "family_inference" not in set(evidence["evidence_scope"])

    alpha = coverage.loc[coverage["accepted_species"].eq("Alpha beta")].iloc[0]
    gamma = coverage.loc[coverage["accepted_species"].eq("Gamma delta")].iloc[0]
    assert bool(alpha["flower_colour_direct"])
    assert bool(alpha["floral_structural_complexity_direct"])
    assert bool(alpha["reproductive_assurance_genus_inferred"])
    assert not bool(alpha["all_three_axes_direct"])
    assert bool(alpha["all_three_axes_resolved"])
    assert gamma["unresolved_axes"] != ""
    assert summary["provider_errors_and_blocked_sources"] == {"europe_pmc": 0, "web_descriptions": 1}
    assert summary["next_continuation_state"]["remaining_provider_tasks"] == 1
    assert not summary["global_fallback_used"]
