from pathlib import Path

import pandas as pd

from island_v2.real_acquisition_ledger import aggregate


def test_actual_collector_columns_map_into_quality_ledger(tmp_path: Path):
    root = tmp_path / "campaign"
    cumulative = root / "cumulative"
    cumulative.mkdir(parents=True)

    pd.DataFrame(
        [
            {"species": "Alpha beta", "status": "completed", "attempts": "1"},
            {"species": "Gamma delta", "status": "completed", "attempts": "1"},
            {"species": "Pending species", "status": "pending", "attempts": "0"},
        ]
    ).to_csv(root / "checkpoint.csv", index=False)
    pd.DataFrame(
        [
            {"species": "Alpha beta", "provider": "gbif", "status": "hit"},
            {"species": "Gamma delta", "provider": "web_descriptions", "status": "hit"},
        ]
    ).to_csv(root / "provider_checkpoint.csv", index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha beta",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "red",
                "source": "gbif_species_description",
                "source_type": "gbif_species_description",
                "source_url": "https://example.test/gbif",
                "raw_description": "Petals red.",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Gamma delta",
                "trait_name": "floral_form",
                "provisional_candidate_value": "tubular",
                "source": "horticultural_site",
                "source_type": "horticultural_site",
                "source_url": "https://example.test/garden",
                "raw_description": "Flowers tubular.",
                "evidence_scope": "species_direct",
            },
            {
                "accepted_species": "Gamma delta",
                "trait_name": "mating_system",
                "provisional_candidate_value": "self_compatible",
                "source": "genus_rule",
                "source_type": "genus_rule",
                "source_url": "",
                "raw_description": "",
                "evidence_scope": "genus_inference",
            },
            {
                "accepted_species": "Gamma delta",
                "trait_name": "flower_primary_color",
                "provisional_candidate_value": "white",
                "source": "family_rule",
                "source_type": "family_rule",
                "source_url": "",
                "raw_description": "",
                "evidence_scope": "family_inference",
            },
        ]
    ).to_csv(cumulative / "v2_reported_candidates.csv", index=False)

    summary = aggregate(root, tmp_path / "out")
    ledger = pd.read_csv(tmp_path / "out" / "all_master_evidence.csv.gz")

    assert summary["total_species_in_shard"] == 3
    assert summary["total_species_processed"] == 2
    assert summary["evidence_records_by_quality"] == {"high": 1, "medium": 1, "low": 1}
    assert summary["direct_coverage"]["flower_colour"] == 1
    assert summary["direct_coverage"]["floral_structural_complexity"] == 1
    assert summary["family_inference_used"] is False
    assert summary["global_fallback_used"] is False
    assert "family_inference" not in set(ledger["evidence_scope"])
    assert set(ledger["matched_source_name"]) == {"Alpha beta", "Gamma delta"}
    assert set(ledger["excerpt"].dropna()) >= {"Petals red.", "Flowers tubular."}
