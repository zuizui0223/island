import pandas as pd

from island_v2.chapter1_pr138_evidence_ladder import build_evidence_ladder


def test_direct_has_priority_and_validated_genus_fills_gaps_only():
    direct = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "resolution_status": "resolved",
                "normalized_value": "open_radial",
                "state_set": '["open_radial"]',
                "selected_quality": "high",
            }
        ]
    )
    low = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "state_set": '["tubular"]',
                "eligible": True,
                "diagnostic_only": False,
                "evidence_scope": "genus_consensus",
                "evidence_quality": "low",
            },
            {
                "accepted_species": "Beta two",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "state_set": '["tubular"]',
                "eligible": True,
                "diagnostic_only": False,
                "evidence_scope": "genus_consensus",
                "evidence_quality": "low",
            },
        ]
    )
    expanded, direct_only, genus_only, _ = build_evidence_ladder(direct, low)
    assert len(expanded) == 2
    alpha = expanded.loc[expanded["accepted_species"].eq("Alpha one")].iloc[0]
    assert alpha["normalized_value"] == "open_radial"
    assert alpha["evidence_origin"] == "direct_species_trait_ledger"
    beta = expanded.loc[expanded["accepted_species"].eq("Beta two")].iloc[0]
    assert beta["evidence_origin"] == "validated_genus_consensus"
    assert len(direct_only) == 1
    assert len(genus_only) == 2
