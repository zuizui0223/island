from __future__ import annotations

import pandas as pd

from island_v2.si_genus_inference import build_combined_genus_inferences


def _config() -> dict:
    return {
        "trait_name": "self_incompatibility",
        "explicit_values": ["SI", "SC"],
        "disqualifying_values": ["mixed_or_variable"],
        "minimum_direct_species": 5,
        "minimum_modal_purity": 1.0,
        "confidence": "low",
        "analysis_role": "sensitivity_only_genus_inference",
    }


def test_multi_source_genus_inference_disqualifies_any_mixed_or_conflicting_genus() -> None:
    master = pd.DataFrame(
        [
            *[(f"Alpha {name}", "Alpha", "FamA") for name in "abcdefg"],
            *[(f"Beta {name}", "Beta", "FamB") for name in "abcdef"],
            ("Gamma a", "Gamma", "FamC"),
            ("Gamma b", "Gamma", "FamC"),
        ],
        columns=["accepted_species", "genus", "family"],
    )
    direct_rows = [
        *[
            (f"Alpha {name}", "Alpha", "FamA", "SC", "source_one")
            for name in "abcdef"
        ],
        *[(f"Beta {name}", "Beta", "FamB", "SI", "source_one") for name in "abcde"],
        ("Beta a", "Beta", "FamB", "mixed_or_variable", "source_two"),
        ("Gamma a", "Gamma", "FamC", "SI", "source_one"),
        ("Gamma a", "Gamma", "FamC", "SC", "source_two"),
    ]
    direct = pd.DataFrame(
        direct_rows,
        columns=[
            "accepted_species",
            "genus",
            "family",
            "standardized_value",
            "source_name",
        ],
    )
    inference, qualifying, report = build_combined_genus_inferences(
        master, direct, _config()
    )
    assert inference["accepted_species"].tolist() == ["Alpha g"]
    assert inference["standardized_value"].tolist() == ["SC"]
    assert inference["candidate_kind"].eq("hierarchical_inference").all()
    assert inference["evidence_scope"].eq("genus_inference").all()
    assert inference["confidence"].eq("low").all()
    assert qualifying["genus"].tolist() == ["Alpha"]
    assert report["n_ambiguous_or_mixed_direct_species"] == 2
    assert report["leave_one_out_n"] == 6
    assert report["leave_one_out_accuracy"] == 1.0
