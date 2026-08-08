from pathlib import Path

import pandas as pd

from analysis.prepare_nhm_adept_bhl_source_package_20260808 import (
    _extract_record,
)
from island_v2.open_web_common import reviewed_source_package_evidence


def test_bhl_record_keeps_multistate_colour_and_separate_structure_traits() -> None:
    record = {
        "_id": "row-1",
        "taxon": "Example plantus",
        "source_id": "123",
        "flower colour": "lavender, white",
        "petal colour": "",
        "corolla colour": "",
        "labellum colour": "",
        "flower symmetry": "zygomorphic",
        "flower shape": "campanulate, tubular",
        "inflorescence arrangement": "panicle, corymb",
        "corolla tube min_ length [cm]": "0.5 cm",
        "corolla tube max_ length [cm]": "1.5 cm",
        "reproduction system": "dioecious",
        "pollination": "wind",
        "reward": "nectar",
    }
    rows = {row["trait_name"]: row for row in _extract_record(record)}

    assert rows["flower_primary_color"]["normalized_value"] == "blue_purple|white"
    assert rows["floral_symmetry"]["normalized_value"] == "zygomorphic"
    assert rows["floral_form"]["normalized_value"] == "bell_campanulate|tubular"
    assert rows["inflorescence_display"]["normalized_value"] == (
        "raceme_spike_panicle|umbel_corymb"
    )
    assert rows["tube_depth_class"]["normalized_value"] == "intermediate"
    assert set(rows) == {
        "flower_primary_color",
        "floral_symmetry",
        "floral_form",
        "inflorescence_display",
        "tube_depth_class",
    }


def test_committed_nhm_bhl_package_passes_the_common_review_gate() -> None:
    root = Path("data/v2/staging/traits/direct_llm_pilot/20260808_nhm_adept_bhl")
    evidence = pd.read_csv(root / "nhm_adept_bhl_evidence.csv.gz", dtype=str).fillna("")
    audit = pd.read_csv(root / "nhm_adept_bhl_manual_audit_200.csv", dtype=str).fillna("")

    selected, scopes, summary = reviewed_source_package_evidence(evidence, audit)

    assert summary["reviewed"] == 200
    assert summary["accepted_correct"] == 200
    assert summary["precision"] == 1.0
    assert summary["cultivar_contamination_rate"] == 0.0
    assert summary["package_gate_passed"]
    assert len(selected) == 27095
    assert selected["source_record_id"].is_unique
    assert selected["source_lineage"].str.startswith("bhl:page:").all()
    assert set(scopes.loc[scopes["production_approved"], "trait_name"]) == {
        "flower_primary_color",
        "floral_symmetry",
        "floral_form",
        "inflorescence_display",
        "tube_depth_class",
    }
    assert not set(selected["trait_name"]).intersection(
        {"sex_system", "pollen_vector_mode", "reward_type"}
    )
