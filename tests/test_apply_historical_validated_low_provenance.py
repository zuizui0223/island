from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).parents[1] / "analysis" / "apply_historical_validated_low_provenance.py"
)
spec = importlib.util.spec_from_file_location(
    "apply_historical_validated_low_provenance", MODULE_PATH
)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_historical_exact_and_semantic_mismatch_are_separated(tmp_path: Path) -> None:
    audit = tmp_path / "audit"
    audit.mkdir()
    _write(
        audit / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv",
        [
            {
                "derived_lineage": 'validated-low:Alpha:flower_primary_color:["red_pink"]',
                "genus": "Alpha",
                "trait_name": "flower_primary_color",
                "inferred_state_set": '["red_pink"]',
                "rule_recovered": "false",
                "support_source_lineages": "[]",
                "support_source_families": "",
                "provenance_rule_tier": "unresolved_historical_secondary",
            },
            {
                "derived_lineage": 'validated-low:Beta:flower_primary_color:["blue_purple","white"]',
                "genus": "Beta",
                "trait_name": "flower_primary_color",
                "inferred_state_set": '["blue_purple","white"]',
                "rule_recovered": "false",
                "support_source_lineages": "[]",
                "support_source_families": "",
                "provenance_rule_tier": "unresolved_historical_secondary",
            },
        ],
    )
    _write(
        audit / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv",
        [
            {
                "accepted_species": "Alpha target",
                "axis": "flower_colour",
                "derived_lineages": 'validated-low:Alpha:flower_primary_color:["red_pink"]',
            },
            {
                "accepted_species": "Beta target",
                "axis": "flower_colour",
                "derived_lineages": 'validated-low:Beta:flower_primary_color:["blue_purple","white"]',
            },
        ],
    )

    wave34 = tmp_path / "wave34"
    wave34.mkdir()
    _write(
        wave34 / "trait_specific_rule_frontier.csv.gz",
        [
            {
                "genus": "Alpha",
                "trait_name": "flower_primary_color",
                "predicted_state_set": '["red_pink"]',
                "n_direct_species": "3",
                "eligible": "true",
                "support_source_lineages": json.dumps(
                    ["url:https://example.org/a", "url:https://example.org/b"]
                ),
            },
            {
                "genus": "Beta",
                "trait_name": "flower_primary_color",
                "predicted_state_set": '["blue_purple"]',
                "n_direct_species": "3",
                "eligible": "true",
                "support_source_lineages": json.dumps(
                    ["baseflor:species:b1", "baseflor:species:b2", "baseflor:species:b3"]
                ),
            },
        ],
    )

    empty_generations: list[tuple[str, Path]] = []
    for label in ("wave35", "wave39", "wave40", "wave41", "wave42", "wave43"):
        root = tmp_path / label
        root.mkdir()
        empty_generations.append((label, root))

    summary = module.apply(audit, [("wave34", wave34)])

    assert summary["exact_rule_semantics_recovered"] == 1
    assert summary["source_provenance_recovered"] == 2
    assert summary["semantic_mismatch_rule_lineages"] == 1
    assert summary["target_cells_with_unresolved_source_provenance"] == 0
    assert summary["target_cells_with_zero_support"] == 0

    unmatched = pd.read_csv(audit / "VALIDATED_LOW_UNMATCHED_RULES.csv", dtype=str).fillna("")
    assert len(unmatched) == 1
    beta = unmatched.iloc[0]
    assert beta["source_provenance_recovered"].lower() == "true"
    assert beta["semantic_mismatch_receipt_state_set"] == '["blue_purple"]'
    assert int(beta["support_source_lineage_count"]) == 3
