from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).parents[1]
    / "analysis"
    / "recover_historical_validated_low_rights_provenance.py"
)
spec = importlib.util.spec_from_file_location("historical_validated_low_provenance", MODULE_PATH)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, object]], columns=None) -> Path:
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False)
    return path


def test_historical_exact_and_semantic_mismatch_are_separated(tmp_path: Path) -> None:
    rights = tmp_path / "rights"
    rights.mkdir()
    exact_lineage = 'validated-low:Alpha:floral_form:["open_radial"]'
    mismatch_lineage = (
        'validated-low:Beta:flower_primary_color:'
        '["blue_purple","green_brown_inconspicuous"]'
    )
    base_rule = {
        "provenance_rule_tier": "unresolved_historical_secondary",
        "rule_recovered": "False",
        "expected_n_direct_species": "0",
        "reconstructed_n_direct_species": "0",
        "direct_species_count_matches": "False",
        "support_source_lineage_count": "0",
        "support_source_lineages": "[]",
        "support_source_family_count": "0",
        "support_source_families": "",
    }
    _write(
        rights / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv",
        [
            {
                **base_rule,
                "derived_lineage": exact_lineage,
                "genus": "Alpha",
                "trait_name": "floral_form",
                "inferred_state_set": '["open_radial"]',
            },
            {
                **base_rule,
                "derived_lineage": mismatch_lineage,
                "genus": "Beta",
                "trait_name": "flower_primary_color",
                "inferred_state_set": '["blue_purple","green_brown_inconspicuous"]',
            },
        ],
    )
    _write(
        rights / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv",
        [
            {
                "accepted_species": "Alpha target",
                "axis": "floral_structural_complexity",
                "derived_lineages": exact_lineage,
            },
            {
                "accepted_species": "Beta target",
                "axis": "flower_colour",
                "derived_lineages": mismatch_lineage,
            },
        ],
    )
    frontier = _write(
        tmp_path / "frontier.csv",
        [
            {
                "genus": "Alpha",
                "trait_name": "floral_form",
                "predicted_state_set": '["open_radial"]',
                "n_direct_species": "3",
                "eligible": "True",
                "support_species": '["Alpha one","Alpha two","Alpha three"]',
                "support_source_lineages": (
                    '["url:https://a.example/1","url:https://a.example/2"]'
                ),
            }
        ],
    )
    sidecar = _write(
        tmp_path / "sidecar.csv",
        [
            {
                "accepted_species": "Beta",
                "trait_names": '["flower_primary_color"]',
                "predicted_state_sets": '["[\\"blue_purple\\"]"]',
                "support_source_lineages": (
                    '["baseflor:species:beta_one","baseflor:species:beta_two",'
                    '"baseflor:species:beta_three"]'
                ),
                "min_direct_species_support": "3",
            }
        ],
    )
    direct_columns = ["accepted_species", "trait_name", "source_lineages"]
    rule_columns = [
        "genus",
        "trait_name",
        "inferred_state_set",
        "n_direct_species",
        "eligible",
    ]
    coverage_columns = [
        "accepted_species",
        "trait_names",
        "source_lineages",
        "quality",
    ]
    direct = _write(tmp_path / "direct.csv", [], direct_columns)
    rules = _write(tmp_path / "rules.csv", [], rule_columns)
    coverage = _write(tmp_path / "coverage.csv", [], coverage_columns)

    summary = module.extend(
        rights,
        None,
        direct,
        frontier,
        sidecar,
        rules,
        coverage,
        direct,
        rules,
        direct,
        rules,
        direct,
        rules,
        direct,
        direct,
        rules,
        direct,
        direct,
        rules,
        direct,
        direct,
    )

    assert summary["exact_rule_generations_recovered"] == 1
    assert summary["source_provenance_recovered_rule_lineages"] == 2
    assert summary["new_historical_exact_rule_recoveries"] == 1
    assert summary["historical_state_set_mismatch_rules"] == 1
    assert summary["unresolved_exact_rule_cells"] == 1
    assert summary["unresolved_source_provenance_cells"] == 0
    assert summary["zero_support_cells"] == 0

    mismatches = pd.read_csv(rights / "VALIDATED_LOW_SEMANTIC_MISMATCHES.csv")
    assert mismatches["genus"].tolist() == ["Beta"]
    cells = pd.read_csv(rights / "VALIDATED_LOW_PROVENANCE_RECOVERY.csv")
    beta = cells.loc[cells["accepted_species"].eq("Beta target")].iloc[0]
    assert not bool(beta["all_rules_exact_rule_recovered"])
    assert bool(beta["all_rules_source_provenance_recovered"])
