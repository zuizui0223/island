from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "recover_validated_low_rights_provenance.py"
spec = importlib.util.spec_from_file_location("recover_validated_low_rights_provenance", MODULE_PATH)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_tiered_rule_provenance_uses_external_congener_and_run329(tmp_path: Path) -> None:
    coverage = _write(
        tmp_path / "coverage.csv",
        [
            {
                "accepted_species": "Alpha target",
                "axis": "flower_colour",
                "quality": "low",
                "source_lineages": 'validated-low:Alpha:flower_primary_color:["red_pink"]',
            },
            {
                "accepted_species": "Beta target",
                "axis": "reproductive_assurance",
                "quality": "low",
                "source_lineages": 'validated-low:Beta:self_incompatibility:["SC"]',
            },
            {
                "accepted_species": "Gamma target",
                "axis": "floral_structural_complexity",
                "quality": "low",
                "source_lineages": 'validated-low:Gamma:floral_form:["open_radial"]',
            },
        ],
    )
    wave53_direct = _write(
        tmp_path / "wave53_direct.csv",
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "flower_primary_color",
                "source_lineages": "url:https://example.org/a",
            }
        ],
    )
    wave53_external = _write(
        tmp_path / "wave53_external.csv",
        [
            {
                "accepted_species": "Alpha two",
                "trait_name": "flower_primary_color",
                "source_lineages": "dataset:dryad.example:row2",
            }
        ],
    )
    wave53_rules = _write(
        tmp_path / "wave53_rules.csv",
        [
            {
                "setting": "current_min3",
                "genus": "Alpha",
                "trait_name": "flower_primary_color",
                "inferred_state_set": '["red_pink"]',
                "n_direct_species": "2",
                "eligible": "true",
                "diagnostic_only": "false",
            },
            {
                "setting": "current_min3",
                "genus": "Beta",
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "n_direct_species": "2",
                "eligible": "false",
                "diagnostic_only": "false",
            },
        ],
    )
    run329_direct = _write(
        tmp_path / "run329_direct.csv",
        [
            {
                "accepted_species": f"Beta {name}",
                "trait_name": "self_incompatibility",
                "source_lineages": f"url:https://example.org/b{name}",
            }
            for name in ("one", "two", "three")
        ],
    )
    run329_rules = _write(
        tmp_path / "run329_rules.csv",
        [
            {
                "setting": "current_min3",
                "genus": "Beta",
                "trait_name": "self_incompatibility",
                "inferred_state_set": '["SC"]',
                "n_direct_species": "3",
                "eligible": "true",
                "diagnostic_only": "false",
            }
        ],
    )

    summary = module.recover(
        coverage,
        wave53_direct,
        wave53_external,
        wave53_rules,
        tmp_path / "out",
        run329_direct_path=run329_direct,
        run329_rules_path=run329_rules,
    )

    assert summary["rule_tier_counts"] == {
        "wave53_current_min3": 1,
        "run329_formal_current_min3": 1,
        "unresolved_historical_secondary": 1,
    }
    assert summary["recovered_rule_lineages"] == 2
    assert summary["rules_with_direct_species_count_mismatch"] == 0
    assert summary["rules_with_zero_support_lineages"] == 0
    assert summary["unresolved_rule_cells"] == 1

    rules = pd.read_csv(tmp_path / "out" / "VALIDATED_LOW_DERIVED_RULE_SUMMARY.csv")
    alpha = rules.loc[rules["genus"].eq("Alpha")].iloc[0]
    assert alpha["reconstructed_n_direct_species"] == 2
    assert alpha["support_source_lineage_count"] == 2
    assert alpha["provenance_rule_tier"] == "wave53_current_min3"

    beta = rules.loc[rules["genus"].eq("Beta")].iloc[0]
    assert beta["reconstructed_n_direct_species"] == 3
    assert beta["support_source_lineage_count"] == 3
    assert beta["provenance_rule_tier"] == "run329_formal_current_min3"

    gamma = rules.loc[rules["genus"].eq("Gamma")].iloc[0]
    assert not bool(gamma["rule_recovered"])
    assert gamma["unmatched_wave53_reason"] == "no_genus_trait_in_wave53"
