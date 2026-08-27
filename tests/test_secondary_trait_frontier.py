from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

import island_v2.secondary_trait_frontier as frontier_module
from island_v2.secondary_trait_frontier import AXES, build_secondary_coverage


@pytest.fixture(autouse=True)
def small_universe(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(frontier_module, "EXPECTED_SPECIES", 4)
    monkeypatch.setattr(frontier_module, "EXPECTED_CELLS", 12)


def _write_fixture(root: Path, *, wrong_axis: bool = False) -> tuple[Path, Path]:
    baseline_dir = root / "baseline"
    frontier_dir = root / "frontier"
    baseline_dir.mkdir()
    frontier_dir.mkdir()
    species = [f"Genus species{i:06d}" for i in range(4)]
    baseline = pd.DataFrame(
        [(name, axis, "", "", "", "", "") for name in species for axis in AXES],
        columns=[
            "accepted_species",
            "axis",
            "trait_composition",
            "trait_names",
            "source_groups",
            "source_lineages",
            "quality",
        ],
    )
    baseline.loc[
        (baseline["accepted_species"] == species[0]) & (baseline["axis"] == "flower_colour"),
        ["trait_composition", "trait_names", "source_groups", "source_lineages", "quality"],
    ] = [
        'flower_primary_color=["white"]',
        "flower_primary_color",
        "fixture_direct",
        "doi:fixture",
        "high",
    ]
    baseline.to_csv(baseline_dir / "secondary_species_axis_coverage.csv.gz", index=False)

    axis = "reproductive_assurance" if wrong_axis else "flower_colour"
    trait = pd.DataFrame(
        [
            {
                "accepted_species": species[1],
                "axis": axis,
                "genus": "Genus",
                "trait_name": "flower_primary_color",
                "predicted_state_set": '["white"]',
                "n_direct_species": "3",
                "n_source_lineages": "3",
                "species_loo_accuracy": "1.0",
                "lineage_loo_accuracy": "1.0",
                "support_source_lineages": '["doi:a", "doi:b", "doi:c"]',
            }
        ]
    )
    trait.to_csv(frontier_dir / "eligible_candidate_species_trait.csv.gz", index=False)
    axis_frame = pd.DataFrame(
        [
            {
                "accepted_species": species[1],
                "axis": axis,
                "genus": "Genus",
                "trait_names": '["flower_primary_color"]',
                "predicted_state_sets": json.dumps(['["white"]']),
                "support_source_lineages": json.dumps(
                    ['["doi:a", "doi:b", "doi:c"]']
                ),
                "quality": "low",
                "analysis_tier": "secondary_robustness",
                "family_inference": "False",
                "global_fallback": "False",
            }
        ]
    )
    axis_frame.to_csv(frontier_dir / "eligible_candidate_species_axis.csv.gz", index=False)
    rules = pd.DataFrame(
        [
            {
                "genus": "Genus",
                "axis": axis,
                "trait_name": "flower_primary_color",
                "eligible": "True",
                "n_direct_species": "3",
                "species_loo_accuracy": "1.0",
                "lineage_loo_accuracy": "1.0",
            }
        ]
    )
    rules.to_csv(frontier_dir / "trait_specific_rule_frontier.csv.gz", index=False)
    pd.DataFrame([{"genus": "Genus", "trait_name": "flower_primary_color"}]).to_csv(
        frontier_dir / "prioritized_direct_acquisition_queue.csv.gz", index=False
    )
    return baseline_dir, frontier_dir


def test_materializes_only_unresolved_trait_specific_low(tmp_path: Path) -> None:
    baseline, frontier = _write_fixture(tmp_path)
    output = tmp_path / "output"
    summary = build_secondary_coverage(
        baseline_dir=baseline, frontier_dir=frontier, output_dir=output
    )

    assert summary["delta"] == {
        "gross_gain_species_axis": 1,
        "loss_species_axis": 0,
        "net_gain_species_axis": 1,
        "by_axis_net_gain": {
            "flower_colour": 1,
            "floral_structural_complexity": 0,
            "reproductive_assurance": 0,
        },
    }
    result = pd.read_csv(output / "wave34_secondary_species_axis_coverage.csv.gz").fillna("")
    original = result.loc[
        (result["accepted_species"] == "Genus species000000")
        & (result["axis"] == "flower_colour")
    ].iloc[0]
    added = result.loc[
        (result["accepted_species"] == "Genus species000001")
        & (result["axis"] == "flower_colour")
    ].iloc[0]
    assert original["quality"] == "high"
    assert original["source_lineages"] == "doi:fixture"
    assert added["quality"] == "low"
    assert added["trait_names"] == "flower_primary_color"
    assert added["source_lineages"] == "doi:a|doi:b|doi:c"


def test_rejects_genus_axis_join_that_drops_trait_identity(tmp_path: Path) -> None:
    baseline, frontier = _write_fixture(tmp_path, wrong_axis=True)
    with pytest.raises(ValueError, match="cross-trait or cross-axis"):
        build_secondary_coverage(
            baseline_dir=baseline,
            frontier_dir=frontier,
            output_dir=tmp_path / "output",
        )
