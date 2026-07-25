from __future__ import annotations

import pandas as pd

from island_v2.incremental_genus_reaudit import incremental_genus_reaudit


def _evidence(species: str, value: str, quality: str = "medium") -> dict[str, str]:
    return {
        "accepted_species": species,
        "trait_name": "floral_symmetry",
        "normalized_value": value,
        "evidence_quality": quality,
        "evidence_scope": "species_direct",
    }


def test_reaudit_emits_only_genus_keys_enabled_by_new_direct_evidence() -> None:
    unresolved = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha target",
                "axis": "floral_structural_complexity",
            },
            {
                "accepted_species": "Beta target",
                "axis": "floral_structural_complexity",
            },
        ]
    )
    baseline = pd.DataFrame(
        [
            _evidence("Alpha one", "actinomorphic"),
            _evidence("Alpha two", "actinomorphic"),
            _evidence("Beta one", "zygomorphic"),
            _evidence("Beta two", "zygomorphic"),
            _evidence("Beta three", "zygomorphic"),
        ]
    )
    added = pd.DataFrame([_evidence("Alpha three", "actinomorphic")])

    incremental, rules, report = incremental_genus_reaudit(
        unresolved=unresolved,
        baseline_evidence=baseline,
        new_evidence=added,
    )

    assert set(zip(incremental["accepted_species"], incremental["axis"])) == {
        ("Alpha target", "floral_structural_complexity")
    }
    assert report["baseline_hypothetical_low_species_axis"] == 1
    assert report["augmented_hypothetical_low_species_axis"] == 2
    assert report["incremental_validated_low_species_axis"] == 1
    assert rules.loc[rules["genus"].eq("Alpha"), "eligible"].all()


def test_reaudit_accepts_integrated_quality_column() -> None:
    unresolved = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha target",
                "axis": "floral_structural_complexity",
            }
        ]
    )
    baseline = pd.DataFrame(
        [
            {
                **_evidence("Alpha one", "actinomorphic"),
                "quality": "medium",
                "evidence_quality": "",
            },
            {
                **_evidence("Alpha two", "actinomorphic"),
                "quality": "medium",
                "evidence_quality": "",
            },
        ]
    )
    added = pd.DataFrame([_evidence("Alpha three", "actinomorphic")])
    incremental, _, report = incremental_genus_reaudit(
        unresolved=unresolved,
        baseline_evidence=baseline,
        new_evidence=added,
    )
    assert len(incremental) == 1
    assert report["baseline_direct_rows"] == 2
