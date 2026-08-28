from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from island_v2 import all_evidence_trait_audit as audit
from island_v2.direct_evidence_exclusions import apply_direct_evidence_exclusions
from island_v2.integrated_trait_coverage import EVIDENCE_COLUMNS

ONTOLOGY = {
    "flower_primary_color": {
        "white",
        "red_pink",
        "yellow_orange",
        "blue_purple",
        "green_brown_inconspicuous",
        "multicolored_variable",
        "other_described",
    },
    "floral_form": {"tubular", "open_radial"},
    "floral_symmetry": {"actinomorphic", "zygomorphic"},
    "self_incompatibility": {"SI", "SC", "mixed_or_variable"},
    "autonomous_selfing_capacity": {"absent", "autonomous"},
    "mating_system": {
        "predominantly_outcrossing",
        "mixed_mating",
        "predominantly_selfing",
        "obligate_selfing",
    },
}


def row(
    species: str,
    trait: str,
    value: str,
    lineage: str,
    *,
    quality: str = "medium",
    source_group: str = "example",
    excerpt: str = "",
) -> dict[str, str]:
    axis = audit.trait_axis(trait)
    return {
        "accepted_species": species,
        "axis": axis,
        "trait_name": trait,
        "normalized_value": value,
        "quality": quality,
        "source_group": source_group,
        "source_provider": source_group,
        "source_url": f"https://example.org/{lineage}",
        "source_record_id": lineage,
        "source_citation": "Example source",
        "source_excerpt": excerpt or f"{trait}={value}",
        "evidence_scope": "species_direct",
        "name_match_method": "accepted_name_exact",
        "source_lineage": lineage,
        "lineage_method": "test",
        "source_run_id": "1",
        "source_artifact": "test",
        "source_file": "test.csv",
        "acceptance_contract": "test",
    }


def test_latest_public_web_loader_preserves_original_artifact_provenance(
    tmp_path: Path,
) -> None:
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "evidence_quality": "medium",
                "evidence_scope": "species_direct",
                "source_provider": "prior.example",
                "source_url": "https://prior.example/alpha",
                "source_record_id": "alpha",
                "source_citation": "Prior source",
                "source_excerpt": "Flowers tubular.",
                "name_match_method": "accepted_name_exact",
                "source_lineage": "page:alpha",
                "inference_rule": "",
                "source_run_id": "prior-run",
                "source_artifact": "prior-artifact",
            }
        ]
    ).to_csv(tmp_path / "broad_web_medium_evidence.csv.gz", index=False)
    manifest = {
        "sources": [
            {
                "source_group": "latest_public_web",
                "run_id": "wrapper-run",
                "artifact_name": "wrapper-artifact",
            }
        ]
    }
    loaded = audit.load_latest_public_web(tmp_path, manifest)
    assert len(loaded) == 1
    assert loaded.iloc[0]["source_run_id"] == "prior-run"
    assert loaded.iloc[0]["source_artifact"] == "prior-artifact"


def test_latest_public_web_loader_coalesces_row_level_quality_columns(
    tmp_path: Path,
) -> None:
    pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "trait_name": "floral_form",
                "normalized_value": "tubular",
                "evidence_quality": "",
                "quality": "medium",
                "evidence_scope": "species_direct",
                "source_lineage": "page:legacy-quality",
            },
            {
                "accepted_species": "Beta two",
                "trait_name": "floral_symmetry",
                "normalized_value": "actinomorphic",
                "evidence_quality": "high",
                "quality": "",
                "evidence_scope": "species_direct",
                "source_lineage": "page:new-evidence-quality",
            },
        ]
    ).to_csv(tmp_path / "broad_web_medium_evidence.csv.gz", index=False)

    loaded = audit.load_latest_public_web(tmp_path, {"sources": []})

    assert len(loaded) == 2
    assert dict(zip(loaded["accepted_species"], loaded["quality"], strict=True)) == {
        "Alpha one": "medium",
        "Beta two": "high",
    }


def test_reviewed_direct_supplement_loader_preserves_contract_and_trait(
    tmp_path: Path,
) -> None:
    record = row(
        "Alpha one",
        "autonomous_selfing_capacity",
        '["autonomous"]',
        "dataset:reviewed",
        quality="high",
    )
    path = tmp_path / "reviewed.csv.gz"
    pd.DataFrame([record]).to_csv(path, index=False)
    loaded = audit.load_reviewed_direct_supplements((path,))

    assert len(loaded) == 1
    assert loaded.iloc[0]["trait_name"] == "autonomous_selfing_capacity"
    assert loaded.iloc[0]["axis"] == "reproductive_assurance"
    assert loaded.iloc[0]["acceptance_contract"] == "test"
    assert loaded.iloc[0]["source_file"] == str(path)


def test_reviewed_direct_supplement_loader_rejects_incomplete_provenance(
    tmp_path: Path,
) -> None:
    record = row("Alpha one", "floral_form", "tubular", "dataset:reviewed")
    record["source_citation"] = ""
    path = tmp_path / "reviewed.csv.gz"
    pd.DataFrame([record]).to_csv(path, index=False)

    with pytest.raises(ValueError, match="incomplete provenance"):
        audit.load_reviewed_direct_supplements((path,))


def test_reviewed_exclusion_prevents_false_conflict_with_corrected_direct_row() -> None:
    wrong = row(
        "Calanthe striata",
        "autonomous_selfing_capacity",
        "autonomous",
        "url:https://europepmc.org/article/AGR/IND500728652",
        quality="high",
    )
    corrected = row(
        "Calanthe striata",
        "autonomous_selfing_capacity",
        "absent",
        "doi:10.3390/horticulturae10101025",
        quality="high",
    )
    exclusions = pd.DataFrame(
        [
            {
                "accepted_species": "Calanthe striata",
                "trait_name": "autonomous_selfing_capacity",
                "normalized_value": "autonomous",
                "source_lineage": "url:https://europepmc.org/article/AGR/IND500728652",
                "reason": "saved excerpt says neither autogamous nor apogamous",
                "reviewer": "reviewer",
                "reviewed_at_utc": "2026-08-13T00:00:00Z",
            }
        ]
    )

    kept, exclusion_audit = apply_direct_evidence_exclusions(
        pd.DataFrame([wrong, corrected], columns=EVIDENCE_COLUMNS),
        exclusions,
    )
    lineages, _ = audit.dedupe_direct_lineages(kept, ONTOLOGY)
    resolved, cell_audit = audit.resolve_direct_cells(lineages)

    assert exclusion_audit.iloc[0]["matched_rows"] == 1
    assert resolved.iloc[0]["normalized_value"] == "absent"
    assert resolved.iloc[0]["resolution_status"] == "resolved"
    assert cell_audit.iloc[0]["classification"] == "single_independent_lineage"


def test_integrated_workflow_reapplies_reviewed_direct_exclusions() -> None:
    workflow = Path(".github/workflows/build-integrated-trait-coverage.yml").read_text(
        encoding="utf-8"
    )
    assert "--direct-evidence-exclusions-csv" in workflow
    assert "direct_evidence_exclusions_20260811.csv" in workflow


def test_lineage_dedup_is_trait_specific_not_axis_specific() -> None:
    evidence = pd.DataFrame(
        [
            row("Alpha one", "floral_form", "tubular", "page:1"),
            row("Alpha one", "floral_symmetry", "actinomorphic", "page:1"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, duplicates = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    assert len(lineages) == 2
    assert duplicates.empty
    assert set(lineages["trait_name"]) == {"floral_form", "floral_symmetry"}


def test_direct_conflicts_are_classified_without_row_order_selection() -> None:
    evidence = pd.DataFrame(
        [
            row("Colour one", "flower_primary_color", "red_pink", "colour:a"),
            row("Colour one", "flower_primary_color", "white", "colour:b"),
            row("Repro one", "self_incompatibility", "SC", "repro:a"),
            row("Repro one", "self_incompatibility", "SI", "repro:b"),
            row("Garden one", "floral_form", "tubular", "garden:wild"),
            row(
                "Garden one",
                "floral_form",
                "open_radial",
                "garden:cultivar",
                excerpt="cultivar selection has open radial flowers",
            ),
            row("Ontology one", "floral_form", "not_an_ontology_value", "bad:a"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    resolved, cell_audit = audit.resolve_direct_cells(lineages)

    colour = resolved.loc[
        resolved["accepted_species"].eq("Colour one")
    ].iloc[0]
    assert colour["classification"] == "true_multistate_variable"
    assert json.loads(colour["state_set"]) == ["red_pink", "white"]

    repro = cell_audit.loc[
        cell_audit["accepted_species"].eq("Repro one")
    ].iloc[0]
    assert repro["classification"] == "unresolved_direct_conflict"
    assert repro["resolution_status"] == "excluded"

    garden = resolved.loc[
        resolved["accepted_species"].eq("Garden one")
    ].iloc[0]
    assert garden["classification"] == "cultivar_contamination"
    assert garden["normalized_value"] == "tubular"

    ontology = cell_audit.loc[
        cell_audit["accepted_species"].eq("Ontology one")
    ].iloc[0]
    assert ontology["classification"] == "source_ontology_mismatch"


def test_explicit_mixed_reproductive_state_resolves_simple_state_conflicts() -> None:
    evidence = pd.DataFrame(
        [
            row(
                "Capsicum pubescens",
                "self_incompatibility",
                "SI",
                "paper:si",
                quality="high",
            ),
            row(
                "Capsicum pubescens",
                "self_incompatibility",
                "mixed_or_variable",
                "dataset:si-sc",
                quality="high",
            ),
            row(
                "Witheringia solanacea",
                "self_incompatibility",
                "SC",
                "paper:sc",
                quality="high",
            ),
            row(
                "Witheringia solanacea",
                "self_incompatibility",
                "SI",
                "paper:si",
                quality="high",
            ),
            row(
                "Witheringia solanacea",
                "self_incompatibility",
                "mixed_or_variable",
                "dataset:si-sc",
                quality="high",
            ),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)

    resolved, cell_audit = audit.resolve_direct_cells(lineages)

    assert set(resolved["accepted_species"]) == {
        "Capsicum pubescens",
        "Witheringia solanacea",
    }
    assert set(resolved["classification"]) == {"true_multistate_variable"}
    assert set(resolved["normalized_value"]) == {"mixed_or_variable"}
    assert set(resolved["state_set"].map(json.loads).map(tuple)) == {
        ("mixed_or_variable",)
    }
    assert set(cell_audit["resolution_status"]) == {"resolved"}


def test_legacy_colour_aliases_remain_multistate_and_internal_conflicts_are_explicit() -> None:
    states, invalid = audit.normalise_state_set(
        "flower_primary_color",
        "cream|violet|black|multicolored/variable",
        ONTOLOGY,
    )
    assert states == (
        "blue_purple",
        "multicolored_variable",
        "other_described",
        "white",
    )
    assert invalid == ()

    evidence = pd.DataFrame(
        [
            row("Delta one", "self_incompatibility", "SC", "same:lineage"),
            row("Delta one", "self_incompatibility", "SI", "same:lineage"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    resolved, cell_audit = audit.resolve_direct_cells(lineages)
    assert resolved.empty
    assert cell_audit.iloc[0]["classification"] == "unresolved_direct_conflict"
    assert cell_audit.iloc[0]["lineage_internal_conflicts"] == 1


def test_thresholds_and_application_remain_trait_specific() -> None:
    evidence_rows = [
        row(f"Beta {index}", "floral_form", value, f"beta:{index}")
        for index, value in enumerate(
            ["tubular", "tubular", "tubular", "open_radial"],
            start=1,
        )
    ]
    evidence_rows += [
        row(f"Gamma {index}", "floral_symmetry", "actinomorphic", f"gamma:{index}")
        for index in range(1, 3)
    ]
    evidence = pd.DataFrame(evidence_rows, columns=EVIDENCE_COLUMNS)
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    cells, _ = audit.resolve_direct_cells(lineages)
    cells["genus"] = cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    old_low = pd.DataFrame(
        columns=["accepted_species", "genus", "axis", "trait_name", "state_set"]
    )
    rules = audit.build_rule_audit(cells, lineages, old_low)

    beta_current = rules.loc[
        rules["setting"].eq("current_min3")
        & rules["genus"].eq("Beta")
        & rules["trait_name"].eq("floral_form")
    ].iloc[0]
    beta_relaxed = rules.loc[
        rules["setting"].eq("relaxed_min3")
        & rules["genus"].eq("Beta")
        & rules["trait_name"].eq("floral_form")
    ].iloc[0]
    assert not bool(beta_current["eligible"])
    assert bool(beta_relaxed["eligible"])

    gamma_min3 = rules.loc[
        rules["setting"].eq("current_min3")
        & rules["genus"].eq("Gamma")
    ].iloc[0]
    gamma_min2 = rules.loc[
        rules["setting"].eq("current_min2_diagnostic")
        & rules["genus"].eq("Gamma")
    ].iloc[0]
    assert not bool(gamma_min3["eligible"])
    assert bool(gamma_min2["eligible"])
    assert bool(gamma_min2["diagnostic_only"])

    master = pd.DataFrame(
        [
            {
                "accepted_species": "Gamma target",
                "genus": "Gamma",
                "n_islands": "3",
                "n_records": "10",
            },
            {
                "accepted_species": "Gamma",
                "genus": "Gamma",
                "n_islands": "3",
                "n_records": "10",
            },
        ]
    )
    inferred = audit.apply_genus_rules(
        master,
        cells,
        rules,
        "current_min2_diagnostic",
    )
    assert len(inferred) == 1
    assert inferred.iloc[0]["accepted_species"] == "Gamma target"
    assert inferred.iloc[0]["trait_name"] == "floral_symmetry"


def test_lineage_loo_removes_every_species_from_the_held_out_lineage() -> None:
    evidence = pd.DataFrame(
        [
            row("Delta one", "self_incompatibility", "SC", "paper:shared"),
            row("Delta two", "self_incompatibility", "SC", "paper:shared"),
            row("Delta three", "self_incompatibility", "SC", "paper:independent"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    cells, _ = audit.resolve_direct_cells(lineages)
    cells["genus"] = cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    old_low = pd.DataFrame(
        columns=["accepted_species", "genus", "axis", "trait_name", "state_set"]
    )

    rules = audit.build_rule_audit(cells, lineages, old_low)
    current = rules.loc[rules["setting"].eq("current_min3")].iloc[0]

    assert current["n_direct_species"] == 3
    assert current["lineage_loo_n"] == 2
    assert current["lineage_loo_correct"] == 2
    assert current["lineage_loo_accuracy"] == 1.0


def test_genus_rules_join_on_accepted_species_genus_not_stale_master_column() -> None:
    master = pd.DataFrame(
        [
            {
                "accepted_species": "Deyeuxia target",
                "genus": "Calamagrostis",
            },
            {
                "accepted_species": "Calamagrostis target",
                "genus": "Deyeuxia",
            },
        ]
    )
    rules = pd.DataFrame(
        [
            {
                "genus": "Calamagrostis",
                "axis": "reproductive_assurance",
                "trait_name": "mating_system",
                "setting": "current_min3",
                "eligible": True,
                "inferred_state_set": '["predominantly_outcrossing"]',
                "inferred_value": "predominantly_outcrossing",
            }
        ]
    )
    direct = pd.DataFrame(columns=["accepted_species", "trait_name"])

    inferred = audit.apply_genus_rules(master, direct, rules, "current_min3")

    assert inferred["accepted_species"].tolist() == ["Calamagrostis target"]
    assert inferred["genus"].tolist() == ["Calamagrostis"]


def test_lineage_loo_does_not_count_same_paper_species_as_independent_validation() -> None:
    evidence = pd.DataFrame(
        [
            row("Epsilon one", "self_incompatibility", "SC", "paper:shared"),
            row("Epsilon two", "self_incompatibility", "SC", "paper:shared"),
            row("Epsilon three", "self_incompatibility", "SI", "paper:independent"),
        ],
        columns=EVIDENCE_COLUMNS,
    )
    lineages, _ = audit.dedupe_direct_lineages(evidence, ONTOLOGY)
    cells, _ = audit.resolve_direct_cells(lineages)
    cells["genus"] = cells["accepted_species"].str.split().str[0]
    lineages["genus"] = lineages["accepted_species"].str.split().str[0]
    old_low = pd.DataFrame(
        columns=["accepted_species", "genus", "axis", "trait_name", "state_set"]
    )

    rules = audit.build_rule_audit(cells, lineages, old_low)
    current = rules.loc[rules["setting"].eq("current_min3")].iloc[0]

    assert current["lineage_loo_n"] == 2
    assert current["lineage_loo_correct"] == 0
    assert current["lineage_loo_accuracy"] == 0.0
    assert not bool(current["eligible"])


def test_species_axis_coverage_uses_quality_precedence_and_exact_denominator() -> None:
    master = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "genus": "Alpha",
                "n_islands": "1",
                "n_records": "1",
            },
            {
                "accepted_species": "Beta two",
                "genus": "Beta",
                "n_islands": "1",
                "n_records": "1",
            },
        ]
    )
    direct = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "axis": "flower_colour",
                "trait_name": "flower_primary_color",
                "state_set": '["white"]',
                "normalized_value": "white",
                "quality": "high",
                "source_groups": "direct",
                "source_lineages": "lineage:a",
                "classification": "single_independent_lineage",
            }
        ]
    )
    low = pd.DataFrame(
        [
            {
                "accepted_species": "Alpha one",
                "axis": "flower_colour",
                "trait_name": "flower_primary_color",
                "state_set": '["red_pink"]',
                "normalized_value": "red_pink",
                "source_lineage": "rule:a",
            },
            {
                "accepted_species": "Beta two",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "state_set": '["SC"]',
                "normalized_value": "SC",
                "source_lineage": "rule:b",
            },
        ]
    )
    coverage = audit.species_axis_coverage(master, direct, low)
    snapshot = audit.coverage_snapshot(coverage)
    assert len(coverage) == 6
    assert snapshot["quality_counts"] == {"high": 1, "medium": 0, "low": 1}
    assert snapshot["unresolved_species_axis"] == 4


def test_acquisition_queue_ranks_cells_unlocked_per_required_acquisition() -> None:
    master = pd.DataFrame(
        [
            {
                "accepted_species": f"Alpha {index}",
                "genus": "Alpha",
                "n_islands": "1",
                "n_records": "1",
            }
            for index in range(100)
        ]
        + [
            {
                "accepted_species": f"Beta {index}",
                "genus": "Beta",
                "n_islands": "1",
                "n_records": "1",
            }
            for index in range(80)
        ]
    )
    coverage = master[["accepted_species"]].copy()
    coverage["axis"] = "reproductive_assurance"
    coverage["quality"] = ""
    rules = pd.DataFrame(
        [
            {
                "setting": "current_min3",
                "genus": "Alpha",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "n_direct_species": 1,
                "dominant_species": 1,
                "required_dominance": 0.95,
                "dominance": 1.0,
                "value_distribution": '{"[\\"SC\\"]":1}',
            },
            {
                "setting": "current_min3",
                "genus": "Beta",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "n_direct_species": 2,
                "dominant_species": 2,
                "required_dominance": 0.95,
                "dominance": 1.0,
                "value_distribution": '{"[\\"SC\\"]":2}',
            },
        ]
    )
    queue = audit.acquisition_queue(master, coverage, rules)
    assert queue.iloc[0]["genus"] == "Beta"
    assert queue.iloc[0]["acquisitions_needed_for_min3_dominance"] == 1
    assert queue.iloc[0]["potential_cells_per_acquisition"] == 80


def test_acquisition_queue_derives_genus_from_accepted_species() -> None:
    master = pd.DataFrame(
        [
            {
                "accepted_species": "Accepted alpha",
                "genus": "StaleSynonymGenus",
                "n_islands": "3",
                "n_records": "10",
            }
        ]
    )
    coverage = pd.DataFrame(
        [
            {
                "accepted_species": "Accepted alpha",
                "axis": "reproductive_assurance",
                "quality": "",
            }
        ]
    )
    rules = pd.DataFrame(
        [
            {
                "setting": "current_min3",
                "genus": "Accepted",
                "axis": "reproductive_assurance",
                "trait_name": "self_incompatibility",
                "n_direct_species": 2,
                "dominant_species": 2,
                "required_dominance": 0.95,
                "dominance": 1.0,
                "value_distribution": '{"[\\"SC\\"]":2}',
            }
        ]
    )

    queue = audit.acquisition_queue(master, coverage, rules)

    assert len(queue) == 1
    assert queue.iloc[0]["genus"] == "Accepted"
    assert queue.iloc[0]["potential_cells_unlocked"] == 1
