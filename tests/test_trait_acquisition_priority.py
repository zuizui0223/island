from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from island_v2 import global_trait_campaign as campaign
from island_v2 import trait_acquisition_priority as priority


def test_repository_trait_priority_contract_has_one_order_and_excludes_guild() -> None:
    acquisition = yaml.safe_load(Path("config/v2_trait_acquisition.yml").read_text("utf-8"))
    cascade = yaml.safe_load(Path("config/trait_fill_cascade.yml").read_text("utf-8"))

    assert tuple(acquisition["priority_order"]) == priority.TARGET_TRAITS
    assert tuple(cascade["target_traits"]) == priority.TARGET_TRAITS
    assert "pollination_functional_guild" not in priority.TARGET_TRAITS
    assert acquisition["per_trait_completion_order"] == [
        "species_direct",
        "qualified_genus_inference",
        "qualified_family_inference",
        "unresolved_no_evidence",
    ]


def _master() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": [
                "Alpha one",
                "Alpha two",
                "Beta one",
                "Delta one",
                "Gamma one",
                "Fern one",
                "Unknown one",
            ],
            "genus": ["Alpha", "Alpha", "Beta", "Delta", "Gamma", "Fern", "Unknown"],
            "family": ["Fam"] * 7,
            "n_islands": [1] * 7,
            "n_records": [1] * 7,
        }
    )


def _geography_inputs() -> dict[str, pd.DataFrame]:
    occurrences = pd.DataFrame(
        {
            "species": [
                "Alpha one",
                "Alpha two",
                "Alpha two",
                "Beta one",
                "Delta one",
                "Gamma one",
                "Unknown one",
            ],
            "island_id": [
                "north",
                "north",
                "south",
                "south",
                "name_trap",
                "unreviewed",
                "tropical",
            ],
        }
    )
    canonical_latitudes = pd.DataFrame(
        {
            "island_id": ["north", "south", "name_trap", "unreviewed", "tropical"],
            "canonical_centroid_lat": ["35", "-20", "30", "36", "10"],
            "geometry_sha256": ["geometry-hash"] * 5,
            "latitude_method": ["equal-area polygon centroid"] * 5,
            "source_url": ["https://example.test/islands"] * 5,
            "source_file_sha256": ["source-hash"] * 5,
        }
    )
    island_registry = pd.DataFrame(
        {
            "island_id": ["north", "south", "name_trap", "unreviewed", "tropical"],
            "source_region_id": [
                "region_north",
                "region_south",
                "temperate_only_in_identifier",
                "region_north",
                "region_tropical",
            ],
            "centroid_lat": ["35", "-20", "30", "36", "10"],
            "identity_review_status": ["accepted"] * 5,
            "identity_source_citation": ["Island coordinates"] * 5,
            "identity_source_url": ["https://example.test/island"] * 5,
        }
    )
    island_assignments = pd.DataFrame(
        {
            "island_id": ["north", "south", "name_trap", "unreviewed", "tropical"],
            "source_region_id": [
                "region_north",
                "region_south",
                "temperate_only_in_identifier",
                "region_north",
                "region_tropical",
            ],
            "assignment_evidence_id": [
                "a_north",
                "a_south",
                "a_trap",
                "a_pending",
                "a_tropical",
            ],
            "review_status": ["accepted", "accepted", "accepted", "pending", "accepted"],
        }
    )
    assignment_evidence = pd.DataFrame(
        {
            "assignment_evidence_id": [
                "a_north",
                "a_south",
                "a_trap",
                "a_pending",
                "a_tropical",
            ],
            "island_id": ["north", "south", "name_trap", "unreviewed", "tropical"],
            "source_citation": ["Assignment evidence"] * 5,
            "source_url": ["https://example.test/assignment"] * 5,
            "review_status": ["accepted", "accepted", "accepted", "pending", "accepted"],
        }
    )
    region_registry = pd.DataFrame(
        {
            "source_region_id": [
                "region_north",
                "region_south",
                "temperate_only_in_identifier",
                "region_tropical",
            ],
            "biogeographic_region": [
                "Palearctic_temperate",
                "temperate",
                "tropical",
                "tropical",
            ],
            "review_status": ["accepted"] * 4,
        }
    )
    return {
        "occurrences": occurrences,
        "canonical_island_latitudes": canonical_latitudes,
        "island_registry": island_registry,
        "island_assignments": island_assignments,
        "island_assignment_evidence": assignment_evidence,
        "source_region_registry": region_registry,
        # Legacy callers may still pass this table.  Its schema and values have
        # no effect on geographic priority classification.
        "source_region_evidence": pd.DataFrame(),
    }


def test_geography_priority_uses_reviewed_geography_not_bombus_fields() -> None:
    audit = priority.build_geography_audit(_master(), **_geography_inputs()).set_index(
        "accepted_species"
    )

    assert audit.loc["Alpha one", "any_reviewed_NH_temperate"] == "true"
    assert (
        audit.loc["Alpha one", "scheduling_region_category"]
        == "reviewed_nh_temperate"
    )
    assert audit.loc["Alpha one", "priority_reviewed_NH_temperate_island_count"] == 1
    assert audit.loc["Alpha one", "priority_assignment_evidence_ids"] == "a_north"
    assert not any("bombus" in column.casefold() for column in audit.columns)
    assert audit.loc["Beta one", "any_reviewed_NH_temperate"] == "false"
    assert audit.loc["Beta one", "any_SH"] == "true"
    assert audit.loc["Beta one", "scheduling_region_category"] == "southern_hemisphere"
    # The source-region identifier contains "temperate", but the reviewed field
    # says tropical; identifiers are never used as geographic evidence.
    assert audit.loc["Delta one", "any_reviewed_NH_temperate"] == "false"
    assert audit.loc["Delta one", "scheduling_region_category"] == "nh_temperate_unreviewed"
    assert audit.loc["Gamma one", "any_reviewed_NH_temperate"] == "false"
    assert audit.loc["Gamma one", "any_NH_temperate_unreviewed"] == "true"
    assert (
        audit.loc["Gamma one", "scheduling_region_category"]
        == "nh_temperate_unreviewed"
    )
    assert audit.loc["Gamma one", "priority_reviewed_NH_temperate_island_count"] == 0
    assert audit.loc["Unknown one", "scheduling_region_category"] == "nh_tropical"
    assert audit.loc["Alpha two", "scheduling_region_category"] == "cross_region"
    assert audit.loc["Alpha two", "crosses_regions"] == "true"
    assert audit.loc["Fern one", "scheduling_region_category"] == "unclassified"


def test_latitude_bands_survive_missing_optional_review_inputs() -> None:
    inputs = _geography_inputs()
    audit = priority.build_geography_audit(
        _master(),
        occurrences=inputs["occurrences"],
        canonical_island_latitudes=inputs["canonical_island_latitudes"],
        island_registry=None,
        island_assignments=None,
        island_assignment_evidence=None,
        source_region_registry=None,
        source_region_evidence=None,
    ).set_index("accepted_species")

    assert audit.loc["Alpha one", "any_reviewed_NH_temperate"] == "unknown"
    assert audit.loc["Alpha one", "scheduling_region_category"] == "nh_temperate_unreviewed"
    assert audit.loc["Unknown one", "scheduling_region_category"] == "nh_tropical"
    assert audit.loc["Beta one", "scheduling_region_category"] == "southern_hemisphere"
    assert audit.loc["Alpha two", "scheduling_region_category"] == "cross_region"
    assert pd.isna(audit.loc["Alpha one", "priority_reviewed_NH_temperate_island_count"])


def test_repository_audit_does_not_gate_latitude_on_optional_review_files(
    tmp_path: Path,
) -> None:
    inputs = _geography_inputs()
    occurrence_csv = tmp_path / "occurrences.csv"
    latitude_csv = tmp_path / "latitudes.csv"
    inputs["occurrences"].to_csv(occurrence_csv, index=False)
    inputs["canonical_island_latitudes"].to_csv(latitude_csv, index=False)

    audit = priority.build_repository_priority_audit(
        _master(),
        master_genus_is_explicit=True,
        occurrence_csv=occurrence_csv,
        canonical_island_latitudes_csv=latitude_csv,
        island_registry_csv=tmp_path / "missing-island-registry.csv",
        island_assignments_csv=tmp_path / "missing-island-assignments.csv",
        island_assignment_evidence_csv=tmp_path / "missing-assignment-evidence.csv",
        source_region_registry_csv=tmp_path / "missing-source-region-registry.csv",
        source_region_evidence_csv=tmp_path / "missing-legacy-region-evidence.csv",
        acquisition_config_path=tmp_path / "missing-acquisition.yml",
        qualified_genus_glob=str(tmp_path / "missing-qualified" / "*.csv"),
        frozen_packet_glob=str(tmp_path / "missing-frozen" / "*.json"),
    ).set_index("accepted_species")

    assert audit.loc["Alpha one", "scheduling_region_category"] == "nh_temperate_unreviewed"
    assert audit.loc["Unknown one", "scheduling_region_category"] == "nh_tropical"
    assert audit.loc["Beta one", "scheduling_region_category"] == "southern_hemisphere"


def test_arctic_cooccurrence_keeps_independent_nh_temperate_priority() -> None:
    inputs = _geography_inputs()
    inputs["occurrences"] = pd.concat(
        [
            inputs["occurrences"],
            pd.DataFrame(
                {
                    "species": ["Alpha one", "Gamma one"],
                    "island_id": ["arctic", "arctic"],
                }
            ),
        ],
        ignore_index=True,
    )
    inputs["canonical_island_latitudes"] = pd.concat(
        [
            inputs["canonical_island_latitudes"],
            pd.DataFrame(
                {
                    "island_id": ["arctic"],
                    "canonical_centroid_lat": ["75"],
                    "geometry_sha256": ["arctic-geometry-hash"],
                    "latitude_method": ["equal-area polygon centroid"],
                    "source_url": ["https://example.test/islands"],
                    "source_file_sha256": ["source-hash"],
                }
            ),
        ],
        ignore_index=True,
    )
    audit = priority.build_geography_audit(_master(), **inputs).set_index(
        "accepted_species"
    )

    assert audit.loc["Alpha one", "scheduling_region_category"] == "unclassified"
    assert audit.loc["Alpha one", "any_reviewed_NH_temperate"] == "true"
    assert audit.loc["Gamma one", "scheduling_region_category"] == "unclassified"
    assert audit.loc["Gamma one", "any_NH_temperate_unreviewed"] == "true"

    ranked, _ = priority.add_priority_sort_keys(audit.reset_index())
    geo_rank = ranked.set_index("accepted_species")["_priority_geo"]
    assert geo_rank.loc["Alpha one"] == 0
    assert geo_rank.loc["Gamma one"] == 1


def _snapshot() -> priority.TraitCoverageSnapshot:
    known = frozenset(
        {"Alpha one", "Alpha two", "Beta one", "Delta one", "Gamma one", "Fern one"}
    )
    return priority.TraitCoverageSnapshot(
        known_species=known,
        applicable_species=known - {"Fern one"},
        source_backed_pairs=frozenset(
            {
                ("Alpha one", "flower_primary_color"),
                ("Beta one", "self_incompatibility"),
            }
        ),
        target_traits=priority.TARGET_TRAITS,
        basis="synthetic accepted direct-cell snapshot",
    )


def test_yield_features_keep_unknown_na_and_count_only_qualified_direct_gaps() -> None:
    master = _master()
    snapshot = _snapshot()
    gaps = priority.build_trait_gap_audit(master, snapshot).set_index("accepted_species")
    assert gaps.loc["Alpha one", "unresolved_trait_cells"] == len(priority.TARGET_TRAITS) - 1
    assert gaps.loc["Alpha two", "unresolved_trait_cells"] == len(priority.TARGET_TRAITS)
    assert gaps.loc["Alpha one", "highest_priority_unresolved_trait"] == "floral_form"
    assert gaps.loc["Alpha one", "highest_priority_unresolved_rank"] == 2
    assert gaps.loc["Alpha two", "highest_priority_unresolved_trait"] == "flower_primary_color"
    assert gaps.loc["Alpha two", "highest_priority_unresolved_rank"] == 1
    assert (
        gaps.loc["Alpha two", "weighted_unresolved_trait_score"]
        > gaps.loc["Alpha one", "weighted_unresolved_trait_score"]
    )
    assert gaps.loc["Fern one", "unresolved_trait_cells_state"] == "not_applicable"
    assert pd.isna(gaps.loc["Fern one", "unresolved_trait_cells"])
    assert gaps.loc["Unknown one", "unresolved_trait_cells_state"] == "unknown"
    assert pd.isna(gaps.loc["Unknown one", "unresolved_trait_cells"])

    ranked, sort_columns = priority.add_priority_sort_keys(gaps.reset_index())
    assert sort_columns[:5] == [
        "_priority_trait_rank_unknown",
        "_priority_trait_rank",
        "_priority_weighted_gap_unknown",
        "_priority_weighted_gap",
        "_priority_unresolved_unknown",
    ]
    ordered = ranked.sort_values(sort_columns)["accepted_species"].tolist()
    assert ordered.index("Alpha two") < ordered.index("Alpha one")

    qualified = pd.DataFrame(
        {
            "accepted_species": ["Alpha one", "Alpha two", "Beta one"],
            "genus": ["Alpha", "Alpha", "Beta"],
            "trait_name": ["self_incompatibility"] * 3,
            "candidate_kind": ["hierarchical_inference"] * 3,
            "evidence_scope": ["genus_inference"] * 3,
            "review_status": ["pending"] * 3,
        }
    )
    rescue = priority.build_genus_rescue_audit(
        master,
        snapshot,
        qualified,
        genus_is_explicit=True,
        input_basis="qualified.csv",
    ).set_index("accepted_species")
    assert rescue.loc["Alpha one", "potential_genus_rescue_cells"] == 1
    assert rescue.loc["Alpha two", "potential_genus_rescue_cells"] == 1
    # Beta already has a direct SI cell, so another direct SI value is not an
    # unresolved focal acquisition and earns no rescue priority.
    assert rescue.loc["Beta one", "potential_genus_rescue_cells"] == 0
    assert rescue.loc["Unknown one", "potential_genus_rescue_state"] == "unknown"
    assert pd.isna(rescue.loc["Unknown one", "potential_genus_rescue_cells"])

    no_name_parse = priority.build_genus_rescue_audit(
        master,
        snapshot,
        qualified,
        genus_is_explicit=False,
        input_basis="qualified.csv",
    )
    assert set(no_name_parse["potential_genus_rescue_state"]) == {"unknown"}


def test_frozen_source_feature_counts_species_direct_only() -> None:
    sources = pd.DataFrame(
        {
            "accepted_species": ["Alpha two", "Alpha two", "Beta one"],
            "evidence_scope": ["species_direct", "species_direct", "genus_inference"],
            "source_chunk_id": ["chunk-1", "chunk-2", "chunk-3"],
            "source_url": ["https://example.test/1", "https://example.test/2", ""],
            "packet_id": ["packet-a", "packet-a", "packet-b"],
        }
    )
    audit = priority.build_frozen_source_audit(
        _master(), sources, input_basis="frozen packets"
    ).set_index("accepted_species")

    assert audit.loc["Alpha two", "frozen_species_direct_source"] == "true"
    assert audit.loc["Alpha two", "frozen_species_direct_source_count"] == 2
    assert audit.loc["Alpha two", "frozen_species_direct_evidence_scope"] == "species_direct"
    assert audit.loc["Beta one", "frozen_species_direct_source"] == "false"
    assert audit.loc["Beta one", "frozen_species_direct_source_count"] == 0

    unknown = priority.build_frozen_source_audit(
        _master(), None, input_basis="unknown:no frozen packet input"
    )
    assert set(unknown["frozen_species_direct_source"]) == {"unknown"}
    assert unknown["frozen_species_direct_source_count"].isna().all()


def test_global_queue_sorts_audit_components_in_declared_order() -> None:
    ledger = pd.DataFrame(
        {
            "accepted_species": ["Geo", "Gap", "Rescue", "Frozen", "Fallback"],
            "genus": ["G"] * 5,
            "family": ["Fam"] * 5,
            "n_islands": [1, 1, 1, 1, 100],
            "n_records": [1, 1, 1, 1, 100],
            "task_status": ["pending"] * 5,
            "scheduling_region_category": [
                "reviewed_nh_temperate",
                "nh_temperate_unreviewed",
                "nh_temperate_unreviewed",
                "nh_temperate_unreviewed",
                "nh_temperate_unreviewed",
            ],
            "unresolved_trait_cells": [1, 14, 13, 13, 13],
            "potential_genus_rescue_cells": [0, 0, 20, 10, 10],
            "frozen_species_direct_source": ["false", "false", "false", "true", "false"],
        }
    )

    batch = campaign.family_balanced_batch(ledger, "task", 5)

    assert batch["accepted_species"].tolist() == [
        "Geo",
        "Gap",
        "Rescue",
        "Frozen",
        "Fallback",
    ]
    assert not any("bombus" in column.casefold() for column in batch.columns)


def test_geographic_sort_uses_independent_nh_temperate_flags() -> None:
    frame = pd.DataFrame(
        {
            "accepted_species": [
                "Reviewed cross",
                "Unreviewed cross",
                "Arctic reviewed mix",
                "Arctic unreviewed mix",
                "Tropical southern cross",
            ],
            "scheduling_region_category": [
                "cross_region",
                "cross_region",
                "unclassified",
                "unclassified",
                "cross_region",
            ],
            "any_reviewed_NH_temperate": ["true", "false", "true", "false", "false"],
            "any_NH_temperate_unreviewed": ["true", "true", "false", "true", "false"],
        }
    )

    ranked, sort_columns = priority.add_priority_sort_keys(frame)
    geo_rank = ranked.set_index("accepted_species")["_priority_geo"].to_dict()

    assert sort_columns == ["_priority_geo"]
    assert geo_rank == {
        "Reviewed cross": 0,
        "Unreviewed cross": 1,
        "Arctic reviewed mix": 0,
        "Arctic unreviewed mix": 1,
        "Tropical southern cross": 4,
    }


def test_lane_specific_sort_ignores_higher_ranked_gaps_that_lane_cannot_fill() -> None:
    frame = pd.DataFrame(
        {
            "accepted_species": ["Color gap", "SI gap"],
            "unresolved_priority_traits": [
                "flower_primary_color",
                "self_incompatibility",
            ],
            "highest_priority_unresolved_rank": [1, 4],
            "weighted_unresolved_trait_score": [15, 12],
            "unresolved_trait_cells": [1, 1],
        }
    )

    ranked, sort_columns = priority.add_priority_sort_keys(
        frame,
        target_traits=["self_incompatibility", "mating_system"],
    )

    assert ranked.sort_values(sort_columns)["accepted_species"].tolist() == [
        "SI gap",
        "Color gap",
    ]


def test_refresh_removes_legacy_bombus_priority_columns() -> None:
    ledger = pd.DataFrame(
        {
            "accepted_species": ["Alpha one"],
            "task_status": ["pending"],
            "any_reviewed_NH_temperate_bombus_applicable": ["true"],
            "priority_nh_temperate_native_bombus": ["true"],
            "priority_nh_temperate_native_bombus_island_count": [1],
            "priority_island_ids": ["north"],
            "priority_source_region_evidence_ids": ["legacy-bombus-evidence"],
        }
    )
    audit = pd.DataFrame(
        {
            "accepted_species": ["Alpha one"],
            "any_reviewed_NH_temperate": ["true"],
            "scheduling_region_category": ["reviewed_nh_temperate"],
        }
    )

    refreshed = priority.merge_priority_audit(ledger, audit)

    assert refreshed.loc[0, "task_status"] == "pending"
    assert refreshed.loc[0, "scheduling_region_category"] == "reviewed_nh_temperate"
    assert not any("bombus" in column.casefold() for column in refreshed.columns)
    assert "priority_island_ids" not in refreshed.columns
    assert "priority_source_region_evidence_ids" not in refreshed.columns
