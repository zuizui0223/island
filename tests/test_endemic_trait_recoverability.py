import pandas as pd

from island_v2.endemic_trait_recoverability import (
    build_minimum_pilot_queue,
    recoverability_headroom,
    strictest_pilot_threshold,
)


def _fixture():
    decisions = pd.DataFrame(
        [
            {
                "regime": "northern_midlatitude",
                "outcome": "generalized_form",
                "n_endemic_status_islands": 35,
                "n_direct_supported_islands": 28,
                "decision": "recoverability_test_before_acquisition",
            },
            {
                "regime": "tropical",
                "outcome": "self_compatibility",
                "n_endemic_status_islands": 40,
                "n_direct_supported_islands": 30,
                "decision": "targeted_trait_acquisition",
            },
        ]
    )
    unsupported = pd.DataFrame(
        [
            {"regime": "northern_midlatitude", "outcome": "generalized_form", "island_id": "i1"},
            {"regime": "northern_midlatitude", "outcome": "generalized_form", "island_id": "i2"},
            {"regime": "northern_midlatitude", "outcome": "generalized_form", "island_id": "i3"},
        ]
    )
    targets = pd.DataFrame(
        [
            {
                "regime": "northern_midlatitude",
                "outcome": "generalized_form",
                "accepted_species": "Alpha one",
                "n_records": 80,
                "priority_rank": 1,
            },
            {
                "regime": "northern_midlatitude",
                "outcome": "generalized_form",
                "accepted_species": "Beta two",
                "n_records": 15,
                "priority_rank": 2,
            },
        ]
    )
    status = pd.DataFrame(
        [
            {
                "island_id": "i1",
                "accepted_species": "Alpha one",
                "origin_status": "native",
                "endemic_status": "endemic",
            },
            {
                "island_id": "i2",
                "accepted_species": "Alpha one",
                "origin_status": "native",
                "endemic_status": "endemic",
            },
            {
                "island_id": "i3",
                "accepted_species": "Beta two",
                "origin_status": "native",
                "endemic_status": "endemic",
            },
            {
                "island_id": "i3",
                "accepted_species": "Introduced three",
                "origin_status": "introduced",
                "endemic_status": "endemic",
            },
        ]
    )
    return decisions, unsupported, targets, status


def test_recoverability_reports_pilot_headroom_without_trait_inference():
    decisions, unsupported, targets, status = _fixture()
    result = recoverability_headroom(
        decisions,
        unsupported,
        targets,
        status,
        thresholds=(100, 50, 10),
        pilot_target=30,
        confirmatory_target=50,
    )

    assert set(result["outcome"]) == {"generalized_form"}
    threshold_100 = result.loc[result["gbif_record_threshold"].eq(100)].iloc[0]
    assert threshold_100["eligible_candidate_species"] == 0
    assert threshold_100["max_support_if_all_candidates_recovered"] == 28
    assert not bool(threshold_100["pilot_reachable"])

    threshold_50 = result.loc[result["gbif_record_threshold"].eq(50)].iloc[0]
    assert threshold_50["eligible_candidate_species"] == 1
    assert threshold_50["unsupported_islands_reached"] == 2
    assert threshold_50["max_support_if_all_candidates_recovered"] == 30
    assert bool(threshold_50["pilot_reachable"])
    assert not bool(threshold_50["confirmatory_reachable"])

    threshold_10 = result.loc[result["gbif_record_threshold"].eq(10)].iloc[0]
    assert threshold_10["unsupported_islands_reached"] == 3
    assert threshold_10["max_support_if_all_candidates_recovered"] == 31


def test_strictest_pilot_threshold_uses_highest_passing_threshold():
    decisions, unsupported, targets, status = _fixture()
    result = recoverability_headroom(
        decisions,
        unsupported,
        targets,
        status,
        thresholds=(100, 50, 10),
    )
    summary = strictest_pilot_threshold(result)

    assert len(summary) == 1
    row = summary.iloc[0]
    assert row["strictest_pilot_reachable_gbif_threshold"] == 50
    assert row["eligible_candidate_species_at_threshold"] == 1
    assert row["unsupported_islands_reached_at_threshold"] == 2
    assert row["max_support_at_threshold"] == 30
    assert bool(row["pilot_reachable_under_tested_thresholds"])


def test_minimum_pilot_queue_stops_at_30_and_keeps_threshold_gate():
    decisions, unsupported, targets, status = _fixture()
    headroom = recoverability_headroom(
        decisions,
        unsupported,
        targets,
        status,
        thresholds=(100, 50, 10),
        pilot_target=30,
    )
    summary = strictest_pilot_threshold(headroom)
    queue, plan = build_minimum_pilot_queue(
        summary,
        unsupported,
        targets,
        status,
        pilot_target=30,
    )

    assert queue["accepted_species"].tolist() == ["Alpha one"]
    assert queue.iloc[0]["selected_gbif_record_threshold"] == 50
    assert queue.iloc[0]["new_islands_added"] == 2
    assert queue.iloc[0]["cumulative_new_islands"] == 2
    assert queue.iloc[0]["required_pilot_gap"] == 2

    row = plan.iloc[0]
    assert row["current_direct_supported_islands"] == 28
    assert row["pilot_gap"] == 2
    assert row["n_selected_species"] == 1
    assert row["n_new_islands_covered"] == 2
    assert bool(row["queue_reaches_pilot"])
