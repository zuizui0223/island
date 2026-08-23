import importlib

import pandas as pd

mod = importlib.import_module("island_v2.endemic_minimum_acquisition_queue")


def test_queue_stops_when_island_gap_is_covered():
    decisions = pd.DataFrame(
        {
            "regime": ["tropical", "northern_midlatitude"],
            "outcome": ["plain_colour", "generalized_form"],
            "gap_to_50": [2, 30],
            "decision": ["targeted_trait_acquisition", "recoverability_test_before_acquisition"],
        }
    )
    unsupported = pd.DataFrame(
        {
            "regime": ["tropical"] * 3,
            "outcome": ["plain_colour"] * 3,
            "island_id": ["i1", "i2", "i3"],
        }
    )
    targets = pd.DataFrame(
        {
            "regime": ["tropical"] * 3,
            "outcome": ["plain_colour"] * 3,
            "accepted_species": ["A a", "B b", "C c"],
            "priority_rank": [1, 2, 3],
            "n_records": [100, 90, 80],
        }
    )
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i2", "i3"],
            "accepted_species": ["A a", "B b", "C c"],
            "origin_status": ["native"] * 3,
            "endemic_status": ["endemic"] * 3,
        }
    )
    queue, summary = mod.build_minimum_queue(decisions, unsupported, targets, status)
    assert queue["accepted_species"].tolist() == ["A a", "B b"]
    assert queue.iloc[-1]["cumulative_new_islands"] == 2
    assert bool(summary.iloc[0]["queue_reaches_target"])
    assert set(queue["outcome"]) == {"plain_colour"}


def test_species_covering_two_islands_counts_both_once():
    decisions = pd.DataFrame(
        {
            "regime": ["tropical"],
            "outcome": ["plain_colour"],
            "gap_to_50": [2],
            "decision": ["targeted_trait_acquisition"],
        }
    )
    unsupported = pd.DataFrame(
        {
            "regime": ["tropical", "tropical"],
            "outcome": ["plain_colour", "plain_colour"],
            "island_id": ["i1", "i2"],
        }
    )
    targets = pd.DataFrame(
        {
            "regime": ["tropical"],
            "outcome": ["plain_colour"],
            "accepted_species": ["A a"],
            "priority_rank": [1],
        }
    )
    status = pd.DataFrame(
        {
            "island_id": ["i1", "i2"],
            "accepted_species": ["A a", "A a"],
            "origin_status": ["native", "native"],
            "endemic_status": ["endemic", "endemic"],
        }
    )
    queue, summary = mod.build_minimum_queue(decisions, unsupported, targets, status)
    assert len(queue) == 1
    assert queue.iloc[0]["new_islands_added"] == 2
    assert bool(summary.iloc[0]["queue_reaches_target"])
