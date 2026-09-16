import pandas as pd

from island_v2.chapter1_status_resolution_process import build_resolution_counts, _audit_config


def test_build_resolution_counts_counts_unique_island_species():
    frame = pd.DataFrame(
        [
            {"island_id": "a", "accepted_species": "sp1", "origin_status": "native"},
            {"island_id": "a", "accepted_species": "sp1", "origin_status": "native"},
            {"island_id": "a", "accepted_species": "sp2", "origin_status": "unresolved"},
            {"island_id": "b", "accepted_species": "sp3", "origin_status": "introduced"},
        ]
    )
    out = build_resolution_counts(frame).set_index("island_id")
    assert int(out.loc["a", "successes"]) == 1
    assert int(out.loc["a", "trials"]) == 2
    assert float(out.loc["a", "resolution_fraction"]) == 0.5
    assert int(out.loc["b", "successes"]) == 1
    assert int(out.loc["b", "trials"]) == 1


def test_audit_config_uses_single_observation_process_outcome():
    cfg = {"model_outcomes": ["x", "y"], "minimum_outcomes_per_vector": 2, "strata": ["all_native"]}
    out = _audit_config(cfg)
    assert out["model_outcomes"] == ["status_resolved"]
    assert out["minimum_outcomes_per_vector"] == 1
    assert out["strata"] == ["all_observed"]
    assert cfg["model_outcomes"] == ["x", "y"]
