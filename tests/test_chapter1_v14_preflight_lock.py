from pathlib import Path
import json


def test_v14_preflight_lock_matches_reordered_architecture():
    lock = json.loads(
        Path("config/chapter1_v14_preflight_result_lock.json").read_text(encoding="utf-8")
    )
    assert lock["status"] == "local_reproduction_pending_ci_promotion"
    assert lock["canonical"] is False
    assert lock["architecture"] == {
        "H1": "seven_response_global_island_syndrome",
        "H2": "selfing_vs_pollinator_facing_floral_decomposition",
        "H3": "independent_global_pollen_limitation",
        "H4": "posthoc_functional_bridge",
    }
    for scope in ("all_analysis", "direct_only"):
        assert all(
            value["joint_q"] < 0.05
            for value in lock["H1"][scope].values()
        )
        assert all(
            value["equal_domain_orientation"] > 0
            for value in lock["H1"][scope].values()
        )
    assert lock["H2"]["named_pollination_architecture"]["robust_named_identity_specific_effect"] is False
    assert lock["claim_ceiling"]["pollen_limitation_mediates_global_syndrome"] is False
