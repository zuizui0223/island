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
        "H4": "tiered_functional_bridge_with_posthoc_discovery_and_support_gated_validation",
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
    raw = lock["H2"]["raw_colour_and_architecture"]
    assert raw["role"] == "raw_state_pollination_syndrome_concordance_not_weighted_score"
    assert (
        raw["colour_conditioned_architecture"]["northern_high_latitude"]
        ["blue_purple_butterfly_form_given_colour"]["direct_q"]
        < 0.05
    )
    assert lock["H4"]["prospective_wild_temporal_replication"]["outcomes_opened"] is False
    assert lock["H4"]["colour_bridge"]["promoted"] is False
    exact = lock["H4"]["exact_H2_score_discovery"]
    assert exact["reproductive_assurance"]["estimate"] < 0
    assert exact["accessibility_generalization"]["estimate"] < 0
    assert exact["accessibility_generalization"]["supplemental_only_two_sided_p"] < 0.05
    reconstruction = lock["H4"]["atomic_reconstruction_sensitivity"]
    assert reconstruction["accessibility_generalization"]["estimate"] < 0
    assert lock["claim_ceiling"]["pollen_limitation_mediates_global_syndrome"] is False
