import pandas as pd

from island_v2.chapter1_trait_first_family_audit import build_trait_first_scores


def test_trait_first_scores_use_historical_family_weights():
    counts = pd.DataFrame([
        {"island_id":"i1","stratum":"all_observed","outcome":"generalized_form","successes":8,"trials":10},
        {"island_id":"i1","stratum":"all_observed","outcome":"actinomorphic_symmetry","successes":4,"trials":10},
        {"island_id":"i1","stratum":"all_observed","outcome":"shallow_open_tube","successes":2,"trials":10},
        {"island_id":"i1","stratum":"all_observed","outcome":"self_compatibility","successes":6,"trials":10},
        {"island_id":"i1","stratum":"all_observed","outcome":"selfing_mating_system","successes":3,"trials":10},
        {"island_id":"i1","stratum":"all_observed","outcome":"autonomous_selfing","successes":9,"trials":10},
    ])
    out=build_trait_first_scores(counts, require_all_components=True).set_index("family")
    expected_access=(1.0*0.8 + 0.75*0.4 + 1.0*0.2)/2.75
    expected_repro=(0.6+0.3+0.9)/3.0
    assert abs(float(out.loc["accessibility_generalization","family_score"])-expected_access) < 1e-12
    assert abs(float(out.loc["reproductive_assurance","family_score"])-expected_repro) < 1e-12
    assert int(out.loc["accessibility_generalization","n_component_traits"]) == 3


def test_available_components_keeps_partial_family_but_complete_rule_drops_it():
    counts = pd.DataFrame([
        {"island_id":"i1","stratum":"all_observed","outcome":"generalized_form","successes":8,"trials":10},
        {"island_id":"i1","stratum":"all_observed","outcome":"actinomorphic_symmetry","successes":4,"trials":10},
    ])
    available=build_trait_first_scores(counts, require_all_components=False)
    complete=build_trait_first_scores(counts, require_all_components=True)
    assert len(available) == 1
    assert available.iloc[0]["family"] == "accessibility_generalization"
    assert complete.empty
