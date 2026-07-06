import pandas as pd

from island_v2.applicability_auto_accept import auto_accept

RULES = {
    "placeholder_locator_markers": ["to be supplied", "to be confirmed"],
    "registered_by": "tester",
    "source_region_rules": {
        "min_independent_evidence_records": 1,
        "require_source_citation": True,
        "require_real_source_locator": False,
        "rules": [
            {
                "rule_id": "afrotropical_no_native_bombus",
                "match_biogeographic_region_regex": "^Afrotropical",
                "require_native_bombus_status": "native_absent",
                "accept_applicability": "structurally_not_applicable",
            },
            {
                "rule_id": "palearctic_temperate_east_asia_native_bombus",
                "match_biogeographic_region_regex": "^Palearctic_temperate",
                "require_native_bombus_status": "native_present",
                "accept_applicability": "applicable",
            },
        ],
    },
    "assignment_rules": {
        "auto_accept": True,
        "require_source_region_accepted": True,
        "require_assignment_citation": True,
        "require_real_assignment_locator": True,
    },
}


def _registry():
    return pd.DataFrame(
        [
            {
                "source_region_id": "cape",
                "biogeographic_region": "Afrotropical_Cape",
                "native_bombus_status": "native_absent",
                "applicability": "structurally_not_applicable",
                "source_region_evidence_id": "EV_cape",
                "review_status": "agent_drafted_pending_review",
                "reviewer": "",
                "decision_date": "2026-07-04",
            },
            {
                "source_region_id": "japan",
                "biogeographic_region": "Palearctic_temperate_East_Asia",
                "native_bombus_status": "native_present",
                "applicability": "applicable",
                "source_region_evidence_id": "EV_japan",
                "review_status": "agent_drafted_pending_review",
                "reviewer": "",
                "decision_date": "2026-07-04",
            },
            {
                "source_region_id": "mystery",
                "biogeographic_region": "Unknown_realm",
                "native_bombus_status": "uncertain",
                "applicability": "unresolved",
                "source_region_evidence_id": "EV_mystery",
                "review_status": "agent_drafted_pending_review",
                "reviewer": "",
                "decision_date": "2026-07-04",
            },
        ]
    )


def _region_evidence():
    return pd.DataFrame(
        [
            {"source_region_evidence_id": "EV_cape", "source_region_id": "cape", "source_citation": "Williams 1998", "source_locator": "pp. 79-152", "review_status": "agent_drafted_pending_review", "reviewer": "", "decision_date": "2026-07-04"},
            {"source_region_evidence_id": "EV_japan", "source_region_id": "japan", "source_citation": "Williams 1998", "source_locator": "pp. 79-152", "review_status": "agent_drafted_pending_review", "reviewer": "", "decision_date": "2026-07-04"},
            {"source_region_evidence_id": "EV_mystery", "source_region_id": "mystery", "source_citation": "", "source_locator": "", "review_status": "agent_drafted_pending_review", "reviewer": "", "decision_date": "2026-07-04"},
        ]
    )


def _assignments():
    return pd.DataFrame(
        [
            {"island_id": "isl_robben", "source_region_id": "cape", "assignment_evidence_id": "AEV_robben", "review_status": "agent_drafted_pending_review", "reviewer": "", "decision_date": "2026-07-04"},
        ]
    )


def _assignment_evidence(locator):
    return pd.DataFrame(
        [
            {"assignment_evidence_id": "AEV_robben", "island_id": "isl_robben", "source_citation": "Boucher & Jarman 1977", "source_locator": locator, "review_status": "agent_drafted_pending_review", "reviewer": "", "decision_date": "2026-07-04"},
        ]
    )


def test_realm_rules_auto_accept_regions_and_their_evidence():
    reg, rev, _asg, _aev, report = auto_accept(
        _registry(), _region_evidence(), _assignments(), _assignment_evidence("to be supplied by human"), RULES, "2026-07-06"
    )
    by = reg.set_index("source_region_id")
    assert by.loc["cape", "review_status"] == "accepted"
    assert by.loc["cape", "reviewer"] == "auto_rule:afrotropical_no_native_bombus"
    assert by.loc["japan", "review_status"] == "accepted"
    # Non-matching realm stays pending (human exception).
    assert by.loc["mystery", "review_status"] == "agent_drafted_pending_review"
    # Linked evidence is accepted so the applicability validator stays satisfied.
    ev = rev.set_index("source_region_evidence_id")
    assert ev.loc["EV_cape", "review_status"] == "accepted"
    assert ev.loc["EV_mystery", "review_status"] == "agent_drafted_pending_review"
    assert report["n_regions_auto_accepted"] == 2


def test_placeholder_assignment_locator_stays_pending():
    _reg, _rev, asg, _aev, report = auto_accept(
        _registry(), _region_evidence(), _assignments(), _assignment_evidence("Locator to be supplied by human reviewer"), RULES, "2026-07-06"
    )
    assert asg.set_index("island_id").loc["isl_robben", "review_status"] == "agent_drafted_pending_review"
    assert report["n_assignments_auto_accepted"] == 0
    assert report["assignments_left_pending"][0]["reason"] == "placeholder assignment locator"


def test_real_assignment_locator_auto_accepts_when_region_accepted():
    _reg, _rev, asg, aev, report = auto_accept(
        _registry(), _region_evidence(), _assignments(), _assignment_evidence("Boucher & Jarman 1977, p. 3, fig. 1"), RULES, "2026-07-06"
    )
    assert asg.set_index("island_id").loc["isl_robben", "review_status"] == "accepted"
    assert aev.set_index("assignment_evidence_id").loc["AEV_robben", "review_status"] == "accepted"
    assert report["n_assignments_auto_accepted"] == 1


def test_already_accepted_rows_are_left_untouched():
    reg = _registry()
    reg.loc[reg["source_region_id"] == "cape", "review_status"] = "accepted"
    reg.loc[reg["source_region_id"] == "cape", "reviewer"] = "human_curator"
    out, _rev, _asg, _aev, _report = auto_accept(
        reg, _region_evidence(), _assignments(), _assignment_evidence("x"), RULES, "2026-07-06"
    )
    # A human-accepted row keeps its reviewer, not overwritten by auto_rule.
    assert out.set_index("source_region_id").loc["cape", "reviewer"] == "human_curator"
