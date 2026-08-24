import importlib

import pandas as pd
import pytest
import typer

mod = importlib.import_module("island_v2.flora_status_support")

FLORA = pd.DataFrame(
    {
        "island_id": ["i1", "i1", "i1", "i1", "i2", "i2", "i2", "i3", "i3"],
        "species": ["A a", "B b", "C c", "D d", "A a", "E e", "F f", "G g", "H h"],
    }
)

STATUS = pd.DataFrame(
    {
        "island_id": ["i1", "i1", "i1", "i2", "i2", "i2", "i3", "i3"],
        "accepted_species": ["A a", "B b", "C c", "A a", "E e", "F f", "G g", "H h"],
        "origin_status": [
            "native", "native", "introduced", "native", "native", "native", "native", "native"
        ],
        "endemic_status": [
            "nonendemic", "endemic", "unresolved", "nonendemic",
            "endemic", "nonendemic", "endemic", "endemic",
        ],
        "status_source": ["GIFT"] * 8,
    }
)

EVIDENCE = pd.DataFrame(
    {
        "accepted_species": ["A a", "B b", "E e", "G g", "A a", "G g"],
        "trait_name": [
            "flower_primary_color", "flower_primary_color", "flower_primary_color",
            "flower_primary_color", "self_incompatibility", "self_incompatibility",
        ],
        "resolution_status": ["resolved"] * 6,
    }
)

COV = pd.DataFrame(
    {
        "island_id": ["i1", "i2", "i3"],
        "analysis_regime": ["northern", "northern", "tropical"],
        "distance_to_continent_km": [10.0, 100.0, 1000.0],
        "spatial_block": ["b1", "b2", "b3"],
    }
)

OUTCOMES = {
    "flower_colour": {"trait_names": ["flower_primary_color"]},
    "self_compatibility": {"trait_names": ["self_incompatibility"]},
}


def test_status_attach_retains_missing_as_unresolved():
    out = mod.attach_floristic_status(FLORA, STATUS).set_index(["island_id", "accepted_species"])
    assert out.loc[("i1", "A a"), "floristic_status"] == "native_nonendemic"
    assert out.loc[("i1", "B b"), "floristic_status"] == "endemic"
    assert out.loc[("i1", "C c"), "floristic_status"] == "introduced"
    assert out.loc[("i1", "D d"), "floristic_status"] == "unresolved"


def test_status_conflict_fails_closed():
    conflict = pd.concat(
        [
            STATUS,
            pd.DataFrame(
                {
                    "island_id": ["i1"],
                    "accepted_species": ["A a"],
                    "origin_status": ["introduced"],
                    "endemic_status": ["unresolved"],
                    "status_source": ["other"],
                }
            ),
        ],
        ignore_index=True,
    )
    out = mod.attach_floristic_status(FLORA, conflict).set_index(["island_id", "accepted_species"])
    assert bool(out.loc[("i1", "A a"), "status_conflict"])
    assert out.loc[("i1", "A a"), "floristic_status"] == "unresolved"


def test_endemism_response_uses_native_denominator_only():
    status = mod.attach_floristic_status(FLORA, STATUS)
    out = mod.island_status_counts(status).set_index("island_id")
    assert out.loc["i1", "n_native_species"] == 2
    assert out.loc["i1", "n_endemic_species"] == 1
    assert out.loc["i1", "endemic_share_of_native"] == pytest.approx(0.5)
    assert out.loc["i1", "n_introduced_species"] == 1


def test_direct_species_by_outcome_uses_final_traits_only():
    result = mod.direct_species_by_outcome(EVIDENCE, OUTCOMES)
    assert result["flower_colour"] == {"A a", "B b", "E e", "G g"}
    assert result["self_compatibility"] == {"A a", "G g"}


def test_direct_support_is_status_stratified_and_not_all_or_nothing():
    status = mod.attach_floristic_status(FLORA, STATUS)
    direct = mod.direct_species_by_outcome(EVIDENCE, OUTCOMES)
    out = mod.build_direct_support(
        status, direct, COV, min_direct_species=1, min_direct_fraction=0.0
    )
    i1_native_colour = out.loc[
        out["island_id"].eq("i1")
        & out["stratum"].eq("all_native")
        & out["outcome"].eq("flower_colour")
    ].iloc[0]
    assert i1_native_colour["n_stratum_species"] == 2
    assert i1_native_colour["n_direct_species"] == 2
    assert bool(i1_native_colour["direct_eligible"])

    i1_native_sc = out.loc[
        out["island_id"].eq("i1")
        & out["stratum"].eq("all_native")
        & out["outcome"].eq("self_compatibility")
    ].iloc[0]
    assert i1_native_sc["n_direct_species"] == 1


def test_support_summary_marks_30_49_as_targeted_zone():
    rows = []
    for index in range(40):
        rows.append(
            {
                "island_id": f"i{index}",
                "stratum": "all_native",
                "outcome": "flower_colour",
                "direct_eligible": True,
                "direct_fraction": 0.6,
                "n_direct_species": 35,
                "analysis_regime": "northern",
                "distance_to_continent_km": float(index + 1),
            }
        )
    summary = mod.summarize_support(
        pd.DataFrame(rows),
        regime_column="analysis_regime",
        isolation_column="distance_to_continent_km",
        pilot_min_islands=30,
        confirmatory_min_islands=50,
    )
    assert summary.iloc[0]["support_class"] == "targeted_acquisition_zone"
    assert summary.iloc[0]["n_direct_eligible_islands"] == 40
    assert summary.iloc[0]["isolation_q95"] > summary.iloc[0]["isolation_q05"]


def test_acquisition_candidates_prioritize_one_record_unlocks():
    status = mod.attach_floristic_status(FLORA, STATUS)
    direct = {"flower_colour": {"A a"}}
    support = mod.build_direct_support(
        status,
        direct,
        COV,
        min_direct_species=2,
        min_direct_fraction=0.0,
    )
    candidates = mod.acquisition_candidates(
        status,
        support,
        direct,
        min_direct_species=2,
        confirmatory_min_islands=50,
        regime_column="analysis_regime",
        isolation_column="distance_to_continent_km",
    )
    assert not candidates.empty
    row = candidates.loc[candidates["accepted_species"].eq("B b")].iloc[0]
    assert row["one_record_unlocks_islands"] == 1


def test_bad_status_values_fail_closed():
    bad = STATUS.copy()
    bad.loc[0, "origin_status"] = "probably_native"
    with pytest.raises(typer.BadParameter):
        mod.collapse_status_ledger(bad)
