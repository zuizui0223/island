import importlib

import pandas as pd

mod = importlib.import_module("island_v2.endemic_trait_support_targets")


def _status(n_islands: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n_islands)],
            "accepted_species": [f"Genus species{i}" for i in range(n_islands)],
            "origin_status": ["native"] * n_islands,
            "endemic_status": ["endemic"] * n_islands,
        }
    )


def _cov(n_islands: int, regime: str = "tropical") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "island_id": [f"i{i}" for i in range(n_islands)],
            "analysis_regime": [regime] * n_islands,
            "log_distance_to_continent_km": [float(i) for i in range(n_islands)],
        }
    )


def _taxa(n_islands: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "accepted_species": [f"Genus species{i}" for i in range(n_islands)],
            "n_records": list(range(n_islands)),
            "n_islands": [1] * n_islands,
        }
    )


def test_status_ceiling_below_50_emits_no_targets():
    status = _status(13)
    evidence = pd.DataFrame(
        columns=["accepted_species", "trait_name", "normalized_value", "resolution_status"]
    )
    summary, _, targets = mod.rank_endemic_targets(
        status, evidence, _taxa(13), _cov(13, "southern_extratropical")
    )
    row = summary.loc[
        (summary.regime == "southern_extratropical")
        & (summary.outcome == "plain_colour")
    ].iloc[0]
    assert row.decision == "status_limited"
    assert targets.empty


def test_targeted_colour_zone_ranks_isolation_edge_before_records():
    status = _status(60)
    # Thirty-one islands already supported for colour; 29 remain unsupported.
    covered = pd.DataFrame(
        {
            "accepted_species": [f"Genus species{i}" for i in range(31)],
            "trait_name": ["flower_primary_color"] * 31,
            "normalized_value": ["white"] * 31,
            "resolution_status": ["resolved"] * 31,
        }
    )
    taxa = _taxa(60)
    # Give a central-isolation candidate more records than edge candidates.
    taxa.loc[taxa.accepted_species.eq("Genus species40"), "n_records"] = 9999
    summary, _, targets = mod.rank_endemic_targets(
        status, covered, taxa, _cov(60), confirmatory_min_islands=50
    )
    row = summary.loc[
        (summary.regime == "tropical") & (summary.outcome == "plain_colour")
    ].iloc[0]
    assert row.decision == "targeted_trait_acquisition"
    colour = targets.loc[
        (targets.regime == "tropical") & (targets.outcome == "plain_colour")
    ]
    assert not colour.empty
    # Isolation-edge candidates sort ahead of non-edge candidates even if a
    # non-edge species has many more occurrence records.
    assert bool(colour.iloc[0].extends_supported_isolation_5_95)


def test_direct_species_is_not_reemitted_as_target():
    status = _status(55)
    evidence = pd.DataFrame(
        {
            "accepted_species": ["Genus species54"],
            "trait_name": ["self_incompatibility"],
            "normalized_value": ["SC"],
            "resolution_status": ["resolved"],
        }
    )
    _, _, targets = mod.rank_endemic_targets(
        status, evidence, _taxa(55), _cov(55), confirmatory_min_islands=50
    )
    si = targets.loc[targets.outcome.eq("self_compatibility")]
    assert "Genus species54" not in set(si.accepted_species)
