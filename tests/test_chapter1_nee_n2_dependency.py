from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_nee_n2_dependency import (
    load_config,
    qualify_n2_support,
    validate_dependency_ledger,
)

CONFIG = Path("config/chapter1_nee_n2_dependency.yml")


def _config() -> dict:
    return load_config(CONFIG)


def _ledger_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "plant_species": "Plantus alpha",
        "accepted_genus": "Plantus",
        "channel_id": "bombus",
        "study_id": "study-1",
        "evidence_design": "D1_channel_exclusion_reproductive_output",
        "evidence_tier": "D1",
        "dependency_estimate": 0.6,
        "dependency_lower": 0.4,
        "dependency_upper": 0.8,
        "source_citation": "Direct exclusion study",
        "source_url": "https://example.org/study-1",
        "review_status": "accepted",
        "evidence_origin": "direct_experimental_reproductive_output",
    }
    row.update(overrides)
    return row


def test_D1_is_primary_and_D3_is_not() -> None:
    ledger = pd.DataFrame(
        [
            _ledger_row(),
            _ledger_row(
                plant_species="Plantus beta",
                study_id="study-2",
                evidence_design="interaction_network_edge",
                evidence_tier="D3",
                dependency_estimate="",
                dependency_lower="",
                dependency_upper="",
                evidence_origin="GloBI_interaction_presence",
            ),
        ]
    )
    primary, receipt = validate_dependency_ledger(ledger, _config())
    assert len(primary) == 1
    assert primary.iloc[0]["evidence_tier"] == "D1"
    assert receipt["n_D3_rows"] == 1
    assert receipt["N2_fitted"] is False


def test_pollination_guild_cannot_be_primary_dependency() -> None:
    ledger = pd.DataFrame([_ledger_row(evidence_origin="LLM_pollination_guild_label")])
    with pytest.raises(ValueError, match="prohibited primary dependency evidence origin"):
        validate_dependency_ledger(ledger, _config())


def test_baker_traits_cannot_be_relabelled_as_dependency() -> None:
    ledger = pd.DataFrame([_ledger_row(evidence_origin="self_compatibility_mating_system")])
    with pytest.raises(ValueError, match="prohibited primary dependency evidence origin"):
        validate_dependency_ledger(ledger, _config())


def test_duplicate_study_species_channel_fails_closed() -> None:
    row = _ledger_row()
    with pytest.raises(ValueError, match="duplicate study"):
        validate_dependency_ledger(pd.DataFrame([row, row]), _config())


def test_dependency_interval_must_contain_estimate() -> None:
    ledger = pd.DataFrame([_ledger_row(dependency_lower=0.7, dependency_estimate=0.6)])
    with pytest.raises(ValueError, match="interval must contain"):
        validate_dependency_ledger(ledger, _config())


def _support_dependency() -> pd.DataFrame:
    rows = []
    for channel_index, channel in enumerate(["bombus", "lepidoptera", "diptera"]):
        for genus_index in range(50):
            rows.append(
                {
                    "source_region_id": f"source-{genus_index % 5}",
                    "accepted_genus": f"Genus_{channel_index}_{genus_index}",
                    "channel_id": channel,
                    "dependency_mean": 0.2 + (genus_index % 6) * 0.1,
                    "dependency_sd": 0.1,
                    "n_species_evidenced": 2,
                    "evidence_tiers": "D1",
                    "primary_dependency_evaluable": True,
                }
            )
    return pd.DataFrame(rows)


def _channel_states() -> pd.DataFrame:
    rows = []
    for channel in ["bombus", "lepidoptera", "diptera"]:
        for index in range(20):
            rows.append(
                {"island_id": f"{channel}-ret-{index}", "channel_id": channel, "channel_state": "retained"}
            )
            rows.append(
                {"island_id": f"{channel}-dis-{index}", "channel_id": channel, "channel_state": "disrupted"}
            )
    return pd.DataFrame(rows)


def test_support_gate_passes_only_from_design_support_not_N2_outcome() -> None:
    report = qualify_n2_support(_support_dependency(), _channel_states(), _config())
    assert report["passed"] is True
    assert report["n_eligible_channels"] == 3
    assert report["n_unique_dependency_resolved_genera_total"] == 150
    assert report["genus_entry_outcome_read"] is False
    assert report["N2_fitted"] is False


def test_D3_cannot_be_marked_primary_evaluable_at_genus_level() -> None:
    dependency = _support_dependency()
    dependency.loc[0, "evidence_tiers"] = "D1|D3"
    with pytest.raises(ValueError, match="D3/non-primary provenance"):
        qualify_n2_support(dependency, _channel_states(), _config())


def test_single_species_genus_cannot_be_exact_point_imputation() -> None:
    dependency = _support_dependency()
    dependency.loc[0, "n_species_evidenced"] = 1
    dependency.loc[0, "dependency_sd"] = 0.0
    with pytest.raises(ValueError, match="single-species genus dependency"):
        qualify_n2_support(dependency, _channel_states(), _config())


def test_support_gate_fails_if_disrupted_support_is_too_small() -> None:
    states = _channel_states()
    states = states.loc[
        ~(
            states["channel_id"].eq("diptera")
            & states["channel_state"].eq("disrupted")
            & ~states["island_id"].isin([f"diptera-dis-{i}" for i in range(10)])
        )
    ].copy()
    report = qualify_n2_support(_support_dependency(), states, _config())
    assert report["passed"] is False
    assert report["n_eligible_channels"] == 2
