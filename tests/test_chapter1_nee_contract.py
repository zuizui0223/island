from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from island_v2.chapter1_nee_contract import (
    EXPECTED_CHANNELS,
    load_contract,
    validate_contract_payload,
)


CONFIG = Path("config/chapter1_nee_double_filter.yml")


def _payload() -> dict[str, object]:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def test_current_nee_contract_is_valid_and_bound_to_frozen_parent() -> None:
    receipt = validate_contract_payload(load_contract(CONFIG))
    assert receipt["contract"] == "chapter1_nee_double_filter_v1"
    assert receipt["parent_contract"] == "chapter1_progressive_analysis_v1"
    assert receipt["frozen_results"] == ["H1", "H2", "H3", "H4"]
    assert tuple(receipt["channels"]) == EXPECTED_CHANNELS
    assert receipt["primary_holdout"] == "leave_one_archipelago_out"
    assert receipt["primary_target"] == "source_available_genus_entry"
    assert receipt["primary_metric"] == "mean_log_loss"


def test_contract_rejects_posthoc_channel_addition() -> None:
    payload = copy.deepcopy(_payload())
    payload["N0_qualification"]["channel_ontology"]["channels"].append("hawkmoth_posthoc")
    with pytest.raises(ValueError, match="channel ontology changed"):
        validate_contract_payload(payload)


def test_contract_rejects_floral_trait_leakage_into_dependency() -> None:
    payload = copy.deepcopy(_payload())
    payload["N0_qualification"]["lineage_dependency"]["prohibited_evidence"].remove(
        "floral_form"
    )
    with pytest.raises(ValueError, match="focal floral architecture leakage"):
        validate_contract_payload(payload)


def test_contract_rejects_reopening_parent_results() -> None:
    payload = copy.deepcopy(_payload())
    payload["parent_freeze"]["reopen_or_rescue_parent_results"] = True
    with pytest.raises(ValueError, match="may not reopen or rescue"):
        validate_contract_payload(payload)


def test_contract_rejects_random_island_split_as_primary() -> None:
    payload = copy.deepcopy(_payload())
    payload["N0_qualification"]["heldout_prediction"]["split"] = "random_island_kfold"
    with pytest.raises(ValueError, match="leave-one-archipelago-out"):
        validate_contract_payload(payload)


def test_contract_rejects_treating_sc_as_channel_dependency() -> None:
    payload = copy.deepcopy(_payload())
    payload["N0_qualification"]["lineage_dependency"]["reproductive_assurance_role"] = (
        "channel_dependency"
    )
    with pytest.raises(ValueError, match="Baker rival"):
        validate_contract_payload(payload)


def test_contract_rejects_n2_rescue_after_failure() -> None:
    payload = copy.deepcopy(_payload())
    payload["N2_channel_dependent_lineage_filter"]["failure_action"] = (
        "try_alternative_channels"
    )
    with pytest.raises(ValueError, match="terminate H5a without rescue"):
        validate_contract_payload(payload)
