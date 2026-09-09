from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from island_v2.chapter1_h5_channel_chain import build_outputs, evaluate_channel_chain


def _config() -> dict:
    path = Path("config/chapter1_h5_channel_chain.yml")
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _row(**overrides) -> dict[str, str]:
    row = {
        "system_id": "sys",
        "island_id": "island-1",
        "channel_id": "large_bee_channel",
        "source_channel_state": "available",
        "source_evidence_id": "source-ref",
        "island_retention_state": "retained",
        "retention_evidence_type": "direct_occurrence",
        "retention_evidence_id": "retention-ref",
        "observation_effort_flower_hours": "10",
        "visit_bouts": "4",
        "no_visit_control_present": "yes",
        "mean_single_visit_conspecific_pollen": "12",
        "mean_no_visit_conspecific_pollen": "2",
        "evidence_scope": "independent_direct",
    }
    row.update({key: str(value) for key, value in overrides.items()})
    return row


def test_retained_positive_visit_chain_computes_effective_service():
    audit = evaluate_channel_chain(pd.DataFrame([_row()]), _config())
    row = audit.iloc[0]
    assert row["source_gate_pass"]
    assert row["retention_gate_pass"]
    assert row["visitation_gate_pass"]
    assert row["single_visit_effectiveness_gate_pass"]
    assert row["visit_bouts_per_flower_hour"] == pytest.approx(0.4)
    assert row["mean_background_adjusted_svd"] == pytest.approx(10.0)
    assert row["effective_pollen_delivery_per_flower_hour"] == pytest.approx(4.0)
    assert row["full_positive_visit_chain_observed"]


def test_adequate_zero_visitation_means_zero_service_not_zero_effectiveness():
    audit = evaluate_channel_chain(
        pd.DataFrame(
            [
                _row(
                    island_retention_state="disrupted",
                    retention_evidence_type="adequate_non_detection",
                    visit_bouts="0",
                    no_visit_control_present="no",
                    mean_single_visit_conspecific_pollen="",
                    mean_no_visit_conspecific_pollen="",
                )
            ]
        ),
        _config(),
    )
    row = audit.iloc[0]
    assert row["visitation_gate_pass"]
    assert not row["single_visit_effectiveness_gate_pass"]
    assert row["single_visit_effectiveness_status"] == "not_applicable_zero_visitation"
    assert row["effective_service_gate_pass"]
    assert row["effective_pollen_delivery_per_flower_hour"] == pytest.approx(0.0)
    assert not row["full_positive_visit_chain_observed"]


def test_positive_visits_without_no_visit_control_withhold_service():
    audit = evaluate_channel_chain(
        pd.DataFrame([_row(no_visit_control_present="no")]),
        _config(),
    )
    row = audit.iloc[0]
    assert row["visitation_gate_pass"]
    assert not row["single_visit_effectiveness_gate_pass"]
    assert not row["effective_service_gate_pass"]
    assert pd.isna(row["effective_pollen_delivery_per_flower_hour"])


def test_structural_absence_cannot_be_relabelled_as_disruption():
    with pytest.raises(ValueError, match="cannot be labelled as island loss"):
        evaluate_channel_chain(
            pd.DataFrame(
                [
                    _row(
                        source_channel_state="structurally_absent",
                        island_retention_state="disrupted",
                    )
                ]
            ),
            _config(),
        )


def test_retained_disrupted_service_contrast_becomes_ready():
    retained = _row(island_id="retained")
    disrupted = _row(
        island_id="disrupted",
        island_retention_state="disrupted",
        retention_evidence_type="adequate_non_detection",
        visit_bouts="0",
        no_visit_control_present="no",
        mean_single_visit_conspecific_pollen="",
        mean_no_visit_conspecific_pollen="",
    )
    audit, readiness, summary = build_outputs(
        pd.DataFrame([retained, disrupted]),
        _config(),
    )
    assert len(audit) == 2
    assert readiness.iloc[0]["retained_disrupted_service_contrast_ready"]
    assert summary["n_channel_contrasts_ready"] == 1
    assert summary["H5_channel_side_status"] == "contrast_ready"
    assert summary["H5_full_status"] == "requires_post_H2_H4_plant_residual_join"


def test_climate_or_syndrome_proxy_cannot_substitute_for_retention_evidence():
    with pytest.raises(ValueError, match="retention_evidence_id"):
        evaluate_channel_chain(
            pd.DataFrame(
                [
                    _row(
                        retention_evidence_type="climate_compatibility_only",
                        retention_evidence_id="",
                    )
                ]
            ),
            _config(),
        )
