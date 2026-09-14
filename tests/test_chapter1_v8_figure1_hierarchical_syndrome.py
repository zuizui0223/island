from __future__ import annotations

import json
from pathlib import Path

import pytest

from island_v2.chapter1_v8_figure1_hierarchical_syndrome import (
    FigureInputError,
    load_inputs,
    render_figure,
)


def _write_sources(tmp_path: Path) -> dict[str, Path]:
    fig3 = tmp_path / "figure3.json"
    fig3.write_text(
        json.dumps(
            {
                "contract": "chapter1_v8_figure3_submission_result_lock_v1",
                "workflow_run_id": 1,
                "artifact_id": 2,
                "artifact_digest": "sha256:figure3",
                "frozen_results": {
                    "support_ladder": "4/4 -> 4/4 -> 0/4",
                    "genus_attenuation_fraction_range": [0.788, 0.859],
                    "conditional_family_to_genus_attenuation_range": [0.706, 0.791],
                },
            }
        )
        + "\n"
    )
    fig4 = tmp_path / "figure4.json"
    fig4.write_text(
        json.dumps(
            {
                "contract": "chapter1_v8_figure4_result_lock_v1",
                "workflow_run_id": 3,
                "artifact_id": 4,
                "artifact_digest": "sha256:figure4",
                "frozen_results": {
                    "v6_palearctic_survival": "99/100",
                    "v6_tropical_survival": "35/75",
                    "geometry_promoted": "0/12",
                    "h5c_interaction_estimate": 0.06495254330943623,
                    "h5c_p_value": 0.4122068222871858,
                    "h5d_qualified": "0/8",
                },
            }
        )
        + "\n"
    )
    n1 = tmp_path / "n1.json"
    n1.write_text(
        json.dumps(
            {
                "contract": "chapter1_nee_n1_result_lock_v1",
                "n1": {
                    "workflow_run_id": 5,
                    "artifact_id": 6,
                    "artifact_digest": "sha256:n1",
                    "gate": {
                        "N1_pass": False,
                        "global_Wald_p_value": 0.6551581817640355,
                        "failure_action": "stop_before_N2_and_keep_frozen_Chapter1",
                    },
                },
            }
        )
        + "\n"
    )
    freeze = tmp_path / "freeze.md"
    freeze.write_text(
        "H4 promotion is 0/16\n"
        "The broad Palearctic primary response survives the finite predeclared MNAR trait-resolution grid\n"
        "Chapter 2 / `izu-core` is reserved for **how and why functionally**\n",
        encoding="utf-8",
    )
    maximal = tmp_path / "maximal.md"
    maximal.write_text("- GloBI breadth: 0/4 promoted.\n", encoding="utf-8")
    return {
        "figure3_lock": fig3,
        "figure4_lock": fig4,
        "n1_lock": n1,
        "submission_freeze": freeze,
        "maximal_integration": maximal,
    }


def test_load_inputs_accepts_frozen_v8_contract(tmp_path: Path) -> None:
    paths = _write_sources(tmp_path)
    values = load_inputs(**paths)
    assert values["support_ladder"] == "4/4 -> 4/4 -> 0/4"
    assert values["genus_attenuation_pct"] == [78.8, 85.9]
    assert values["geometry_promoted"] == "0/12"
    assert values["n1_pass"] is False
    assert values["n1_failure_action"] == "stop_before_N2_and_keep_frozen_Chapter1"


def test_load_inputs_rejects_promoted_n1(tmp_path: Path) -> None:
    paths = _write_sources(tmp_path)
    payload = json.loads(paths["n1_lock"].read_text())
    payload["n1"]["gate"]["N1_pass"] = True
    paths["n1_lock"].write_text(json.dumps(payload) + "\n")
    with pytest.raises(FigureInputError, match="failed N1 gate"):
        load_inputs(**paths)


def test_render_figure_writes_submission_surfaces(tmp_path: Path) -> None:
    paths = _write_sources(tmp_path)
    output_dir = tmp_path / "rendered"
    manifest = render_figure(output_dir=output_dir, **paths)
    assert manifest["new_biological_models_fitted"] is False
    assert manifest["new_p_values_generated"] is False
    assert manifest["n2_opened"] is False
    assert manifest["frozen_results"]["h4_promoted"] == "0/16"
    for suffix in ("png", "svg", "pdf"):
        assert (output_dir / f"chapter1_v8_figure1_hierarchical_syndrome.{suffix}").is_file()
    assert (output_dir / "chapter1_v8_figure1_manifest.json").is_file()
