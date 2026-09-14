from pathlib import Path

import pandas as pd

from island_v2.chapter1_v8_figure2_branching import (
    REGIME_CONTEXTS,
    RESPONSES,
    SOURCE_MODES,
    STRATA,
    load_inputs,
    render_figure,
)


def _write_fixture(root: Path) -> None:
    for short, offset in (("all", 0.0), ("direct", 0.01)):
        rows = []
        for stratum in STRATA:
            for context in REGIME_CONTEXTS:
                for response in RESPONSES:
                    sign = -1.0 if context == "tropical" and response == "accessibility_generalization" else 1.0
                    rows.append(
                        {
                            "context_layer": "analysis_regime",
                            "axis_set": "universal_plant_response",
                            "stratum": stratum,
                            "support_tier": "confirmatory",
                            "context": context,
                            "syndrome": response,
                            "distance_slope": sign * (0.05 + offset),
                            "cluster_robust_se": 0.01,
                            "p_value": 0.01,
                            "n_islands": 100,
                            "q_axis_family": 0.02,
                            "axis_supported": True,
                        }
                    )
            for response in RESPONSES:
                rows.append(
                    {
                        "context_layer": "biogeographic_realm",
                        "axis_set": "universal_plant_response",
                        "stratum": stratum,
                        "support_tier": "confirmatory",
                        "context": "Palearctic",
                        "syndrome": response,
                        "distance_slope": 0.07 + offset,
                        "cluster_robust_se": 0.012,
                        "p_value": 0.01,
                        "n_islands": 120,
                        "q_axis_family": 0.02,
                        "axis_supported": True,
                    }
                )
        path = root / f"branching/{short}/global_branch_distance_slopes.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(path, index=False)

    pathway_rows = []
    for mode_index, mode in enumerate(SOURCE_MODES):
        for stratum in STRATA:
            pathway_rows.append(
                {
                    "context_layer": "biogeographic_realm",
                    "source_mode": mode,
                    "stratum": stratum,
                    "context": "Palearctic",
                    "support_tier": "confirmatory",
                    "threshold": 50,
                    "model": "attraction_conditional_on_selfing_core",
                    "status": "fit",
                    "n_unique_islands": 120,
                    "n_clusters": 20,
                    "distance_estimate": 0.08 + 0.005 * mode_index,
                    "distance_se": 0.02,
                    "distance_p": 0.01,
                    "selfing_core_estimate": 0.01,
                    "selfing_core_se": 0.01,
                    "selfing_core_p": 0.2,
                    "distance_q": 0.02,
                }
            )
    path = root / "source-adjusted/source_adjusted_pathway_realm.csv"
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(pathway_rows).to_csv(path, index=False)


def test_load_inputs_fail_closed_design(tmp_path: Path) -> None:
    _write_fixture(tmp_path)
    tables = load_inputs(tmp_path)
    assert len(tables["regime"]) == 16
    assert len(tables["palearctic"]) == 8
    assert len(tables["pathway"]) == 8


def test_render_v8_figure2(tmp_path: Path) -> None:
    artifact = tmp_path / "artifact"
    output = tmp_path / "output"
    _write_fixture(artifact)
    manifest = render_figure(
        artifact_root=artifact,
        output_dir=output,
        source_run_id=1,
        source_artifact_id=2,
        source_artifact_digest="sha256:test",
    )
    assert manifest["contract"] == "chapter1_v8_figure2_biogeographic_branching_v1"
    assert manifest["new_biological_models_fitted"] is False
    assert manifest["n_regime_primary_rows"] == 16
    assert (output / "chapter1_v8_figure2_biogeographic_branching.png").is_file()
    assert (output / "chapter1_v8_figure2_biogeographic_branching.svg").is_file()
    assert (output / "chapter1_v8_figure2_biogeographic_branching.pdf").is_file()
