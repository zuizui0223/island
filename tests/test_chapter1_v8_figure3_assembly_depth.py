from pathlib import Path

import pandas as pd

from island_v2.chapter1_v8_figure3_assembly_depth import (
    AXES,
    CANONICAL_SOURCE_MODE,
    SCOPES,
    SOURCE_MODES,
    STAGES,
    STRATA,
    load_inputs,
    render_figure,
)


def _write_fixture(effect_root: Path, pr142_root: Path) -> None:
    attenuation_rows = []
    for scope in SCOPES:
        for stratum in STRATA:
            for index, mode in enumerate(SOURCE_MODES):
                observed = 0.05 + 0.001 * index
                family = observed * 0.72
                genus = observed * 0.18
                attenuation_rows.append(
                    {
                        "evidence_scope": scope,
                        "source_mode": mode,
                        "stratum": stratum,
                        "context": "Palearctic",
                        "observed_score_vector_norm": observed,
                        "after_family_residual_vector_norm": family,
                        "after_genus_residual_vector_norm": genus,
                        "family_attenuation_fraction": 0.28,
                        "genus_attenuation_fraction": 0.82,
                        "conditional_genus_attenuation_fraction": 0.75,
                    }
                )
    effect_root.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(attenuation_rows).to_csv(
        effect_root / "taxonomic_vector_attenuation.csv", index=False
    )

    effects_rows = []
    for scope in SCOPES:
        for stratum in STRATA:
            for stage_index, stage in enumerate(STAGES):
                for axis in AXES:
                    value = [0.04, 0.025, 0.006][stage_index]
                    effects_rows.append(
                        {
                            "evidence_scope": scope,
                            "source_mode": CANONICAL_SOURCE_MODE,
                            "stratum": stratum,
                            "context": "Palearctic",
                            "stage": stage,
                            "axis": axis,
                            "distance_slope": value,
                            "ci_low": value - 0.01,
                            "ci_high": value + 0.01,
                        }
                    )
    pd.DataFrame(effects_rows).to_csv(effect_root / "taxonomic_effects_long.csv", index=False)

    for scope in SCOPES:
        rows = []
        for stratum in STRATA:
            rows.append(
                {
                    "context_layer": "biogeographic_realm",
                    "stratum": stratum,
                    "context": "Palearctic",
                    "n_source_modes_expected": 4,
                    "observed_supported_modes": 4,
                    "after_family_supported_modes": 4,
                    "after_genus_supported_modes": 0,
                    "V2_taxonomic_depth_classification": "compatible_with_genus_level_assembly_beyond_family",
                }
            )
        path = pr142_root / f"taxonomic-depth/{scope}/classification.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(rows).to_csv(path, index=False)


def test_load_frozen_assembly_depth_design(tmp_path: Path) -> None:
    effect = tmp_path / "effect"
    primary = tmp_path / "primary"
    _write_fixture(effect, primary)
    tables = load_inputs(effect, primary)
    assert len(tables["attenuation"]) == 16
    assert len(tables["effects"]) == 24
    assert len(tables["classification"]) == 4


def test_render_v8_figure3(tmp_path: Path) -> None:
    effect = tmp_path / "effect"
    primary = tmp_path / "primary"
    output = tmp_path / "output"
    _write_fixture(effect, primary)
    manifest = render_figure(
        effect_root=effect,
        pr142_root=primary,
        output_dir=output,
        effect_run_id=1,
        effect_artifact_id=2,
        effect_artifact_digest="sha256:effect",
        primary_run_id=3,
        primary_artifact_id=4,
        primary_artifact_digest="sha256:primary",
    )
    assert manifest["contract"] == "chapter1_v8_figure3_assembly_depth_v1"
    assert manifest["frozen_support_ladder"] == "4/4 -> 4/4 -> 0/4"
    assert manifest["new_biological_models_fitted"] is False
    assert (output / "chapter1_v8_figure3_assembly_depth.png").is_file()
    assert (output / "chapter1_v8_figure3_assembly_depth.svg").is_file()
    assert (output / "chapter1_v8_figure3_assembly_depth.pdf").is_file()
