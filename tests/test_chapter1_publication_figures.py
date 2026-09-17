from __future__ import annotations

import json
import os
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest

from island_v2.chapter1_publication_data import (
    PublicationInputs,
    build_publication_source_data,
    publication_text_is_clean,
)
from island_v2.chapter1_publication_figures import PUBLICATION_LABELS, render_publication_bundle


@pytest.fixture(scope="module")
def input_root() -> Path:
    value = os.environ.get("CH1_PUBLICATION_INPUT_ROOT")
    if not value:
        pytest.skip("CH1_PUBLICATION_INPUT_ROOT is required for artifact integration tests")
    root = Path(value)
    assert root.exists()
    return root


@pytest.fixture(scope="module")
def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def publication_inputs(input_root: Path, repo_root: Path) -> PublicationInputs:
    return PublicationInputs(
        isolation_dir=input_root / "isolation",
        atomic_dir=input_root / "atomic",
        glopl_overlap_dir=input_root / "glopl_overlap",
        glopl_global_dir=input_root / "glopl_global",
        glopl_shape_dir=input_root / "glopl_shape",
        functional_bridge_dir=input_root / "functional_bridge",
        unified_lock=repo_root / "config/chapter1_v13_unified_island_syndrome_result_lock.json",
        functional_lock=repo_root / "config/chapter1_v13_functional_bridge_result_lock.json",
    )


def test_publication_labels_exclude_internal_repository_jargon() -> None:
    for text in PUBLICATION_LABELS:
        assert publication_text_is_clean(text), text


def test_source_data_reproduce_frozen_support_counts(
    publication_inputs: PublicationInputs, tmp_path: Path
) -> None:
    paths = build_publication_source_data(publication_inputs, tmp_path)

    map_islands = pd.read_csv(paths["figure1_islands"])
    glopl_sites = pd.read_csv(paths["figure1_glopl_sites"])
    tested = pd.read_csv(paths["extended_data_table6"])
    atomic = pd.read_csv(paths["figure2_atomic"])
    omnibus = pd.read_csv(paths["figure2_omnibus"])

    assert map_islands["island_id"].nunique() == 8265
    assert map_islands.loc[map_islands["trait_informed"], "island_id"].nunique() == 4453
    assert glopl_sites["site_key"].nunique() == 1248
    assert glopl_sites.loc[glopl_sites["is_frozen_island_site"], "site_key"].nunique() == 197
    assert tested["island_id"].nunique() == 37
    assert tested["trait_informed"].all()

    primary = omnibus.loc[omnibus["evidence_scope"].eq("all_analysis")]
    direct = omnibus.loc[omnibus["evidence_scope"].eq("direct_only")]
    assert int(primary["n_unique_islands"].sum()) == 4334
    assert int(direct["n_unique_islands"].sum()) == 4256

    assert len(atomic) == 48
    assert set(atomic["evidence_scope"]) == {"all_analysis", "direct_only"}
    assert atomic.groupby(["evidence_scope", "context"]).size().eq(6).all()


def test_figure3_source_data_match_frozen_locks(
    publication_inputs: PublicationInputs, tmp_path: Path
) -> None:
    paths = build_publication_source_data(publication_inputs, tmp_path)
    figure3 = pd.read_csv(paths["figure3"])

    global_row = figure3.loc[figure3["row_id"].eq("global_distance")].iloc[0]
    assert global_row["estimate"] == pytest.approx(0.07936837202478536)
    assert global_row["se"] == pytest.approx(0.037734195473258715)
    assert global_row["two_sided_p"] == pytest.approx(0.035434833922249484)

    traits = figure3.set_index("row_id")
    assert traits.loc["trait_autonomous_selfing_primary", "estimate"] == pytest.approx(
        -0.446724, abs=1e-6
    )
    assert traits.loc["trait_actinomorphic_symmetry_primary", "estimate"] == pytest.approx(
        -0.381194, abs=1e-6
    )
    assert traits.loc["trait_generalized_form_primary", "estimate"] == pytest.approx(
        -0.184104, abs=1e-6
    )
    assert traits.loc["trait_self_compatibility_primary", "estimate"] == pytest.approx(
        -0.117727, abs=1e-6
    )


def test_publication_bundle_contains_main_and_extended_outputs(
    publication_inputs: PublicationInputs, tmp_path: Path
) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "bundle"
    build_publication_source_data(publication_inputs, source_dir)
    outputs = render_publication_bundle(source_dir, output_dir)

    expected_stems = {
        "figure1_global_coverage",
        "figure2_recurrent_plant_response",
        "figure3_pollen_limitation_bridge",
        "extended_data_figure1_atomic_support",
        "extended_data_figure2_sampling_overlap",
        "extended_data_figure3_glopl_sensitivity",
        "extended_data_figure4_functional_bridge_sensitivity",
    }
    for stem in expected_stems:
        for suffix in (".pdf", ".svg", ".png"):
            path = output_dir / ("main" if stem.startswith("figure") else "extended_data") / f"{stem}{suffix}"
            assert path in outputs
            assert path.exists()
            assert path.stat().st_size > 1000

    manifest = json.loads((source_dir / "publication_figure_manifest.json").read_text())
    assert manifest["counts"]["frozen_islands"] == 8265
    assert manifest["counts"]["trait_informed_islands"] == 4453
    assert manifest["counts"]["glopl_sites"] == 1248
    assert manifest["counts"]["glopl_tested_islands"] == 37


def test_svg_publication_text_has_no_internal_labels(
    publication_inputs: PublicationInputs, tmp_path: Path
) -> None:
    source_dir = tmp_path / "source"
    output_dir = tmp_path / "bundle"
    build_publication_source_data(publication_inputs, source_dir)
    render_publication_bundle(source_dir, output_dir)

    for svg in output_dir.rglob("*.svg"):
        root = ET.parse(svg).getroot()
        visible_text = " ".join(
            element.text or ""
            for element in root.iter()
            if element.tag.rsplit("}", 1)[-1] == "text"
        ).lower()
        assert "v13" not in visible_text
        assert "workflow" not in visible_text
        assert "artifact" not in visible_text
        assert "main figure" not in visible_text
