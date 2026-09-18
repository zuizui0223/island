from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

from island_v2.chapter1_h4_pollimcrop_transportability import (
    admitted_species_for_outcome_read,
    aggregate_analysis_cells,
    exact_binomial,
    fit_two_way_clustered_trait,
    load_preflight,
    read_outcome_rows,
    require_support_before_outcome_read,
    run_transportability,
)


MAPPING = yaml.safe_load(
    Path("config/chapter1_h4_pollimcrop_response_mapping_v1.yml").read_text(
        encoding="utf-8"
    )
)


DOWNLOADER_SCRIPT = Path("scripts/download_chapter1_h4_pollimcrop_locked_dataset.py")
DOWNLOADER_SPEC = importlib.util.spec_from_file_location(
    "h4_pollimcrop_locked_downloader",
    DOWNLOADER_SCRIPT,
)
assert DOWNLOADER_SPEC is not None and DOWNLOADER_SPEC.loader is not None
DOWNLOADER = importlib.util.module_from_spec(DOWNLOADER_SPEC)
DOWNLOADER_SPEC.loader.exec_module(DOWNLOADER)


def _cells(effect: float = -0.6, n_publications: int = 24) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for i in range(n_publications):
        article = f"study_{i:02d}"
        continent = "Europe" if i % 2 else "Asia"
        for state in (0, 1):
            species = f"Species {i:02d}_{state}"
            for replicate in range(2):
                rows.append(
                    {
                        "article_code": article,
                        "species": species,
                        "accepted_species": species,
                        "continent": continent,
                        "country": f"Country{i % 4}",
                        "locality": f"Locality{i}",
                        "experiment_year": 2010 + (i % 10),
                        "supplement_type": "hand_cross",
                        "scale": "flower",
                        "crop_part": "fruit",
                        "PL_effectsize": (
                            0.5
                            + 0.03 * (i % 10)
                            + effect * state
                            + 0.015 * (((i + 2 * state) % 5) - 2)
                            + 0.002 * replicate
                        ),
                        "autonomous_selfing": state,
                        "accessibility_generalization_score": float(state),
                    }
                )
    raw = pd.DataFrame(rows)
    return aggregate_analysis_cells(raw)


def test_exact_binomial_is_strict() -> None:
    assert exact_binomial("Malus domestica") == "Malus domestica"
    assert exact_binomial("Malus domestica Borkh.") == ""
    assert exact_binomial("malus domestica") == ""
    assert exact_binomial("Malus") == ""


def _committed_support_lock() -> dict:
    return {
        "contract": "chapter1_h4_pollimcrop_transportability_preflight_v1",
        "lock_contract": (
            "chapter1_h4_pollimcrop_transportability_preflight_result_lock_v1"
        ),
        "outcomes_read": False,
        "support": {
            "H4a_reproductive_assurance": {"evaluable": True},
            "H4b_accessibility_generalization": {"evaluable": False},
        },
        "decision": {
            "admitted_hypotheses": ["H4a_reproductive_assurance"],
            "support_failed_hypotheses": ["H4b_accessibility_generalization"],
            "outcome_extraction_authorized": True,
        },
    }


def test_load_preflight_requires_committed_support_result_lock(tmp_path: Path) -> None:
    value = _committed_support_lock()
    value.pop("lock_contract")
    path = tmp_path / "raw_preflight.json"
    path.write_text(json.dumps(value), encoding="utf-8")
    with pytest.raises(Exception, match="requires the committed support-result lock"):
        load_preflight(path)


def test_load_preflight_accepts_consistent_committed_support_lock(tmp_path: Path) -> None:
    path = tmp_path / "support_lock.json"
    path.write_text(json.dumps(_committed_support_lock()), encoding="utf-8")
    got = load_preflight(path)
    assert got["decision"]["admitted_hypotheses"] == ["H4a_reproductive_assurance"]


def test_downloader_requires_committed_support_result_lock() -> None:
    value = _committed_support_lock()
    value.pop("lock_contract")
    with pytest.raises(Exception, match="committed support-result lock"):
        DOWNLOADER.validate_support_lock(value)


def test_downloader_rejects_inconsistent_admission_list() -> None:
    value = _committed_support_lock()
    value["decision"]["admitted_hypotheses"] = ["H4b_accessibility_generalization"]
    with pytest.raises(Exception, match="admission list is inconsistent"):
        DOWNLOADER.validate_support_lock(value)


def test_no_support_refuses_before_outcome_read() -> None:
    preflight = {
        "support": {
            "H4a_reproductive_assurance": {"evaluable": False},
            "H4b_accessibility_generalization": {"evaluable": False},
        }
    }
    with pytest.raises(Exception, match="outcome values must remain unread"):
        require_support_before_outcome_read(preflight)


def test_outcome_reader_uses_frozen_semicolon_comma_locale(tmp_path: Path) -> None:
    path = tmp_path / "pollimcrop.csv"
    path.write_text(
        "article_code;species;continent;country;locality;experiment_year;"
        "supplement_type;scale;crop_part;PL_effectsize\n"
        "study_a;Malus domestica;Europe;Spain;site;2020;H;flower;fruit;0,42\n"
        "study_b;Prunus avium;Europe;Spain;site;2021;H;flower;fruit;9,99\n",
        encoding="utf-8-sig",
    )
    preflight = {
        "support": {
            "H4a_reproductive_assurance": {"evaluable": True},
            "H4b_accessibility_generalization": {"evaluable": False},
        },
        "figshare": {"csv_sha256": hashlib.sha256(path.read_bytes()).hexdigest()},
        "source_format": {
            "delimiter": ";",
            "decimal_mark": ",",
            "encoding": "utf-8-sig",
        },
    }
    rows = read_outcome_rows(
        path,
        preflight,
        allowed_species={"Malus domestica"},
    )
    assert rows["accepted_species"].tolist() == ["Malus domestica"]
    assert rows["PL_effectsize"].tolist() == pytest.approx([0.42])


def test_two_way_cluster_fit_recovers_negative_trait_effect() -> None:
    cells = _cells(effect=-0.6)
    result = fit_two_way_clustered_trait(cells, "autonomous_selfing")
    assert result["evaluable"] is True
    assert result["estimate"] < -0.5
    assert result["one_sided_negative_p"] < 0.025
    assert result["n_publications"] == 24
    assert result["n_species"] == 48


def test_publication_weights_recomputed_after_predictor_filtering() -> None:
    cells = _cells(effect=-0.5, n_publications=12)
    # Delete one state from half the publications for H4b only.
    mask = (
        cells["article_code"].isin([f"study_{i:02d}" for i in range(6)])
        & cells["accessibility_generalization_score"].eq(1.0)
    )
    cells.loc[mask, "accessibility_generalization_score"] = np.nan
    result = fit_two_way_clustered_trait(
        cells,
        "accessibility_generalization_score",
    )
    assert result["evaluable"] is True
    assert result["estimate"] < 0


def test_run_transportability_does_not_use_failed_hypothesis() -> None:
    cells = _cells(effect=-0.6)
    preflight = {
        "support": {
            "H4a_reproductive_assurance": {"evaluable": True},
            "H4b_accessibility_generalization": {"evaluable": False},
        }
    }
    result = run_transportability(cells, preflight, MAPPING)
    assert result["primary_results"]["H4a_reproductive_assurance"]["evaluable"]
    assert result["primary_results"]["H4b_accessibility_generalization"] == {
        "evaluable": False,
        "reason": "frozen_preflight_support_gate_failed",
        "outcomes_used": False,
    }


def test_outcome_species_union_uses_only_supported_hypotheses() -> None:
    traits = pd.DataFrame(
        {
            "accepted_species": ["Species one", "Species two", "Species three"],
            "autonomous_selfing": [1, np.nan, 0],
            "generalized_form": [np.nan, 1, 0],
            "actinomorphic_symmetry": [np.nan, 1, 0],
            "shallow_open_tube": [np.nan, 1, 0],
        }
    )
    parent = yaml.safe_load(
        Path("config/chapter1_h4_prospective_temporal_replication_v1.yml").read_text(
            encoding="utf-8"
        )
    )
    preflight = {
        "support": {
            "H4a_reproductive_assurance": {"evaluable": True},
            "H4b_accessibility_generalization": {"evaluable": False},
        }
    }
    allowed = admitted_species_for_outcome_read(traits, parent, preflight)
    assert allowed == {"Species one", "Species three"}
    assert "Species two" not in allowed
