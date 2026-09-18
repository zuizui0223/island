from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest
import yaml


SCRIPT = Path("scripts/preflight_chapter1_h4_pollimcrop_transportability.py")
SPEC = importlib.util.spec_from_file_location("h4_pollimcrop_preflight", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CONFIG = yaml.safe_load(
    Path("config/chapter1_h4_pollimcrop_transportability_preflight_v1.yml").read_text(
        encoding="utf-8"
    )
)


def test_metadata_reader_does_not_read_outcome_columns() -> None:
    csv_bytes = (
        "line;article_code;species;country;PL_effectsize;SUP_FS_mean;OBS_FS_mean\n"
        "1;Study_2020_10.1/x;Malus domestica;Spain;0,9;0,8;0,2\n"
    ).encode()
    frame, schema = MODULE.metadata_only_frame(csv_bytes, CONFIG)
    assert "PL_effectsize" in schema
    assert "PL_effectsize" not in frame.columns
    assert "SUP_FS_mean" not in frame.columns
    assert frame.loc[0, "species"] == "Malus domestica"


def test_exact_binomial_is_strict() -> None:
    assert MODULE.exact_binomial("Malus domestica") == "Malus domestica"
    assert MODULE.exact_binomial("Malus domestica Borkh.") == ""
    assert MODULE.exact_binomial("malus domestica") == ""
    assert MODULE.exact_binomial("Malus") == ""

def test_wrong_dataset_scope_is_rejected_before_support() -> None:
    metadata = pd.DataFrame(
        {
            "article_code": ["study_a"],
            "species": ["Malus domestica"],
            "country": ["Spain"],
        }
    )
    with pytest.raises(Exception, match="metadata scope mismatch"):
        MODULE.validate_expected_scope(metadata, CONFIG)


def test_response_mapping_is_frozen_and_uses_same_csv_locale() -> None:
    mapping = MODULE.load_response_mapping(
        Path("config/chapter1_h4_pollimcrop_response_mapping_v1.yml")
    )
    assert mapping["status"] == "frozen_before_row_level_PolLimCrop_outcome_read"
    assert mapping["source_method_basis"]["response_column"] == "PL_effectsize"
    assert mapping["source_method_basis"]["csv_format"] == CONFIG["source"]["csv_format"]

