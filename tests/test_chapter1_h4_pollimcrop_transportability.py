from __future__ import annotations

import io
from pathlib import Path

import pandas as pd
import yaml

from scripts.preflight_chapter1_h4_pollimcrop_transportability import (
    exact_binomial,
    metadata_only_frame,
)


CONFIG = yaml.safe_load(
    Path("config/chapter1_h4_pollimcrop_transportability_preflight_v1.yml").read_text(
        encoding="utf-8"
    )
)


def test_metadata_reader_does_not_read_outcome_columns() -> None:
    csv_bytes = (
        "line,article_code,species,country,PL_effectsize,SUP_FS_mean,OBS_FS_mean\n"
        "1,Study_2020_10.1/x,Malus domestica,Spain,0.9,0.8,0.2\n"
    ).encode()
    frame, schema = metadata_only_frame(csv_bytes, CONFIG)
    assert "PL_effectsize" in schema
    assert "PL_effectsize" not in frame.columns
    assert "SUP_FS_mean" not in frame.columns
    assert frame.loc[0, "species"] == "Malus domestica"


def test_exact_binomial_is_strict() -> None:
    assert exact_binomial("Malus domestica") == "Malus domestica"
    assert exact_binomial("Malus domestica Borkh.") == ""
    assert exact_binomial("malus domestica") == ""
    assert exact_binomial("Malus") == ""
