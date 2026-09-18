from __future__ import annotations

import importlib.util
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import yaml


SCRIPT = Path("scripts/screen_chapter1_h4_prospective_methods.py")
SPEC = importlib.util.spec_from_file_location("h4_methods_screen", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)

CONFIG = yaml.safe_load(
    Path("config/chapter1_h4_prospective_methods_preflight_v1.yml").read_text(
        encoding="utf-8"
    )
)


def test_methods_extraction_drops_nested_results() -> None:
    root = ET.fromstring(
        """
        <article>
          <body>
            <sec>
              <title>Materials and methods</title>
              <p>We hand pollinated flowers in a natural population.</p>
              <sec>
                <title>Study site</title>
                <p>The field site was 35.20 N, 139.40 E.</p>
              </sec>
              <sec>
                <title>Results</title>
                <p>Supplemented flowers produced twice as many seeds.</p>
              </sec>
            </sec>
          </body>
        </article>
        """
    )
    texts = MODULE._section_texts(root, CONFIG)
    joined = "\n".join(texts)
    assert "hand pollinated" in joined
    assert "35.20 N" in joined
    assert "twice as many seeds" not in joined


def test_coordinate_parser_handles_decimal_cardinal_pair() -> None:
    coordinates = MODULE._coordinates("Field site: 35.20 N, 139.40 E.")
    assert coordinates == [{"lat": 35.2, "lon": 139.4}]


def test_species_matching_is_exact_to_frozen_table() -> None:
    lookup = MODULE._species_lookup(
        pd.DataFrame(
            {
                "accepted_species": [
                    "Campanula punctata",
                    "Bombus ardens",
                ]
            }
        )
    )
    assert MODULE._matched_species(
        "We studied Campanula punctata and Campanula unknownensis.",
        lookup,
    ) == ["Campanula punctata"]


def test_config_forbids_methods_and_results_text_outputs() -> None:
    outputs = CONFIG["outputs"]
    assert outputs["methods_text_materialized"] is False
    assert outputs["results_text_materialized"] is False
    assert outputs["effect_sizes_materialized"] is False
