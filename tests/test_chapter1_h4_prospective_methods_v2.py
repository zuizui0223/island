from __future__ import annotations

import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path("scripts").resolve()))
import screen_chapter1_h4_prospective_methods_v2 as v2  # noqa: E402


def _config() -> dict:
    return {
        "section_policy": {
            "methods_include_regex": [
                "method",
                "study site",
                "study species",
                "study system",
                "experimental",
            ],
            "hard_exclude_regex": ["result", "discussion", "conclusion"],
            "target_species_section_regex": [
                "study species",
                "focal species",
                "study system",
            ],
        },
        "design_flags": {
            "supplementation_terms": ["pollen supplement", "hand pollin"],
            "natural_control_terms": ["open pollin", "natural pollin"],
            "field_terms": ["field site", "natural population"],
            "exclusion_terms": ["greenhouse", "crop"],
        },
    }


def test_title_species_has_priority_over_broad_methods_species() -> None:
    xml = ET.fromstring(
        """
        <article><body>
          <sec><title>Methods</title>
            <p>Vegetation around the field site included Quercus crassipes.</p>
            <sec><title>Study species</title>
              <p>We studied Salvia elegans in natural populations.</p>
            </sec>
          </sec>
        </body></article>
        """
    )
    sections = v2._selected_sections(xml, _config())
    lookup = {
        "salvia elegans": "Salvia elegans",
        "quercus crassipes": "Quercus crassipes",
    }
    target, source, broad = v2._target_species(
        article_title="Pollen limitation in Salvia elegans",
        sections=sections,
        lookup=lookup,
        config=_config(),
    )
    assert target == ["Salvia elegans"]
    assert source == "article_title"
    assert "Quercus crassipes" in broad


def test_target_species_section_used_when_title_has_no_binomial() -> None:
    xml = ET.fromstring(
        """
        <article><body>
          <sec><title>Methods</title>
            <sec><title>Study species</title>
              <p>The focal taxon was Phacelia secunda.</p>
            </sec>
          </sec>
        </body></article>
        """
    )
    sections = v2._selected_sections(xml, _config())
    target, source, _ = v2._target_species(
        article_title="Long-lived flowers and pollen limitation",
        sections=sections,
        lookup={"phacelia secunda": "Phacelia secunda"},
        config=_config(),
    )
    assert target == ["Phacelia secunda"]
    assert source == "target_species_section"


def test_nested_results_are_not_included_in_methods_hash_text() -> None:
    xml = ET.fromstring(
        """
        <article><body>
          <sec><title>Methods</title>
            <p>Open pollination and hand pollination were performed at a field site.</p>
            <sec><title>Results</title>
              <p>SECRET_RESULT_VALUE 0.91.</p>
            </sec>
          </sec>
        </body></article>
        """
    )
    sections = v2._selected_sections(xml, _config())
    joined = "\n".join(text for _, text in sections)
    assert "hand pollination" in joined
    assert "SECRET_RESULT_VALUE" not in joined


def test_normalise_frame_does_not_expand_sampling_frame() -> None:
    frame = pd.DataFrame(
        {
            "source_database": ["OpenAlex", "Crossref"],
            "source_record_id": ["a", "b"],
            "doi": ["10.1/A", ""],
            "title": ["A", "B"],
            "publication_date": ["2020-01-01", "2021-01-01"],
            "publication_year": ["2020", "2021"],
            "query_family": ["pollen limitation", "pollen supplementation"],
        }
    )
    out = v2._normalise_frame(frame)
    assert len(out) == 2
    assert out.loc[0, "doi"] == "10.1/a"


def test_outcome_exposure_exclusions_are_mechanical() -> None:
    frame = pd.DataFrame(
        {
            "source_database": ["OpenAlex", "OpenAlex"],
            "source_record_id": ["a", "b"],
            "doi": ["10.1002/ece3.6884", "10.0000/keep"],
            "title": [
                "Pollinator dependence but no pollen limitation for eight plants occurring north of the Arctic Circle",
                "Unexposed study",
            ],
            "publication_date": ["2020-01-01", "2021-01-01"],
            "publication_year": ["2020", "2021"],
            "query_family": ["pollen limitation", "pollen limitation"],
        }
    )
    frame = v2._normalise_frame(frame)
    lock = {
        "exclusion_keys": {
            "doi": ["10.1002/ece3.6884"],
            "title_exact": [],
        }
    }
    kept, removed = v2._apply_exposure_exclusions(frame, lock)
    assert removed == 1
    assert kept["doi"].tolist() == ["10.0000/keep"]
