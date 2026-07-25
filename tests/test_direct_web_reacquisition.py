from __future__ import annotations

import pandas as pd

from island_v2.direct_web_reacquisition import (
    _direct_evidence,
    _shared_full_text_sources,
    _wfo_aliases_from_payloads,
    combine_direct_web_evidence,
    select_reproductive_targets,
)
from island_v2.v1_category_search import extract_evidence_from_sources


def test_select_targets_requires_two_filled_axes_and_missing_reproduction() -> None:
    coverage = pd.DataFrame(
        [
            {"accepted_species": species, "axis": axis, "after_quality": quality}
            for species, qualities in {
                "Plantus alpha": ("high", "medium", ""),
                "Plantus beta": ("low", "", ""),
                "Plantus gamma": ("high", "medium", "low"),
                "Genus": ("high", "medium", ""),
            }.items()
            for axis, quality in zip(
                (
                    "flower_colour",
                    "floral_structural_complexity",
                    "reproductive_assurance",
                ),
                qualities,
                strict=True,
            )
        ]
    )

    targets = select_reproductive_targets(coverage, max_species=10)

    assert targets["accepted_species"].tolist() == ["Plantus alpha"]
    assert targets.loc[0, "existing_quality_score"] == 5


def test_full_text_candidates_become_medium_species_direct_evidence() -> None:
    full_text = pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alpha",
                "source_text": (
                    "Plantus alpha has a mixed mating system and is self-compatible."
                ),
                "source_url": (
                    "https://www.ebi.ac.uk/europepmc/webservices/rest/"
                    "PMC123/fullTextXML"
                ),
                "source_citation": "Author. DOI:10.1000/example",
                "source_type": "europe_pmc_full_text_xml",
                "evidence_scope": "species_direct",
                "pmcid": "PMC123",
                "source_record_id": "PMC123:passage",
                "source_lineage": "doi:10.1000/example",
            }
        ]
    )

    extracted = extract_evidence_from_sources(_shared_full_text_sources(full_text))
    evidence = _direct_evidence(extracted, full_text)

    assert set(evidence["trait_name"]) == {"mating_system", "self_incompatibility"}
    assert set(evidence["evidence_quality"]) == {"medium"}
    assert set(evidence["evidence_scope"]) == {"species_direct"}
    assert set(evidence["axis"]) == {"reproductive_assurance"}
    assert set(evidence["source_lineage"]) == {"doi:10.1000/example"}


def test_wfo_alias_parser_keeps_species_rank_and_skips_infraspecific_names() -> None:
    match_payload = {
        "match": {"wfo_id": "wfo-accepted"},
        "parsedName": {"canonical_form": "Plantus alpha"},
    }
    graph_payload = {
        "data": {
            "taxonNameById": {
                "currentPreferredUsage": {
                    "id": "wfo-accepted-2026-06",
                    "hasName": {"id": "wfo-accepted"},
                    "hasSynonym": [
                        {
                            "id": "wfo-synonym",
                            "fullNameStringPlain": "Oldgenus alpha Author",
                        },
                        {
                            "id": "wfo-variety",
                            "fullNameStringPlain": "Othera beta var. minor Author",
                        },
                    ],
                }
            }
        }
    }

    aliases = _wfo_aliases_from_payloads(
        "Plantus alpha",
        match_payload,
        graph_payload,
        max_aliases=20,
    )

    assert [row["search_name"] for row in aliases] == ["Oldgenus alpha"]
    assert aliases[0]["name_resolution_lineage"].endswith(":wfo-synonym")


def test_combine_direct_web_evidence_preserves_denominator_and_precedence(
    tmp_path,
) -> None:
    coverage_path = tmp_path / "coverage.csv"
    evidence_path = tmp_path / "evidence.csv"
    axes = [
        "flower_colour",
        "floral_structural_complexity",
        "reproductive_assurance",
    ]
    pd.DataFrame(
        [
            {
                "accepted_species": species,
                "axis": axis,
                "after_quality": quality,
            }
            for species, qualities in {
                "Plantus alpha": ("high", "medium", ""),
                "Plantus beta": ("", "", ""),
            }.items()
            for axis, quality in zip(axes, qualities, strict=True)
        ]
    ).to_csv(coverage_path, index=False)
    pd.DataFrame(
        [
            {
                "accepted_species": "Plantus alpha",
                "axis": "reproductive_assurance",
                "evidence_quality": "medium",
                "source_lineage": "doi:alpha",
            },
            {
                "accepted_species": "Plantus beta",
                "axis": "floral_structural_complexity",
                "evidence_quality": "medium",
                "source_lineage": "doi:beta",
            },
        ]
    ).to_csv(evidence_path, index=False)

    summary = combine_direct_web_evidence(
        coverage_csv=coverage_path,
        evidence_csvs=[evidence_path],
        output_dir=tmp_path / "combined",
    )

    assert summary["denominator"]["species_axis"] == 6
    assert summary["net_change"]["pure_gain_species_axis"] == 2
    assert summary["after"]["quality_counts"] == {
        "high": 1,
        "medium": 3,
        "low": 0,
    }
    assert summary["after"]["species_by_filled_axis_count"] == {
        "0": 0,
        "1": 1,
        "2": 0,
        "3": 1,
    }
