from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).parents[1] / "analysis" / "recover_citation_hash_rights_provenance.py"
)
spec = importlib.util.spec_from_file_location(
    "recover_citation_hash_rights_provenance", MODULE_PATH
)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_recovers_multiple_provenance_tables_and_excludes_unresolved_map_rows(
    tmp_path: Path,
) -> None:
    coverage = _write(
        tmp_path / "coverage.csv",
        [
            {
                "source_lineages": (
                    "citation:abc|citation:giftref|citation:0123456789abcdef01234567"
                )
            },
            {
                "source_lineages": (
                    "citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497"
                )
            },
            {"source_lineages": "citation:keil1975:pectis-cylindrica-autogamy"},
        ],
    )
    broad = _write(
        tmp_path / "broad.csv",
        [
            {
                "source_lineage": "citation:abc",
                "source_provider": "moeller_etal_2017_dryad",
                "source_url": "https://doi.org/10.5061/dryad.577q1",
                "source_citation": "Moeller et al. 2017",
            }
        ],
    )
    reviewed_lineages = _write(
        tmp_path / "reviewed_lineages.csv",
        [
            {
                "source_lineage": "citation:giftref",
                "source_provider": "gift_v3_2_direct",
                "source_url": "https://gift.uni-goettingen.de/api/example",
                "source_citation": "Example underlying flora",
            }
        ],
    )

    summary = module.recover(
        coverage, (broad, reviewed_lineages), tmp_path / "out"
    )
    assert summary["distinct_citation_lineages"] == 5
    assert summary["source_provenance_recovered"] == 4
    assert summary["unresolved_citation_lineages"] == 1
    assert summary["unresolved_hex24_lineages"] == 1
    assert summary["unresolved_hex20_lineages"] == 0
    assert summary["unresolved_named_lineages"] == 0
    assert summary["rights_granted"] is False
    assert summary["scientific_database_modified"] is False

    mapping = pd.read_csv(
        tmp_path / "out" / "CITATION_LINEAGE_FAMILY_MAP.csv", dtype=str
    )
    lookup = dict(
        zip(mapping["source_lineage"], mapping["source_family"], strict=True)
    )
    assert lookup["citation:abc"] == "dataset:dryad"
    assert lookup["citation:giftref"] == "publication_citation:giftref"
    assert (
        lookup[
            "citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497"
        ]
        == "publication:chung_etal_2000_heredity"
    )
    assert (
        lookup["citation:keil1975:pectis-cylindrica-autogamy"]
        == "publication:keil_1975_pectis_cylindrica"
    )
    assert "citation:0123456789abcdef01234567" not in lookup

    audit = pd.read_csv(
        tmp_path / "out" / "CITATION_PROVENANCE_AUDIT.csv", dtype=str
    ).fillna("")
    chung = audit.loc[
        audit["source_lineage"].eq(
            "citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497"
        )
    ].iloc[0]
    assert "10.1046/j.1365-2540.2000.00781.x" in chung["source_urls"]


def test_conflicting_reviewed_family_mapping_fails_closed(tmp_path: Path) -> None:
    coverage = _write(
        tmp_path / "coverage.csv", [{"source_lineages": "citation:abc"}]
    )
    first = _write(
        tmp_path / "first.csv",
        [
            {
                "source_lineage": "citation:abc",
                "source_provider": "moeller_etal_2017_dryad",
                "source_url": "https://doi.org/10.5061/dryad.577q1",
                "source_citation": "Moeller et al. 2017",
            }
        ],
    )
    second = _write(
        tmp_path / "second.csv",
        [
            {
                "source_lineage": "citation:abc",
                "source_provider": "flora_of_australia_official",
                "source_url": "https://profiles.ala.org.au/opus/foa",
                "source_citation": "Different source",
            }
        ],
    )

    try:
        module.recover(coverage, (first, second), tmp_path / "out")
    except ValueError as exc:
        assert "conflicting source families" in str(exc)
    else:
        raise AssertionError("conflicting lineage-family mapping should fail closed")
