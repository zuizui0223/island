from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd


MODULE_PATH = Path(__file__).parents[1] / "analysis" / "recover_citation_hash_rights_provenance.py"
spec = importlib.util.spec_from_file_location("recover_citation_hash_rights_provenance", MODULE_PATH)
assert spec is not None and spec.loader is not None
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _write(path: Path, rows: list[dict[str, object]]) -> Path:
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def test_recovers_citation_families_and_manual_receipts(tmp_path: Path) -> None:
    coverage = _write(
        tmp_path / "coverage.csv",
        [
            {"source_lineages": "citation:abc"},
            {"source_lineages": "citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497"},
            {"source_lineages": "citation:Fischer-Rahelivololona-2007-Adansonia-29-269-315"},
        ],
    )
    reviewed = _write(
        tmp_path / "reviewed.csv",
        [
            {
                "source_lineage": "citation:abc",
                "source_provider": "moeller_etal_2017_dryad",
                "source_url": "https://doi.org/10.5061/dryad.577q1",
                "source_citation": "Moeller et al. 2017",
            }
        ],
    )

    summary = module.recover(coverage, reviewed, tmp_path / "out")
    assert summary["distinct_citation_lineages"] == 3
    assert summary["source_provenance_recovered"] == 3
    assert summary["unresolved_citation_lineages"] == 0
    assert summary["rights_granted"] is False
    assert summary["scientific_database_modified"] is False

    mapping = pd.read_csv(tmp_path / "out" / "CITATION_LINEAGE_FAMILY_MAP.csv", dtype=str)
    lookup = dict(zip(mapping["source_lineage"], mapping["source_family"], strict=True))
    assert lookup["citation:abc"] == "dataset:dryad"
    assert lookup["citation:Chung-Chung-Oh-Epperson-2000-Heredity-85-490-497"] == "publication:chung_etal_2000_heredity"
    assert lookup["citation:Fischer-Rahelivololona-2007-Adansonia-29-269-315"] == "publication:fischer_rahelivololona_2007_adansonia"
