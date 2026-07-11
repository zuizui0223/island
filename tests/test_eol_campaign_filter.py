from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from island_v2.eol_campaign_filter import filter_eol_candidates_to_campaign


def test_filter_eol_candidates_to_campaign_preserves_exclusion_audit(tmp_path: Path) -> None:
    ledger_path = tmp_path / "campaign_ledger.csv.gz"
    candidate_path = tmp_path / "trait_candidates.csv.gz"
    output_dir = tmp_path / "output"

    pd.DataFrame(
        [
            {"accepted_species": "Campanula microdonta", "family": "Campanulaceae"},
            {"accepted_species": "Poa annua", "family": "Poaceae"},
        ]
    ).to_csv(ledger_path, index=False, compression="gzip")

    pd.DataFrame(
        [
            {
                "accepted_species": "Campanula microdonta",
                "trait_name": "flower_primary_color",
                "standardized_value": "white",
                "source_record_id": "eol-1",
            },
            {
                "accepted_species": "Poa annua",
                "trait_name": "pollen_vector_mode",
                "standardized_value": "abiotic_wind",
                "source_record_id": "eol-2",
            },
            {
                "accepted_species": "Cycas revoluta",
                "trait_name": "pollen_vector_mode",
                "standardized_value": "biotic",
                "source_record_id": "eol-3",
            },
        ]
    ).to_csv(candidate_path, index=False, compression="gzip")

    report = filter_eol_candidates_to_campaign(
        campaign_ledger=ledger_path,
        candidate_csv=candidate_path,
        output_dir=output_dir,
    )

    retained = pd.read_csv(report["retained_candidates_csv"], dtype=str).fillna("")
    excluded = pd.read_csv(report["excluded_candidates_csv"], dtype=str).fillna("")
    coverage = pd.read_csv(report["coverage_by_trait_csv"], dtype=str).fillna("")
    summary = json.loads(Path(report["summary_json"]).read_text(encoding="utf-8"))

    assert set(retained["accepted_species"]) == {"Campanula microdonta", "Poa annua"}
    assert excluded["accepted_species"].tolist() == ["Cycas revoluta"]
    assert summary["n_campaign_species"] == 2
    assert summary["n_retained_candidate_rows"] == 2
    assert summary["n_excluded_candidate_rows"] == 1

    flower = coverage.loc[coverage["trait_name"].eq("flower_primary_color")].iloc[0]
    assert int(flower["n_species_with_direct_evidence"]) == 1
    assert float(flower["pct_campaign_species_covered"]) == 50.0
